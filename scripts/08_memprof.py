"""Memory profile of a command and all its child processes (Linux /proc; no extra dependencies).

    python scripts/08_memprof.py OUT_PREFIX [--interval S] -- CMD [ARGS...]

Every S seconds (default 1) it walks the process tree and reads /proc/<pid>/smaps_rollup:
  pss_gb      sum of proportional set size: shared pages split between the processes that map them, so forked
              workers sharing the parent's memory are counted once. This is the footprint to plan around.
  rss_sum_gb  sum of RSS: counts shared pages once per process (what Slurm's MaxRSS reports; an overestimate)
  rss_max_gb  largest single process RSS (what /usr/bin/time -v reports; an underestimate for multi-process tools)
Writes OUT_PREFIX.mem.tsv (timeline) and OUT_PREFIX.mem.json (peaks, wall time, exit code, cgroup memory.peak if
readable; memory.peak also counts page cache, so it is an upper bound). Exits with the command's exit code.
"""
import argparse
import json
import os
import subprocess
import sys
import time


def children(pid):
    out = []
    try:
        for tid in os.listdir(f"/proc/{pid}/task"):
            with open(f"/proc/{pid}/task/{tid}/children") as fh:
                out += [int(c) for c in fh.read().split()]
    except OSError:
        pass
    return out


def tree(root):
    seen, stack = [], [root]
    while stack:
        p = stack.pop()
        seen.append(p)
        stack += children(p)
    return seen


def rollup(pid):
    v = {}
    try:
        with open(f"/proc/{pid}/smaps_rollup") as fh:
            for line in fh:
                k, _, rest = line.partition(":")
                if k in ("Rss", "Pss"):
                    v[k] = int(rest.split()[0])
    except OSError:
        return None
    return v


def cgroup_peak():
    try:
        rel = open("/proc/self/cgroup").read().strip().split("::", 1)[1]
        d = "/sys/fs/cgroup" + rel
        # walk up to the Slurm job cgroup (the step cgroup may not expose memory.peak)
        while d.startswith("/sys/fs/cgroup") and len(d) > len("/sys/fs/cgroup"):
            f = os.path.join(d, "memory.peak")
            if os.path.exists(f) and ("job_" in d or "step_" in d):
                return d, int(open(f).read()) / 2**30
            d = os.path.dirname(d)
    except (OSError, IndexError, ValueError):
        pass
    return None, None


def main():
    if "--" not in sys.argv:
        sys.exit(__doc__)
    i = sys.argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--interval", type=float, default=1.0)
    a = ap.parse_args(sys.argv[1:i])
    cmd = sys.argv[i + 1:]
    t0 = time.time()
    proc = subprocess.Popen(cmd)
    peak = {"pss_gb": 0.0, "rss_sum_gb": 0.0, "rss_max_gb": 0.0, "nproc": 0}
    peak_t = {}
    with open(f"{a.out}.mem.tsv", "w") as tsv:
        tsv.write("t_s\tnproc\tpss_gb\trss_sum_gb\trss_max_gb\n")
        while proc.poll() is None:
            rs = [r for r in (rollup(p) for p in tree(proc.pid)) if r]
            if rs:
                row = {"nproc": len(rs),
                       "pss_gb": sum(r.get("Pss", 0) for r in rs) / 2**20,
                       "rss_sum_gb": sum(r.get("Rss", 0) for r in rs) / 2**20,
                       "rss_max_gb": max(r.get("Rss", 0) for r in rs) / 2**20}
                t = time.time() - t0
                tsv.write(f"{t:.1f}\t{row['nproc']}\t{row['pss_gb']:.3f}\t{row['rss_sum_gb']:.3f}\t{row['rss_max_gb']:.3f}\n")
                tsv.flush()
                for k, v in row.items():
                    if v > peak[k]:
                        peak[k], peak_t[k] = v, round(t, 1)
            time.sleep(a.interval)
    cg_dir, cg_peak = cgroup_peak()
    summary = {"cmd": cmd, "exit_code": proc.returncode, "wall_s": round(time.time() - t0, 1),
               **{f"peak_{k}": round(v, 3) for k, v in peak.items()},
               **{f"t_peak_{k}": v for k, v in peak_t.items()},
               "cgroup_memory_peak_gb": round(cg_peak, 3) if cg_peak else None, "cgroup": cg_dir}
    json.dump(summary, open(f"{a.out}.mem.json", "w"), indent=1)
    print(json.dumps(summary), file=sys.stderr)
    sys.exit(proc.returncode)


if __name__ == "__main__":
    main()

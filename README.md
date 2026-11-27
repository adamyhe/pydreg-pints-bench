# pydreg-pints-bench

A benchmark of two callers of transcribed regulatory elements on nascent-transcription data:

* [pydreg](https://pypi.org/project/pydreg/) 1.0.0, a Python port of dREG;
* [PINTS](https://pypi.org/project/pyPINTS/) (pyPINTS) 1.3.0.

The benchmark covers 15 libraries in 5 cell types (K562, GM12878, HCT116, A549, Jurkat):

* 3′ run-on assays: PRO-seq, GRO-seq, ChRO-seq, leChRO-seq;
* TSS assays: PRO-cap, GRO-cap, CAGE.

Both tools are scored with one symmetric, element-level protocol against several families of labels:

* the PINTS paper's labels;
* the dREG paper's test sets;
* H3K27ac/DHS labels;
* PINTS calls on a matched TSS library;
* accessible but TSS-silent DHSs.

Intervals come from a paired 5 Mb block bootstrap, with sensitivity checks for GC-matched negatives and for CpG-island
vs non-CpG-island positives.

## Results

* **[RESULTS_SUMMARY.md](RESULTS_SUMMARY.md)** has the consolidated findings and the run-time table.
* **[docs/pydreg_vs_pints.html](docs/pydreg_vs_pints.html)** is the same results as an interactive page with
  per-library charts. It is self-contained and can be served from any static host.
* **[results/summary/report.md](results/summary/report.md)** has the full per-library tables: every label set,
  condition and metric, with intervals.
* [REVIEW.md](REVIEW.md) reviews the PINTS paper's comparison. [BENCHMARK.md](BENCHMARK.md) gives the design and
  label definitions. RESULTS_PHASE1/2.md are earlier write-ups, superseded by the summary.

In short:

* **3′ run-on data:** pydreg is better at strict false-positive rates, especially on enhancer-like
  (non-CpG-island) elements.
* **TSS assays:** PINTS is better.
* **Transcribed vs untranscribed DHSs:** the tools are close, with PINTS modestly ahead on enhancer-like sites.
* **Speed and memory:** PINTS is faster overall. pydreg needs a GPU but far less memory.

## Layout

```
config/libraries.tsv     libraries: phase, cell, assay, genome, dREG role, PINTS exp type, bigWig sources
config/site.env.example  site settings (copy to config/site.env)
scripts/00*_fetch*.sh    downloads, run on a login node (internet needed)
scripts/01*_build_labels*.py, 03_evaluate.py, 04_summarize.py, 05_bootstrap.py, 06_report.py, 07_gc_diagnostic.py
slurm/                   one sbatch script per stage, sub.sh (submit wrapper), submit_*.sh (dependency chains)
src/ppbench/             library code: labels, calls, evaluate, bootstrap, gc, genome
tests/                   unit tests
results/                 eval/ (point metrics), bootstrap/*/ci.tsv (intervals), summary/ (tables and report)
```

Not in the repo, rebuilt by the scripts:

* `data/`: libraries, references and built labels;
* `results/calls/`: caller output, about 11 GB;
* `results/bootstrap/*/element_scores.tsv.gz`.

## Requirements

* [uv](https://docs.astral.sh/uv/) and Python ≥ 3.12. Run `uv sync` to install pydreg[gpu], pyPINTS and the
  analysis stack into `.venv/`.
* An NVIDIA GPU with CUDA for pydreg's CuPy backend. The original runs used one Tesla P100 per library.
* About 64 GB of RAM for PINTS on the deepest libraries. It peaked at 58 GB on K562 PRO-seq G1.
* kentUtils on `PATH` or at `KENT_BIN`: `bigWigToBedGraph`, `bedGraphToBigWig`, `bigWigMerge`, `bigWigInfo`,
  `twoBitToFa`, `liftOver`.
* Slurm, for the provided job scripts. Each stage is also a plain Python script you can run directly.
* `curl` and `gzip`.

## Configuration

```bash
cp config/site.env.example config/site.env   # then edit
```

Every setting is optional. Python reads environment variables first, then `config/site.env`.

| Setting | Default | Used for |
|---|---|---|
| `PPB_SLURM_CLUSTER`, `PPB_SLURM_ACCOUNT` | Slurm defaults | `--clusters` / `--account` on every job |
| `PPB_SBATCH_ARGS` | none | extra sbatch flags for every job, e.g. `--partition=cpu` |
| `PPB_GPU_SBATCH_ARGS` | none | extra flags for GPU jobs, e.g. `--partition=gpu --gres=gpu:a100:1` |
| `PPB_SLURM_LOGS` | `logs/slurm/` | job logs |
| `PPB_SCRATCH_ROOT` | `/workdir/$USER` | node-local job scratch (falls back to `$TMPDIR`) |
| `KENT_BIN` | `PATH`, then `/programs/kentUtils/bin` | kentUtils binaries |
| `DREG_K562_HMM_BED`, `GM12878_UW_DHS_BED` | none (public fallbacks) | the two non-public inputs, below |

The sbatch headers carry the partitions and GPU type of the original cluster (Cornell BioHPC cbsubscb). Flags in
`PPB_SBATCH_ARGS` and `PPB_GPU_SBATCH_ARGS`, or flags passed to `slurm/sub.sh`, override them.

### Non-public inputs

Two dREG-paper files were read from Danko lab storage. They are not publicly downloadable.

| Setting | File | What it feeds | Fallback when unset |
|---|---|---|---|
| `DREG_K562_HMM_BED` | `hg19.k562.new_hmm2b.post2.bed`: K562 GRO-cap HMM TREs | A manifest check that the dREG negative-rebuild recipe reproduces the released K562 negatives; also part of the TSS signal defining K562 TSS-silent DHS negatives | The check is skipped; K562 TSS-silent negatives use GRO-cap signal only |
| `GM12878_UW_DHS_BED` | `uw.merge.narrowPeak.bed`: GM12878 UW DNase peaks, merged | Every GM12878 DHS-based label: rebuilt dREG negatives, H3K27ac/DHS labels, quiet negatives, TSS-silent DHSs, and the restriction of HiDRA negatives to accessible sites | Public UCSC UW DNase GM12878 replicates 1 + 2 (`wgEncodeUwDnaseGm12878PkRep{1,2}`), concatenated |

With the fallbacks, the GM12878 and K562 TSS-silent label sets differ slightly from ours, so the numbers won't match
`results/` exactly. All other inputs are public: ENCODE, GEO, UCSC, FANTOM5, ChIP-Atlas, the dREG data hub, and the
PINTS paper's supplementary tables.

## Running

Each phase is fetch → labels and callers → evaluate → bootstrap → report. The `submit_*.sh` scripts chain the
Slurm stages with dependencies.

```bash
uv sync
cp config/site.env.example config/site.env            # edit for your site

# phase 1: K562 + GM12878, hg19
bash scripts/00_fetch_data.sh 2>&1 | tee logs/00_fetch_data.log
bash slurm/sub.sh slurm/00_gpu_smoke.sbatch            # check CuPy sees the GPU
bash slurm/submit_phase1.sh

# phase 2: HCT116, A549, K562 PRO-cap, hg38
bash scripts/00b_fetch_phase2.sh
bash slurm/submit_phase2.sh

# phase 3: Jurkat, hg19
bash scripts/00c_fetch_jurkat.sh
bash slurm/submit_phase3.sh

# re-evaluate and bootstrap all 15 libraries (TSS-silent negatives, GC and CpG conditions), then the report
bash slurm/submit_rerun_p12.sh
```

* **One job:** `bash slurm/sub.sh [sbatch flags] slurm/<stage>.sbatch`. Array jobs select libraries by phase:
  `--export=ALL,PPB_PHASE=1|2|3|all`.
* **Outputs:**
  * caller output: `results/calls/<lib>/{pydreg,pints}/`;
  * point metrics: `results/eval/<lib>/`;
  * intervals: `results/bootstrap/<lib>/ci.tsv`;
  * report: `results/summary/report.md`.
* **GC diagnostic:** `.venv/bin/python scripts/07_gc_diagnostic.py`.

## Tests

```bash
uv run pytest tests
```

These cover:

* bootstrap weights reproducing the point metrics;
* equivalence with the PINTS-paper evaluation rule;
* GC matching;
* condition handling.

Some tests use the built label tables under `data/`.

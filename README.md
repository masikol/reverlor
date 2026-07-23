# REVERLOR

**RE**peat **VER**ification using **LO**ng **R**eads

A bioinformatics tool that finds inexact interspersed repeats in genomic sequences via minimap2 self-alignment, then verifies whether those repeats are resolved (spanned) by long reads from a BAM alignment file. Repeats are classified as *resolved* or *unresolved* based on a configurable read-through (read spanning) threshold.

**Status:** v0.0.a (alpha) -- 2026-07-23

---

## Table of Contents

- [Motivation](#motivation)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Usage](#usage)
  - [BAM File Preparation](#bam-file-preparation)
  - [Full Pipeline (reverlor.py)](#full-pipeline-reverlorpy)
  - [Find Only (reverlor_find.py)](#find-only-reverlor_findpy)
  - [Verify Only (reverlor_verify.py)](#verify-only-reverlor_verifypy)
- [Options Reference](#options-reference)
- [Output Files](#output-files)
- [Examples](#examples)

---

## Motivation

1. Accurate reconstruction of repeat regions (e.g., rRNA clusters, IS elements, prophages) is necessary for complete and contiguous genome assemblies.

2. Even long reads (Oxford Nanopore, PacBio) sometimes fail to span a repeat end-to-end.

3. A program that would find such repeats and verify them (check how many reads span them end-to-end) would be useful.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/reverlor.git
cd reverlor
```

### 2. Create and activate a virtual environment

```bash
python3 -m venv reverlor_venv
source reverlor_venv/bin/activate
```

### 3. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 4. Test installation

```bash
pip install -r requirements_with_tests.txt
python3 -m pytest tests
```

---

## Quick Start

```bash
# Full pipeline: find repeats + verify with long reads
python3 reverlor/reverlor.py genome.fasta reads.bam output_dir/

# With verbose output
python3 reverlor/reverlor.py -vv genome.fasta reads.bam output_dir/

# With custom parameters
python3 reverlor/reverlor.py \
    --min-repeat-len 127 \
    --span 10 \
    --threads 4 \
    genome.fasta reads.bam output_dir/
```

---

## Pipeline Overview

### Step 1: Repeat Finding (`reverlor_find`)

1. Self-align the reference FASTA using minimap2 (via `mappy` Python package).
2. Extract alignment coordinates as raw BED entries.
3. Merge overlapping/adjacent repeats using `pybedtools`.
4. Filter by minimum repeat length.
5. Output: `repeats.bed`.

### Step 2: Repeat Verification (`reverlor_verify`)

1. Read repeat regions from the BED file created at Step 1.
2. For each repeat, compute **“shoulder” coordinates** (flanking regions upstream and downstream).
3. Use `pysam` to find long reads spanning from the upstream shoulder position to the downstream one.
4. If fewer than `--span` reads span the repeat, mark it as **unresolved**.
5. Output: `unresolved_repeats.bed`.

---

## Usage

### BAM File Preparation

For best performance, we strongly recommend to pass **sorted and indexed** BAM files to `reverlor.py`.

Here is a quick example on how to create such file:

```bash
minimap2 -a genome.fasta reads.fastq.gz \
    | samtools view -F 4 \
    | samtools sort -O BAM -o reads.sorted.bam

samtools index reads.sorted.bam
```

### Full Pipeline (`reverlor.py`)

Runs both repeat finding and verification in sequence.

```bash
python3 reverlor/reverlor.py [options] <fasta> <input_bam> <output_dir>
```

**Positional arguments:**

| Argument | Description |
|----------|-------------|
| `fasta` | Path to input FASTA file (reference genome) |
| `input_bam` | Path to input BAM file (sorted, indexed long reads) |
| `output_dir` | Path to output directory |

---

### Find Only (`reverlor_find.py`)

Finds repeats without verification.

```bash
python3 reverlor/reverlor_find.py [options] <input_fasta> <output_dir>
```

**Positional arguments:**

| Argument | Description |
|----------|-------------|
| `input_fasta` | Path to input FASTA file |
| `output_dir` | Path to output directory |

---

### Verify Only (`reverlor_verify.py`)

Verifies pre-existing repeat predictions against a **sorted and indexed** BAM file.

```bash
python3 reverlor/reverlor_verify.py [options] <input_fasta> <input_bed> <input_bam> <output_dir>
```

**Positional arguments:**

| Argument | Description |
|----------|-------------|
| `input_fasta` | Path to input FASTA file |
| `input_bed` | Path to input BED file (repeat predictions) |
| `input_bam` | Path to **sorted and indexed** input BAM file |
| `output_dir` | Path to output directory |

---

## Options Reference

### Common Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-V, --version` | flag | -- | Show version and exit |
| `-v, --verbose` | count | 0 | Verbosity level: `-v` (warnings), `-vv` (info), `-vvv` (debug) |
| `-t, --threads` | int | 1 | Number of CPU threads (passed to minimap2) |
| `--keep-tmp` | flag | False | Keep temporary files |
| `--tmpdir` | str | system | Temporary directory path |

### Repeat Finding Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-l/--min-repeat-len` | int | 200 | Minimum repeat length to report |
| `-i/--min-repeat-interval` | int | 100 | Minimum interval between repeats (shorter ones get merged) |
| `-k/--minimap-k` | int | 19 | minimap2 k-mer length |
| `-w/--minimap-w` | int | 19 | minimap2 minimizer window size |
| `-m/--minimap-m` | int | 65 | minimap2 min chain score |
| `--minimap-x` | choice | None | minimap2 preset: `map-ont`, `lr:hq`, `map-hifi`, `map-pb`, `map-iclr`, `asm5`, `asm10`, `asm20` |

### Repeat Verification Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-s/--span` | int | 5 | Minimum number os spanning reads to consider a repeat resolved |
| `-u/--shoulder-len` | int | 200 | Shoulder length (bp) for coordinate checking |
| `--samtools-f` | int | [] | samtools `-f` flag (includes segments having these flags, repeatable) |
| `--samtools-F` | int | [256] | samtools `-F` flag (excludes segments having these flags, repeatable) |

#### Filter reads using SAM flags (`--samtools-f, --samtools-F`)

- These flags are applied at Step 2.3 (see the [Pipeline Overview](#pipeline-overview) above).
- Default `-F 256` excludes secondary alignments from verification.

- Add `--samtools-F 2048` to also exclude supplementary alignments:
```
reverlor/reverlor.py \
    --samtools-F 256 \
    --samtools-F 2048 \
    genome.fasta \
    reads.sorted.bam \
    output/dir/
```

- Add `--samtools-f 16` to only include reverse-complementedly mapped reads.
```
reverlor/reverlor.py \
    --samtools-F 256 \
    --samtools-f 16 \
    genome.fasta \
    reads.sorted.bam \
    output/dir/
```

---

## Output Files

Reverlor outputs BED files.

### Full Pipeline (`reverlor.py`)

| File | Description |
|------|-------------|
| `repeats.bed` | Detected repeat regions |
| `unresolved_repeats.bed` | Repeats with insufficient spanning coverage |

### Find Only (`reverlor_find.py`)

| File | Description |
|------|-------------|
| `repeats.bed` | Detected repeat regions |

### Verify Only (`reverlor_verify.py`)

| File | Description |
|------|-------------|
| `unresolved_reapeats.bed` | Repeats with insufficient spanning coverage |


---

## Examples

### Example 1: Basic usage on a bacterial genome

```bash
python3 reverlor/reverlor.py \
    genome.fasta \
    reads.sorted.bam \
    output_dir/
```

### Example 2: Exclude secondary and supplementary alignments

```bash
python3 reverlor/reverlor.py \
    --samtools-F 256 \
    --samtools-F 2048 \
    genome.fasta \
    reads.sorted.bam \
    output_dir/
```

### Example 3: Find repeats only, no verification (no BAM required)

```bash
python3 reverlor/reverlor_find.py \
    genome.fasta \
    output_dir/
```

### Example 4: Verify existing repeat predictions (`repeats.bed`)

```bash
python3 reverlor/reverlor_verify.py \
    --span 10 \
    --shoulder-len 300 \
    genome.fasta \
    repeats.bed \
    reads.sorted.bam \
    output_dir/
```

### Example 5: Debug mode with temporary files

```bash
python3 reverlor/reverlor.py \
    -vvv \
    --keep-tmp \
    --tmpdir /tmp/reverlor_debug \
    genome.fasta \
    reads.sorted.bam \
    output_dir/
```

### Example 6: Full pipeline with verbose messages

```bash
python3 reverlor/reverlor.py \
    --min-repeat-len 200 \
    --min-repeat-interval 50 \
    --span 5 \
    --shoulder-len 200 \
    --samtools-F 256 \
    -t 8 \
    -vv \
    genome.fasta \
    reads.sorted.bam \
    output_full/
```

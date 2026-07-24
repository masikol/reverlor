# An Archived Step of Reverlor Development: Compare Repeat Finder Programs

Here is the code to compare different repeat finder programs:

1. minimap2 self-alignment.
2. RepeatScout.
3. phRAIDER.
4. TotalRepeats.
5. RepeatModeler.
6. GenericRepeatFinder.
7. REPrise.

## Contents

1. `reverlor_test_version/` an old development version of reverlor. There is no `reverlor_verify` in it. It supports multiple repeat finder programs though (`--finder` option).

2. `test_repeat_detection/`: there is a pipeline that tests accuracy of repeat finders.

## Running the pipeline

### 1. Install dependencies.

1. minimap2 version 2.30-r1287: [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2).
2. RepeatScout version 1.0.7: [https://github.com/Dfam-consortium/RepeatScout](https://github.com/Dfam-consortium/RepeatScout).
3. phRAIDER version 2.0 [https://github.com/karroje/phRAIDER](https://github.com/karroje/phRAIDER).
4. TotalRepeats.jar: `7435616` commit on `master` branch: [https://github.com/rkalendar/TotalRepeats/tree/743561689d2a199dd91495fe7b83e45b2ee5d82b](https://github.com/rkalendar/TotalRepeats/tree/743561689d2a199dd91495fe7b83e45b2ee5d82b).
5. RepeatModeler version 2.0.9 [https://github.com/Dfam-consortium/RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler).
6. GenericRepeatFinder version 1.0.2 [https://github.com/bioinfolabmu/GenericRepeatFinder](https://github.com/bioinfolabmu/GenericRepeatFinder).
7. REPrise: `6c1d829` commit on `master` branch: [https://github.com/hmdlab/REPrise/tree/6c1d829af2e88c079215da056dbe1b16a9a96967](https://github.com/hmdlab/REPrise/tree/6c1d829af2e88c079215da056dbe1b16a9a96967).
8. bedtools version v2.31.1: [https://bedtools.readthedocs.io/en/latest/content/installation.html](https://bedtools.readthedocs.io/en/latest/content/installation.html).
9. Python (3.12.3) packages:
```bash
pip install \
    polars==1.40.1 \
    Mutation-Simulator==3.0.2 \
    seaborn==0.13.2 \
    biopython==1.84 \
    numpy==1.26.4 \
    matplotlib==3.10.0
```

### 2. Set parameters, paths and dependencies.

Parameters, paths and dependencies are (quite awkwardly, I admit!) defined in two files:

1. `test_repeat_detection/pipeline.sh`;
2. `test_repeat_detection/mock_repeats_settings.py`.


#### 2.1 Parameters in `pipeline.sh`

There is a “Configuration” block in `pipeline.sh`:
```bash
# >>> Configuration >>>

FINDERS=(
    minimap2
    repeat-scout
    phraider
    total-repeats
    repeat-modeler
    grf
    reprise
)

MUTATION_TYPE='SNP' # choice
# MUTATION_TYPE='SNP_INS' # choice
# MUTATION_TYPE='SNP_DEL' # choice
# MUTATION_TYPE='SNP_INDEL' # choice

N_DESIRED_REPEAT_COPIES=3

N_THREADS=6

REPO_DIR='/mnt/data/Max/repos/reverlor'
WORKDIR='/mnt/data/Max/max_disser/work_dirs/reverlor'

# <<< Configuration <<<
```

Please set each of these varibales as you think is reasonable for you.

`FINDERS` is an array of available repeat finder programs to test. Feel free to comment out those you don’t intend to test.
`MUTATION_TYPE` is type of mutations performed by the `mutation-simulator` program. Leave only the desired one uncommented.
`N_DESIRED_REPEAT_COPIES`: number of repeat copies the script `insert_mock_repeats.py` will insert into the .
`REPO_DIR` is a directory that contains the directory `cmp_repeat_finders/`.

`WORKDIR` is an output directory for the pipeline. It must contain two input fasta files:

1. `data/Mycoplasma_mycoides_JCVI-syn3.0.fasta`. Can be downloaded here: [CP014940.1](https://www.ncbi.nlm.nih.gov/nuccore/CP014940.1).
2. `data/pUC18.fasta` (prodided: `cmp_repeat_finders/data/pUC18.fasta`).

Together, these two files represent an input mock genome into which mock repeats will be inserted by the pipeline.

So if your `WORKDIR` is `/some/path` then it must contain `/some/path/data/Mycoplasma_mycoides_JCVI-syn3.0.fasta` and `/some/path/data/pUC18.fasta`.

#### 2.1 Parameters in `pipeline.sh`

`RATE_FROM`, `RATE_TO`, and `RATE_STEP` define rates of single nucleotide replacements and indels performed by Metation-Simulator.

`N_REPEAT_REPLICATES`: is the number of test replicates. For example, if `N_REPEAT_REPLICATES` is 10, then `N_DESIRED_REPEAT_COPIES` mock repeats will be inserted in random coordinates into 10 replicates of the input mock genome: for each finder program and for each mutation rate.

All other settings seem very obvious. Feel free to ask anyway.

### 3. Actually run the pipeline

```bash
bash test_repeat_detection/pipeline.sh
```

### 4. See the results

Feel free to play around in `assess_pipeline_result/explore_accuracy.ipynb`

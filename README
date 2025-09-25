# Reciprocal Best Matching (RBM)

A computational pipeline called Reciprocal Best Matching (RBM), which is designed for evaluating protein complex structure predictions with unknown stoichiometry during CASP16 assessment.


## Directory Structure

- `RBM/` - Implementation of the Reciprocal Best Matching algorithms.
- `main.py` - Main script for running the pipeline.
- `input/` - Directory for example input structures and data.


## Overview

We developed and used RBM in our assessment of CASP16 Phase 0 challenge, in which stoichiometry information is not provided to the predictors. This pipeline provides an unbiased evaluation method that:
- Penalizes both over-prediction and under-prediction of subunits
- Maintains strong correlation with established CASP metrics
- Integrates with existing CASP assessment pipelines in regular challenges
- Can be extended to any scores that assess a binary protein pair.

## Installation

Our pipeline requires the following software:
- DockQ
- lDDT
- TMscore

To install DockQ, please refer to:

`https://github.com/wallnerlab/DockQ`

To install lddt, please refer to:

`https://anaconda.org/bioconda/lddt`

To install TMscore, you can consider installing it from source:

`https://github.com/pylelab/USalign`



## Usage

Basic usage:

```bash
python main.py --input_dir <input_directory> \
               --target <target_name> \
               --name <structure_name> \
               --output_dir <output_directory> \
               [optional arguments]
```

### Required Arguments

- `--input_dir`: Directory containing input data.
- `--target`: Target name provided by CASP organizers.
- `--name`: Name of the input structure used in the assessment stage.
- `--output_dir`: Directory for output files.

### Optional Arguments

- `--qs_cutoff`: Cutoff for QS_best (default: 10)
- `--dockq_path`: Path to DockQ executable (default: "DockQ")
- `--lddt_path`: Path to lDDT executable (default: "lddt")
- `--tmscore_path`: Path to TMscore executable (default: "TMscore")
- `--scores`: Choose 1-6 scores from: ICS, IPS, QS_best, DockQ, lDDT, TMscore (default: ICS, IPS, QS_best, DockQ, lDDT, TMscore)
- `--n_cpu`: Number of CPU cores to use (some software are computationally intensive, so multi-processing is necessary) (default: 48)
- `--antibody`: Flag for antibody structure analysis
- `--chainAs`: Chain IDs for the antibody component
- `--chainBs`: Chain IDs for the antigen component
- `--rbm_version`: RBM scoring version ('min', 'all', or 'average', default: 'min'. For more details, please refer to our paper)
- `--interface_weight`: Interface weighting method ('log2', 'log10', or 'linear', default: 'log10'. For more details, please refer to our paper)




## Output

The pipeline generates:
- Interface quality scores (ICS, IPS, QS_best)
- Structural similarity metrics (DockQ, lDDT, TM-score)
- Per-interface evaluation results
- Per-model evaluation results


## Citation

If you found our method useful, please cite our paper:

If you found our CASP16 assessment results useful for your research, please consider citing our CASP16 assessment paper:

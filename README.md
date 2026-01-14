# Reciprocal Best Matching (RBM)

A computational pipeline called Reciprocal Best Matching (RBM), which is designed for evaluating protein complex structure predictions with unknown stoichiometry during CASP16 assessment.


## Overview

We developed and used RBM in our assessment of CASP16 Phase 0 challenge, in which stoichiometry information is not provided to the predictors. This pipeline provides an unbiased evaluation method that:
- Penalizes both over-prediction and under-prediction of subunits
- Maintains strong correlation with established CASP metrics
- Integrates with existing CASP assessment pipelines in regular challenges
- Can be extended to any scores that assess a binary protein pair.

## Installation

Our pipeline requires the following software:
- [DockQ](https://github.com/wallnerlab/DockQ)
- [lDDT](https://anaconda.org/bioconda/lddt)
- [TMscore](https://github.com/pylelab/USalign)

```
pip install DockQ
conda install -c conda-forge -c bioconda lddt tmalign
```

## Input Example

To execute the pipeline, target structure, model structures, a list of models, and the OST raw results need to be provided. The directory structure should be as follows:


```
input/
└── H0208/                                    # Target directory (target name)
    ├── H0208.pdb                             # Reference structure (1 file)
    ├── H0208.txt                             # Model list file (in our case, this file is from the prediction center)
    ├── model/                                # Model structures directory
    │   ├── H0208TS014_1                      # Model file (in PDB format)
    │   ├── H0208TS014_2                      ......
    │   ├── H0208TS014_3                      ......
    │   ├── ...                               
    │   └── H0208TS331_5                      ......
    └── ost/                                  # OpenStructure results
        ├── H0208TS014_1.json                 # OST comparison result for model 1
        ├── H0208TS014_2.json                 ......
        ├── H0208TS014_3.json                 ......
        ├── H0208TS022_2.json                 ......
        ├── ...                               
        └── H0208TS331_5.json                 ......
```


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


## Output Example

The pipeline generates:
- Interface quality scores (ICS, IPS, QS_best)
- Structural similarity metrics (DockQ, lDDT, TM-score)
- Per-interface evaluation results
- Per-model evaluation results

```
output/
└── H0208/
    ├── IPS/                              # Interface Patch Similarity scores
    │   ├── H0208TS014_1.result              # Per-interface IPS for each model
    │   └── ...
    ├── ICS/                              # Interface Contact Similarity scores
    │   ├── H0208TS014_1.result              # Per-interface ICS for each model
    │   └── ...
    ├── QS_best/                          # QS-score (best match)
    │   └── H0208.result         # All models combined
    ├── DockQ/                            # DockQ scores
    │   └── H0208.result         # All models combined
    ├── lDDT/                             # Local Distance Difference Test
    │   └── H0208.result         # All models combined
    ├── TMscore/                          # TM-score (structure alignment)
    │   └── H0208.result         # All models combined
    ├── per_interface_scores/             # Consolidated per-interface scores
    │   └── H0208.result         # Best scores for each interface
    └── per_model_scores_v3/              # Final model rankings (v3, recommended)
        └── H0208.result
```

## Citation

If you found our method useful, please cite our paper:

If you found our CASP16 assessment results useful for your research, please consider citing our CASP16 assessment paper:

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

```bash
# Create and activate a conda environment
conda create -n rbm python=3.9
conda activate rbm

# Install required packages

pip install DockQ
conda install -c conda-forge -c bioconda lddt tmalign
```

## Input Requirements

The pipeline now uses a **file-based approach** where you specify individual file paths instead of directory structures. You need to provide:

1. **Reference PDB file**: The target/reference structure
2. **Model PDB file**: The predicted model structure
3. **OST JSON file**: OpenStructure comparison results containing chain mappings and contacts

## Usage

Basic usage:

```bash
python main.py \
  --reference_pdb <path/to/reference.pdb> \
  --model_pdb <path/to/model.pdb> \
  --ost_json <path/to/ost_comparison.json> \
  --output_dir <output_directory> \
  --target_name <target_name>
```

### Example

```bash
python main.py \
  --reference_pdb ./data/H0208/H0208.pdb \
  --model_pdb ./data/H0208/model/H0208TS014_1 \
  --ost_json ./data/H0208/ost/H0208TS014_1.json \
  --output_dir ./results \
  --target_name H0208 \
  --model_name H0208TS014_1
```

### Required Arguments

- `--reference_pdb`: Path to reference/target PDB file
- `--model_pdb`: Path to model PDB file
- `--ost_json`: Path to OpenStructure (OST) JSON file with chain mappings and contacts
- `--output_dir`: Directory for output files
- `--target_name`: Target name for organizing output files (e.g., H0208)

### Optional Arguments

- `--model_name`: Model name to use in output files (default: basename of model_pdb)
- `--dockq_path`: Path to DockQ executable (default: "DockQ")
- `--lddt_path`: Path to lDDT executable (default: "lddt")
- `--tmscore_path`: Path to TMscore executable (default: "/home2/s439906/software/USalign/TMscore")
- `--scores`: Choose 1-6 scores from: ICS, IPS, QS_best, DockQ, lDDT, TMscore (default: all)
- `--antibody`: Flag for antibody structure analysis
- `--chainAs`: Chain IDs for the antibody component (for antibody mode)
- `--chainBs`: Chain IDs for the antigen component (for antibody mode)
- `--rbm_version`: RBM scoring version ('min', 'all', or 'average', default: 'min'. For more details, please refer to our paper)
- `--interface_weight`: Interface weighting method ('log2', 'log10', or 'linear', default: 'log10'. For more details, please refer to our paper)
- `--keep_tmp`: Keep temporary interface files in interface_tmp directory (default: False, will remove after completion)


## Output

Results are saved in `output_dir/model_name/`:

- **Per-interface scores**: `model_name.interface_scores` (tab-separated file with IPS, ICS, QS, DockQ, lDDT, TM-score for each interface)
- **Final RBM score**: `model_name.model_scores` (final model evaluation)
- **Individual metric files**: `model_name.ips`, `model_name.ics`, `model_name.qs`, `model_name.dockq`, `model_name.lddt`, `model_name.tm`

Temporary interface files are stored in `interface_tmp/` and automatically removed after completion (use `--keep_tmp` to preserve them).

## Citation

If you found our method useful, please cite our paper:

If you found our CASP16 assessment results useful for your research, please consider citing our CASP16 assessment paper:

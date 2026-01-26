# Reciprocal Best Matching (RBM)

Reciprocal Best Matching (RBM) is a computational pipeline called designed for evaluating protein complex structure predictions with unknown stoichiometry during CASP16 assessment.


## Overview

We developed and used RBM in our assessment of CASP16 Phase 0 challenge, in which stoichiometry information is not provided to the predictors. This pipeline provides an unbiased evaluation method that:
- Penalizes both over-prediction and under-prediction of subunits
- Maintains strong correlation with established CASP scores
- Integrates with existing CASP assessment pipelines in regular challenges
- Can be extended to any scores that assess a binary protein pair.

## Installation

Our pipeline requires the following software:
- [DockQ](https://github.com/wallnerlab/DockQ)
- [lDDT](https://anaconda.org/bioconda/lddt)
- [TMscore](https://github.com/pylelab/USalign)

Installation instructions:

```bash
# Create and activate a conda environment
conda create -n rbm python=3.9
conda activate rbm

# Install required packages
pip install numpy
pip install DockQ
conda install -c conda-forge -c bioconda lddt tmalign
```

**Important**: If you have not added DockQ, lDDT, and TMscore to your PATH, you must provide their full paths using `--dockq_path`, `--lddt_path`, and `--tmscore_path` arguments when running the pipeline.

## Input Requirements

The pipeline now uses a **file-based approach** where you specify individual file paths instead of directory structures. You need to provide:

1. **Reference PDB file**: The target(reference) structure
2. **Model PDB file**: The model structure
3. **OST JSON file**: OpenStructure comparison results containing `chem_groups`, `chem_mapping`, and `reference_contacts`, and `model_contacts`

Note: this software needs to pre-compute chain mappings and contacts using [OpenStructure](https://git.scicore.unibas.ch/schwede/openstructure) (OST). The best way to do this is to use the Docker image files provided by OST developers. An example command to pre-compute chain mappings and contacts is:

```bash
apptainer run \
  -B <path/to/data>:/work/data \
  -B <path/to/output>:/work/output \
  /path/to/openstructure.sif \
  compare-structures \
    -m /work/data/<model_name>.pdb \
    -mf pdb \
    -r /work/data/<reference_name>.pdb \
    -rf pdb \
    -o /work/output/<model_name>.json \
    --ics \
    --ips \
```

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
  --model_pdb ./input_example/H0208/model/H0208TS014_1 \
  --ost_json ./input_example/H0208/ost/H0208TS014_1.json \
  --output_dir ./output_example \
  --target_name H0208 \
  --model_name H0208TS014_1 \
  --dockq_path /path/to/DockQ \
  --lddt_path /path/to/lddt \
  --tmscore_path /path/to/TMscore
```

Note: The tool paths (`--dockq_path`, `--lddt_path`, `--tmscore_path`) are only needed if these tools are not in your system PATH.

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
- `--interface_weight`: Interface weighting method ('log10' or 'linear', default: 'log10'. For more details, please refer to our paper)
- `--keep_tmp`: Keep temporary interface files in interface_tmp directory (default: False, will remove after completion)


## Output

Results are saved in `output_dir/model_name/`:

- **Per-interface scores**: `model_name.interface_scores` (tab-separated file with IPS, ICS, QS, DockQ, lDDT, TM-score for each interface)
- **Final RBM score**: `model_name.model_scores` (final model evaluation)

- **Individual score files**: `model_name.ips`, `model_name.ics`, `model_name.qs`, `model_name.dockq`, `model_name.lddt`, `model_name.tm`

<details>
<summary>Click to view details of each file.</summary>

Below are descriptions of each individual score file produced in your output directory.

### IPS (Interface Patch Similarity) Scores:

**Format**: Tab-separated with columns: `model_name`, `category`, `chainpair`, `matchpair`, `size1`, `size2`, `score`

**Example output**:
```
H0217TS014_2	reference	C:G	H:A	21	26	0.0
H0217TS014_2	reference	C:G	H:G	21	26	0.8039
H0217TS014_2	reference	C:G	B:A	21	26	0.8039
H0217TS014_2	reference	C:G	B:G	21	26	0.0
```

### ICS (Interface Contact Similarity) Scores:

**Format**: Tab-separated with columns: `model_name`, `category`, `chainpair`, `matchpair`, `size1`, `size2`, `score`

**Example output**:
```
H0217TS014_2	reference	C:G	H:A	21	26	0.0
H0217TS014_2	reference	C:G	H:G	21	26	0.8276
H0217TS014_2	reference	C:G	B:A	21	26	0.8276
H0217TS014_2	reference	C:G	B:G	21	26	0.0
```

### QS_best:

**Format**: Tab-separated with columns: `model_name`, `category`, `chainpair`, `matchpair`, `score`

**Example output**:
```
H0217TS014_2	reference	J:L	I:F	0.0
H0217TS014_2	reference	J:L	I:L	0.0
H0217TS014_2	reference	J:L	C:F	0.0
H0217TS014_2	reference	J:L	C:L	0.0
```

### DockQ:

**Format**: Tab-separated with columns: `model_name`, `category`, `chainpair`, `matchpair`, `score`

**Example output**:
```
H0217TS014_2	reference	A:H	E:A	0.854
H0217TS014_2	reference	H:A	A:E	0.854
H0217TS014_2	prediction	D:H	B:A	0.801
H0217TS014_2	prediction	H:D	A:B	0.801
```

### lDDT (Local Distance Difference Test) Scores:

**Format**: Tab-separated with columns: `model_name`, `category`, `chainpair`, `matchpair`, `score`

**Example output**:
```
H0217TS014_2	reference	A:H	E:A	0.8967
H0217TS014_2	reference	H:A	A:E	0.8967
H0217TS014_2	prediction	D:H	B:A	0.887
H0217TS014_2	prediction	H:D	A:B	0.887
```

### TM-score:

**Format**: Tab-separated with columns: `model_name`, `category`, `chainpair`, `matchpair`, `score`

**Example output**:
```
H0217TS014_2	reference	A:H	E:A	0.9822
H0217TS014_2	reference	H:A	A:E	0.9822
H0217TS014_2	prediction	D:H	B:A	0.9828
H0217TS014_2	prediction	H:D	A:B	0.9828
```

### interface_scores: Aggregated Per-Interface Scores

**Format**: Tab-separated with header row: `model`, `category`, `chainpair`, `matchpair`, `weight`, `ips`, `ics`, `qsbest`, `dockq`, `lddt`, `tm`

**Example output**:
```
model	category	chainpair	matchpair	weight	ips	ics	qsbest	dockq	lddt	tm
H0217TS014_2	REF	C:G	H:G	23.5	0.8039	0.8276	0.8714	0.784	0.8879	0.982
H0217TS014_2	REF	E:H	D:G,J:A	3.5	0.5	0.5	0.7246	0.706	0.8899	0.9805
H0217TS014_2	REF	F:L	na	14.0	0.0	0.0	0.0	0.0	0.0	0.0
```

### model_scores: Final Model-Level Scores

**Format**: Tab-separated with header row: `model`, `ips`, `ics`, `qsbest`, `dockq`, `lddt`, `tm`

**Example output**:
```
model	ips	ics	qsbest	dockq	lddt	tm
H0217TS014_2	0.609	0.5868	0.6064	0.59	0.6806	0.7587
```



</details>

Temporary interface files are stored in `interface_tmp/` and automatically removed after completion (use `--keep_tmp` to preserve them).

## Citation

If you found our method useful, please cite our paper:

If you found our CASP16 assessment results useful for your research, please consider citing our CASP16 assessment paper:

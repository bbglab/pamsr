# PAMSR

## Overview

**PAMSR** is a pipeline designed to investigate the statistical constraints underlying the detection of trinucleotide mutational signatures. It tests the detectability of a given signature across different levels of mutational activity.

## Workflow

The pipeline consists of **two sequential stages**:

### Signature Refitting

The first stage performs a mutational signature refitting on an input mutation dataset using [SigProfilerAssignment](https://github.com/SigProfilerSuite/SigProfilerAssignment).

Following the refitting step, a signature of interest can be selected as the target signature for the subsequent power analysis.

### Synthetic Data Generation and Power Analysis

The second stage generates synthetic mutation-count data based on the characteristics of the input dataset. The selected signature is introduced at controlled levels of injection, allowing to study the behaviour of the signature reconstruction across different signature activites and sample sizes. Signature reconstruction is studied using [mSigAct](https://github.com/steverozen/mSigAct).

By systematically varying factors such as the sample size and the signature injection, PAMSR characterize the statistical power and limitations of mutational signature refitting methods.

## Installation

PAMSR is distributed as a Nextflow pipeline through this repository. A conda recipe can be found in the `assets` folder, so one can easily prepare an enviorment with nextflow and run the pipeline:

```bash
# Clone the github repository
git clone https://github.com/bbglab/pamsr.git
# Prepare the conda enviorment
cd pamsr
conda env create -f assets/conda_recipes/nf_env.yml
```

## Requirements

PAMSR requires the following components:

- Nextflow.
- Singularity runtime used to execute the PAMSR components.
- Containers are automatically pulled from Docker Hub:
  - docker.io/gomdomingoa/msigact:v0.1.0: R-based statistical analysis component, containing mSigAct.
  - docker.io/gomdomingoa/pamsr-py:v0.1.0: python-based data processing component, containing various libraries for data processing.
  - docker.io/ferriolcalvet/sigprofiler_assignment:1.1.3: python-based analysis component, containing SigProfilerAssignment.

## Usage

### Signature Reffiting

For this stage, the input consists in a `samplesheet` that leads to mutational matrices and a `metadata` file that lists the `sample_id` and `group` of the samples that will be used. Both can be inputed as `.csv` or `.tsv` files (which must be specified in the corresponding parameter).

The samplesheet must contain, at least, the columns `input_path` and `delimiter`. It should look as follows:

`samplesheet.csv`:

```csv
input_path,delimiter,...
/home/...../input/mutational_matrices/a_data.csv,csv,...
/home/...../input/mutational_matrices/b_data.csv,csv,...
/home/...../input/mutational_matrices/c_data.tsv,tsv,...
```
The samplesheet must contain, at least, the column `sample_id`, whoose values must match the name of the sample that is indicated in the first row of the mutational matrix. It should look as follows:

`metadata.csv`:

```csv
sample_id,group,...
P105A,lung,...
P105p,lung,...
P105f,lung,...
P105C,lung,...
P105B,lung,...
```

Then, it is recommended to choose a `project_name` (by default will be set to `test`), and the pipeline will be ready to use.

```bash
nextflow run main.nf \
    --project_name "test" \
    --input_mode "mutational_matrix" \
    --analysis_mode "signature_profiler_assignment" \
    --metadata_ip_mm "assets/input_preparation/mutational_matrices/metadata_ip_mm.csv" \
    --samplesheet_ip_mm "assets/input_preparation/mutational_matrices/samplesheet_ip_mm.csv" \
    --metadata_delim "csv" \
    --samplesheet_delim "csv" \
    --genome_assembly "GRCh38" \
    -with-singularity
```

### Synthetic Data Generation and Power Analysis

For this stage, one should have in mind which mutational signature will be the focus of the study. The input can be specified as before, but for this step is higly recommended to try the `tuned_reconstruction` mode. This allows the user to reconstruct the input data from the mutational matrices derived from the results of the signature reffiting, removing the contribution of a given signature (generally, the one that is under study).

This new type of input consists in a `samplesheet` that leads to an `activity` matrix and a `signatures_catalog` used by the signature refitting tool. The `metadata` file consists in a file that lists the `sample_id` and `group` of the samples that will be used. Again, both files be inputed as `.csv` or `.tsv` files (which must be specified in the corresponding parameter).

The samplesheet must contain, at least, the columns `input_path` and `delimiter`. It should look as follows:

`samplesheet.csv`:

```csv
input_path,delimiter,type
/home/.../input/data/signatures_activities.txt,tsv,activities
/home/.../input/data/solution_signatures.txt,tsv,signatures_catalog
```
The samplesheet must contain, at least, the column `sample_id`, whoose values must match the name of the sample that is indicated in the first row of the mutational matrix. It should look as follows:

`metadata.csv`:

```csv
sample_id,...
P105A,lung,...
P105p,lung,...
P105f,lung,...
P105C,lung,...
P105B,lung,...
```

Then, it is recommended to choose a `project_name` (by default will be set to `test`), and the pipeline will be ready to use.

```bash
nextflow run main.nf \
    --project_name "test" \
    --input_mode "tuned_reconstruction" \
    --analysis_mode "power_analysis" \
    --metadata_ip_tr "assets/input_preparation/tuned_reconstruction/metadata_ip_tr.csv" \
    --samplesheet_ip_tr "assets/input_preparation/tuned_reconstruction/metadata_ip_tr.csv" \
    --metadata_delim "csv" \
    --samplesheet_delim "csv" \
    --genome_assembly "GRCh38" \
    --target_signature_eliminate_reconstruction "SBS31" \
    --target_signature_injection "SBS31" \
    --target_signature_pa "SBS31" \
    -with-singularity
```

## Credits

pamser was originally written by @gomdomingoa and supervised by @koszulordie.

## Citations

An list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## Documentation

---
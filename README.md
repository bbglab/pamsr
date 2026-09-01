# PAMSR

Power Analysis for Mutational Signature Reconstruction

## Overview

**PAMSR** is a pipeline designed to investigate the statistical constraints underlying the detection of trinucleotide mutational signatures.

...

## Workflow

The pipeline consists of **two sequential stages**:

### 1. Signature Refitting

The first stage performs mutational signature refitting on an input mutation dataset using [SigProfilerAssignment](https://github.com/SigProfilerSuite/SigProfilerAssignment).

Following the refitting step, a signature of interest can be selected as the target signature for the subsequent power analysis.

### 2. Synthetic Data Generation and Power Analysis

The second stage generates synthetic mutation-count data based on the characteristics of the input dataset. The selected signature is introduced at controlled levels of injection, allowing to study the behaviour of the signature reconstruction across different signature activites and sample sizes. Signature reconstruction is studying using [mSigAct](https://github.com/steverozen/mSigAct).

By systematically varying factors such as the sample size and the signature injection, PAMSR can be used to characterize the statistical power and limitations of mutational signature refitting methods.

## Installation 
PAMSR is distributed as a Nextflow pipeline through this repository.

## Requirements

PAMSR requires the following components:

- Nextflow.
- Singularity runtime used to execute the PAMSR components.
- Containers are automatically pulled from Docker Hub:
  - pamsr-r: R-based statistical analysis component, containing mSigAct.
  - pamsr-py: Python-based data processing component, containing SigProfilerAssignment.

## Quick Start

...

## Input

...

## Output

...

## Parameters

...

## Test Run

...

## Reproducibility

...

## Credits

pamser was originally written by Alberto Domingo Gómez and supervised by @koszulordie.

## Citations

An list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## Documentation

Find the documentation ([link to docs](https://github.com/bbglab/deepCSA/tree/main/docs)).

We are working to provide the biggest possible detail on the [usage](docs/usage.md) and explanation of the rationale and [tools](docs/tools.md), but this is still in progress.

## License
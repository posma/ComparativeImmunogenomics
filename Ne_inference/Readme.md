# Ne_inference: Effective Population Size Inference Using PSMC

This directory is part of the [ComparativeImmunogenomics](https://github.com/posma/ComparativeImmunogenomics) project. It provides tools for simulating sequencing reads, running PSMC, and analyzing inferred demographic histories of primates and other mammals.

---

## Contents

- `simulate_uniform_reads.py` – Simulates uniform Illumina reads from genome assemblies at specified coverage.
- `Snakefile` - Snakemake rules for data preparation and PSMC analysis for multiple dipoid genomes
- `parse_PSMC_output.py` – Parses PSMC output, computes summary statistics, and outputs tables.

---
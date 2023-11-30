# European Green Woodpecker Genome Assembly

## Overview

This repository contains scripts and tools for the assembly and analysis of the genome of the European Green Woodpecker (*Picus viridis*). The assembly process involves various steps, including chromosome naming, visualization using Circos, calculation of pi diversity, population genetics statistics, and analysis of the site frequency spectrum (SFS).

## Repository Structure

- **Circos:** This directory contains scripts and files related to the visualization of the genome using Circos.
  - `circos_picus_colaptes.R`: Circos script for visualizing the genome with related species (*Colaptes auratus*).
  - `circos_picus_gallus.R`: Circos script for visualizing the genome with the reference model organism, the Chicken (*Gallus gallus*).
  - `pviridis_chromosomes_names.csv`: CSV file containing chromosome names for *Picus viridis*.
  - `rename_picus_chr_v2.py`: Python script for renaming *Picus viridis* chromosomes.
  - `Rename_picus_contigs.ipynb`: Jupyter notebook for renaming *Picus viridis* contigs.
  - `rename_tablev2.csv`: CSV file containing renaming information.
  - `run_mummer.sh`: Shell script for running MUMmer, a system for aligning whole genome sequences.

- **pi_diversity.sh:** Shell script for calculating pi diversity, a measure of nucleotide diversity.

- **popgen_stats_picus.R:** R script for computing population genetics statistics.

- **SFS:** This directory contains scripts related to the analysis of the site frequency spectrum (SFS).
  - `customgraphics.py`: Python script for custom graphics related to SFS.
  - `SFS_picus.py`: Python script for computing the SFS for *Picus viridis*.
  - `sfs_tools.py`: Python script containing tools for SFS analysis.

## Components description

1. **Genome Comparison:**
   - Execute `run_mummer.sh` for genome comparison using MUMmer. Used as an input for sytneny analysis and Circos plotting.

2. **Chromosome Naming and Visualization:**
   - Use `rename_picus_chr_v2.py` and `Rename_picus_contigs.ipynb` for renaming chromosomes and contigs accordingly to the genome of *Gallus gallus*.
   - Run Circos scripts (`circos_picus_colaptes.R` and `circos_picus_gallus.R`) for visualization.

3. **Population Genetics Analysis:**
   - Run `pi_diversity.sh` for calculating Watterson's Pi diversity using vcftools.
   - Execute `popgen_stats_picus.R` for additional population genetics statistics.

4. **Site Frequency Spectrum Analysis:**
   - Utilize scripts in the `SFS` directory for analyzing the site frequency spectrum.

## License

This project is licensed under the MIT License - see the LICENSE file for details.


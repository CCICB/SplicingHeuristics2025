# SplicingHeuristics2025

Supporting code for figures in our 2025 Splicing Heuristics Manuscript (Sullivan et al.), _"Data-driven insights to inform splice-altering variant assessment."_

This repository provides the code and processing pipelines used to generate the figures and analyses for the manuscript. The primary goal of this project is to explore and assess splice-altering variant effects using data-driven heuristics. The variant data referenced in this repository can be linked to IDs in [SpliceVarDB](https://www.splicevardb.org), though explicit variant information is not provided here.

---

## Repository Structure

The repository is organized into several directories, each with a specific role in the analysis pipeline:

1. `exons/`
Contains the scripts and processing logic for working with Gencode files to extract data about known exons and introns. 

2. `branchpoints/`
Processes the experimentally identified branchpoints from Mercer et al. Note hg19 genome build.

3. `requirements/`
Includes the scripts to compute splicing requirements and create sequence logos.

4. `heuristics/`
Houses the code for generating data-driven heuristics for splice-altering variant evaluation and spliceogenicity.

5. `novel/`
Contains the code to apply splicing requirements to novel variants and the plot for pseudoexon inclusion mechanisms.

6. `outcomes/`
Plotting the splicing transcript outcomes based on the location of the splice-altering variant.


# SARS-CoV-2_BW_lineage 
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![DOI](https://zenodo.org/badge/566477400.svg)](https://zenodo.org/badge/latestdoi/566477400)

## Description: 
This code represents the implementation of our methodology for preprocessing the Genomic DB and make the respective analyses. 
With this code, we generated the images and results of the manuscript [1] accepted on Infection and available at https://doi.org/10.21203/rs.3.rs-2285898/v1.  

# Structure of the code (Pseudocode): 
Manipulation of .fasta files: (bash commands)
1) Importing dataset
2) Linearize fasta file with multiple sequences
3) Remove those with incomplete collection date
4) Create of a comparison table with ATCG contents using [2] [seqtk](https://github.com/lh3/seqtk)
5) Determine quality of sequences (those with N % < 10)
6) Extract all high-quality sequences 
7) Compute a sequence alignment using [3] [nextalign](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextalign-cli.html) application 
8) Linearize MAFFT output fasta file
9) Calling [4][R](https://www.r-project.org/) scripts for data engineering 
 * preFilter_raw_metatable.R
 * Filter_preTable.R
 * Analyze_MexSARS_indels.R
10) Create longidutinal pdfs with [5] [pdftk](https://www.pdflabs.com/t/pdftk/)
11) Finally get the consensus per variant in each category
 * Mutations_table.R
 * Get_consensus.R
 * Get_mutations.R	

## Table of Contents:
This repository contains different types of formats. We have different Scripts 
and xlsx files which contain data until 7th November
- The /bash/main_commands.sh file: contains the general pipeline followed in this work
- The /code/preFilter_raw_metatable.R File: apply data cleaning and preprocessing
- The /code/Filter_preTable.R file: contains data engineering process
- The /code/Analyze_MexSARS_indels.R file: creates a table to analyze nt mutations (substitutions) per positions using a linearized fasta file
- The /code/Mutations_table.R file: explore all the alignments and creates table comparing of every genome by position
- The /code/Get_consensus.R file: parse a table containing SARS-CoV-2 nt positions as rows and each nucleotide (including gaps and Ns)
- The /code/Get_mutations.R file: parse a table containing SARS-CoV-2 nt positions as rows and each nucleotide (including gaps and Ns) creates a 
summary table of all  possible observed mutations per position
- The /trees/BW_BA.5.6.2_with_UShER_refs_Ancestral_states.nexus file: contains the phylogenetic reconstruction of BW.x and BA.5.6.2 
genomes wordwide, with reference sequences (shown as Supplementary Figure 3). Inner nodes show the ancestral state reconstruction of 
specific mutation events. Reference genomes were obtained using UShER
- The /trees/PastML_geographical_inference.zip file: contains the putative origin of the BW clade with a maximum likelihood ancestral 
region reconstruction made with PastML

## How to run the scripts

To run our code:
- Download the full content of the directory containing this README file.
- Make sure you have installed R (we developed the code under version: R.4.2.2 (Innocent and Trusting; released on 2022-10-31), so we suggest this or subsequent versions). 
Also, please verify the corresponding programs / dependencies mention in the section below

We used the following version of the programs/libraries :
- Nextstrain v.2.3.0
- Nextalign v2.11.0 
- UShER v0.6.2
- TreeTime v0.9.5
- IQTree v2.1.2
- PastML v1.9.34
- Seqtk-1.3 (r106)
- See Credits for other programs mention before

## Credits
To use our original or adapted codes, please cite our work https://doi.org/10.21203/rs.3.rs-2285898/v1 as [1], see reference at the end of this document.

The R script codes present in this directory has been written by Rodrigo García-López, Xaira Rivera-Gutierrez and Mauricio Rosales-Rivera.
We gratefully acknowledge the Consorcio Mexicano de Vigilancia Genómica (CoViGen-Mex) and the authors from the originating laboratories responsible for obtaining the specimens and the submitting laboratories from which genetic sequence data were generated and shared via the GISAID initiative included in Supplementary Table 1, as well as the Unidad de Investigación Médica de Yucatán of the Instituto Mexicano del Seguro Social, which collected most of the BW lineage and BA.5.6 samples used in this study. We appreciate the computer assistance provided by Jerome Verleyen, and Juan Manuel Hurtado. The project LANCAD-UNAM-DGTIC-396 of the Dirección General de Cómputo y Tecnologías de la Información (DGTIC-UNAM) provided supercomputing resources in MIZTLI.

The first version of this work appeared on Research Square in November 2022 as preprint: https://doi.org/10.21203/rs.3.rs-2285898/v1. 

### FULL CITATION
SARS-CoV-2 BW lineage, a fast-growing Omicron variant from southeast Mexico bearing relevant escape mutations. García-López R*, Rivera-Gutiérrez X, Rosales-Rivera M, Taboada B*, Zárate S, Muñoz-Medina JE, Roche B, Herrera-Estrella A, Gómez-Gil B, Sanchez-Flores A, Arias CF. Infection. 2023 Apr 14. doi: 10.1007/s15010-023-02034-7.

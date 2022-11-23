# SARS-CoV-2_BW.1

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Description: 
This code represents the implementation of our methodology for preprocessing the Genomic DB. 
With this code, we generated the images and results of the manuscript [1] submitted on Infection and available at https://doi.org/10.21203/rs.3.rs-2285898/v1. 

# Structure of the code (Pseudocode): 
Manipulation of .fasta files: (bash commands)
1) Importing dataset
2) Linearize fasta file with multiple sequences
3) Remove those with incomplete collection date
4) Create of a comparison table with ATCG contents using [2] [seqtk](https://github.com/lh3/seqtk)
5) Determine quality of sequences (those with N % < 10)
6) Extract all high-quality sequences 
7) Compute a sequence alignment using [3] [MAFFT](https://academic.oup.com/mbe/article/30/4/772/1073398) application 
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
## How to run the scripts

To run our code:
- Download the full content of the directory containing this README file.
- Make sure you have installed R (we developed the code under version: R.4.2.2 (Innocent and Trusting; released on 2022-10-31), so we suggest this or subsequent versions). 
Also, please verify the corresponding programs / dependencies mention in the section below

We used the following version of the programs/libraries :
- MAFFT v.7.0.5
- Seqtk-1.3 (r106)
- See Credits for other programs mention before

## Credits
To use our original or adapted codes, please cite our work https://doi.org/10.21203/rs.3.rs-2285898/v1 as [1], see reference at the end of this document.

The R script codes present in this directory has been written by Rodrigo García-López, Xaira Rivera-Gutierrez and Mauricio Rosales-Rivera.
We gratefully acknowledge the Consorcio Mexicano de Vigilancia Genómica (CoViGen-Mex) and the authors from the originating laboratories responsible for obtaining the specimens and the submitting laboratories from which genetic sequence data were generated and shared via the GISAID initiative included in Supplementary Table 1, as well as the Unidad de Investigación Médica de Yucatán of the Instituto Mexicano del Seguro Social, which collected most of the BW.1 and BA.5.6 samples used in this study. We appreciate the computer assistance provided by Jerome Verleyen, and Juan Manuel Hurtado. The project LANCAD-UNAM-DGTIC-396 of the Dirección General de Cómputo y Tecnologías de la Información (DGTIC-UNAM) provided supercomputing resources in MIZTLI.

The first version of this work appeared on Research Square in November 2022 as preprint: https://doi.org/10.21203/rs.3.rs-2285898/v1. 

### Aditional Cites


### FULL CITE HERE
Rodrigo García-López, Xaira Rivera-Gutiérrez, Mauricio Rosales-Rivera et al. SARS-CoV-2 BW.1, a fast-growing Omicron variant from southeast Mexico bearing relevant escape mutations, 21 November 2022, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2285898/v1]
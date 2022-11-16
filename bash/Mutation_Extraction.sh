# Started: 2022-11-27
# by Xaira Rivera-Gutiérrez for Dr. Arias laboratory at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license

#Comparing weekly tables versus consensus tables

# All available .tsv weekly and consensus tables from BQ.1, BA.5.6.2 and BW.1:

# Obtain a nucleotide position of interest for each week in order to observe weekly nucleotide changes:
# In this example, I'll be using nucleotide position 22942 in weekly descriptions, but it works with any position:

for f in 22W*.tsv ; do echo $f; cat $f | awk -F "\t" '$1 == 22942 {print $0}' ; done > 22942Pos.tsv

# Compare weekly consensus with an overall consensus in order to detect nucleotide substitutions for each position of the consensus:
# In this example, I'll be using all files that are named like '22*_cons.tsv', which are a consensus of the nucleotides per week, and the file 'BW.1_cons-Mex.tsv', which is a consensus of all the sequences
# As a result, you will get a tsv file that shows nucleotides prevalences and percentages in both files only if they are different.

for f in 22*_cons.tsv; do echo $f; paste $f BW.1_cons.tsv | awk -F "\t" '$8 != $18 {print $0}' | grep -v '-' > Compare_$f; done

# Extract nucleotide substitutions per week:
# Here we will work with the same files that the previous example
# As a result, we will get each week nucleotides changes (previous nucleotide, nucleotide position and new nucleotide), without prevalences.
for f in 22*_cons.tsv; do echo $f; paste $f BW.1_cons.tsv | awk -F "\t" '$8 != $18 {print $0}' | grep -v '-' | awk '{print $8 $1 $18}' > Muts_$f; done

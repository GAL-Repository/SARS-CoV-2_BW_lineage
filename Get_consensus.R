# Started: 2022-10-27
# by Rodrigo García-López for Prof. Arias's Viral Analysis Group at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# This script is intended to parse a table containing SARS-CoV-2 nt positions as rows and each nucleotide (including gaps and Ns) that were produced with script Analyze_MexSARS_indels.R
# Test in R
# infile <- "/home/rod/Documents/01_Projects/SARS/Furtivo/BW1/02_Analyses/Lineage_World/BW.1.tsv"
# prefix <- "test.tsv"

### LOAD INPUT AND FUNCTIONS ###
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) { # at least, two arguments should be included: <infile> <prefix_output>
  stop("A minimum of 2 arguments are mandatory: cat table.tsv|Rscript Get_consensus.R <infile> <prefix_output>", call.=FALSE)
}
infile <- as.character(args[1]) # Get the name of the input file (full path is optional)
prefix <- as.character(args[2])  # Get a string handle to create output names (full path is optional)

point_consensus <- function(vect){
	pos_max <- which.max(vect) # Get the position of the largest element. If more than 1 has the same value,  this will take the first one only
	tot <- sum(vect) # Get the sum (all observed items)
	cons <- names(pos_max) # Get the actual Nt
	value_max <- vect[pos_max]
	perc <- round(value_max*100/tot,2)
	out <- cbind(cons,value_max,perc)
	names(out) <- c("Consensus", "Cons=obs", "Cons %")
	return(out)
}

### MAIN ###
print(infile)
df <- read.table(infile, sep="\t",header=T, skip=0, comment.char='',quote="",fill=FALSE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
sub <- df[,c("-", "A", "T", "C", "G", "N")] # Get a subset with the nucleotide contents (or gaps if that was the case)
out <- as.data.frame(t(apply(sub,1, point_consensus)))
write.table(cbind(sub,out),paste0(prefix,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

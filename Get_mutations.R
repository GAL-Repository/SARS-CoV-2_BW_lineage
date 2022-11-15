# Started: 2022-11-07
# by Rodrigo García-López for Prof. Arias's Viral Analysis Group at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# This script is intended to parse a table containing SARS-CoV-2 nt positions as rows and each nucleotide (including gaps and Ns) that were produced with script Analyze_MexSARS_indels.R. The objective is getting a list of mutations per table even if position is repeated.
# Test in R
# infile <- "/home/rod/Documents/01_Projects/SARS/Furtivo/BW1/2022-11-07/02_Analyses/Lineage_World/BW.1.tsv"
# prefix <- "test.tsv"

### LOAD INPUT AND FUNCTIONS ###
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) { # at least, two arguments should be included: <infile> <prefix_output>
  stop("A minimum of 4 arguments are mandatory: Rscript Get_mutations.R <infile> <prefix_output>", call.=FALSE)
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

all_mutations <- function(vect){ # This gets each row, which should have columns: "-", "A", "T", "C", "G", "N", "Ref", "Ref=obs %"
# 	if(vect["Ref=obs %"]<100){
		target <- vect[!names(vect)==as.character(vect["Ref"])][1:5] # Subset nt positions that are different from the reference
		nt <- names(target)[which(target > 0)] # Extract positions having any item other than the ref
		out <- t(sapply(nt, function(x){o <- c(paste0(vect["Ref"],rownames(vect),x),target[x],round(target[x]*100/sum(vect[1:6]),2))})) # code nt changes and count total obs (absolute and relative)
		colnames(out) <- c("Mut(Nt)","Obs", "Rel"); rownames(out) <- NULL
		return(out)
# 	}else{return(c(0,0,0))}
}

### MAIN ###
print(infile)
df <- read.table(infile, sep="\t",header=T, skip=0, comment.char='',quote="",fill=FALSE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
sub <- df[,c("-", "A", "T", "C", "G", "N", "Ref", "Ref=obs %")] # Get a subset with the nucleotide contents (or gaps if that was the case)
# out <- as.data.frame(t(apply(sub,1, point_consensus))) # Get the consensus

outtab <- as.data.frame(matrix(0, ncol=3, nrow=2))
colnames(outtab) <- c("Mut(Nt)","Obs", "Rel")
for(i in 1:nrow(sub)){ # Get all mutations in a table
	if(sub[i,"Ref=obs %"]<100){
		outtab <- rbind(outtab, all_mutations(sub[i,]))
	}
}
outtab <- outtab[3:nrow(outtab),] # Remove the two first columns (no longer used)

write.table(as.matrix(outtab),paste0(prefix,".tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)





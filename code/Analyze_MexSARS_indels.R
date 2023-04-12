# Updated 2023-02-15: Data was updated. Adjusted script for new sequences added for review (added all until 2022-11-30)
# Updated 2023-02-12: All analyses were updated.
# Updated 2022-10-26: For use with BQ.1 and BW.1 variants (wordwide and in Mexico)
# Started 2022-08-31 by Rodrigo Garcia-Lopez for Prof. Arias laboratory at IBt, UNAM
# Under GNU GPLv3 license
# This script takes a linearized and modified fasta file and uses it for creating a table to analyze nt mutations (substitutions) per positions
# Read the table
df <- read.table(unz("analysis_seq_input.zip", "analysis_seq_input.tsv"), header=F, quote="", sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
meta <- read.table("01_data/AllVariants-metadata.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="") # Load metadata with id folio for xref. IDs should be unique
for_removal <- meta[is.na(meta[,3]),2] # Get IDs for items that should be removed because they have no date
# for_removal <- c(for_removal, "EPI_ISL_14262716")# manually append a known issue
# for_removal <- unlist(sapply(for_removal,function(x){grep(x,df[,1])})) # and get the matching indices
# df <- df[-for_removal,]
dim(df)
# [1] 13883     2

remove <- which(table(df[,1])>1) # Get which are repeated
remove <- which(!is.na(remove[df[,1]]))
if(length(remove)>1){df <- df[-remove,]} # to remove them (we cannot determine which is correct)
dim(df)
# [1] 13883     2
refnum <- grep("Wuhan",df[,1]) # UPDATE 2023-02-13: Locate wuhan reference (it is not the first sequence now)
df <- rbind(df[refnum,],df[-refnum,]) # UPDATE 2023-02-13: and make it the first sequence again
nam <- sub("^>","",df[,1]) # We can now save the names of all genomes
df <- df[,-1] # and their sequences, we no longer need anything else
ref <- unlist(strsplit(df[1], split = "")) # Extract the first item, the Wuhan reference
save.image("checkpoint1.Rdata")
# load("checkpoint1.Rdata")

# Repeat, but this time create new matrix with postions as rows, genomes as columns
all <- matrix("-", nrow=length(ref), ncol=length(df))
for(i in 1:length(df)){
	print(i)
	temp <- unlist(strsplit(df[i], split = ""))
	for(j in 1:length(temp)){
		all[j,i] <- temp[j]
	}
}
colnames(all) <- nam
rownames(all) <- 1:nrow(all)
rm(nam)
rm(df)
save.image("checkpoint2.Rdata")

# load("checkpoint2.Rdata")
write.table(all,"Nt_table.tsv", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE) # and export it as a tsv file
# IMPORTANT NOTE: I'd rather use a larger number of rows than columns, so I'll stick with the second method .
#### ANALYSIS ####
# Now, we must determine which are valid positions according to the original genome (first row)
# colnames(all)[1]
# [1] "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"
# table(all[,1])
#    -    A    C    G    T
# 3226 8954 5492 5863 9594

# First, create a map of positions with any nt present (not gaps, "-")
ref_map <- ref <- all[,1] # Start two vectors, one for the actual seq (ref) and one for the indices (ref_map)
ref_name <- colnames(all)[1] # This is just for backup
gap <- ok <- 1 # Start indices at 1
for(i in 1:length(ref)){ # Go position by position in the genome
	if(ref[i]!="-"){ # if no gap is found
		ref_map[i] <- ok # add the index
		ok <- ok+1 # and advance 1 step
		gap <- 1 # Also, reset gap count if it ended
	}else{ # If a gap is found
		ref_map[i] <- paste(ok-1,"i",gap,sep="") # Append the gap count, which also has the position of the last non-gap position
		gap <- gap + 1 # and advance the gap count
	}
}
rownames(all) <- ref_map # replace names with the new reference-based map as it was updated
# Next, evaluate insertions
dir.create("02_Analyses", showWarnings=FALSE) # Create an output folder
pos_ins <- ref_map[grep("i",ref_map)] # add a list of positions with insertions (nt not found in the reference but gaps)
pos_only_in_ref <- ref_map[grep("i",ref_map,invert=TRUE)] # this holds the contrary, those not in insertions.
save.image("checkpoint2.Rdata")
write.table(cbind("Position"=ref_map, "Ins"=grepl("i",ref_map)),"02_Analyses/Insertions_by_position_2023-02-15.tsv", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE) # and export it as a table
# ins_ref <- which(all[1,]=="-") # Detect the actual insertions

### WHOLE TABLE EXPLORATION ###
# by position
std_ambiguities <- function(vect){ # This functions gets a character vector (i.e. one item per genome per position) and ignores rare ambiguities, returning a standard table
	vect <- vect[-1] # The first item is the reference, thus, we remove it
	in_tab <- table(vect)
	valid <- c("-","A","T","C","G","N") # These will be acceptable characters
	out <- sapply(valid, function(x){in_tab[x]}) # read whichever items are present
	names(out) <- valid # Set names to valid items
	out[is.na(out)] <- 0 # Treat NAs as 0s
	out["N"] <- out["N"]+sum(in_tab)-sum(out) # Some items may have other ambiguities, treat them as Ns
	return(out)
}
Nperc <- function(vect){ # This function gets the character vector for a genome (row) or position (column) and creates a table
	matrix(0, ncol=6)
	tab <- std_ambiguities(vect) #

}
position_variation_Nt <- t(sapply(rownames(all),function(x){print(x);std_ambiguities(all[x,])}))
rownames(position_variation_Nt) <- rownames(all)
ref_obs <- sapply(1:nrow(position_variation_Nt),function(x){position_variation_Nt[x,ref[x]]})

write.table(cbind(position_variation_Nt,"Ref"=ref, "Ref=obs"=ref_obs),"02_Analyses/All_genomes_position_composition_2023-02-15.tsv", sep="\t", quote=FALSE, col.names=NA, row.names=TRUE) # and export it as a table

### SUBSET TABLE EXPLORATION ###
load("checkpoint2.Rdata")
# ### FUNCTIONS ###
find_items <- function(queryStrings, dataf, colname){ # IMPORTANT. This is used to avoid partial index matching. Use it to hash a table and use rownames to match a vector. colname may be numeric
	searchVect <- dataf[,colname]
	names(searchVect) <- rownames(dataf)
	out <- searchVect[unlist(queryStrings)]
	return(out)
}
strip_more <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	string <- gsub("\\+","",string) # remove + characters
	string <- gsub("\\.","p",string) # replace "." with p
	string <- gsub("\\(","_",string) # and remove "()", replace with _
	string <- gsub("\\)","_",string) #
	return(string)
}
group_map <- function(vect) { # Get a mutation table (columns = genomes; rows = mutations), and a vector of the same rowsize. Then use the vector to define how groups should be splitted and build a list of df with the table observations per group (a map of positions)
	group_names <- levels(as.factor(vect))
	group_positions <- sapply(group_names,function(x){grep(paste0("^",strip_more(x),"$"),strip_more(vect))})
	return(group_positions)
}
std_ambiguities <- function(vect){ # This functions gets a character vector (i.e. one item per genome per position) and ignores rare ambiguities, returning a standard table. The only difference with the one above is that we must not ignore the first item
# 	vect <- vect[-1] # The first item is the reference, thus, we remove it
	in_tab <- table(vect)
	valid <- c("-","A","T","C","G","N") # These will be acceptable characters
	out <- sapply(valid, function(x){in_tab[x]}) # read whichever items are present
	names(out) <- valid # Set names to valid items
	out[is.na(out)] <- 0 # Treat NAs as 0s
	out["N"] <- out["N"]+sum(in_tab)-sum(out) # Some items may have other ambiguities, treat them as Ns
	return(out)
}
plot_curves <- function(inmat,item_name,group,int) { # get the output table from function std_ambiguities, then create a plot depicting different graphs
	y=100
	ymax <- as.numeric(rownames(inmat)[nrow(inmat)])
	pdf(paste0("02_Analyses/",group,"/",item_name,".pdf"),width=18, height=10)
	oripar <- par(no.readonly=TRUE) # save default params
	par(oma = c(1.5, 0, 1, 6)) # This is just a creative fix to plot the legend outside
	par(mar=c(5.1, 5.1, 4.1, 2.1))
	plot(las=2, 1, type = "n", xlim = c(0, ymax), ylim = c(0, y), frame.plot=FALSE, xlab="Genomic position (Nt)", ylab="Prevalence (%)", xaxt='n', yaxt='n', cex.lab=1.5, main=paste0(group,": ",item_name," | N = ",int), cex.main=2)
	axis(1, at=c(seq(0,ymax,2500),ymax),cex=1.5)
	axis(2, las=1, at=seq(0,y,5),cex.axis=1.5)
	abline(h=seq(5,100,5), col="lightgray", lty=3)
	lines(100-as.numeric(inmat[,"Ref=obs %"]), lwd=2, col="aquamarine2")
	lines(inmat[,"N %"], lwd=2, col="violetred", lty=1)
	lines(inmat[,"- %"], lwd=2, col="royalblue", lty=2)
	Genome <- diff(c(1, 266, 13468, 21556, 21563, 25385, 25393, 26221, 26245, 26473, 26523, 27192, 27202, 27388, 27394, 27756, 27888, 27894, 28260, 28274, 29534, 29558, 29675, ymax)) # This vector holds each interval (all ORFs). Endings were adjusted with +1 for calculations
	cols = c("01 5-UTR"="cornflowerblue", "02 ORF1a"="coral1","03 ORF1b"="turquoise3", "Sep1"=NA, "04 S"="chartreuse2", "Sep2"=NA, "05 ORF3a"="purple", "Sep3"=NA, "06 E"="firebrick", "Sep4"=NA, "07 M"="gold2", "Sep5"=NA, "08 ORF6"="hotpink", "Sep6"=NA, "09 ORF7a"="forestgreen", "10 ORF7b"="darkblue", "Sep7"=NA, "11 ORF8"="darkslategray", "Sep8"=NA, "12 N"="darkorange2", "Sep9"=NA, "14 ORF10"="brown2", "15 3-UTR"="mediumpurple")
	sub_cols <- cols[!is.na(cols)]
	par(oma = c(0, 0, 1, 6)) # Adjust for plotting the genes
	barplot(las=1,cbind(Genome), horiz = TRUE, beside = FALSE, col=cols, border=NA, xaxt='n',add=T)
	par(fig = c(0, 1, 0, 1), oma = c(1.5, 0, 1, 0.5), mar = c(8.1, 0, 7.1, 0), new = TRUE)
	plot(0,0,type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab="", ylab="")
	legend("topright",legend=c("Non-reference","N %","Gap %"), lty=c(1,1,2), col=c("aquamarine2","violetred","royalblue"), title="Lines", lwd=3, bg="white")
	legend("right",legend=names(sub_cols), pch=15, col=sub_cols, pt.cex=1.5, title="Genomic positions", bg="white")
	par(oripar) # reset parameters
	dev.off()
	return(NULL)
}
plot_curves_png <- function(inmat,item_name,group,int) { # get the output table from function std_ambiguities, then create a plot depicting different graphs
	y=100
	ymax <- as.numeric(rownames(inmat)[nrow(inmat)])
	png(paste0("02_Analyses/",group,"/",item_name,".png"),width=18, height=10, units="in", res=100)
	oripar <- par(no.readonly=TRUE) # save default params
	par(oma = c(1.5, 0, 1, 6)) # This is just a creative fix to plot the legend outside
	par(mar=c(5.1, 5.1, 4.1, 2.1))
	plot(las=2, 1, type = "n", xlim = c(0, ymax), ylim = c(0, y), frame.plot=FALSE, xlab="Genomic position (Nt)", ylab="Prevalence (%)", xaxt='n', yaxt='n', cex.lab=1.5, main=paste0(group,": ",item_name," | N = ",int), cex.main=2)
	axis(1, at=c(seq(0,ymax,2500),ymax),cex=1.5)
	axis(2, las=1, at=seq(0,y,5),cex.axis=1.5)
	abline(h=seq(5,100,5), col="lightgray", lty=3)
	lines(100-as.numeric(inmat[,"Ref=obs %"]), lwd=2, col="aquamarine2")
	lines(inmat[,"N %"], lwd=2, col="violetred", lty=1)
	lines(inmat[,"- %"], lwd=2, col="royalblue", lty=2)
	Genome <- diff(c(1, 266, 13468, 21556, 21563, 25385, 25393, 26221, 26245, 26473, 26523, 27192, 27202, 27388, 27394, 27756, 27888, 27894, 28260, 28274, 29534, 29558, 29675, ymax)) # This vector holds each interval (all ORFs). Endings were adjusted with +1 for calculations
	cols = c("01 5-UTR"="cornflowerblue", "02 ORF1a"="coral1","03 ORF1b"="turquoise3", "Sep1"=NA, "04 S"="chartreuse2", "Sep2"=NA, "05 ORF3a"="purple", "Sep3"=NA, "06 E"="firebrick", "Sep4"=NA, "07 M"="gold2", "Sep5"=NA, "08 ORF6"="hotpink", "Sep6"=NA, "09 ORF7a"="forestgreen", "10 ORF7b"="darkblue", "Sep7"=NA, "11 ORF8"="darkslategray", "Sep8"=NA, "12 N"="darkorange2", "Sep9"=NA, "14 ORF10"="brown2", "15 3-UTR"="mediumpurple")
	sub_cols <- cols[!is.na(cols)]
	par(oma = c(0, 0, 1, 6)) # Adjust for plotting the genes
	barplot(las=1,cbind(Genome), horiz = TRUE, beside = FALSE, col=cols, border=NA, xaxt='n',add=T)
	par(fig = c(0, 1, 0, 1), oma = c(1.5, 0, 1, 0.5), mar = c(8.1, 0, 7.1, 0), new = TRUE)
	plot(0,0,type = "n", bty = "n", xaxt = "n", yaxt = "n",xlab="", ylab="")
	legend("topright",legend=c("Non-reference","N %","Gap %"), lty=c(1,1,2), col=c("aquamarine2","violetred","royalblue"), title="Lines", lwd=3, bg="white")
	legend("right",legend=names(sub_cols), pch=15, col=sub_cols, pt.cex=1.5, title="Genomic positions", bg="white")
	par(oripar) # reset parameters
	dev.off()
	return(NULL)
}
non_ins_graphs <- function(map,group) {
	dir.create(paste0("02_Analyses/",group), showWarnings=FALSE)
	top <- max(unlist(lapply(map,length)))
	for(i in names(map)){
		ref_only_non_ins <- ref[as.numeric(names(pos_only_in_ref))]
		names(ref_only_non_ins) <- 1:length(ref_only_non_ins)
		print(i)
		sub_genomes <- map[[i]] # get the index of the subset (e.g. Month)
		position_variation_temp <- t(sapply(as.numeric(names(pos_only_in_ref)),function(x){std_ambiguities(all[x,sub_genomes])})) # go through all the selected genomes and count the total times each nucleotide is observed
		tot <- length(sub_genomes)
		temp <- pos_only_in_ref; names(temp) <- NULL
		rownames(position_variation_temp) <- temp # append the actual position names
		ref_obs <- sapply(as.numeric(rownames(position_variation_temp)),function(x){position_variation_temp[x,ref_only_non_ins[x]]})
		position_variation_temp <- cbind("Position"=rownames(position_variation_temp),position_variation_temp,"Ref"=ref_only_non_ins[pos_only_in_ref], "Ref=obs"=ref_obs, position_variation_temp*100/tot,"Ref=obs %"=ref_obs*100/tot)
		colnames(position_variation_temp)[10:15] <- paste(colnames(position_variation_temp)[10:15], "%")
		write.table(position_variation_temp,paste0("02_Analyses/",group,"/",i,".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE) # and export it as a table
# 		print("ok")
		plot_curves(position_variation_temp, i, group, tot)
# 		plot_curves_png(position_variation_temp, i, group, tot)
	}
}
save.image("checkpoint2.Rdata")

# ### LOAD METADATA ###
# load("checkpoint2.Rdata")
meta <- read.table("01_data/AllVariants-metadata.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=FALSE, fill=T, quote="", row.names=NULL) # Load metadata with id folio for xref. IDs should be unique
rownames(meta) <- meta[,2] # use gisaid ids as names
dim(meta)
# [1] 14162    17
# For both tables to match, we can get use the GISAID ID as names
old_genomeNames <- colnames(all)
colnames(all) <- sapply(old_genomeNames, function(x){strsplit(x,"\\|")[[1]][2]})  # get the gisaid id
# test <- cbind(colnames(all)[2:ncol(all)],meta[colnames(all)[2:ncol(all)],1])
test <- cbind(colnames(all),meta[colnames(all),2])
# UPDATE 2023-02-10: The follwing items have no traceable metadata (Wuhan-Hu-1's missing metadata is expected)
# test[is.na(test[,2]),1]
#                                          MN908947 (Wuhan-Hu-1/2019)
#                                                                  NA
# hCoV-19/Brazil/DF-LACENDF-530033378/2022|EPI_ISL_16293675|2022-11-2
#                                                  "EPI_ISL_16293675"
#               hCoV-19/Greece/345325/2022|EPI_ISL_16383297|2022-11-4
#                                                  "EPI_ISL_16383297"
items_noMetadata <- all[,is.na(test[,2])] # keep them in backup object
all <- all[,!is.na(test[,2])] # and filter them from the main dataframe
# dim(all)
# [1] 29903 13880
meta <- meta[colnames(all),] # UPDATE 2022-02-10: subset and order the rest of the metadata (this should match the all df in length)
# dim(meta)
# [1] 13880    17
all_Months <- find_items(colnames(all),meta, "Month") # Now, months
map_Months <- group_map(all_Months)
save.image("checkpoint2.Rdata")

# ### ACTUAl ANALYSES ###
# Here, we'll get comparisons per type of variant, location and time
# first, for BQ.1
load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BQ.1",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BQ.1-World"
non_ins_graphs(map_Weeks, "Week_BQ.1-World")

# next, BA.5.6.2
load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BA.5.6.2",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BA.5.6.2-World"
non_ins_graphs(map_Weeks, "Week_BA.5.6.2-World")

# Then, BW.1
load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BW.1",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1-World"
non_ins_graphs(map_Weeks, "Week_BW.1-World")

# Finally, BW.1.1
load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BW.1.1",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1.1-World"
non_ins_graphs(map_Weeks, "Week_BW.1.1-World")

# next, BA.5.6.2 in Mexico
load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BA.5.6.2",]
meta <- meta[meta["Region_L2"]=="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BA.5.6.2-Mexico"
non_ins_graphs(map_Weeks, "Week_BA.5.6.2-Mexico")

# next, BA.5.6.2 for World without Mexico
load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BA.5.6.2",]
meta <- meta[meta["Region_L2"]!="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BA.5.6.2-Except_Mexico"
non_ins_graphs(map_Weeks, "Week_BA.5.6.2-Except_Mexico")

load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BW.1",]
meta <- meta[meta["Region_L2"]!="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1-Except_Mexico"
non_ins_graphs(map_Weeks, "Week_BW.1-Except_Mexico")

load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BW.1",]
meta <- meta[meta["Region_L2"]=="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1-Mexico"
non_ins_graphs(map_Weeks, "Week_BW.1-Mexico")

load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BW.1.1",]
meta <- meta[meta["Region_L2"]!="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1.1-Except_Mexico"
non_ins_graphs(map_Weeks, "Week_BW.1.1-Except_Mexico")

load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BW.1.1",]
meta <- meta[meta["Region_L2"]=="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1.1-Mexico"
non_ins_graphs(map_Weeks, "Week_BW.1.1-Mexico")

load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BQ.1",]
meta <- meta[meta["Region_L2"]!="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BQ.1-Except_Mexico"
non_ins_graphs(map_Weeks, "Week_BQ.1-Except_Mexico")

load("checkpoint2.Rdata")
# select the desired items (subset the table)
meta <- meta[meta["Lineage"]=="BQ.1",]
meta <- meta[meta["Region_L2"]=="Mexico",]
all <- all[,rownames(meta)]
# Analyze weeks
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BQ.1-Mexico"
non_ins_graphs(map_Weeks, "Week_BQ.1-Mexico")

# Next, compare by variant
load("checkpoint2.Rdata")
# Analyze lineages
all_Lineage <- find_items(colnames(all),meta, "Lineage") # now, lineages
map_Lineage <- group_map(all_Lineage)
group <- "Lineage_World"
non_ins_graphs(map_Lineage, "Lineage_World")

# repeat, only for Mexico
load("checkpoint2.Rdata")
meta <- meta[meta["Region_L2"]=="Mexico",]
all <- all[,rownames(meta)]
# Analyze lineages
all_Lineage <- find_items(colnames(all),meta, "Lineage") # now, lineages
map_Lineage <- group_map(all_Lineage)
group <- "Lineage_Mexico"
non_ins_graphs(map_Lineage, "Lineage_Mexico")

# repeat, World without Mexico
load("checkpoint2.Rdata")
meta <- meta[meta["Region_L2"]!="Mexico",]
all <- all[,rownames(meta)]
# Analyze lineages
all_Lineage <- find_items(colnames(all),meta, "Lineage") # now, lineages
map_Lineage <- group_map(all_Lineage)
group <- "Lineage_World_except_Mexico"
non_ins_graphs(map_Lineage, "Lineage_World_except_Mexico")

# repeat, only for Yucatan
load("checkpoint2.Rdata")
# meta <- meta[meta["Lineage"]=="BQ.1",]
meta <- meta[meta["Region_L2"]=="Mexico",]
meta <- meta[meta["Region_L3"]=="Yucatan",]
all <- all[,rownames(meta)]
# Analyze lineages
all_Lineage <- find_items(colnames(all),meta, "Lineage") # now, lineages
map_Lineage <- group_map(all_Lineage)
group <- "Lineage_Yucatan"
non_ins_graphs(map_Lineage, "Lineage_Yucatan")

# repeat, only for Yucatan
load("checkpoint2.Rdata")
meta <- meta[meta["Lineage"]=="BQ.1",]
meta <- meta[meta["Region_L2"]=="Mexico",]
meta <- meta[meta["Region_L3"]=="Yucatan",]
all <- all[,rownames(meta)]
# Analyze lineages
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BQ.1-Yucatan"
non_ins_graphs(map_Weeks, "Week_BQ.1-Yucatan")

load("checkpoint2.Rdata")
meta <- meta[meta["Lineage"]=="BA.5.6.2",]
meta <- meta[meta["Region_L2"]=="Mexico",]
meta <- meta[meta["Region_L3"]=="Yucatan",]
all <- all[,rownames(meta)]
# Analyze lineages
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BA.5.6.2-Yucatan"
non_ins_graphs(map_Weeks, "Week_BA.5.6.2-Yucatan")

load("checkpoint2.Rdata")
meta <- meta[meta["Lineage"]=="BW.1",]
meta <- meta[meta["Region_L2"]=="Mexico",]
meta <- meta[meta["Region_L3"]=="Yucatan",]
all <- all[,rownames(meta)]
# Analyze lineages
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1-Yucatan"
non_ins_graphs(map_Weeks, "Week_BW.1-Yucatan")

load("checkpoint2.Rdata")
meta <- meta[meta["Lineage"]=="BW.1.1",]
meta <- meta[meta["Region_L2"]=="Mexico",]
meta <- meta[meta["Region_L3"]=="Yucatan",]
all <- all[,rownames(meta)]
# Analyze lineages
all_Weeks <- find_items(colnames(all),meta, "Week") # now, weeks
map_Weeks <- group_map(all_Weeks)
group <- "Week_BW.1.1-Yucatan"
non_ins_graphs(map_Weeks, "Week_BW.1.1-Yucatan")


# UPDATE 2022-11-16: There is a rare clade of BA.5.6.2 sequences from the world that are phylogenetically closer to Mexican BA.5.6.2 sequences
load("checkpoint2.Rdata")
meta2 <- read.table("01_data/Rare_BA.5.6.2_World_cluster.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="")
sub <- meta2[,2]
sub <- c(sub, "EPI_ISL_15686813")
all <- all[,sub] # subset only those in the clade
# Analyze lineages
all_Lineage <- find_items(colnames(all),meta, "Lineage") # now, weeks
map_Lineage <- group_map(all_Lineage)
group <- "Clade_World_BA.5.6.2_close_to_Mex"
non_ins_graphs(map_Lineage, "Clade_World_BA.5.6.2_close_to_Mex")

load("checkpoint2.Rdata")
meta2 <- read.table("01_data/Rare_BA.5.6.2_World_cluster2_mini.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="")
sub <- meta2[,2]
sub <- c(sub, "EPI_ISL_15686813")
all <- all[,sub] # subset only those in the clade
# Analyze lineages
all_Lineage <- find_items(colnames(all),meta, "Lineage") # now, weeks
map_Lineage <- group_map(all_Lineage)
group <- "Clade_World_BA.5.6.2_close_to_Mex_minicluster"
non_ins_graphs(map_Lineage, "Clade_World_BA.5.6.2_close_to_Mex_minicluster")

load("checkpoint2.Rdata")
meta2 <- read.table("01_data/Rare_BA.5.6.2_World_cluster3_Brazil.tsv", header=T, sep ='\t',stringsAsFactors = FALSE, check.names=F, fill=T, quote="")
sub <- meta2[,2]
sub <- c(sub, "EPI_ISL_15686813")
all <- all[,sub] # subset only those in the clade
# Analyze lineages
all_Lineage <- find_items(colnames(all),meta, "Lineage") # now, weeks
map_Lineage <- group_map(all_Lineage)
group <- "Clade_World_BA.5.6.2_close_to_Mex_Brazil2seqs"
non_ins_graphs(map_Lineage, "Clade_World_BA.5.6.2_close_to_Mex_Brazil2seqs")

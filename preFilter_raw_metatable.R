# Started: 2022-11-08
# by Rodrigo García-López for Prof. Arias's Viral Analysis Group at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# This script is intended to take a subset of the GISAID metadata files and exctract locations
# tsv <- "01_data/gisaid_hcov-19_2022_11_meta.tsv" # This had bad coding for new lines
tsv <- "01_data/gisaid_hcov-19_MexCov2_Template.tsv"
df <- read.table(tsv, header=T, quote="", sep="\t", check.names=FALSE, stringsAsFactors=FALSE, fill=TRUE)
# df <- read.table("01_data/2022_04_27_DB_MexCov2.tsv",header=T, sep='\t', skip=0, comment.char='',fill=T, check.names=FALSE, stringsAsFactors = FALSE)
# df <- read.table(unz("01_data/17Enero2022_Metadata_BDMexCov2_FINAL.zip", "17Enero2022_Metadata_BDMexCov2_FINAL.tsv"), header=T, quote="\"", sep="\t", check.names=FALSE, fill=TRUE)
print("Starting dimensions:")
dim(df)
# [1] 1018   15
colnames(df)
#  [1] "name"               "accession_id"       "type"
#  [4] "passage_details"    "collection_date"    "submission_date"
#  [7] "location"           "covv_add_location"  "host"
# [10] "covv_add_host_info" "gender"             "age"
# [13] "patient_status"     "specimen"           "Lineage"
target <- c("name","accession_id","collection_date", "submission_date", "location", "gender", "age","Lineage")
df <- df[,target]
nam <- c("Virus name","Accession ID", "Collection date", "Submission date", "Location", "Gender","Patient age","Lineage")
colnames(df) <- nam

# for_removal <- c("Host", "Additional location information", "Sampling strategy", "Last vaccinated", "Passage", "Specimen", "AA Substitutions", "Clade", "Additional host information", "Patient status") # DEPRECATED
# df <- df[,-sapply(for_removal, function(x) {grep(x, colnames(df))})] # remove items
dim(df)
# [1] 7404    8
### Fix and append columns ###
# Remove rare punctuations or accents
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	return(string)
}
names(df) <- strip(names(df))
df[,"Region_L1"] <- sapply(df[,"Location"], function(x){spl <- strsplit(x, " / ");spl[[1]][1]})
df[,"Region_L2"] <- sapply(df[,"Location"], function(x){spl <- strsplit(x, " / ");spl[[1]][2]})
df[,"Region_L3"] <- sapply(df[,"Location"], function(x){spl <- strsplit(x, " / ");spl[[1]][3]})

for_removal <- c("Location")
df <- df[,-sapply(for_removal, function(x) {grep(x, colnames(df))})] # remove items

print("Ending dimensions:")
dim(df)
# [1] 4846   10
write.table(df, "01_data/preFiltered.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

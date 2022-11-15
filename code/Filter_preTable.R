# Started: 2022-10-26
# by Rodrigo García-López for Prof. Arias's Viromics Group at IBt, UNAM, Cuernavaca, Mexico.
# Under GNU GPLv3 license
# This script is intended for usage with tables produced by script preFilter_raw_metatable.R. It fixes some known issues with metadata from gisaid, and adds some other fields such as Week, Month, age group, etc.
df <- read.table("01_data/preFiltered.tsv", header=T, quote="", sep="\t", check.names=FALSE)
# df <- read.table("01_data/Xai_tabla_filtrada_Delta_fixed5.tsv",header=T, sep='\t', skip=0, comment.char='',fill=FALSE, check.names=FALSE, stringsAsFactors = FALSE)
print("Starting total items:")
dim(df)
# [1] 4846   10
# Add a category of simplified lineages to group BW.1 and BA.5.6.2
table(df[,"Lineage"])
# BA.5.6.2     BQ.1     BW.1
#      850     3957       39
df[,"Simp_Lineage"] <- df[,"Lineage"]
df[df[,"Simp_Lineage"] == "BA.5.6.2","Simp_Lineage"] <- "BA.5.6.2.x"
df[df[,"Simp_Lineage"] == "BW.1","Simp_Lineage"] <- "BA.5.6.2.x"
table(df[,"Simp_Lineage"])
# BA.5.6.2.x       BQ.1
#        889       3957

# Fix dates
df[,"Collection date"] <-  as.Date(df[,"Collection date"],"%d/%m/%Y")
df[,"Submission date"] <- as.Date(df[,"Submission date"],"%d/%m/%Y")
# Remove all those newer than the last 16 nov date
print("Last date in the whole table:")
test_date <- max(as.Date(df[,"Collection date"]));test_date
# [1] "2022-08-24"
# Fix gender
df[df[,"Gender"] == "71","Gender"] <- "unknown"
df[df[,"Gender"] == "O","Gender"] <- "unknown"
df[,"Gender"] <- sub("[Ff]emale","F",df[,"Gender"])
df[,"Gender"] <- sub("[Mm]ale","M",df[,"Gender"])
df[,"Gender"] <- sub("[Uu]nknown","u",df[,"Gender"])
# Fix age
df[, "Patient age"] <- round(as.numeric(df[, "Patient age"])) # remove all non-numeric and round to integers
df[is.na(df[, "Patient age"]), "Patient age"] <- "u"
# DEPRECATED: START
# # # # Fix identifiers
# # # df[,"ID Folio"] <- sub("Sin información","u",df[,"ID Folio"])
# # # df[,"SINAVE ID"] <- sub("Sin información","u",df[,"SINAVE ID"])
# # # df[,"SINOLAVE ID"] <- sub("Sin información","u",df[,"SINOLAVE ID"])
# DEPRECATED: END
# Fix type of patient

# Append the month
df[,"Month"] <- format(as.Date(df[,"Collection date"]),"%Y-%m(%b)")
pdf("01_data/Monthly_genomes.pdf",width=7)
	par(oma=c(5,1,1,1))
	barplot(las=2,table(df[,"Month"]), border=F, col="cornflowerblue", ylab="Total Genomes")
dev.off()
# Append week
	# To do this, first create a calendar (week 1 starts on sunday of the first complete week in the year and reaches up to 52)
	week <- data.frame("week"=c(paste0("20W",sprintf('%0.2d', rep(1:52,each=7))),paste0("21W",sprintf('%0.2d', rep(1:52,each=7))),paste0("22W",sprintf('%0.2d', rep(1:52,each=7)))))
	# Now rename the rows to use them as a hash (dictionary)
	rownames(week) <- seq(as.Date("2020/01/05"), as.Date("2022/12/31"), by="day")[1:nrow(week)]
	# write.table(week, "week_calendar.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	# Now used the newly created dictionary to append the corresponding week
	df[,"Week"] <- week[as.character(df[,"Collection date"]),1]
pdf("01_data/Weekly_genomes.pdf",width=7)
	par(oma=c(5,1,1,1))
	barplot(las=2,table(df[,"Week"]), border=F, col="chartreuse3", ylab="Total Genomes")
dev.off()
# Next, append age categories
	# First create dictionary
	age1 <- data.frame("Age"=as.character(rep(paste0(sprintf('%0.2d',1:14),":",sapply(seq(0,130,10),function(x) paste0(x,"-",x+9))),each=10)), stringsAsFactors = FALSE)
	rownames(age1) <- 0:(nrow(age1)-1) # and use age1s as rownames
	age1[as.numeric(rownames(age1))>=100,1] <- "11:100+" # This was later added to collate all 100 an over
	df[,"Age_range_by10"] <- age1[as.character(df[,"Patient age"]),1] # Now append it to the main dataframe
	# Repeat with smaller bins
	age2 <- data.frame("Age"=rep(paste0(sprintf('%0.2d',1:28),":",sapply(seq(0,135,5),function(x) paste0(x,"-",x+4))),each=5), stringsAsFactors = FALSE)
	rownames(age2) <- 0:(nrow(age2)-1) # and use age2s as rownames
	age2[as.numeric(rownames(age2))>=100,1] <- "21:100+" # This was later added to collate all 100 an over
	df[,"Age_range_by5"] <- age2[as.character(df[,"Patient age"]),1] # Now append it to the main dataframe
	# Repeat with a new classification
	age3 <- data.frame("age3"=c(rep("1:0-12", 13), rep("2:13-25",13),rep("3:26-45",20), rep("4:46-60",15), rep("5:61-74",14), rep("6:75+",55)), stringsAsFactors = FALSE)
	rownames(age3) <- 0:(nrow(age3)-1) # and use age3s as rownames
	df[,"Age_range_manual"] <- age3[as.character(df[,"Patient age"]),1] # Now append it to the main dataframe
	# Added a new age scheme matching vaccination
	age4 <- data.frame("age4"=c(rep("1:0-17", 18), rep("2:18-29",12),rep("3:30-39",10), rep("4:40-49",10), rep("5:50-59",10), rep("6:60+",100)), stringsAsFactors = FALSE)
	rownames(age4) <- 0:(nrow(age4)-1) # and use age3s as rownames
	df[,"Age_vac"] <- age4[as.character(df[,"Patient age"]),1] #
pdf("01_data/age5_genomes.pdf")
	barplot(las=2,table(df[,"Age_range_by10"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_range_by5"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_range_manual"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_vac"]),border=F, col="coral1")
dev.off()
	fix_age <- function(vect){ # Ages may be missing. Fix all and convert them to character.
		vect[vect=="u"]=NA # change "u"s to NAs to prevent warnings
		vect <- sprintf('%0.3d',as.numeric(vect))
		vect[vect=="NA"]="u" # Revert NAs (they are now text) to "u"s
		return(vect)
	}
	df[,"Patient age"] <- fix_age(df[,"Patient age"])
print("Resulting total items:")
dim(df)
# [1] 4846   17
write.table(df, "01_data/AllVariants-metadata.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

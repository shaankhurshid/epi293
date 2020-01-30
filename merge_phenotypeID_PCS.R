
###################################
## This code is used to generate phenotype file for SNPTest and ensure sample order the same as imputed data
## You may need to do necessary modification for your own project
###################################

# Dependencies
library(data.table)
library(stringr)

#########################################################################################################
# PART 1: WITH REFERENCE
#########################################################################################################

################################### Process imputed data
# Read in imputed samples
sam.w <- read.table("~/Documents/MPH/EPI293/final/final_QC_with.imputed.samples")[-c(1,2),]

# Extract the IDs
id.w = as.numeric(substring(sam.w[,1],1,as.numeric(regexpr("_",sam.w[,1]))-1))

# Add the IDs back to sam.w
sam.w <- cbind(sam.w,id.w)

# Set DT for later merge
setDT(sam.w)

# Read in original PLINK phenotype file
phe = fread("~/Documents/MPH/EPI293/final/Quantitative_trait_PLINK_format.txt")
colnames(phe) <- c('FID','IID','trait')

# Read in PC file
mds <- read.csv("~/Documents/MPH/EPI293/final/mds.csv")
colnames(mds) <- c('FID','IID','PC1','PC2')

# Set keys and merge
phe <- merge(phe,mds,by='IID')

# Cleanup
phe <- phe[,FID.y := NULL]
setnames(phe,'FID.x','FID')

# Reorder combined phenotype/PC file by index obtained from imputation file
ordered <- phe[id.w]

# Rejoin index for reference
ordered <- cbind(ordered,id.w)

# Recapitulate V3 (missing) from imputed dataset
setDF(ordered)
ordered <- merge(ordered,sam.w[,c('V3','id.w')],by='id.w',sort=FALSE)

# Recapitulate IDs from imputed dataset
ordered <- cbind(ordered,sam.w[,c('V1','V2')])

# Format for SNPTest
# Pick the things SNPTest likes
ordered <- ordered[,c('V1','V2','V3','trait','PC1','PC2')]

# Has to be a matrix...ugh
ordered <- as.matrix(ordered)
# Cleaning empty spaces?
trait_clean <- str_replace_all(ordered[,4],' ','')
PC1_clean <- str_replace_all(ordered[,5],' ','')
PC2_clean <- str_replace_all(ordered[,6],' ','')

# Combine back
ordered <- cbind(ordered,trait_clean,PC1_clean,PC2_clean)
ordered <- ordered[,c("V1","V2","V3","trait_clean","PC1_clean","PC2_clean")]

typ = matrix(c("0","0","0","P","C","C"),nrow=1)
colnames(typ) = colnames(ordered) = c("ID_1","ID_2","missing","trait",'PC1','PC2')

output <- rbind(typ,ordered)

# Write out
write.table(output, "~/Documents/MPH/EPI293/final/Quantitative_trait_SNPTEST_format_with_MDS.txt", col.names=T, row.names=F, quote=F, sep=" ")

#########################################################################################################
# PART 2: WITHOUT REFERENCE
#########################################################################################################

################################### Process imputed data
# Read in imputed samples
sam.w <- read.table("~/Documents/MPH/EPI293/final/final_QC_without.imputed.samples")[-c(1,2),]

# Extract the IDs
id.w = as.numeric(substring(sam.w[,1],1,as.numeric(regexpr("_",sam.w[,1]))-1))

# Add the IDs back to sam.w
sam.w <- cbind(sam.w,id.w)

# Set DT for later merge
setDT(sam.w)

# Read in original PLINK phenotype file
phe = fread("~/Documents/MPH/EPI293/final/Quantitative_trait_PLINK_format.txt")
colnames(phe) <- c('FID','IID','trait')

# Read in PC file
mds <- read.csv("~/Documents/MPH/EPI293/final/mds.csv")
colnames(mds) <- c('FID','IID','PC1','PC2')

# Set keys and merge
phe <- merge(phe,mds,by='IID')

# Cleanup
phe <- phe[,FID.y := NULL]
setnames(phe,'FID.x','FID')

# Reorder combined phenotype/PC file by index obtained from imputation file
ordered <- phe[id.w]

# Rejoin index for reference
ordered <- cbind(ordered,id.w)

# Recapitulate V3 (missing) from imputed dataset
setDF(ordered)
ordered <- merge(ordered,sam.w[,c('V3','id.w')],by='id.w',sort=FALSE)

# Recapitulate IDs from imputed dataset
ordered <- cbind(ordered,sam.w[,c('V1','V2')])

# Format for SNPTest
# Pick the things SNPTest likes
ordered <- ordered[,c('V1','V2','V3','trait','PC1','PC2')]

# Has to be a matrix...ugh
ordered <- as.matrix(ordered)
# Cleaning empty spaces?
trait_clean <- str_replace_all(ordered[,4],' ','')
PC1_clean <- str_replace_all(ordered[,5],' ','')
PC2_clean <- str_replace_all(ordered[,6],' ','')

# Combine back
ordered <- cbind(ordered,trait_clean,PC1_clean,PC2_clean)
ordered <- ordered[,c("V1","V2","V3","trait_clean","PC1_clean","PC2_clean")]

typ = matrix(c("0","0","0","P","C","C"),nrow=1)
colnames(typ) = colnames(ordered) = c("ID_1","ID_2","missing","trait",'PC1','PC2')

output <- rbind(typ,ordered)

# Write out
write.table(output, "~/Documents/MPH/EPI293/final/Quantitative_trait_SNPTEST_format_without_MDS.txt", col.names=T, row.names=F, quote=F, sep=" ")

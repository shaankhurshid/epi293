# This is a script to perform GWAS using gentype
# This is part of HSPH EPI 293

# Dependencies
library(data.table)

# Original GWAS
gwas <- fread("~/Documents/MPH/EPI293/plink.assoc.linear")

# QC of original data
imiss <- fread("~/Documents/MPH/EPI293/plink.imiss")
lmiss <- fread("~/Documents/MPH/EPI293/plink.lmiss")
freq <- fread("~/Documents/MPH/EPI293/plink.frq")
hwe <- fread("~/Documents/MPH/EPI293/plink.hwe")

bad_snps <- fread("~/Documents/MPH/EPI293/bad_SNP_uniq.txt")

# Summary data
## Frequency missing
mean(imiss$F_MISS); sd(imiss$F_MISS)  ## per individual
mean(lmiss$F_MISS); sd(imiss$F_MISS)  ## per SNP

# Sequential QC
all_snps <- hwe$SNP

## Missingness
c(nrow(lmiss[F_MISS > 0.05]),nrow(lmiss[F_MISS > 0.05])/nrow(lmiss))  ## N below 5%
missing <- lmiss[F_MISS > 0.05]
qc_1 <- all_snps[!(all_snps %in% missing$SNP)]

## HWE
c(nrow(hwe[P < 1e-6]),nrow(hwe[P < 1e-6])/nrow(hwe))  ## HWE below p=0.05
hwe <- hwe[P < 1e-6]
qc_2 <- qc_1[!(qc_1 %in% hwe$SNP)]

## MAF
c(nrow(freq[MAF < 0.02]),nrow(freq[MAF < 0.02])/nrow(freq))  ## N below 5%
maf <- freq[MAF < 0.02]
qc_3 <- qc_2[!(qc_2 %in% maf$SNP)]

## QQ plot
obs <- -log10(gwas$P)
exp <- -log10(rank(gwas$P)/(nrow(gwas)+1))

png("~/Documents/MPH/EPI293/qq_plot.png",width=600, height=500)

plot(exp,obs,main="QQ plot",xlab="Expected -log p-value",ylab="Observed -log p-value")
lines(c(0,10),c(0,10),col="red")
text(1,3.5,'lambda = 1.06')

dev.off()

## Genomic Control pamameter lambda
qchisq(median(gwas$P,na.rm=T),df=1,lower.tail=F)/0.455

#Manhattan Plot
png("~/Documents/MPH/EPI293/mh_plot.png",width=600, height=450)

plot(gwas$BP/1e6,-log10(gwas$P),xlab="position",ylab="-log p-value")

dev.off()

#Near the top hit
png("~/Documents/MPH/EPI293/mh_zoomed.png",width=600, height=450)
plot(gwas$BP/1e6,-log10(gwas$P),xlim=c(163.5e5/1e6,166e5/1e6),xlab="position",ylab="-log p-value",pch=19)
segments(163e5/1e6,8,166.9e5/1e6,8,col='black',lty=2)
dev.off()

#SNPs with GW significance
sig_original <- gwas[-log10(gwas$P) > 8]
  
##############################################################
##############################################################

## read in association analysis result using imputed dosage

w <- read.table("trait_lab3_QC_with.imputed.txt",header=T)[,-c(2,3,7:19,25)]
wo <- read.table("trait_lab3_QC_without.imputed.txt",header=T)[,-c(2,3,7:19,25)]
setDT(w); setDT(wo)

## read information file for imputed dosage 

iw <- read.table("lab3_QC_with.imputed.info",header=T)[,c(1:3,7)]
iwo <- read.table("lab3_QC_without.imputed.info",header=T)[,c(1:3,7)]
setDT(iw); setDT(wo)

## add prediction R2 to association result
## Add SNP name column to association analysis result
w[,SNP := paste0(alternate_ids,":",position)]
wo[,SNP := paste0(alternate_ids,":",position)]

## Merge
setkey(w,'SNP'); setkey(wo,'SNP')
w <- w[iw,nomatch=0]
wo <- wo[iwo,nomatch=0]

## Remove p-value = NA
w <- w[!is.na(frequentist_add_pvalue)]
wo <- wo[!is.na(frequentist_add_pvalue)]

# Plot imputation quality (Rsq)
png("~/Documents/MPH/EPI293/rsq_hist_w.png",width=600, height=450)

hist(w$Rsq,main='',xlab='')
title(main="Imputation Quality (with reference panel)",xlab='R-squared')

dev.off()

png("~/Documents/MPH/EPI293/rsq_hist_wo.png",width=600, height=450)

hist(wo$Rsq,main='',xlab='')
title(main="Imputation Quality (without reference panel)",xlab='R-squared')

dev.off()

## filter out SNPs with small Rsq value < 0.3
w_poor_imp <- w[Rsq <= 0.3]
w <- w[Rsq > 0.3]

wo_poor_imp <- wo[Rsq <= 0.3]
wo <- wo[Rsq > 0.3]

## plot and compare GWAS results based on imputed SNPs and genotype SNPs
png("~/Documents/MPH/EPI293/gwas_compare.png",width=600, height=450)
col<-c("black","#1f78b4","#fccde5")

plot(0,0,xlab="Position",ylab="-log p-value",xlim=c(163.5e5/1e6,166.7e5/1e6),
     ylim=c(0,max(-log10(w[,"frequentist_add_pvalue"]))),pch=1)

points(w$position/1e6,-log10(w$frequentist_add_pvalue),col=col[2],pch=19,cex=0.8)
points(wo$position/1e6,-log10(wo$frequentist_add_pvalue),col=col[3],pch=19,cex=0.8)
points(gwas$BP/1e6,-log10(gwas$P),pch=19,col=col[1],cex=0.8)

legend(16.565,10.6,legend=c("No imputation","Imputation (w/ reference)","Imputation (w/o reference)"),
       col=col,pch=19,bty='n',x.intersp=0.5,y.intersp=1)

segments(163e5/1e6,8,166.9e5/1e6,8,col='black',lty=2)

dev.off()

# SNPs with GW significance
## No reference
sig_wo <- wo[-log10(wo$frequentist_add_pvalue) > 8][order(frequentist_add_pvalue)]

## With reference
sig_w <- w[-log10(w$frequentist_add_pvalue) > 8][order(frequentist_add_pvalue)]




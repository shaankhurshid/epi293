
## plot the PC value computed by EIGENSTRAT

eigen <- read.table("~/Documents/MPH/EPI293/final/GWA.pca.evec",header=F,as.is=T)
head(eigen[,1:10])

png('~/Documents/MPH/EPI293/final/final_eigen.png',height=500,width=600)
plot(eigen$V2,eigen$V3,main="EIGENSTRAT plot",xlab="PC1",ylab="PC2")
dev.off()

## plink.raw is generated using the following command
## plink --ped lab2.ped --map lab2.map --recodeA

## MDS plot of lab2 data

a=read.table("~/Documents/MPH/EPI293/final/plink.raw",header=T,as.is=T)
head(a[,1:10])

row.names(a)=paste(a$FID,a$IID,sep=".")

data=a[,7:7839]

# FYI. R code to generate principal component vector, not working with missing genotype
#pc=prcomp(data)

# for genotype with missing data, we could use MDS
# generate MDS vector
library(MASS)
d=dist(data)
mds = isoMDS(d,k=10)
p=as.data.frame(mds$points)
p$names=rownames(p)

png('~/Documents/MPH/EPI293/final/final_mds.png',height=500,width=600)
plot(p$V1,p$V2,main="MDS plot",xlab="MDS1",ylab="MDS2")
dev.off()

## Compare PCA and MDS 
png('~/Documents/MPH/EPI293/final/final_compare.png',height=500,width=700)
layout(matrix(1:6,nrow=2))
plot(p$V1,eigen$V2)
plot(p$V2,eigen$V3)
plot(p$V3,eigen$V4)
plot(p$V4,eigen$V5)
plot(p$V5,eigen$V6)
plot(p$V6,eigen$V7)
dev.off()

## output the MDS factors to text files, which can be used as input for PLINK in order to correct for population stratification in GWAS
mds <- cbind(a$FID,a$IID,p$V1,p$V2)
write.table(mds,"~/Documents/MPH/EPI293/final/MDS_covariates_for_PLINK.txt",col.name=F,row.name=F,quote=F,sep="\t")

## (1) MDS plot for simulated population samples

a=read.table("~/Documents/MPH/EPI293/output_sequence_1.txt",header=F,as.is=T,colClasses="character")

m = matrix(unlist(strsplit(a$V2,split="")),nrow=40,byrow=TRUE)

rownames(m)=a$V1

d=dist(m)+1e-5
mds = isoMDS(d,k=2)
p=as.data.frame(mds$points)
p$names=rownames(p)

png("~/Documents/MPH/EPI293/mds_sequence_1.png",height=500,width=600)
plot(p$V1,p$V2,main="MDS plot simulation 1",xlab="MDS1",ylab="MDS2")
points(p[a$V1=="POP1:",]$V1,p[a$V1=="POP1:",]$V2,col="red",pch=19)
dev.off()

## (2) MDS plot for simulated population samples

a=read.table("~/Documents/MPH/EPI293/output_sequence_2.txt",header=F,as.is=T,colClasses="character")

m = matrix(unlist(strsplit(a$V2,split="")),nrow=400,byrow=TRUE)

d=dist(m)+1e-5
mds = isoMDS(d,k=2)
p=as.data.frame(mds$points)
p$names=rownames(p)

png("~/Documents/MPH/EPI293/mds_sequence_2.png",height=500,width=600)
plot(p$V1,p$V2,main="MDS plot simulation 2",xlab="MDS1",ylab="MDS2")
points(p[a$V1=="POP1:",]$V1,p[a$V1=="POP1:",]$V2,col="red",pch=19)
dev.off()

## (3) MDS plot for simulated population samples

a=read.table("~/Documents/MPH/EPI293/output_sequence_3.txt",header=F,as.is=T,colClasses="character")

m = matrix(unlist(strsplit(a$V2,split="")),nrow=400,byrow=TRUE)

d=dist(m)+1e-5
mds = isoMDS(d,k=2)
p=as.data.frame(mds$points)
p$names=rownames(p)

png("~/Documents/MPH/EPI293/mds_sequence_3.png",height=500,width=600)
plot(p$V1,p$V2,main="MDS plot simulation 3",xlab="MDS1",ylab="MDS2")
points(p[a$V1=="POP1:",]$V1,p[a$V1=="POP1:",]$V2,col="red",pch=19)
dev.off()

## (4) MDS plot for simulated population samples

a=read.table("~/Documents/MPH/EPI293/output_sequence_4.txt",header=F,as.is=T,colClasses="character")

m = matrix(unlist(strsplit(a$V2,split="")),nrow=400,byrow=TRUE)

d=dist(m)+1e-5
mds = isoMDS(d,k=2)
p=as.data.frame(mds$points)
p$names=rownames(p)

png("~/Documents/MPH/EPI293/mds_sequence_4.png",height=500,width=600)
plot(p$V1,p$V2,main="MDS plot simulation 4",xlab="MDS1",ylab="MDS2")
points(p[a$V1=="POP1:",]$V1,p[a$V1=="POP1:",]$V2,col="red",pch=19)
dev.off()

## (5) MDS plot for simulated population samples

a=read.table("~/Documents/MPH/EPI293/output_sequence_5.txt",header=F,as.is=T,colClasses="character")

m = matrix(unlist(strsplit(a$V2,split="")),nrow=400,byrow=TRUE)

d=dist(m)+1e-5
mds = isoMDS(d,k=2)
p=as.data.frame(mds$points)
p$names=rownames(p)

png("~/Documents/MPH/EPI293/mds_sequence_5.png",height=500,width=600)
plot(p$V1,p$V2,main="MDS plot simulation 5",xlab="MDS1",ylab="MDS2")
points(p[a$V1=="POP1:",]$V1,p[a$V1=="POP1:",]$V2,col="red",pch=19)
dev.off()

############################################
## Genomic Control parameter under the NULL

disease1 = rbinom(200,1,0.3)
disease2 = rbinom(200,1,0.1)
summary(disease1)
summary(disease2)
disease = c(disease1,disease2)

f=function(v){
	fisher.test(disease,v)$p.value
}

## data from simulation 5 ##
a=read.table("~/Documents/MPH/EPI293/output_sequence_5.txt",header=F,as.is=T,colClasses="character")

m = (matrix(as.numeric(unlist(strsplit(a$V2,split=""))),nrow=400,byrow=TRUE))

# all SNP
pval = apply(m,2,f)
summary(pval)
hist(pval)

# common SNP
maf = apply(m,2,mean,na.rm=T)
maf[maf>0.5] = 1-maf[maf>0.5]
summary(maf)

## QQ plot
png("~/Documents/MPH/EPI293/qq_sim5_all.png",height=500,width=600)
obs = -log10(pval)
exp = -log10(rank(pval)/(length(pval)+1))
plot(exp,obs,main="QQ plot",xlab="Expected -log p",ylab="Observed -log p")
lines(c(0,10),c(0,10),col="red")
dev.off()

## Genomic Control pamameter lambda
lambda_all <- qchisq(median(pval,na.rm=T),df=1,lower.tail=F)/0.455

## Focus on common SNPs
pval1 = pval[maf>0.05]
summary(pval1)
hist(pval1)

## QQ plot

png("~/Documents/MPH/EPI293/qq_sim5_common.png",height=500,width=600)
obs = -log10(pval1)
exp = -log10(rank(pval1)/(length(pval1)+1))
plot(exp,obs,main="QQ plot",xlab="Expected -log p",ylab="Observed -log p")
lines(c(0,10),c(0,10),col="red")
dev.off()

## Genomic Control pamameter lambda
lambda_common <- qchisq(median(pval1,na.rm=T),df=1,lower.tail=F)/0.455



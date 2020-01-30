# Format SNPTEST output for LDSC
setDT(w)

a <- w

w <- w[,c('SNP','alternate_ids','position','alleleA','alleleB','frequentist_add_beta_1',
          'frequentist_add_se_1','frequentist_add_pvalue','frequentist_add_info')]

w[,SNP := str_replace(SNP,'20:','')]
colnames(w) <- c('snpid','hg20chr','bp','a1','a2','beta','se','pval','info')

a <- as.matrix(w)

write.table(a,file='~/Documents/MPH/EPI293/final/imputed_with_for_ldsc2.txt',sep='\t',quote=FALSE)

write.table(gwas_adjusted,file='~/Documents/MPH/EPI293/final/gwas_adjusted_for_ldsc.txt',sep='\t',quote=FALSE)

bim <- read.table('~/Documents/MPH/EPI293/final/final_QC.bim')
bim <- bim[bim$V2 %in% gwas_adjusted$SNP,]
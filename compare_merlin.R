# Compare output of untransformed versus inverse normal transformed eQTL association

# Dependencies
library(data.table)

# load eQTL tables
untransformed <- fread(file='~/Documents/MPH/EPI293/assoc_noInv-fastassoc-chr17.tbl')
transformed <- fread(file='~/Documents/MPH/EPI293/assoc-fastassoc-chr17.tbl')

# Ratio of heritability estimates and p-values
## Join
setkey(untransformed,'SNP','TRAIT'); setkey(transformed,'SNP','TRAIT')
combined <- transformed[untransformed,nomatch=0]

## Average ratio of transformed/untransformed
combined[,':='(p_ratio = PVALUE/i.PVALUE,
               h_ratio = H2/i.H2)]

summary(combined$p_ratio)
summary(combined$h_ratio)
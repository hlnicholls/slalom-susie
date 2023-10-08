setwd('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/')
library(data.table)
library(dplyr)

############################################
# 1. Getting R2 column for sumstats tsv file
############################################
16_50724938_G_A
#LD for all variants in a locus calculated by plink in LD_with_lead_for_all_SNPs_per_locus.ipynb
ld <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/ld/APOE_subset_for_ld_calculation.ld')

#Filter for only LD with lead SNP
filtered_df <- ld %>% filter(SNP_A == '19:44908822:C:T' | SNP_B == '19:44908822:C:T')

# Get list of all associated variants from sumstats tsv.gz with : separators
sumstat <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/APOE_LDL_locus_sumstat_flip_check.txt.gz')
snps <- select(sumstat, ID)
snps$ID <- gsub("_(\\d+)_([A-Z])_([A-Z])", ":\\1:\\2:\\3", snps$ID)

merged <- merge(ld, snps, by.x = "SNP_B", by.y = "ID")

colnames(merged)[1] <- 'ID'
merged$ID <- gsub(":", "_", merged$ID)
merged <- select(merged, 'ID', 'R2')
df <- merge(sumstat, merged, by='ID', all.x=TRUE)

############################
# 2. Calculate 't_dentist_s'
############################

n_samples_constant <- 94595
df$r <- (n_samples_constant * sum(df$R2, na.rm = TRUE)) / (n_samples_constant * sum(!is.na(df$R2)))

lead <- df[df$ID == "19_44908822_C_T", ]
lead_z <- lead$beta / lead$se

df$t_dentist_s <- ((df$beta / df$se) - df$r* lead_z)^2 / (1 - df$r^2)
df$dentist_outlier <- as.integer(df$t_dentist_s < 1e-4 & df$R2 > 0.6)
fwrite(df, 'APOE_LDL_locus_sumstat_with_dentist.txt.gz', sep = '\t')

############################

# Original Python code for dentist in SLALOM:

#df["r"] = np.nansum(n_samples * ld, axis=1) / np.nansum(n_samples * ~np.isnan(ld), axis=1)

#df["t_dentist_s"] = ((df.beta / [df.se](http://df.se/)) - df.r * lead_z) ** 2 / (1 - df.r ** 2)

#Outliers identified by threshold:
#P DENTIST-S < 10-4 and r2 > 0.6

############################

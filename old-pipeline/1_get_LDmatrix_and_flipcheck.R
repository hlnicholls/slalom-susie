library(data.table)
library(dplyr)
library(tidyverse)
library(reticulate)

# Input variables for script: plink file path and target name, target_chrom, target_pos
# Also col numbers need manually checking at line 56

gwas_file_path <- "/Users/hn9/Documents/Analysis/FM-comparison/gwas-examples/NOD2-Crohns/28067908-GCST004132-EFO_0000384.h.tsv.gz"
target <- "NOD2"
target_chrom <- "16"
target_pos <- 50724938
start_pos <- target_pos - 500000
end_pos <- target_pos + 500000

# Define the name of the conda environment with bioconda plink installed
# or run from line 21 with path to plink installation
conda_env <- "fm_env"

command <- paste("/Users/hn9/anaconda3/envs/fm_env/bin/plink",
                 "--bfile", paste0("/Users/hn9/Documents/Analysis/FM-comparison/ukb_v3_downsampled10k/ukb_v3_chr", target_chrom, ".downsampled10k"),
                 "--allow-extra-chr --recode A",
                 "--chr", target_chrom,
                 "--from-bp", start_pos,
                 "--to-bp", end_pos,
                 "--maf 0.001",
                 "--out", paste0(target, "_locus_UKBB.txt"))

use_condaenv(conda_env)

system(command, intern = TRUE)

########################################################################################

# Calculate LD correlation 
data <- fread(paste0(target, '_locus_UKBB.txt.raw'))
data <- data[, -c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
df <- cor(data, method = "pearson", use="pairwise.complete.obs")
#fwrite(df, paste0(target, '_locus_UKBB_ld.txt.gz'))

########################################################################################

# Read in whole GWAS summary statistics
sumstat <- fread(file = gwas_file_path, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE, stringsAsFactors = FALSE)

sumstat <- sumstat[sumstat$hm_chrom == target_chrom & sumstat$hm_pos >= start_pos & sumstat$hm_pos <= end_pos, ]
column_with_nas <- "hm_chrom"
sumstat <- sumstat[complete.cases(sumstat [, column_with_nas]), ]

sumstat$z<-sumstat$beta/sumstat$standard_error

########################################################################################
# Col numbers need to be manually changed depending on their actual position:
########################################################################################
colnames(sumstat)
colnames(sumstat)[c(1, 2, 3, 4, 5, 6, 20, 19, 17)] <- c('ID','rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'p', 'se')
sumstat <- sumstat[c('ID', 'rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'p', 'beta', 'se', 'z')]

##############################################################################################
## Viewing if there is allele concordance between LD matrix and sum stats ##
df1_transpose <- t(data)
df1_transpose <- data.frame(Index = rownames(df1_transpose), df1_transpose, row.names = NULL)
colnames(df1_transpose)[1] <- 'SNP'

# Getting position to merge with sum stats
df1_transpose <- df1_transpose %>%
  mutate(chrom = str_extract(SNP, "^\\d+"),
         position = str_extract(SNP, "(?<=:)\\d+(?=:|$)"))
df1_transpose <- df1_transpose[,c('SNP', 'position')]

#Re-making SNP column into ID column to match sum stats ID (using "_" only)
df1_transpose <- df1_transpose %>%
  separate(SNP, into = c("chromosome", "Col2", "ref_LD",  "alt_LD", "alt_LD2"), sep = "[:,_]", remove = FALSE)
df1_transpose$ID <- paste(df1_transpose$chromosome, df1_transpose$Col2, df1_transpose$ref_LD, df1_transpose$alt_LD, sep = "_")
#concordance_test <- merge(sumstat, df1_transpose, by='ID')

sumstat$ID <- paste(sumstat$chromosome, sumstat$position, sumstat$allele1, sumstat$allele2, sep = "_")
concordance_test <- merge(sumstat, df1_transpose)

concordance_test <- concordance_test[,c('ID','allele1', 'allele2', 'SNP')]
colnames(concordance_test)[4] <- 'UKBB_LD_raw_SNP'

##############################################################################################

## Keeping only overlapping SNPs in both datasets ##
sumstat_filtered <- subset(sumstat, ID %in% concordance_test$ID)

col_indices <- which(colnames(df) %in% concordance_test$UKBB_LD_raw_SNP)
row_indices <- which(rownames(df) %in% concordance_test$UKBB_LD_raw_SNP)
matrix_filtered <- df[row_indices, col_indices]
colnames(matrix_filtered) <- NULL
rownames(matrix_filtered) <- NULL

sumstat_flip <- sumstat_filtered 

concordance_test2 <- concordance_test %>%
  separate(UKBB_LD_raw_SNP, into = c("Col1", "position", "allele1", "allele2", "alt2"), sep = "[:,_]")

sumstat_flip <- sumstat_flip[order(sumstat_flip$ID), ]
concordance_test2 <- concordance_test2[order(concordance_test2$ID), ]

sumstat_flip$z <- ifelse(sumstat_flip$allele1 != concordance_test2$allele1 |
                           sumstat_flip$allele2 != concordance_test2$allele2,
                         -sumstat_flip$z, sumstat_flip$z)

# Check for nans or inf values
any_na <- any(is.na(matrix_filtered))
print(any_na) 
any_inf <- any(sapply(matrix_filtered, function(x) any(is.infinite(x))))
print(any_inf) 

any_na <- any(is.na(sumstat_flip))
print(any_na) 
any_inf <- any(sapply(sumstat_flip, function(x) any(is.infinite(x))))
print(any_inf) 

setwd("/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci")
fwrite(sumstat_flip, paste0(target, '_locus_sumstat_flip_check.txt.gz'), sep = '\t')
fwrite(matrix_filtered, paste0(target, '_locus_ukbb_ld.txt.gz'), sep = '\t', row.names=FALSE, col.names = FALSE)


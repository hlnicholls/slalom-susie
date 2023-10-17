library(susieR)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(reticulate)
library(foreach)
library(doParallel)
set.seed(1)

#Check-list to run
#All path directories are correct
#All column names and SNP IDs match input data
#PLINK directory is setup (line 107)
#check sample size input for susie

#need list of lead snps, lead snp ID expected to be in file name of LD matrix and sumstats
lead_snps <- fread('./susie-r-pipeline/loci_list.txt')
lead_snps$target <- paste(lead_snps$chrom, lead_snps$pos, lead_snps$ref, lead_snps$alt, sep = "_")

# set threshold for DENTIST outlier detection
nlog10p_dentist_s_threshold <- 1e-4
r2_threshold <- 0.6

# Define the number of workers (adjust the value based on available CPU cores)
num_workers <- 16

cl <- makeCluster(num_workers)
registerDoParallel(cl)

foreach(i = 1:nrow(lead_snps), .packages = c("data.table", "dplyr", "stringr", "tidyr", "reticulate", "susieR")) %dopar% {
  target <- as.character(lead_snps[i, "target"])
  target_chrom <- as.character(lead_snps[i, "chrom"])

  ######################################################################################################
  # Read in sumstats and LD for locus (using lead SNP ID for target file name)
  
  setwd("./susie-r-pipeline/loci")
  sumstat <- fread(paste0(target, '_locus_sumstat_flip_check.txt.gz'),
                   sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE,
                   stringsAsFactors = FALSE)
  
  ld <- fread(paste0(target, '_locus_ukbb_ld.txt.gz'),
              sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE,
              stringsAsFactors = FALSE)
  setwd("./susie-r-pipeline/intermediate-results")
  ######################################################################################################
  ## Viewing if there is allele concordance between LD matrix and sum stats ##
  df1_transpose <- t(ld) # ld dataframe needs SNP IDs in columns
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
  
  #sumstat$ID <- paste(sumstat$chromosome, sumstat$position, sumstat$ref, sumstat$alt, sep = "_")
  concordance_test <- merge(sumstat, df1_transpose)
  
  concordance_test <- concordance_test[,c('ID','ref', 'alt', 'SNP')]
  colnames(concordance_test)[4] <- 'UKBB_LD_raw_SNP'
  
  sumstat_filtered <- subset(sumstat, ID %in% concordance_test$ID)
  
  rownames(ld) <- colnames(ld)
  col_indices <- which(colnames(ld) %in% concordance_test$UKBB_LD_raw_SNP)
  row_indices <- which(rownames(ld) %in% concordance_test$UKBB_LD_raw_SNP)
  matrix_filtered <- ld[row_indices, col_indices]
  colnames(matrix_filtered) <- NULL
  rownames(matrix_filtered) <- NULL
  
  ######################################################################################################
  ## Allele flip check ##
  sumstat_flip <- sumstat_filtered 
  
  concordance_test2 <- concordance_test %>%
    separate(UKBB_LD_raw_SNP, into = c("Col1", "position", "ref", "alt", "alt2"), sep = "[:,_]")
  
  sumstat_flip <- sumstat_flip[order(sumstat_flip$ID), ]
  concordance_test2 <- concordance_test2[order(concordance_test2$ID), ]
  
  sumstat_flip$Z <- ifelse(sumstat_flip$ref != concordance_test2$ref |
                             sumstat_flip$alt != concordance_test2$alt,
                           -sumstat_flip$Z, sumstat_flip$Z)
  
  fwrite(sumstat_flip, paste0(target, '_locus_sumstat_flip_check.txt.gz'), sep = '\t')
  fwrite(matrix_filtered, paste0(target, '_locus_ukbb_ld.txt.gz'), sep = '\t', row.names=FALSE, col.names = FALSE)
  
  cat(paste(target, 'sumstat shape:', dim(sumstat_flip), 'ld matrix shape:', dim(sumstat_flip)))
  ######################################################################################################
  
  ######################################################################################################
  ## Calculate DENTIST ##
  ############################
  # 1. Get LD for all variants with lead SNP
  ############################
  # need an anaconda environment with PLINK installed:
  conda_env <- "finemapping_env"
  target_ID <-gsub("_", ":", target)
  command <- paste("/Users/hn9/anaconda3/envs/finemapping_env/bin/plink", #change directory to PLINK here
                   "--bfile", paste0("/Users/hn9/Documents/Analysis/FM-comparison/ukb_v3_downsampled10k/ukb_v3_chr", target_chrom, ".downsampled10k"),
                   "--allow-extra-chr",
                   "--r2",
                   "--ld-snp", target_ID,
                   "--ld-window-kb 1000",
                   "--ld-window 99999",
                   "--ld-window-r2 0",
                   "--out", paste0(target, "_locus_UKBB"),
                   "--silent")
  
  use_condaenv(conda_env)
  
  system(command, intern = TRUE)
  
  lead_ld <- fread(paste0(target, '_locus_UKBB.ld'),
                   sep = " ", header = TRUE, check.names = FALSE, data.table = FALSE,
                   stringsAsFactors = FALSE)
  
  #Filter for only LD with lead SNP
  filtered_df <-  lead_ld %>% filter(SNP_A == target_ID | SNP_B == target_ID)
  
  # Get list of all associated variants from sumstats tsv.gz with : separators
  sumstat <- sumstat_flip
  snps <- select(sumstat, ID)
  snps$ID <- gsub("_(\\d+)_([A-Z])_([A-Z])", ":\\1:\\2:\\3", snps$ID)
  
  merged <- merge(lead_ld, snps, by.x = "SNP_B", by.y = "ID")
  
  colnames(merged)[1] <- 'ID'
  merged$ID <- gsub(":", "_", merged$ID)
  merged <- select(merged, 'ID', 'R2')
  df <- merge(sumstat, merged, by='ID', all.x=TRUE)
  
  ####################################################
  # 2. Calculate 't_dentist_s' and 'nlog10p_dentist_s'
  ####################################################
  
  n_samples_constant <- 500000
  df$r <- (n_samples_constant * sum(df$R2, na.rm = TRUE)) / (n_samples_constant * sum(!is.na(df$R2)))
  
  lead <- df[df$ID == target, ]
  lead_z <- lead$Z
  
  df$t_dentist_s <- (df$Z - df$r*lead_z)^2 / (1 - df$r^2)
  df$nlog10p_dentist_s <- pchisq(df$t_dentist_s, df = 1, lower.tail = FALSE) / (-log10(10))
  n_dentist_s_outlier <- sum(df$R2 > r2_threshold & df$nlog10p_dentist_s > nlog10p_dentist_s_threshold)
  # other threshold option in SLALOM paper:
  #n_dentist_s_outlier <- sum(df$R2 > df$t_dentist_s)
  
  print(target)
  cat('Number of DENTIST outliers detected:', n_dentist_s_outlier, '\n')
  df$dentist_outlier <- ifelse(df$R2 > r2_threshold & df$nlog10p_dentist_s > nlog10p_dentist_s_threshold, 1, 0)
  fwrite(df, paste0(target, '_locus_sumstat_flipcheck_with_dentist.txt.gz'), sep = '\t')
  cat(paste(target, 'sumstat with dentist shape:', dim(df), 'ld matrix shape:', dim(matrix_filtered)))
  ######################################################################################################
  
  ######################################################################################################
  ## SuSiE Fine-mapping ##
  z_vector <- unlist(df$Z)
  ld_matrix <- as.matrix(matrix_filtered)
  fitted <- susie_rss(z = z_vector, R = ld_matrix, n = 500000,  L = 10)
  
  sets <- susie_get_cs(fitted,
                       X = df,
                       coverage = 0.9,
                       min_abs_corr = 0.1)
  
  # Get credible sets
  process_credible_set <- function(cs, z_vector, pip) {
    result <- data.frame(
      position = cs,
      z_vector = z_vector[cs],
      PIP = pip[cs]
    )
    result <- result[order(result$z_vector, decreasing = TRUE),]
    return(result)
  }
  
  all_results <- lapply(fitted$sets$cs, process_credible_set, z_vector = z_vector, pip = fitted$pip)
  combined_results <- bind_rows(all_results, .id = "CredibleSet")
  
  combined_results <- combined_results %>%
    rowwise() %>%
    mutate(
      ID = df$ID[position],
      PValue = df$pval[position],
      Z = df$Z[position],
      dentist_outlier = df$dentist_outlier[position]
    )
  
  setwd("./susie-r-pipeline/results")
  fwrite(
    x = combined_results,
    file = paste0(target, '_susie_locus_ukbb.txt.gz'),
    sep = "\t",
    quote = FALSE,
    na = "NA",
    row.names = FALSE,
    col.names = TRUE,
    compress = "gzip"
  )
  
}

stopCluster(cl)

setwd("~/Documents/slalom-susie/FM-benchmarking")
library(CARMA)
library(data.table)
library(dplyr)
library(susieR)
library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

#######################################################################################################
# 1. No changes to D3
#######################################################################################################
# CARMA D3 Fine-mapping

d3_df <- as.data.frame(coloc_test_data$D3)

sumstats_d3 <- select(d3_df, snp, position, beta, varbeta, N, sdY,type,MAF)
rownames(sumstats_d3) <- as.character(1:nrow(sumstats_d3))
sumstats_d3$Z <- D3$beta/sqrt(D3$varbeta)
ld_d3 <- subset(d3_df, select = -c(snp, position, beta, varbeta, N, sdY,type,MAF))

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats_d3$Z
ld.list[[1]]<-as.matrix(ld_d3)
lambda.list[[1]]<-1

CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                                 outlier.switch=T)

carma_res_original = sumstats_d3 %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    carma_res_original$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}

carma_res_original$Outlier <- ifelse(row.names(carma_res_original) %in% CARMA.results[[1]]$Outlier$Index, 1, 0)
fwrite(carma_res_original, 'carma_res_original.csv')
#######################################################################################################
# SuSiE D3 Fine-mapping

#S3=runsusie(D3)

z_vector <- unlist(sumstats_d3$Z)
ld_matrix <- as.matrix(ld_d3)
SSr=susie_rss(z = z_vector,R=ld_matrix,n=1000)
sets <- susie_get_cs(SSr,
                     coverage = 0.9,
                     min_abs_corr = 0.1)

#susie_plot(fitted, y="PIP")

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

all_results <- lapply(SSr$sets$cs, process_credible_set, z_vector = z_vector, pip = SSr$pip)
combined_results <- bind_rows(all_results, .id = "CredibleSet")
susie_res_original <- merge(sumstats_d3, combined_results, by='position', all.x=TRUE)
fwrite(susie_res_original, 'susie_res_original.csv')
#######################################################################################################

#######################################################################################################
# 2. Changing sign of lead SNP z-score
#######################################################################################################
# CARMA D3 Fine-mapping

d3_df <- as.data.frame(coloc_test_data$D3)

sumstats_d3 <- select(d3_df, snp, position, beta, varbeta, N, sdY,type,MAF)
rownames(sumstats_d3) <- as.character(1:nrow(sumstats_d3))
sumstats_d3$Z <- D3$beta/sqrt(D3$varbeta)
lead_snp_row <- which.max(abs(sumstats_d3$Z))
sumstats_d3$Z[lead_snp_row] <- -sumstats_d3$Z[lead_snp_row]
sumstats_d3$ExpectedOutlier <- 0
sumstats_d3$ExpectedOutlier[lead_snp_row] <- 1
ld_d3 <- subset(d3_df, select = -c(snp, position, beta, varbeta, N, sdY,type,MAF))

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats_d3$Z
ld.list[[1]]<-as.matrix(ld_d3)
lambda.list[[1]]<-1

CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                     outlier.switch=T)

carma_res_lead = sumstats_d3 %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    carma_res_lead$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}

carma_res_lead$Outlier <- ifelse(row.names(carma_res_lead) %in% CARMA.results[[1]]$Outlier$Index, 1, 0)
fwrite(carma_res_lead, 'carma_res_lead.csv')
#######################################################################################################
# SuSiE D3 Fine-mapping

#S3=runsusie(D3)

z_vector <- unlist(sumstats_d3$Z)
ld_matrix <- as.matrix(ld_d3)
SSr=susie_rss(z = z_vector,R=ld_matrix,n=1000)
sets <- susie_get_cs(SSr,
                     coverage = 0.9,
                     min_abs_corr = 0.1)

#susie_plot(fitted, y="PIP")

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

all_results <- lapply(SSr$sets$cs, process_credible_set, z_vector = z_vector, pip = SSr$pip)
combined_results <- bind_rows(all_results, .id = "CredibleSet")
susie_res_lead <- merge(sumstats_d3, combined_results, by='position', all.x=TRUE)
fwrite(susie_res_lead, 'susie_res_lead.csv')
#######################################################################################################

#######################################################################################################
# 3. Changing sign of second most significant SNP z-score
#######################################################################################################
# CARMA D3 Fine-mapping

d3_df <- as.data.frame(coloc_test_data$D3)

sumstats_d3 <- select(d3_df, snp, position, beta, varbeta, N, sdY,type,MAF)
rownames(sumstats_d3) <- as.character(1:nrow(sumstats_d3))
sumstats_d3$Z <- D3$beta/sqrt(D3$varbeta)
second_largest_z <- sort(abs(sumstats_d3$Z), decreasing = TRUE)[2]
second_largest_row <- which(abs(sumstats_d3$Z) == second_largest_z)
sumstats_d3$Z[second_largest_row] <- -sumstats_d3$Z[second_largest_row]
sumstats_d3$ExpectedOutlier <- 0
sumstats_d3$ExpectedOutlier[second_largest_row] <- 1
ld_d3 <- subset(d3_df, select = -c(snp, position, beta, varbeta, N, sdY,type,MAF))

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats_d3$Z
ld.list[[1]]<-as.matrix(ld_d3)
lambda.list[[1]]<-1

CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                     outlier.switch=T)

carma_res_second_snp = sumstats_d3 %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    carma_res_second_snp$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}

carma_res_second_snp$Outlier <- ifelse(row.names(carma_res_second_snp) %in% CARMA.results[[1]]$Outlier$Index, 1, 0)
fwrite(carma_res_second_snp, 'carma_res_second_snp.csv')
#######################################################################################################
# SuSiE D3 Fine-mapping

#S3=runsusie(D3)

z_vector <- unlist(sumstats_d3$Z)
ld_matrix <- as.matrix(ld_d3)
SSr=susie_rss(z = z_vector,R=ld_matrix,n=1000)
sets <- susie_get_cs(SSr,
                     coverage = 0.9,
                     min_abs_corr = 0.1)

#susie_plot(fitted, y="PIP")

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

all_results <- lapply(SSr$sets$cs, process_credible_set, z_vector = z_vector, pip = SSr$pip)
combined_results <- bind_rows(all_results, .id = "CredibleSet")
susie_res_second_snp <- merge(sumstats_d3, combined_results, by='position', all.x=TRUE)
fwrite(susie_res_second_snp, 'susie_res_second_snp.csv')
#######################################################################################################

#######################################################################################################
# 4. Adding a copy of an independent SNP to create perfect linkage
# Duplicating lead variant (SNP with highest z-score)
#######################################################################################################
# CARMA D3 Fine-mapping

d3_df <- as.data.frame(coloc_test_data$D3)

sumstats_d3 <- select(d3_df, snp, position, beta, varbeta, N, sdY,type,MAF)
rownames(sumstats_d3) <- as.character(1:nrow(sumstats_d3))
sumstats_d3$Z <- D3$beta/sqrt(D3$varbeta)
lead_snp_row <- which.max(abs(sumstats_d3$Z))
duplicated_row <- sumstats_d3[lead_snp_row, , drop = FALSE]
sumstats_d3 <- rbind(sumstats_d3, duplicated_row)
sumstats_d3$DuplicatedLead <- 0
duplicated_rows <- which(abs(sumstats_d3$Z) == abs(sumstats_d3$Z[lead_snp_row]))
sumstats_d3$DuplicatedLead[duplicated_rows] <- 1

ld_d3 <- subset(d3_df, select = -c(snp, position, beta, varbeta, N, sdY,type,MAF))
ld_d3_dup <- ld_d3
lead_snp_col_index <- which(colnames(ld_d3_dup) == "LD.s105")
new_col <- ld_d3_dup[, lead_snp_col_index, drop = FALSE]
ld_d3_dup <- cbind(ld_d3_dup, new_col)
colnames(ld_d3_dup)[ncol(ld_d3_dup)] <- "LD.s105_DUP"
new_row <- ld_d3_dup[lead_snp_col_index, , drop = FALSE]
ld_d3_dup <- rbind(ld_d3_dup, new_row)
rownames(ld_d3_dup)[nrow(ld_d3_dup)] <- "LD.s105_DUP"
ld_d3_dup["LD.s105_DUP", "LD.s105_DUP"] = 1

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats_d3$Z
ld.list[[1]]<-as.matrix(ld_d3_dup)
lambda.list[[1]]<-1

CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,
                     outlier.switch=T)

carma_res_copied_snp = sumstats_d3 %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
    carma_res_copied_snp$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  }
}

carma_res_copied_snp$Outlier <- ifelse(row.names(carma_res_copied_snp) %in% CARMA.results[[1]]$Outlier$Index, 1, 0)
fwrite(carma_res_copied_snp, 'carma_res_copied_snp.csv')
#######################################################################################################
# SuSiE D3 Fine-mapping

#S3=runsusie(D3)

z_vector <- unlist(sumstats_d3$Z)
ld_matrix <- as.matrix(ld_d3_dup)
SSr=susie_rss(z = z_vector,R=ld_matrix,n=1000)
sets <- susie_get_cs(SSr,
                     coverage = 0.9,
                     min_abs_corr = 0.1)

#susie_plot(fitted, y="PIP")

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

all_results <- lapply(SSr$sets$cs, process_credible_set, z_vector = z_vector, pip = SSr$pip)
combined_results <- bind_rows(all_results, .id = "CredibleSet")
susie_res_copied_snp <- merge(sumstats_d3, combined_results, by='position', all.x=TRUE)
fwrite(susie_res_copied_snp, 'susie_res_copied_snp.csv')
#######################################################################################################

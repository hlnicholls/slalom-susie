library(CARMA)
library(data.table)
library(dplyr)
##############################################################################################################################
sumstats <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/APOE_LDL_locus_sumstat_with_dentist.txt.gz')
ld <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/APOE_LDL_locus_ukbb_ld.txt.gz')

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats$z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = T)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
APOE_time_T <- difftime(end_time, start_time, units = "secs")

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = F)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
APOE_time_F <- difftime(end_time, start_time, units = "secs")
##############################################################################################################################
##############################################################################################################################
sumstats <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/CLCN6_locus_sumstat_with_dentist.txt.gz')
ld <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/CLCN6_locus_ukbb_ld.txt.gz')

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats$z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = T)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
CLCN6_time_T <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = F)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
CLCN6_time_F <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")
##############################################################################################################################

##############################################################################################################################
sumstats <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/NOD2_locus_sumstat_with_dentist.txt.gz')
ld <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/NOD2_locus_ukbb_ld.txt.gz')

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats$z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = T)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
NOD2_time_T <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = F)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
NOD2_time_F <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")
##############################################################################################################################
##############################################################################################################################
sumstats <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/PCSK9_locus_sumstat_with_dentist.txt.gz')
ld <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/PCSK9_locus_ukbb_ld.txt.gz')

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats$z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = T)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
PCSK9_time_T <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = F)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
PCSK9_time_F <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")
##############################################################################################################################
##############################################################################################################################
sumstats <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/PTK2B_locus_sumstat_with_dentist.txt.gz')
ld <- fread('/Users/hn9/Documents/GitHub/fine-mapping-inf/susieinf/loci/PTK2B_locus_ukbb_ld.txt.gz')

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstats$z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = T)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
PTK2B_time_T <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")

start_time <- Sys.time()
time_taken <- system.time({
  CARMA.results <- CARMA(z.list, ld.list, lambda.list = lambda.list, outlier.switch = F)
})
end_time <- Sys.time()
print(paste("CPU Time taken: ", time_taken['elapsed'], " seconds"))
print(paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds"))
PTK2B_time_F <- paste("Wall-clock Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds")
##############################################################################################################################
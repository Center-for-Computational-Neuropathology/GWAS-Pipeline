#Load required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(qqman))


#replication
data <- fread("/sc/arion/projects/tauomics/Shrishtee/plinkfile_PSP/plink_files_PSP/Replication_Upenn/PSP_LR_withsex_all_chrom.assoc.logistic.ADD_gwas_without_NA")

#Removing Non-autopsy confirmed
data <- fread("/sc/arion/projects/tauomics/Shrishtee/plinkfile_PSP/plink_files_PSP/Association_results_new/PSP_notautopsy_Confirmed_all_chrom.assoc.LR.ADD_gwas_without_NA")

#qc
head(data)
range(data$P)
range(data$STAT)
# data_filtered <- data[abs(data$STAT) <= 5, ]

data <- data[data$P != 1, ]
data <- data[data$P != 0, ]
range(data$P)
data <- mutate(data, P = -log10(P))
range(data$P)

head(data)

colnames(data)[4] <- "Ref_allele"

#FREQUENCY
freq <- fread("/sc/arion/projects/tauomics/Shrishtee/plinkfile_PSP/plink_files_PSP/plink.frq.cc")
head(freq)
common_cols <- intersect(colnames(data), colnames(freq))
data_merge <- merge(data, freq, by = common_cols)

# Print the updated data
head(data_merge)

#Order by chrom and position
data_merge <- data_merge[order(data_merge$CHR, data_merge$BP), ]
head(data_merge)

########## MAKE A DATAFRAME WITH ONLY THE FILTERED 7M SNPs
data_7M_kurt <- fread("/sc/arion/projects/tauomics/Shrishtee/plinkfile_PSP/COJO/Kurt_summ_stats.ma")
head(data_7M_kurt)
filtered_data_7M <- data_merge[data_merge$SNP %in% data_7M_kurt$SNP, ]
range(filtered_data_7M$P)
range(data_merge$P)

####################

#Convert SNPs to variant IDs for locus zoom
rm(freq)
# snp_parts <- strsplit(data_merge$SNP, ":", fixed = TRUE)
# data_merge$SNP <- sapply(snp_parts, function(x) paste0(x[1], ":", x[2], "_", x[3], "/", x[4]))
# rm(snp_parts)
# head(data_merge)
# 
# data_merge <- data_merge %>% select(SNP, P, OR, SE, MAF_A)
# head(data_merge)
# 
# fwrite(data_merge, 
#        file = "/sc/arion/projects/tauomics/Shrishtee/plinkfile_PSP/plink_files_PSP/Association_results_new/PSP_March26_8703_samples_with_sex_16M.tsv", 
#        sep = "\t")


### prep 7M for locus zoom
head(filtered_data_7M)
filtered_data_7M <- filtered_data_7M[order(filtered_data_7M$CHR, filtered_data_7M$BP), ]
head(filtered_data_7M)
range(filtered_data_7M$P)

snp_parts <- strsplit(filtered_data_7M$SNP, ":", fixed = TRUE)
filtered_data_7M$SNP <- sapply(snp_parts, function(x) paste0(x[1], ":", x[2], "_", x[3], "/", x[4]))
rm(snp_parts)
head(filtered_data_7M)

fwrite(filtered_data_7M, 
       file = "/sc/arion/projects/tauomics/Shrishtee/plinkfile_PSP/plink_files_PSP/Association_results_new/PSP_GWAS_8703_samples_with_sex_7M_withoutNonautopsyconfirmed_All_columns.tsv",
       sep = "\t")

filtered_data_7M <- filtered_data_7M %>% select(SNP, P, OR, SE, MAF_A)
head(filtered_data_7M)
range(filtered_data_7M$P)


fwrite(filtered_data_7M, 
       file = "/sc/arion/projects/tauomics/Shrishtee/plinkfile_PSP/plink_files_PSP/Association_results_new/PSP_March26_8703_samples_with_sex_7M_non_autopsy_status.tsv", 
       sep = "\t")



# data$A1 <- sapply(strsplit(data$SNP, ":"), function(x) x[3])
# data$A2 <- sapply(strsplit(data$SNP, ":"), function(x) x[4])
# 
# data <- na.omit(data)
# 
# head(data)
# 
# data <- within(data, {
#   temp <- A1
#   A1[Ref_allele != temp] <- A2[Ref_allele != temp]
#   A2[Ref_allele != temp] <- temp[Ref_allele != temp]
#   rm(temp)  # Remove temporary variable
# })
# 
# data <- data %>% select(-A1, -A2)




# #Select columns to plot a Manhattan
# gwasResults <- data %>% dplyr::select(SNP, CHR, BP, P)
# # colnames(gwasResults) <- c("SNP", "CHR", "BP", "P")
# 
# #Check the frequency of SNPs on each chromosome
# as.data.frame(table(gwasResults$CHR))
# gwasResults$CHR <- as.numeric(gsub("[^0-9]", "", gwasResults$CHR))
# 
# # Plot Manhattan plot
# manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 110), cex = 0.6, cex.axis = 0.9,
#           col = c("blue4", "orange3"), suggestiveline = TRUE, genomewideline = TRUE)
# 
# 
# # Create a Manhattan plot
# ggplot(data, aes(x = CHR, y = NEG_LOG10_P)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = rainbow(length(unique(data$CHR)))) +
#   labs(title = "Manhattan Plot",
#        x = "Genomic Position (BP)",
#        y = "-log10(p-value)",
#        color = "Chromosome") +
#   theme_minimal() +
#   theme(legend.position = "bottom")

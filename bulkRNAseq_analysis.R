## BULK RNASEQ analysis
# use R 4.1.1 for deseq analysis
library(DESeq2)
library(tidyverse)

#Reading data
# gene names in "Genes" column, and each sample as a single column
counts_data <- read.csv("bulkRNAseq_counts.txt", sep = "\t", 
                        header = T, row.names = "Genes")
# calculate # reads for each sample
read_depths <- colSums(counts_data)
#print(read_depths)
plot(read_depths) # varying depths for the samples

sample_names <- colnames(counts_data)

# inconsistent number of samples and columns in count matrix
library(GEOquery)
geo_data <- getGEO("GSE176078")
# sample IDs in geo_data$GSE176078_series_matrix.txt.gz$title
# disease subtype in geo_data$GSE176078_series_matrix.txt.gz$`clinical_subtype:ch1`
# entries 20 and 23 are absent
# creating a dataframe for sample name and disease subtype
sample_conditions <- data.frame(matrix(nrow = 24, ncol = 2))
colnames(sample_conditions) <- c("sample_name", "tumor_subtype")
sample_ids <- geo_data$GSE176078_series_matrix.txt.gz$title[-c(20,23)]
cancer_subtypes <- geo_data$GSE176078_series_matrix.txt.gz$`clinical_subtype:ch1`[-c(20,23)]
sample_conditions$sample_name <- sample_ids
sample_conditions$tumor_subtype <- cancer_subtypes
row.names(sample_conditions) <- sample_conditions$sample_name
sample_conditions <- sample_conditions[-c(1)]

# row names to be in same order as that of columns
sample_cleaned <- data.frame(matrix(nrow=24, ncol=2))
colnames(sample_cleaned) <- c("sample_name", "cancer_subtype")
i = 1
for (colid in colnames(counts_data)) {
  subtype_cancer <- sample_conditions[row.names(sample_conditions)==colid,]
  print(subtype_cancer)
  sample_cleaned[i, "cancer_subtype"] <- subtype_cancer
  sample_cleaned[i, "sample_name"] <- colid
  
  i = i +1
}

row.names(sample_cleaned) <- sample_cleaned$sample_name
sample_cleaned <- sample_cleaned[-c(1)]
# creating DESEQ2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData = sample_cleaned,
                              design = ~cancer_subtype)
# Differential gene expression analysis
dds <- DESeq(dds)
res <- results(dds)

#filtering based on adj p
results_df <- data.frame(res)
res_df <- results_df %>% drop_na(padj)
row.names(res_df[res_df$padj < 0.05, ])

#plotting volcano plot
plot(results_df$log2FoldChange, -log10(results_df$padj),
            xlab = "Log2FC", ylab = "-log10(Adjsuted p value")
abline(h = 5, col = "red")
abline(v = c(-2,2), col = "red")








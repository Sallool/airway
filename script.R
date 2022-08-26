# script to perform differential gene expression analysis using DESeq2 package
# setwd()


BiocManager::install("airway")
browseVignettes("airway")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

data("airway")
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]

sample_info$dex <- gsub("trt", "treated", sample_info$dex )

names(sample_info) <- c("cellLine", "dexamethasone")

write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)


countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# Step 1: preparing count data ................

# read in counts data
countsData <- read.csv("counts_data.csv")
head(countsData)

# read in sample info
colData <- read.csv("sample_info.csv")

# making sure the row names in colData matches to column names in counts_data
all(colnames(countsData) %in% rownames(colData))


# are they in the same order
all(colnames(countsData) == rownames(colData))


# Step 2: construct a DESeqDataSet object ..............

dds <- DESeqDataSetFromMatrix(countData = countsData,
                       colData = colData,
                       design = ~ dexamethasone)
dds

# pre-filtering: removing rows with low gene counts
# keeping row that have at least 10 reads total

keep <- rowSums(counts(dds)) >= 10
keep

dds <- dds[keep, ]
dds

# set the factor level

dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# Step 3: Run DESeq ..................

dds <- DESeq(dds)

res <- results(dds)
res

# Explore results ...............

summary(res)

"to adjust p-value"

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)


# contrasts
resultsNames(dds)

# MA plot
plotMA(res)









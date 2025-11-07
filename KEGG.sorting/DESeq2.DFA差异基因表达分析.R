library(DESeq2)
library(GUniFrac)
metadata <- read.csv("Design_76.csv", header = T, row.names = 1)
otutab1 <- read.csv("81.gene.abundance.csv", header = T, row.names = 1)
otutab2 <- otutab1[,rownames(metadata)]
otutab <- otutab2[rowSums(otutab2)>0, ]
#taxonomy <- read.csv("taxa20_2%.csv", header = T, row.names = 1)
#countData <- as(otutab, "matrix") #convert the countData object to a matrix
# this resemble as.matrrix()
head(countData)
# Build the DESeq2 Object
#construct the data object from matrix of counts and the metadata table
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metadata,
                              design = ~ Compartment)
# filter the data
dds <- dds[rowSums(counts(dds)) > 1,]
dds
dds <- estimateSizeFactors(dds) # Normalize the Count Data
sizeFactors(dds)
# Estimate the Dispersion
dds<- estimateDispersions(dds)

dds <- DESeq(dds)
res <- results(dds)
res <- results(dds, contrast = c("Compartment", "Root", "Bulk"))
resOrdered <- res[order(res$padj), ]
summary(res)
table(res$padj < 0.05)   # 查看padj小于0.05的个数
res[order(res$padj < 0.05),]
res1 <- res[which(res$padj < 0.05),]
sum(res$padj < 0.05, na.rm = TRUE)
write.csv(as.data.frame(res1), file="Rotation_Root_differential.csv")
write.csv(as.data.frame(res), file="Rotation_Root_all.csv")
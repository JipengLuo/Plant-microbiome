BiocManager::install('GUniFrac')
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

### 替换ASV成简单名称
diff <- read.csv("NPKM-D2 to D0.csv", header = T, row.names = 1)
ID <- read.csv("ID.csv", header = T, row.names = 1)
id.new <- ID[rownames(diff), ]
write.csv(id.new, file = "NPKM-D2 to D0_new_table.csv")
********************************************************************

a <- read.csv("CK-D1 to D0.csv", header = T, row.names = 1)
b <- read.csv("CK-D2 to D0.csv", header = T, row.names = 1)
c <- read.csv("NPK-D1 to D0.csv", header = T, row.names = 1)
d <- read.csv("NPK-D2 to D0.csv", header = T, row.names = 1)
e <- read.csv("NPKM-D1 to D0.csv", header = T, row.names = 1)
f <- read.csv("NPKM-D2 to D0.csv", header = T, row.names = 1)
g <- rbind(a, b, c, d, e, f)
write.csv(g, file = "Differential abundant ASVs.csv")
#DESeq()函数同时完成filterring和Dispersion estimation的分析
# dds$group <- relevel(dds$group, "NonSmoker") 该两种形式都可行
dds$group <- factor(dds$group, levels = c("Root", "Soil"))

dds <- DESeq(dds)
#Extract the Results Table
res <- results(dds)
res
# We extract some results of interest
mcols(res, use.names=TRUE)
mcols(dds,use.names=TRUE)[1:4,1:4]
substr(names(mcols(dds)),1,10)
head(assays(dds)[["mu"]])
head(dispersions(dds))
head(mcols(dds)$dispersion)
sizeFactors(dds)
head(coef(dds))
###################################
# Compare Differential Abundance Between Groups Using Contrast
res <- results(dds, contrast = c("group", "Smoker", "NonSmoker") )
res
# Adjust p-Values Using FDR
# how many OTUs that have a p-value below and greater than 0.01
sum(res$pvalue < 0.01, na.rm=TRUE )
table(is.na(res$pvalue))
# checke the No. of p lower than 0.1
table(res$padj < 0.1) # or use: table(res[,"padj"] < 0.1)
res_Sig <- res[which(res$padj < 0.1 ),]
head(res_Sig[order(res_Sig$log2FoldChange),])
tail(res_Sig[order( res_Sig$log2FoldChange ),])
# Export the results
write.table(as.data.frame(res_Sig), file="res_Sig.xls", sep="\t", quote = F,col.names = NA)


# Diagnose and Improve the Testing Results
plotMA(res)
# Diagnostic Plot Using the plotDispEsts()
plotDispEsts(dds, ylim = c(1e-2, 1e3))

# Clustering with Heatmap
rld <- rlog(dds) #trransformation
vst <- varianceStabilizingTransformation(dds)

par(mfrow = c(1, 3))
plot(log2(1 + counts(dds, normalized=TRUE)[,1:2] ), main="Ordinary log2 transformation",
     col="#00000020", pch=20, cex=0.3 )
plot(assay(rld)[,1:2], main="regularized-logarithm transformation
(rlog)", col="#00000020", pch=20, cex=0.3 )
plot(assay(vst)[,1:2], main="Variance stabilizing transformations
(VST).", col="#00000020", pch=20, cex=0.3 )
head(assay(rld))[,1-3]

# Clustering using heatmap
library("gplots")
library("RColorBrewer")
library("genefilter")
library(SummarizedExperiment)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ),10)
heatmap.2(assay(rld)[topVarGenes,], scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))


************using p.adjust() to correct for FDR*******

res[,"padj"] <- p.adjust(res_fdr$pval, method = "BH")



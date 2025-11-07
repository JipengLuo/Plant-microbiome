#载入所有所需R包
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)
egg <- read.table("eggnog.emapper.annotations", sep="\t", header=T, quote = "",stringsAsFactors = FALSE)
# 选取前100行数据
# egg100 <- egg[1:100,]

gene_ids <- egg$query  #提取id列
#有的基因没有注释到会显示为""，需使用逻辑值索引去除未注释到的
eggnog_lines_with_go <- egg$GOs!= ""
eggnog_lines_with_go
#将一个GeneId对应多个GOId的宽数据格式转换位长数据格式
eggnog_annoations_go <- str_split(egg[eggnog_lines_with_go,]$GOs, ",")

gene2go <- data.frame(
           gene = rep(gene_ids[eggnog_lines_with_go], times=sapply(eggnog_annoations_go,length)), 
           term = unlist(eggnog_annoations_go))
gene2go$term[gene2go$term=='-'] <- NA
gene2go <- na.omit(gene2go)
term2gene1 <- gene2go[, c(2, 1)]

library(clusterProfiler)
gene_list <- gene2go$gene[1:80]  #选择geneID和GO名称
term2gene <- gene2go[ ,c(2,1)]

df <- enricher(gene = gene_list,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               TERM2GENE = term2gene)
head(df)
barplot(df)
dotplot(df)

# y轴的标签通常是GO term (就是文字的那个)而不是GO id。clusterProfiler包同样
# 提供了函数对ID和term互相转换。go2term(), go2ont()
df <- as.data.frame(df)
head(df)
dim(df)
df1 <- go2term(df$ID)
dim(df1)
head(df1)
df$term <- df1$Term
df2 <- go2ont(df$ID)
dim(df2)
head(df2)
df$Ont <- df2$Ontology
head(df)
df3 <- df %>%
       dplyr::select(c("term","Ont","pvalue"))
head(df3)

library(ggplot2)
ggplot(df3, aes(x=term,y= -log10(pvalue))) +
  geom_col(aes(fill=Ont)) +
  coord_flip() + labs(x="") +
  theme_bw()



##############KEGG富集分析#####################################################
library(stringr)
library(dplyr)
library(clusterProfiler)
egg <- read.table("eggnog.emapper.annotations", sep="\t", header=T, quote = "")
colnames(egg)
egg[egg==""] <- NA
gene2ko <- egg %>%
           dplyr::select(GID = query, Ko = KEGG_ko) %>% na.omit()
head(gene2ko)
head(gene2go)

eggnog_lines_with_ko <- egg$KEGG_Ko!= ""
eggnog_lines_with_ko
#将一个GeneId对应多个GOId的宽数据格式转换位长数据格式
eggnog_annoations_ko <- str_split(egg$KEGG_ko, ",")

gene2ko <- data.frame(
  gene = rep(gene_ids, times=sapply(eggnog_annoations_ko, length)), 
  term = unlist(eggnog_annoations_ko))

gene2ko[,2] <- gsub("ko:", "", gene2ko[,2])  #将'ko:'替换成''，空的
head(gene2ko)
gene2ko$term[gene2ko$term=="-"] <- NA
gene2ko <- na.omit(gene2ko)
names(gene2ko) <- c("gene", "Ko")

#kegg_info.RData这个文件里有pathway2name这个对象
load(file = "kegg_info.RData")
#ko2pathway <- as.data.frame(ko2pathway)
#write.table(pathway2name,file="pathway2name.txt",sep="\t",row.names=F,quote=F)
pathway2gene <- gene2ko %>% 
  left_join(ko2pathway, by = "Ko") %>%
   dplyr::select(Pathway, gene) %>% na.omit()

str(gene2ko)
str(ko2pathway)
head(ko2pathway)
pathway2name
df <- enricher(gene=gene_list,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             TERM2GENE = pathway2gene,
             TERM2NAME = pathway2name)
dotplot(df)
barplot(df)





#进一步处理，将间接注释补全，将GOid翻译为GOterm与GOontology；为直接注释补充为间接注释
term2gene <- buildGOmap(term2gene1)
#将GoId转换为GoTerm
go2term <- go2term(term2gene$GO)
#将GoId转换为GoOnt
go2ont <- go2ont(term2gene$GO)

#读取gene列表，感兴趣的基因存为gene1.txt文件，存为一列
gene1 <- read.table("gene1.txt", header = FALSE, stringsAsFactors = FALSE)
gene1 <- gene1$V1[1:nrow(gene1)]

#使用enricher函数进行富集分析，这里设置pvalueCutoff = 1, qvalueCutoff = 1以展示所有结果
df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff = 1)
#dotplot 气泡图
#横轴为GeneRatio，代表该GO term下富集到的基因个数占列表基因总数的比例
#纵轴为富集到的GO Terms的描述信息，showCategory指定展示的GO Terms的个数
p1 <- dotplot(df, showCategory = 15 ,title = "barplot for enricher")
p1

#分别重新设置点的颜色，点的大小，Y轴标签的换行
p2 <- p1 + scale_color_continuous(low = "purple", high = "green") + scale_size(range = c(5, 15)) + scale_y_discrete(labels = function(y) str_wrap(y, width = 20))
p2     

#cnetplot 关系网络图
#barplot和dotplot都只显示了最显著的GO terms，cnetplot可显示哪些基因参与了这些terms
#灰色的点代表基因，黄色的点代表富集到的GO terms
#如果一个基因位于一个GO Terms下，则将该基因与GO连线
#黄色节点的大小对应富集到的基因个数，默认showCategory设置top5富集到的GO terms
#设置node_label为："category", "gene", "all"（默认）, "none"来控制标签显示
p3 <- cnetplot(df, node_label = "all", showCategory = 6)
p3

#设置circular = TRUE展示为环形
p4 <- cnetplot(df, circular = TRUE, colorEdge = TRUE, node_label = "category", showCategory = 6)
p4

#读取基因表达量文件，第一列为geneid，第二例为表达量，然后按clusterProfiler包的要求处理
d <- read.table("biaoda_files.txt", header = FALSE, sep = "\t")
geneList <- d[,2]
names(geneList) <- as.character(d[,1])
geneList <- sort(geneList, decreasing = TRUE)

#设置基因标签颜色与表达量相关 foldChange=geneList，并调整基因点的颜色、大小
p5 <- cnetplot(df, node_label = "all", showCategory = 6, foldChange = geneList, circular = TRUE, colorEdge = TRUE) +scale_color_continuous(low = "green", high = "red") + scale_size(range = c(5, 15))
p5

p6 <- emapplot(df2,pie="count", pie_scale=1, layout="kk")
p6


#读入多个基因集，第二个
gene2 <- read.table("0.txt", header = FALSE)
gene2 <- gene2$V1[1:nrow(gene2)]
#第三个
gene3 <- read.table("50.txt", header = FALSE)
gene3 <- gene3$V1[1:nrow(gene3)]
#将多个基因集合并为一个list
gene_cluster = list(gene1 = gene1, gene2 = gene2, gene3 = gene3)

df2 <- compareCluster(gene_cluster, fun = 'enricher',TERM2GENE = term2gene, TERM2NAME = go2term,pvalueCutoff = 1, qvalueCutoff = 1)

p7 <- dotplot(df3,  showCategory = 5) + scale_y_discrete(labels = function(y) str_wrap(y, width = 50)) +  scale_size(range = c(3, 10)) + scale_color_continuous(low = "purple", high = "green")
p6








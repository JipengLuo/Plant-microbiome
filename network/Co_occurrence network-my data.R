library(igraph)
library(dplyr)
library(Hmisc)

## 读入OTU/ASV表格，列为样本，行为物种
design <- read.csv("Design_74_理化和坐标.csv", header = T, row.names = 1)
otu.spe <- read.csv('RPKM.KO.abundance_9353_74samples.csv', header = T,row.names = 1)

## 选择Irrigated and Unirrigated samples
irr.samp <- otu.spe[which(grepl("irrigated", colnames(otu.spe)))] # 选择列名中含有'irrigated'的列名
unirr.samp01 <- otu.spe[which(grepl("drought", colnames(otu.spe)))]  # 选择列名中含有'drought'的列名
dim(unirr.samp01)

# 应该先从unirrigated farms中随机选择28个样品, 然后再合并, 选择至少在一半样品中出现的OTUs
# 随机选择30列样品作图
set.seed(123)  # 设置随机种子以保证结果可复现
unirr.samp <- unirr.samp01[, sample(ncol(unirr.samp01), 28)]  # replace = F, 不放回取样

# 合并数据
otu_df <- cbind(irr.samp, unirr.samp) %>% t()
dim(otu_df)

otu.percent <- otu_df/rowSums(otu_df)  # Transfer into percentage in each sample
# 表示每个OTU占所有样本中总相对丰度的比例，也可以理解为该OTU的总体平均丰度（跨所有样本）
# 方法一
otu.dorm <- otu.percent[, colSums(otu.percent)/sum(otu.percent)>0.0001]  # Remove rare otu with total RA<0.01%
dim(otu.dorm)
# 方法二;  #筛选出在任意样本(至少一个样品)中相对丰度>0.0001(即0.01%)的otu
otu.dorm <- otu.percent[ ,colSums(otu.percent>0.0001) > 0]
dim(otu.dorm)
sum(otu.dorm)/sum(otu.percent)  # Check coverage

otu.present <- otu.dorm
otu.present[otu.dorm !=0] <- 1 #Remove rare OTUs that occurring in at least 4 samples.
otu.abundance <- otu.dorm[rowSums(otu.present) >= 28,]
otuz <- t(otu.abundance)
dim(otuz)
#raresamples <- colnames(otu)[colnames(otu) %in% colnames(otuz)] # 注意行和列的转换
otu_rare01 <- otu.spe[rownames(otuz), colnames(otuz)]
dim(otu_rare01)
otu_rare <- log(otu_rare01 + 1)
write.csv(otu_rare01, file = "KO.Species.count_56samples_0.01%_occurr.in.half.samp.csv")

## 定义一些颜色
col_g <- "#C1C1C1"
cols <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423",
          "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")

# 需要对丰度表格样品名重命名，为irrigated1,irrigated2... drought1, drought2...
trt_id <- c("irrigated", "drought") ## 定义样本的关键词，然后从样本名抓取处理的样本
split_otu <- lapply(apply(sapply(trt_id, function(x){grep(x, colnames(otu_rare))}), 2, 
                          FUN= function(x){otu_rare[,x]}), function(x){x[(which(rowSums(x) !=0)),]})
# function(x){x[(which(rowSums(x) !=0)),]}  # 选择行和非零的行
# function(x){x[-(which(rowSums(x)==0)),]}  # 移除行和为 0 的行，保留剩余的行（即非全零行）

# 1）. 获取分组索引
# indices <- sapply(trt_id, function(x) grep(x, colnames(otu_rare)))
# print(indices)
# 2.提取子矩阵
# sub_matrices <- apply(indices, 2, function(x) otu_rare[, x])
# print(lapply(sub_matrices, dim))  # 检查每个子矩阵的维度
# print(lapply(sub_matrices, function(x) sum(x > 0)))  # 检查非零值数量

# 3) 保留矩阵中 非全零行（即行和不为 0 的行）
# result <- lapply(sub_matrices, function(x) x[(which(rowSums(x) != 0)), ])
# print(lapply(result, dim))  # 检查结果维度,都是0； 换个方法，保留rowSums(x) != 0, 选择行和非零的行

## 有此处主要聚焦展示绘图方法，构建网络时没有对输入数据进行筛选
## 此外，笔者建议在样本量较小的情况下，推荐采用MENA的方法构建网络
g <- lapply(split_otu, function(x){
 #occor <- WGCNA::corAndPvalue(t(x)/colSums(x), method = 'pearson')
  occor <- WGCNA::corAndPvalue(t(x),method = 'spearman')
  mtadj <- multtest::mt.rawp2adjp(unlist(occor$p), proc='BH')
  adpcor <- mtadj$adjp[order(mtadj$index),2]
  occor.p <- matrix(adpcor, dim(t(x)/colSums(x))[2])
  ## R value
  occor.r <- occor$cor
  diag(occor.r) <- 0
  occor.r[occor.p>0.01|abs(occor.r)<0.95] = 0  # 先用 RMT 获得最佳的r值和p值; Species: 0.75; KO: 0.7972
  occor.r[is.na(occor.r)]=0
  g <-  graph_from_adjacency_matrix(occor.r, weighted = TRUE, mode = 'undirected')
  # 删除自相关
  g <- simplify(g)
  # 删除孤立节点
  g <- delete_vertices(g, which(degree(g)==0) )
  return(g)
})

# 也可以用MENA的方法构建保存成如下形式接着跑代码，R中计算基于RMT的cutoff值太慢了（不推荐）
saveRDS(g, file = 'Network_KOs_filtered.(Count).rda')

g <- readRDS("Network_KOs_filtered.(Count).rda")


## 提取第一个网络演示
g1 <- g[[1]]
# plot(g[[1]])

## 设置网络的weight，为计算模块性做准备
E(g1)$correlation <- E(g1)$weight
range(E(g1)$weight)

E(g1)$weight <- abs(E(g1)$weight)  # 相关系数取绝对值
# The proportion of positive and negative edges
sum(E(g1)$weight > 0)
sum(E(g1)$weight < 0)


set.seed(007)
V(g1)$modularity <- membership(cluster_fast_greedy(g1))

V(g1)$label <- V(g1)$name  # node对应的gene name or species name
V(g1)$label <- NA
modu_sort <- V(g1)$modularity %>% table() %>% sort(decreasing = T) #统计每个模块node的数量
top_num <- 15  # 提取前18个modules
modu_name <- names(modu_sort[1:15])  # name of modules
modu_cols <- cols[1:length(modu_name)]  # color for each module
names(modu_cols) <- modu_name
V(g1)$color <- V(g1)$modularity
V(g1)$color[!(V(g1)$color %in% modu_name)] <- col_g
V(g1)$color[(V(g1)$color %in% modu_name)] <- modu_cols[match(V(g1)$color[(V(g1)$color %in% modu_name)],modu_name)]
V(g1)$frame.color <- V(g1)$color

E(g1)$color <- col_g
for ( i in modu_name){
  col_edge <- cols[which(modu_name==i)]
  otu_same_modu <- V(g1)$name[which(V(g1)$modularity==i)]
  E(g1)$color[(data.frame(as_edgelist(g1))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g1))$X2 %in% otu_same_modu)] <- col_edge
}

write_graph(g1, file = "Irrigated farms_KOs_filtered.(Count).graphml", format = "graphml")

# Importing data
#g <- readRDS("network_KOs.rda")

# 计算 layout
sub_net_layout <- layout_with_fr(g1, niter=999, grid = 'nogrid')
## 可视化并输出
par(font.main = 4)
plot(g1, layout=sub_net_layout, edge.color = E(g1)$color, vertex.size=1.8)
title(main = paste0('Nodes=', length(V(g1)$name),', ','Edges=', nrow(data.frame(as_edgelist(g1)))))

# Exporting the figure; 导出图片
pdf(paste0("Irrigated farms_Count_KOs.pdf"), encoding="MacRoman", width=16, height=8)
par(font.main=4)
plot(g1,layout=sub_net_layout, edge.color = E(g1)$color, vertex.size=1.8)
title(main = paste0('Nodes=',length(V(g1)$name),', ','Edges=',nrow(data.frame(as_edgelist(g1)))))
# 添加图例，显示颜色代表的模块
legend("right", legend = modu_name, col = modu_cols, pch = 14, bty = "n", title = "Modules")
dev.off()


# Network topological parameters
length(E(g1))  # number of Edge
length(V(g1))  # number of Vector/Node
mean(igraph::degree(g1))  # average.degree
diameter(g1, directed=F, weights=NA) # network diameter
vertex_connectivity(g1)    # 节点连通性(vertex_connectivity)
edge_connectivity(g1)      # 边连通性 (edge_connectivity), also called 'adhension'
mean_distance(g1, directed=F)   # average path length
#mean_distance(g1)       # average path length
# The density of a graph is the ratio of the number of edges and the number of possible edges
edge_density(g1, loops=FALSE)  # connectance or density of graph
edge_density(simplify(g1), loops=FALSE)   # density, this is also right, but different
transitivity(g1)     # 'average' 计算全局Avg.clustering.coeffi; 'local'计算每个节点的局部聚类系数; 'localundirected'无向图

# The vertex and edge betweenness are defined by the number of geodesics (shortest paths) going through a vertex or an edge
betweenness(g1, directed = F, weights = NULL) # node/vertex betweenness
edge_betweenness(g1, e = E(g1), directed = F, weights = NULL, cutoff = -1)  # edge betweenness

# Centralization is a method for creating a graph level centralization measure from the centrality scores of the vertices
centr_degree(g1)$centralization  # Centralize a graph according to the degrees of vertices
centr_clo(g1, mode="all")$centralization  # Centralize a graph according to the closeness of vertices
centr_eigen(g1, directed=FALSE)$centralization  # Centralize a graph according to the eigenvector centrality of vertices
### Betweenness centrality
centr_betw(g1, directed=F, normalized=T)$centralization # The graph level centrality index
centr_betw(g1, directed=F, normalized=T)$res  # The node-level centrality scores

centr_clo(g1, mode = "all")      # closeness centrality of network
centr_clo(g1, mode = "all")$res  # closeness centrality of each node
centralization.closeness(g1, mode = "all", normalized = T) # norm = T,结果除以总节点数;Centralize a graph according to the closeness of vertices
centr_eigen(g1, directed=FALSE)$centralization

# Modularity
fc = cluster_fast_greedy(g1, weights=NULL)
modularity = modularity(g1, membership(fc))
modularity

# The most centralized graph according to eigenvector centrality
centr_eigen(g1)$centralization

# Calculate centralization from pre-computed scores
deg <- degree(g1)
tmax <- centr_degree_tmax(g1, loops=FALSE)
centralize(deg1, tmax)

# The most centralized graph according to eigenvector centrality
g0 <- graph( c(2,1), n=10, dir=FALSE )
g1 <- make_star(10, mode="undirected")
centr_eigen(g0)$centralization
centr_eigen(g1)$centralization



#***********************************************************
a <- degree(g1)  # the number of edges each node has
write.csv(a, file = "node.degree.csv")

## Closeness measures how many steps is required to access every other vertex from a given vertex
b <- closeness(g1, mode="all", weights=NULL, normalized = T)  # Closeness centrality of vertices
c <- betweenness(g1, directed = TRUE, weights = NULL) # node betweenness

d <- data.frame(degree = a, 
                closeness = b,
                betweenness = c)

write.csv(d, file = "Degree,betweenness,closness_Unirrigated_Taxon.csv")


# 绘制边权重分布图
# 加载 ggplot2 包
library(ggplot2)
## 提取个网络演示
g1 <- g[[1]]
E(g1)$correlation <- E(g1)$weight
range(E(g1)$weight)

rvalue <- as.matrix(E(g1)$weight)
cor_vals <- rvalue[upper.tri(rvalue, diag = FALSE)]
hist(rvalue, breaks = 50, xlim = c(-0.7,1),
     main = "Distribution of Correlations", 
     xlab = "Correlation Coefficient", 
     col = "skyblue", border = "white")


hist(cor_vals, breaks = 50, xlim = c(-1,1),
     main = "Distribution of Correlations in irrigated farms", 
     xlab = "Correlation Coefficient", 
     col = "gray", border = "white", prob = TRUE)  # prob = TRUE 转成密度图

lines(density(cor_vals), col = "blue", lwd = 2)
abline(v = 0, col = "red", lty = 2)


# 提取边权重
# edge_weights <- E(g1)$weight
# 或者，如果已经赋值给了 correlation 属性
edge_weights <- E(g1)$weight
# 将边权重转换为数据框，方便 ggplot 使用
weight_df <- data.frame(weight = edge_weights)

library(dplyr)
weight_df %>% group_by(weight) %>% summarize(n())

# 1). 绘制直方图, hist plot
theme_set(theme_bw())
ggplot(weight_df, aes(x = weight)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(x = "Edge Weight (Correlation)", 
       y = "Frequency", 
       title = "Distribution of Edge Weights in Co-occurrence Network") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13)) + 
  scale_x_continuous(limits = c(-1,1), breaks = seq(-1,1,0.2)) + ylim(0,12000)

ggsave("Edge_weights_hist_unirrigated.pdf", width = , height = )


# 2). 绘制密度图 density plot
ggplot(weight_df, aes(x = weight)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(x = "Edge Weight (Correlation)", 
       y = "Density", 
       title = "Density Plot of Edge Weights in Co-occurrence Network") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  ) + 
  scale_x_continuous(limits = c(-1,1), breaks = seq(-1,1,0.2)) + ylim(0,12000)

ggsave("Edge_weights_density.plot_unirrigated.pdf", width = , height = )



library(igraph)

# 假设实际网络有 n 个节点和 m 条边
n <- vcount(g1)  # node
m <- ecount(g1)  # edge

# 1. 生成与实际网络相同节点数和边数的随机网络
random_network <- sample_gnm(n, m, directed = FALSE)

actual_clustering <- transitivity(g1, type = "global")
random_clustering <- transitivity(random_network, type = "global")

print(paste("实际网络聚类系数:", actual_clustering))
print(paste("随机网络聚类系数:", random_clustering))

# 2. 度分布分析
actual_degree_dist <- degree_distribution(g1)
random_degree_dist <- degree_distribution(random_network)

# 绘制实际网络和随机网络的度分布
# Exporting the figure; 导出图片
pdf(paste0("Irrigated actual network and random network.KO.pdf"), encoding="MacRoman", width=8, height=6)
par(font.main=5)
plot(actual_degree_dist, type="l", col="blue", ylim = c(0,0.1), xlab="Degree", ylab="Frequency", main="Degree Distribution")
lines(random_degree_dist, col = "red")
legend("topright", legend = c("Irrigated Network.KO", "Random Network"), col = c("blue", "red"), lty = 1)
dev.off()

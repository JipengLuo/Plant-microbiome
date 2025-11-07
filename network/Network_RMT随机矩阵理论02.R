# Using WGCNA to calculate correlations between each pair of species
library(WGCNA)
BiocManager::install("GO.db")
otu_rare11 <- read.csv('RPKM.KO.abundance_9353_74samples.csv', header = T,row.names = 1,stringsAsFactors = F)
otu_rare11 <- log(otu_rare22 + 1)
## 选择irrigated and unirrigated samples
irr.samp <- otu_rare11[which(grepl("irrigated", colnames(otu_rare11)))] # 选择列名中含有'irrigated'的列名
unirr.samp01 <- otu_rare11[which(grepl("drought", colnames(otu_rare11)))]  # 选择列名中含有'drought'的列名
dim(unirr.samp01)

# 随机选择30列样品作图
set.seed(123)  # 设置随机种子以保证结果可复现
unirr.samp <- unirr.samp01[, sample(ncol(unirr.samp01), 28)]  # replace = F, 不放回取样
# 合并数据
# otu_rare <- cbind(irr.samp, unirr.samp) 
# dim(otu_rare)

# Irrigated farms
cor.mat <- corAndPvalue(t(irr.samp), method = 'spearman')  # obtain a large list contains several elements
range(cor.mat$cor)

# 直方图观察 correlations 分布；将矩阵转为向量形式（只保留上三角，排除对角线）
rvalue <- cor.mat$cor
cor_vals <- rvalue[upper.tri(rvalue, diag = FALSE)]
hist(cor_vals, breaks = 50, xlim = c(-0.7,1),
     main = "Distribution of Correlations", 
     xlab = "Correlation Coefficient", 
     col = "skyblue", border = "white")


hist(cor_vals, breaks = 50, xlim = c(-1,1),
     main = "Distribution of Correlations in irrigated farms", 
     xlab = "Correlation Coefficient", 
     col = "gray", border = "white", prob = TRUE)  # prob = TRUE 转成密度图

lines(density(cor_vals), col = "blue", lwd = 2)
abline(v = 0, col = "red", lty = 2)


#### or use function cor() to calculate P-value
#cor.mat <- cor(otu.abundance)
### using deconvolution method to infer direct correlation dependencies
### 利用去卷积方法除去间接相互作用
source ("decon.R") 
cor.decon <- decon(abs(cor.mat$cor))
plot(cor.decon, cor.mat$cor, xlim = c(0.65,1), main = "Correlation Plot")

library(vegan)
mantel(data.frame(cor.decon), data.frame(cor.mat$cor))  # p<0.001

## Using RMT identify the appropriate similarity of the thresholds to infer networks
cor.decon <- cor.mat$cor
range(cor.decon)

require(RMThreshold)
rm.get.threshold(cor.decon)
# Choose the threshold value based on the shift of the fitting curves from random to exponent distribution
# Alternatively, we can set range of threshold
rm.get.threshold(cor.decon, interval = c(0.35, 0.95)) # irrigated: 0.77; unirrigated: 0.746 

###adjust P-value
cor.p <- matrix(p.adjust(cor.mat$p), nrow = dim(cor.mat$p)[1])

##cut-off matrix; cor.decon <- decon(abs(cor.mat$cor))
cor.cut <- cor.decon #cor.cut is the correlation dataset
cor.cut[abs(cor.decon)<0.75|cor.p>0.01] <- 0

require(igraph)
cor.net <- graph_from_adjacency_matrix(cor.cut, diag = F, weighted = TRUE,mode = "undirected")
cor.net.clean <- induced_subgraph(cor.net, degree(cor.net)>0) #select the node > 0 edge
# 将igraph weight属性赋值到igraph.weight
igraph.weight = E(cor.net.clean)$weight
sum(igraph.weight > 0)  # number of postive correlation
sum(igraph.weight < 0)  # number of negative correlation
# set edge color，postive correlation 设定为blue, negative correlation设定为red
E.color = igraph.weight
E.color = ifelse(E.color>0, "blue",ifelse(E.color<0, "red","grey"))
E(cor.net.clean)$color = as.character(E.color)

#ly <- layout_in_circle(cor.net.clean)
ly <- layout_on_sphere(cor.net.clean)
ly <- layout.fruchterman.reingold(cor.net.clean)
ly <- layout_with_fr(cor.net.clean)
plot.igraph(cor.net.clean, layout = ly, vertex.size=8, vertex.label=NA,
            edge.curved=F, vertex.color="forestgreen")

*******************按照网络模块绘图*************************************************
### plot modularity; 按照模块绘图
# 2. 按照模块设置节点的颜色;添加节点以及边的颜色
# 选取包含节点数量前18个模块赋予不同的颜色，剩余模块赋予灰色
col_g <- "grey71"
cols <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423",
          "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")

cfg <- cluster_fast_greedy(as.undirected(cor.net.clean))
plot(cfg, as.undirected(cor.net.clean), vertex.label= NA, vertex.size=5)
length(cfg)  ## number of communities
membership(cfg)  ## community membership for each node
modularity(cfg) # how modular the graph partitioning is
###  对于每个community, 可以通过如下方式得到其子图
nodes <- V(cor.net.clean)[cfg$membership == 1]
g <- induced_subgraph(cor.net.clean, nodes)
plot(g, layout = layout_in_circle)


set.seed(007)
V(cor.net.clean)$modularity <- membership(cluster_fast_greedy(cor.net.clean))

V(cor.net.clean)$label <- V(cor.net.clean)$name
V(cor.net.clean)$label <- NA
modu_sort <- V(cor.net.clean)$modularity %>% table() %>% sort(decreasing = T)
top_num <- 15
modu_name <- names(modu_sort[1:15])
modu_cols <- cols[1:length(modu_name)]
names(modu_cols) <- modu_name
V(cor.net.clean)$color <- V(cor.net.clean)$modularity
V(cor.net.clean)$color[!(V(cor.net.clean)$color %in% modu_name)] <- col_g
V(cor.net.clean)$color[(V(cor.net.clean)$color %in% modu_name)] <- modu_cols[match(V(cor.net.clean)$color[(V(cor.net.clean)$color %in% modu_name)],modu_name)]
V(cor.net.clean)$frame.color <- V(cor.net.clean)$color

cfg1 <- cluster_fast_greedy(as.undirected(cor.net.clean))
plot(cfg1, as.undirected(cor.net.clean), vertex.label= NA, vertex.size=5)

# We can also plot the communities without relying on their built-in plot:
V(cor.net.clean)$community <- cfg$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(cor.net.clean, vertex.color=colrs[V(cor.net.clean)$community],vertex.label= NA, vertex.size=5)


igraph.col = taxa[V(cor.net.clean)$name,]
levels(igraph.col$phylum) = c("forestgreen","red","skyblue","yellow","black","pink","gray","cyan")
V(cor.net.clean)$color = as.character(igraph.col$phylum)
#igraph.size =  otu.dorm[,V(cor.net.clean)$name]
#igraph.size1 =log((V(cor.net.clean)$abundance)*100)
#V(cor.net.clean)$abundance = igraph.size
#plot.igraph(cor.net.clean, layout = ly, vertex.size=8,vertex.label=NA, edge.curved=F, vertex.color=levels(igraph.col$phylum))

taxa <- read.csv("Generalist_specialist_taxa-终稿.csv",header = T,row.names =1, stringsAsFactors = FALSE)
## In igraph, Vertex(V) mean nodes.
V(cor.net.clean)$Domain <- taxa[V(cor.net.clean)$name, 1]
V(cor.net.clean)$Phylum <- c(taxa[V(cor.net.clean)$name,2])
V(cor.net.clean)$Class <- c(taxa[V(cor.net.clean)$name,3])
V(cor.net.clean)$Order <- c(taxa[V(cor.net.clean)$name,4])
V(cor.net.clean)$Family <- c(taxa[V(cor.net.clean)$name,5]) 
V(cor.net.clean)$Genus <- c(taxa[V(cor.net.clean)$name,6])
V(cor.net.clean)$sign <- c(taxa[V(cor.net.clean)$name,7])
V(cor.net.clean)$abundance <- colSums(otu.dorm)[V(cor.net.clean)$name]/10 
############################################################
### Only remain the positive edges
positive.cut <- cor.cut[upper.tri(cor.cut)]
positive <- cor.mat$cor[upper.tri(cor.mat$cor)][positive.cut!=0]
positive[positive>0] <- "Positive"
E(cor.net.clean)$positve <- positive

write.graph(cor.net.clean,file = "Bulk soil_samples.graphml",format = "graphml")

#############explore the positive and negative correlations


###########  Calculating the keystone taxa
library(ggplot2)
a <- V(cor.net.clean)$Family
b <- degree(cor.net.clean)  # mean(igraph::degree(cor.net.clean))
c <- closeness(cor.net.clean, mode = "all", normalized = T)
d <- centralization.betweenness(cor.net.clean)$res  # The node-level centrality scores

taxadg1 <- data.frame(taxa=a, Average_degree=b)
theme_set(theme_bw())
p1 <- ggplot(taxadg1, aes(x=taxa, y=Average_degree))+
       geom_bar(stat="identity", fill = "grey31") + 
       theme(axis.title = element_text(size=16),
             axis.text.x = element_text(size=16, hjust=1, vjust=1,colour = "black"),
             axis.text.y = element_text(size=16),
             panel.grid = element_blank()) +  
             coord_flip() + 
  scale_y_continuous(breaks = seq(0,1000,200)) +
### scale_y_continuous(expand = c(0,0))//这个可以去掉与X轴间隙
       #scale_y_continuous(expand = c(0, 0)) +###这个可以去掉与Y轴间隙
       geom_hline(yintercept=100, linetype=3, color="blue")
p1 

##########   Closeness 
taxadg2 <- data.frame(taxa=a, Closeness_centrality=c)
theme_set(theme_bw())
p2 <- ggplot(taxadg2,aes(x=taxa, y=Closeness_centrality))+
  geom_bar(stat="identity", fill = "grey31") + 
  theme(axis.title = element_text(size=16,colour = "black"),
        axis.text.x = element_text(size=16, hjust=1, vjust=1),
        axis.text.y = element_text(size=16),
        panel.grid = element_blank()) +
        coord_flip() #+ scale_y_continuous(expand = c(0, 0))
p2

##########   betweenness 
taxadg3 <- data.frame(taxa=a, Betweenness_centrality=d)
theme_set(theme_bw())
p3 <- ggplot(taxadg3, aes(x=taxa, y=Betweenness_centrality))+
  geom_bar(stat="identity", fill = "grey31") + 
  theme(axis.title = element_text(size=16,colour = "black"),
        axis.text.x = element_text(size=16, hjust=1, vjust=1),
        axis.text.y = element_text(size=16),
        panel.grid = element_blank()) +
        coord_flip() +
        scale_y_continuous(breaks = seq(0,55000,10000)) + 
        #scale_y_continuous(expand = c(0, 0)) +
        geom_hline(yintercept=5000, linetype=3, color="blue")
p3

library(ggpubr)
p13 <- ggarrange(p1, p2, p3, labels = c("A","B","C"), nrow = 1, align = "h") #common.legend= TRUE, legend= "bottom"
p13

ggsave("Keystone_taxa_Family_45_samples.pdf", width = 21, height = 11)



*********************计算网络拓扑学参数******************************************
#### Calculating topology parameters
length(E(cor.net.clean))  #edges numbers
length(V(cor.net.clean))  #nodes numbers
mean(igraph::degree(cor.net.clean))   # average.degree
a <- degree(cor.net.clean)
write.csv(a, file = "node.degree.csv")
diameter(cor.net.clean, directed=F, weights=NA) # or delete weights
mean_distance(cor.net.clean, directed=F)   # average path length
#average.path.length = average.path.length(cor.net.clean) 
transitivity(cor.net.clean)  # clustering.coefficient


### Betweenness centrality;一下四种结果相同
centr_betw(cor.net.clean, directed=F, normalized=T)
### The node-level centrality scores
centr_betw(cor.net.clean, directed=F, normalized=T)$res
### The graph level centrality index
centr_betw(cor.net.clean, directed=F, normalized=T)$centralization
centr_degree(cor.net.clean)$res

centralization.closeness(cor.net.clean, mode = "all", normalized = T) # norm = T,结果除以总节点数
centr_clo(cor.net.clean, mode = "all")$res

### Compute eigenvector centrality scores
eigen_centrality(cor.net.clean,directed = FALSE,scale = TRUE,weights = NULL)$vector

centr_eigen(cor.net.clean, directed = FALSE)$centralization
### Assortativity
assortativity_degree(cor.net.clean, directed=F)

### Modularity;   read.graph(cor.net.clean,file = "test.graphml",format = "graphml")
fc = cluster_fast_greedy(cor.net.clean, weights =NULL)
modularity = modularity(cor.net.clean, membership(fc))
modularity
## Closeness measures how many steps is required to access every other vertex from a given vertex
b <- closeness(cor.net.clean, mode="all", weights=NA)  # Closeness centrality of vertices
as.matrix(b)
c <- cbind(as.matrix(a), as.matrix(b))
write.csv(c, file = "Root44444.csv")
centr_clo(cor.net.clean, mode="all", normalized=T)  # Centralize a graph according to the closeness of vertices

#平均加权度（average weighted degree），仅适用于含权网络
#average_weight_degree <- mean(strength(igraph))
### 幂律分布
d <- degree(cor.net.clean, mode="in")
fit1 <- fit_power_law(d+1, 10)
fit2 <- fit_power_law(d+1, 10, implementation="R.mle")

connectance = edge_density(cor.net.clean,loops=FALSE)
centralization.betweenness = centralization.betweenness(cor.net.clean)$centralization
centralization.degree = centralization.degree(cor.net.clean)$centralization
connect.neighborhood(cor.net.clean, mode = "all")
edge.connectivity = edge_connectivity(cor.net.clean)
hubs <- hub_score(cor.net.clean, weights=NA)$vector

### Betweenness
betweenness(cor.net.clean, directed=T, weights=NA)
edge_betweenness(cor.net.clean, directed=T, weights=NA)
centr_betw(cor.net.clean, directed=T, normalized=T)
### Modularity;   read.graph(cor.net.clean,file = "test.graphml",format = "graphml")
fc = cluster_fast_greedy(cor.net.clean, weights =NULL)
modularity = modularity(cor.net.clean, membership(fc))
modularity

fc = cluster_edge_betweenness(cor.net,weights =NULL)
modularity = modularity(cor.net,membership(fc))

fc = cluster_walktrap(cor.net.clean,weights =NULL)
modularity = modularity(cor.net.clean,membership(fc))

ly <- layout_with_kk(cor.net.clean)
plot.igraph(cor.net.clean,
            layout = ly, 
            vertex.size = log(100000*V(cor.net.clean)$abundance), 
            vertex.label.font = 3,
            vertex.label.cex = .8,
            vertex.color = factor(c(taxa[V(cor.net.clean)$name,1])),
            vertex.label = c(taxa[V(cor.net.clean)$name,4]))


#######采用Erdos –Rényi 构建random network
g <- erdos.renyi.game(1067, 8235, type = c("gnm"),
                      directed = FALSE,
                      loops = FALSE)
degree_distribution(g)
transitivity(g) # clustering.coefficient
average.path.length(g) 

fc = cluster_fast_greedy(g, weights =NULL)
modularity = modularity(g, membership(fc))
modularity

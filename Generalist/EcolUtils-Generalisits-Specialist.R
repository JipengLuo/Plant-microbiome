# 安装EcoLUtils包
# devtools::install_github("GuillemSalazar/EcolUtils")
library(vegan)
library(spaa)
library(EcolUtils)
# Specialist/Generalist classification of OTUs based on niche width and permutation algorithms
# ASVs in row, samples in column
comm.tab1 <- read.csv("Species.count_74.csv", row.names=1, header=TRUE, check.names = FALSE)  #, comment.char="@"
comm.tab1[1:5,1:5]
colnames(comm.tab1)
design <- read.csv("Design_74_理化和坐标.csv", row.names=1, header=TRUE)
comm.tab <- comm.tab1[which(grepl("irrigated", colnames(comm.tab1)))]  #Select unirrigated columns

# Or  方法二
library(dplyr)
comm.tab <- comm.tab1 %>% select(matches("^irrigated")) # ^通配符
dim(comm.tab)

# comm.tab <- comm.tab1[,rownames(design)]
# write.csv(comm.tab, file = "RPKM.KOs.count_74.csv")
comm.tab <- comm.tab[rowSums(comm.tab)>0,]  # 保留在所有样品中序列数之和>100的KOs
comm.tab <- t(comm.tab)
#comm.tab <- comm.tab[ ,which(colSums(comm.tab)>50)]
dim(comm.tab)

# 指定随机置换（permutation）的次数，用于生成空模型（null model）以比较实际数据与随机数据的生态位宽度
res <- spec.gen(comm.tab, niche.width.method = "levins", perm.method = "quasiswap", n= 1000, probs = c(0.025, 0.975))
table(res$sign)
# 对于非零的最大值，输出的数组元素通常会是 1（因为归一化后值在 [0, 1] 范围内，向上取整后都变成 1）
comm.tab.bin <- ceiling(comm.tab/max(comm.tab))

# res1 <- read.csv("Root_generalists.csv", row.names=1,header=TRUE)
plot(colSums(comm.tab), colSums(comm.tab.bin)/dim(comm.tab.bin)[1],
     col = res$sign, pch = 19, log="x", xlab = "Abundance", ylab = "Occurrence")
legend("bottomright", levels(res$sign), col= 1:3, pch= 19, inset= 0.01, cex= 1)
# width = 8 cm, height = 7.5 cm

write.csv(res, file = "Irrigated_Generalists_Species_count.csv")



#*************************************************
# 导入数据，作图
res <- read.csv("Irrigated_Generalists_KO.count.csv", header = T, row.names = 1)

# 自定义各组的点颜色
unique(res$sign)  # 查看 res$sign 的唯一值
# 定义颜色向量，长度与因子水平相同
colors <- c("#001f5f", "#007f2c","#db133b")[as.numeric(res$sign)]
# 绘制散点图
plot(log2(colSums(comm.tab)), colSums(comm.tab.bin)/dim(comm.tab.bin)[1],
     col=colors, pch=19, cex=1, cex.axis=1.3, cex.lab=1.3, 
     log="x", xlab="Log2 (Relative abundance)", ylab="Occurrence")
# 添加图例; topleft,
legend("bottomright", legend=levels(res$sign), col=c("#001f5f", "#007f2c","#db133b"), pch=19, cex=1)

# 或者用 ifelse 函数来定义颜色向量
# 创建一个与 res$sign 对应的颜色向量。例如：
# 根据 res$sign 的值分配颜色
colors <- ifelse(res$sign == "GENERALIST", "#001f5f",  # 正值用蓝色
                 ifelse(res$sign == "SPECIALIST", "#db133b",  # 负值用红色
                        "#007f2c"))  # 零值用灰色


root <- read.csv("Root_generalists.csv", row.names=1,header=TRUE)
new.root <- comm.tab1[rownames(root), ]
rhizo <- read.csv("Rhizosphere_generalists.csv", row.names=1,header=TRUE)
new.rhizo <- comm.tab1[rownames(rhizo), ]
new <- rbind(new.root, new.rhizo)

write.csv(new, file = "All_generalists_ZOTUtable.csv")




## Rarefaction of a community matrix with permutations
library(vegan)
data(varespec)
rrarefy.perm(varespec*100)

# Pairwise comparisons for Permutational Multivariate Analysis of Variance Using Distance Matrices
library(vegan)
data(dune)
data(dune.env)
adonis.pair(vegdist(dune),dune.env$Management)




### Classification of OTU's seasonality
comm.tab<-read.table("OTUtable_Salazar_etal_2015_Molecol.txt",sep="\t",row.names=1,header=TRUE,comment.char="@")
comm.tab<-t(comm.tab[,1:60])
comm.tab<-comm.tab[,which(colSums(comm.tab)>0)]
res<-seasonality.test(comm.tab,n=10)

# Computation of the abundance-weighted mean of an environmental variable for all OTUs
# in a community and statistical comparison to randomized communities
niche.val(comm.tab, dune.env, n = 1000, probs = c(0.025, 0.975))

# Computation of the range of an environmental variable for all OTUs 
# in a community and statistical comparison to randomized communities.
niche.range(comm.tab, env.var, n = 1000, probs = c(0.025, 0.975))


### Split moving-window analysis based on multivariate community data and permutations.
library(vegan)
data("varespec")
data("varechem")
tmp<-smwda(vegdist(varespec),varechem$N)
plot(tmp$windows$env.var.mean,tmp$windows$stat.real.zscore,col=tmp$windows$sign,type="b",pch=19)



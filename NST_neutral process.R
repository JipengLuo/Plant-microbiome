#install.packages("NST")
library(NST)

data(tda)
comm=tda$comm
bray=beta.g(comm,dist.method="bray")
bray.3col=dist.3col(bray)
group=tda$group

tnst=tNST(comm=comm, group=group, rand=20,
          output.rand=TRUE, nworker=1)
#检验各组ST、NST的分布情况及各组ST、NST差异的显著性。
nst.bt=nst.boot(nst.result=tnst, group=NULL, rand=99,
                trace=TRUE, two.tail=FALSE, out.detail=FALSE,
                between.group=FALSE, nworker=1)
#ST和NST组间进行Permutational multivariate ANOVA
nst.pova=nst.panova(nst.result=tnst, rand=99)

#可视化
group_box(tnst$index.pair.grp, group = tnst$index.pair.grp$group, 
          mode = 3,p_value2 = T)



#**************************** My data ******************************
# 1. NST --> community assembly processes
comm1 <- read.csv("Gene.count.81_9818-终稿.csv", header = T, row.names = 1)
design <- read.csv("Design_74_理化和坐标.csv", header = T, row.names = 1)
design$rowname <- rownames(design)

comm2 <- t(comm1[,rownames(design)])
bray = beta.g(comm2, dist.method="bray")
bray.3col = dist.3col(bray)
group = as.data.frame(design$Irrigation.adj01)
rownames(group) <- rownames(design)

tnst = tNST(comm=comm2, group=group, rand=30, output.rand=TRUE, nworker=1)

write.csv(tnst$index.pair.grp, file = "NST_Gene.count_9818.csv")
tnst$index.pair


#检验各组ST、NST的分布情况及各组ST、NST差异的显著性。
nst.bt=nst.boot(nst.result=tnst, group=NULL, rand=99,
                trace=TRUE, two.tail=FALSE, out.detail=FALSE,
                between.group=FALSE, nworker=1)

#ST和NST组间进行Permutational multivariate ANOVA
nst.pova=nst.panova(nst.result=tnst, rand=99)

#可视化
group_box(tnst$index.pair.grp, group = tnst$index.pair.grp$group, 
          mode = 3,p_value2 = T)





# 2. Sloan Neutral Community Model

library(devtools)
library(remotes)
#install_github("Russel88/MicEco")
library(MicEco)

otutab <- read.csv("Abundance.rpkm.9818.csv", header = T, row.names = 1)
design <- read.csv("Design_74_理化和坐标.csv", header = T, row.names = 1)
design.irr <- subset(design, Irrigation.adj01 == "irrigated")
design.unirr <- subset(design, Irrigation.adj01 == "non-irrigated")

otu.irr <- otutab[,rownames(design.irr)]
otu.unirr <- otutab[,rownames(design.unirr)]

neutral.fit(otu.irr)    #使用的是最大似然估计拟合模型，R2计算方法也不同
neutral.fit(otu.unirr)  #使用的是最大似然估计拟合模型，R2计算方法也不同

OTU.irr <- otu.irr[rowSums(otu.irr)>50,]
OTU.unirr <- otu.unirr[rowSums(otu.unirr)>50,]

#install.packages("pctax")
library(pctax)
data(OTU.irr)
ncm_res <- ncm(OTU.irr)
plot(ncm_res)   # width = 9.5 cm, height = 7 cm


data(OTU.unirr)
ncm_res <- ncm(OTU.unirr)
plot(ncm_res)   # width = 9.5 cm, height = 7 cm

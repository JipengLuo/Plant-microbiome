# 画sloan图
require(phyloseq)
#remotes::install_github("DanielSprockett/reltools")
require(reltools)
require(ggplot2)
library(minpack.lm)
library(Hmisc)
library(stats4)
library(vegan)

spp1 <- read.csv("Gene.count.81_9818-终稿.csv", header = T, row.names = 1)
#spp <- decostand(spp1, "normalize")  # ,"standardize"
taxa <- read.csv("gene7ko_gene_ID对照_97707.csv", row.names=1, header=T)
spp <- spp1[rownames(spp1) %in% rownames(taxa),]

design <- read.csv("Design_74_理化和坐标.csv", header = T, row.names = 1)
design.irr <- subset(design, Irrigation.adj == "irrigated")
design.unirr <- subset(design, Irrigation.adj == "non-irrigated")
spp.irr1 <- spp1[,rownames(design.irr)]
spp.unirr1 <- spp1[,rownames(design.unirr)]

spp.irr <- spp.irr1[rowSums(spp.irr1) > 20,]
spp.unirr <- spp.unirr1[rowSums(spp.unirr1) > 20,]

set.seed(123)  # 设置随机种子以保证结果可复现
unirr.random <- spp.unirr1[, sample(ncol(spp.unirr1), 28)]  # replace = F, 不放回取样

# 共有的KOs
# common.row <- intersect(rownames(spp.irr1), rownames(spp.unirr1))
# spp2 <- spp1[common.row,]
# write.csv(spp2, file = "unigene.count.63936.csv")


# 计算中性模型中获得的各个参数
spp.out <- fit_sncm(t(spp.unirr), pool=NULL, taxon=NULL)   ### [ ,1:129]
otu <- spp.out$predictions
fit <- spp.out$fitstats
m <- fit[1]   # m=0.097
r2 <- fit[4]  # r2=0.252
# write.csv(spp.out$predictions, file = "predictions.csv")
# otu <- read.csv(file = "predictions.csv", row.names=1,header=T)

library(MicEco)
library(ggplot2)
### otu, an OTU-table with taxa as columns and samples as rows
res = neutral.fit(t(spp.unirr)) 

# 判断点的位置，分成三组，红-绿-黄
require(plyr)
out2 = mutate(otu, group=NA)
out2$group[otu[,2]<otu[,4]]="lower"    # 低于下界
out2$group[otu[,2]>otu[,5]]="upper"    # 高于上界
out2$group[(otu[,2]>=otu[,4])&(otu[,2]<=otu[,5])]="middle"  # 中间

out2$group <- factor(out2$group, levels = c("lower", "middle", "upper"))
mycols<-c("orange", "gray60", "#008B8B")

theme_set(theme_bw())
p1 = ggplot() +
     geom_point(data = out2, aes(x=log(p), y=freq, color= group),size = 1.7)+
        scale_colour_manual(values = mycols) +
        annotate("text",x= -9, y=0.20,label=paste("m = ", round(m,3), sep=''), size=5)+
        annotate("text",x= -9, y=0.25,label=paste("R2 = ", round(r2,3), sep=''), size=5)
p1

p2 <- p1 +
  geom_line(data = otu, aes(x=log(p),y=freq.pred),size = 1,linetype = 1, color ="blue") +
  geom_line(data = otu, aes(x=log(p),y=pred.lwr),size = 0.8,linetype = 2, color ="blue") +
  geom_line(data = otu, aes(x=log(p),y=pred.upr),size = 0.8,linetype = 2, color ="blue") +
  xlab("log10(mean relative abundance)") + ylab("Occurrence frequency in irrigated")
p2

plot_theme = theme(panel.background=element_blank(),
                   #panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=15),
                   legend.position="none",    # right
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=15),
                   text=element_text(family="sans", size=15))

p3 = p2 + plot_theme
p3

ggsave("Unirrigated farms (46 samples).pdf", width = 6.5, height = 4.9)

############# 计算三组OTU比例的饼图
dim(out2)
low = nrow(out2[out2[,7]== "lower",])
med = nrow(out2[out2[,7]== "middle",])
high = nrow(out2[out2[,7]== "upper",])

type <- c('med','high','low')
nums <- c(med, high, low)
df <- data.frame(type = type, nums = nums)
label_value <- paste('', round(df$nums/sum(df$nums) * 100, 1), '%', sep = '')
label_value
label <- paste(df$type, label_value, sep = ' ')
label
p6 <- ggplot(data = df, aes(x = 1, y = nums, fill = type)) + 
      geom_bar(stat = 'identity', position = 'stack', width = 0.5)+
      scale_fill_manual(name='', labels=c(label[2], label[3], label[1]),
                        values=c("DarkCyan","orange","gray60")) +
    coord_polar(theta = 'y') +
    labs(x = '', y = '', title = '') +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(legend.position = "right")
 
p6

p7 = p6 + theme(panel.background=element_blank(),
               panel.grid=element_blank(),
               legend.background=element_blank(),
               legend.key=element_blank(),
               legend.text= element_text(size=18),
               plot.title = element_text(hjust = 0.5)) + 
               labs(title="Unirrigated farms (28 samples)")    #Generalist
p7 
 
ggsave("Unirrigated farms (46 samples)_pie.chart.pdf", width = 8, height = 5)


#*************** Unirrigated farms ***************************************************
# 计算中性模型中获得的各个参数
spp.out <- fit_sncm(t(spp.irr), pool=NULL, taxon=NULL)   ### [ ,1:129]
otu <- spp.out$predictions
fit <- spp.out$fitstats
m <- fit[1]     # m=0.108
r2 <- fit[4]    # r2=0.109
#write.csv(spp.out$predictions, file = "predictions.csv")
#otu <- read.csv(file = "predictions.csv", row.names=1,header=T)

library(MicEco)
library(ggplot2)
### otu, an OTU-table with taxa as columns and samples as rows
res = neutral.fit(t(spp.irr)) 

#判断点的位置，分成三组，红-绿-黄
require(plyr)
out2 = mutate(otu, group=NA)
out2$group[otu[,2]<otu[,4]]="lower"  ##低于下界
out2$group[otu[,2]>otu[,5]]="upper"  ##高于上界
out2$group[(otu[,2]>=otu[,4])&(otu[,2]<=otu[,5])]="middle"  ##中间

out2$group <- factor(out2$group, levels = c("lower", "middle", "upper"))
mycols<-c("orange", "gray60", "#008B8B")

theme_set(theme_bw())
p4 = ggplot() +
  geom_point(data = out2, aes(x=log(p), y=freq, color= group), size = 1.7)+
  scale_colour_manual(values = mycols) +
  annotate("text",x= -8, y=0.15,label=paste("R2 = ", round(r2,3), sep=''), size=5) +
  annotate("text",x= -8, y=0.20,label=paste("m = ", round(m,3), sep=''), size=5)
 
p4

p5 <- p4 +
  geom_line(data = otu, aes(x=log(p),y=freq.pred),size = 1,linetype = 1, color ="blue") +
  geom_line(data = otu, aes(x=log(p),y=pred.lwr),size = 0.8,linetype = 2, color ="blue") +
  geom_line(data = otu, aes(x=log(p),y=pred.upr),size = 0.8,linetype = 2, color ="blue") +
  xlab("log10(mean relative abundance)") + ylab("Occurrence frequency in unirrigated")
p5

plot_theme = theme(panel.background=element_blank(),
                   #panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=15),
                   legend.position="none",    # right
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=15),
                   text=element_text(family="sans", size=15))

p6 = p5 + plot_theme
p6

ggsave("Irrigated farms (28 samples).pdf", width = 6.5, height = 4.9)

############# 计算三组OTU比例的饼图
dim(out2)
low = nrow(out2[out2[,7]== "lower",])
med = nrow(out2[out2[,7]== "middle",])
high = nrow(out2[out2[,7]== "upper",])

type <- c('med','high','low')
nums <- c(med, high, low)
df <- data.frame(type = type, nums = nums)
label_value <- paste('', round(df$nums/sum(df$nums) * 100, 1), '%', sep = '')
label_value
label <- paste(df$type, label_value, sep = ' ')
label

p7 <- ggplot(data = df, aes(x = 1, y = nums, fill = type)) + 
  geom_bar(stat = 'identity', position = 'stack', width = 0.5)+
  scale_fill_manual(name='', labels=c(label[2], label[3], label[1]),
                    values=c("DarkCyan","orange","gray60")) +
  coord_polar(theta = 'y') +
  labs(x = '', y = '', title = '') +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(legend.position = "right")

p7

p8 = p7 + theme(panel.background=element_blank(),
                panel.grid=element_blank(),
                legend.background=element_blank(),
                legend.key=element_blank(),
                legend.text= element_text(size=18),
                plot.title = element_text(hjust = 0.5)) + 
  labs(title="Irrigated farms")  

p8 

ggsave("Irrigated farms (28 samples)_pie.chart.pdf", width = 8, height = 5)


#  合并图表
library(cowplot)
p11 <- plot_grid(p3, p6, nrow = 1, align = "h")
p11

# 12cm, 4.5 cm
ggsave("Neutral_model_KOs.pdf", width = 12, height = 5)




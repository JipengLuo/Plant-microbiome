## 计算 niche width
library(spaa)

### Irrigated fields
genera <- read.csv("Irrigated_Generalists_KO.count.csv", header = T, row.names = 1)
table(genera$sign)   # GENERALIST    NON SIGNIFICANT    SPECIALIST 
                     #    1476            3715             4104 
zotu <- read.csv("RPKM.KOs.count_74.csv", header = T, row.names = 1)
relative <- zotu/colSums(zotu) # Relative abundance
gener <- zotu[rownames(genera), ]
gener.rel <- relative[rownames(genera), ]

## 列为物种名，行为样品名
otu <- t(gener)
width <- niche.width(otu, method = c("levins")) # the average of niche width from all taxa occurring in all samples

### 计算log转化的mean relative abundance
mean <- log10(apply(gener.rel, 1, mean))
head(mean)

aa <- rbind(width, mean)
rownames(aa) <- c("niche", "abundance")
bb <- data.frame(t(aa))
cc <- cbind(bb, genera)
cc <- cc %>%
  dplyr::arrange(dplyr::case_match(sign,
                                   "NON SIGNIFICANT" ~ 1,
                                   "SPECIALIST" ~ 2,
                                   "GENERALIST" ~ 3))  # GENERALIST 最后绘制

# 作图, distribution of generalists and specialists
library(ggplot2)
color <- c("GENERALIST"="#4169B2", "NON SIGNIFICANT"="#cccccc", "SPECIALIST"="#00bfc4")  # #00bfc4, #0c86f4

theme_set(theme_bw())
p1 <- ggplot(cc, aes(x= abundance, y= niche, color= sign)) + 
  geom_point(size= 2, alpha= 0.7) +
  scale_color_manual(values = color) +
  theme(axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "NA") +
  labs(x= "Mean relative abundance (log10)", y= "Niche width (levins)")

p1

ggsave("Niche_width vs abundance_KO_Generalists.pdf", p1, width = 7, height = 4.6)




#*****************************************************************************
### Unirrigated fields
genera <- read.csv("Unirrigated_Generalists_KO.csv", header = T, row.names = 1)
table(genera$sign)   # GENERALIST    NON SIGNIFICANT    SPECIALIST 
                     #    866            4770             3509 
zotu <- read.csv("RPKM.KOs.count_74.csv", header = T, row.names = 1)
relative <- zotu/colSums(zotu) # Relative abundance
gener <- zotu[rownames(genera), ]
gener.rel <- relative[rownames(genera), ]

## 列为物种名，行为样品名
otu <- t(gener)
width <- niche.width(otu, method = c("levins")) # the average of niche width from all taxa occurring in all samples

### 计算log转化的mean relative abundance
mean <- log10(apply(gener.rel, 1, mean))
head(mean)

aa <- rbind(width, mean)
rownames(aa) <- c("niche", "abundance")
bb <- data.frame(t(aa))
cc <- cbind(bb, genera)
cc <- cc %>%
  dplyr::arrange(dplyr::case_match(sign,
                                   "NON SIGNIFICANT" ~ 1,
                                   "SPECIALIST" ~ 2,
                                   "GENERALIST" ~ 3))  # GENERALIST 最后绘制

# 作图, distribution of generalists and specialists
library(ggplot2)
color <- c("GENERALIST"="#BB2BA0", "NON SIGNIFICANT"="#cccccc", "SPECIALIST"="#00bfc4")

theme_set(theme_bw())
p2 <- ggplot(cc, aes(x= abundance, y= niche, color= sign)) + 
  geom_point(size= 2, alpha= 0.7) +
  scale_color_manual(values = color) +
  theme(axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.position = "NA") +
  labs(x= "Mean relative abundance (log10)", y= "Niche width (levins)")

p2
#ggsave("Niche_width vs abundance_KO_Generalists.pdf", p2, width = 7, height = 4.6)


library(ggpubr)
p11 <- ggarrange(p1, p2, labels = c("(a)","(b)"), nrow = 1, ncol = 2, align = "hv",
                 font.label = list(size = 24, face = "bold", color = "black"), 
                 widths = c(1, 1),     # 每个图宽度
                 common.legend = FALSE) 
p11

ggsave("Distribution of KOs_Irrigated_Non-irrigated.pdf", p11, width = 11, height = 4.5)








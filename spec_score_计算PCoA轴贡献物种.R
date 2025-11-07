package.list <- c('vegan', 'ade4', 'viridis', 'gplots', 'indicspecies')
for (package in package.list) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

design <- read.csv("Design_74_理化和坐标.csv", header = T, row.names = 1)
design.irr <- subset(design, Irrigation.adj == "irrigated")
design.unirr <- subset(design, Irrigation.adj == "non-irrigated")
fish01 <- read.csv("Species.count.csv", header = T, row.names = 1) # colname, species
fish02 <- t(fish01[,rownames(design.unirr)])

fish <- fish02[,rowSums(fish02)>100]
dim(fish)

# Calculate Bray-Curtis 
fish.db <- vegdist(fish, method = "bray")

# Define Order of Sites
order <- rev(attr(fish.db, "Labels"))  

# Plot Heatmap
levelplot(as.matrix(fish.db)[, order], aspect = "iso", col.regions = inferno, 
          xlab = "Doubs Site", ylab = "Doubs Site", scales = list(cex = 0.5), 
          main = "Bray-Curtis Distance")

fish.ward <- hclust(fish.db, method = "ward.D2")

# Plot Cluster
par(mar = c(1, 5, 2, 2) + 0.1)
plot(fish.ward, main = "Doubs River Fish: Ward's Clustering", 
     ylab = "Squared Bray-Curtis Distance")

# PCoA analysis
fish.pcoa <- cmdscale(fish.db, eig = TRUE, k = 3) 

explainvar1 <- round(fish.pcoa$eig[1] / sum(fish.pcoa$eig), 3) * 100
explainvar2 <- round(fish.pcoa$eig[2] / sum(fish.pcoa$eig), 3) * 100
explainvar3 <- round(fish.pcoa$eig[3] / sum(fish.pcoa$eig), 3) * 100
sum.eig <- sum(explainvar1, explainvar2, explainvar3)
sum.eig

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)

# Initiate Plot
plot(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2], #ylim = c(-0.2, 0.7),
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, 
     cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
points(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2],
       pch = 19, cex = 2, bg = "gray", col = "blue") # col = factor(design$Irrigation.adj)
text(fish.pcoa$points[ ,1], fish.pcoa$points[ ,2], 
     labels = row.names(fish.pcoa$points))


# First we calculate the relative abundances of each species at each site
fishREL <- fish   # colname, species
for(i in 1:nrow(fish)){
  fishREL[i, ] = fish[i, ] / sum(fish[i, ])
} 

#library(BiodiversityR)
# Some problems in 2023 with `BiodiversityR`
# Mostly (only) for species scores in PCoA.
# Alternative is to load function from source file
source("spec.scores.function.R")
# Now, we use this information to calculate and add species scores
fish.pcoa <- add.spec.scores(fish.pcoa, fishREL, method = "pcoa.scores")
plot(fish.pcoa$points[, 1], fish.pcoa$points[, 2], 
     xlab = "PCoA Axis 1", ylab = "PCoA Axis 2", col = "blue",
     main = "PCoA Ordination", type = "p", pch = 16) # ,col = "blue"
text(fish.pcoa$cproj[ ,1], fish.pcoa$cproj[ ,2], 
     labels = row.names(fish.pcoa$cproj), col = "black")

# A more quantitative way of identifying influential species involves determining the 
# correlation coefficient of each species along the PCoA axes
spe.corr <- add.spec.scores(fish.pcoa, fishREL, method = "cor.scores")$cproj
corrcut <- 0.7       # user defined cutoff
imp.spp <- spe.corr[abs(spe.corr[, 1]) >= corrcut | abs(spe.corr[, 2]) >= corrcut, ]

library(dplyr)
imp.spp
imp.spp01 <- na.omit(imp.spp)

write.csv(imp.spp01, file = "PCoA轴关联物种-spec.scores分析_Unirrigated.csv")

# Permutation Test for Species Abundances Across Axes
fit <- envfit(fish.pcoa, fishREL, perm = 999)
fit
# 导出结果
r.value <- as.data.frame(fit$vectors$r)
pval <- as.data.frame(fit$vectors$pvals)
r.value$variable <- rownames(r.value)
stats <- data.frame(r.value = r.value, pvalue = pval)
write.csv(stats, file = "add.species.score.csv")

#******************************Indicator species******************************
# Alternatively, we may ask about the habitat preferences of each species
design <- read.csv("Design_74_理化和坐标.csv", header = T, row.names = 1)
# design.irr <- subset(design, Irrigation.adj == "irrigated")
# design.unirr <- subset(design, Irrigation.adj == "non-irrigated")
fish01 <- read.csv("Species.count.csv", header = T, row.names = 1) # colname, species
fish02 <- t(fish01[,rownames(design)])
fish <- fish02[,colSums(fish02)>50]
dim(fish)
plant <- read.csv("NDVI_paritial correlation.偏相关分析数据.csv", header = T, row.names = 1)

fish <- t(fish01[,rownames(design)])
fish <- fish[,colSums(fish)>50]
# Run PERMANOVA with adonis function
adonis2(fish ~ plant$NDVI, method = "bray", permutations = 999)

# To discern how individual species relate to groups of sites using Indicator species analysis
indval <- multipatt(fish, cluster = plant$NDVI, func = "IndVal.g", control = how(nperm=999))
summary(indval)

# we may ask about the habitat preferences of each species
fish.rel <- decostand(fish, method = "total")
phi <- multipatt(fish.rel, cluster = quality, func = "r.g", control = how(nperm=999))
summary(phi)




library(dplyr)

As <- read.delim("AsCycDB_results2.txt",header = F)
As.map <- read.delim("asgene.map", sep = " ", header = F)
names(As.map) <- c("id", "gene", "source")
As2 <- As[,1:2][which(As$V12>=60),]
names(As2) <- c("orf","id")
As3 <- left_join(As2, As.map, by = c("id"))
As4 <- As3[!is.na(As3$gene),]
write.csv(As4, "AsgeneDB_anno.csv", row.names = F)


M <- read.delim("MCyc_results2.txt",header = F)
M.map <- read.delim("id2gene.map.txt", sep = "\t", header = F)
names(M.map) <- c("id", "gene", "source")
M2 <- M[,1:2][which(M$V12>=60),]
names(M2) <- c("orf","id")
M3 <- left_join(M2, M.map, by = c("id"))
M4 <- M3[!is.na(M3$gene),]
write.csv(M4, "MCyc_anno.csv", row.names = F)

N <- read.delim("NCyc_results2.txt",header = F)
N.map <- read.delim("id2gene.map.2019Jul.txt", sep = "\t", header = F)
names(N.map) <- c("id", "gene")
N2 <- N[,1:2][which(N$V12>=60),]
names(N2) <- c("orf","id")
N3 <- left_join(N2, N.map, by = c("id"))
N4 <- N3[!is.na(N3$gene),]
write.csv(N4, "NCyc_anno.csv", row.names = F)

P <- read.delim("PCyc_results2.txt",header = F)
P.map <- read.delim("id2genemap.txt", sep = "\t", header = F)
names(P.map) <- c("id", "gene", "source")
P2 <- P[,1:2][which(P$V12>=60),]
names(P2) <- c("orf","id")
P3 <- left_join(P2, P.map, by = c("id"))
P4 <- P3[!is.na(P3$gene),]
write.csv(P4, "PCyc_anno.csv", row.names = F)

S <- read.delim("SCyc_results2.txt",header = F)
S.map <- read.delim("id2gene.2020Mar.map", sep = "\t", header = F)
names(S.map) <- c("id", "gene", "source")
S2 <- S[,1:2][which(S$V12>=60),]
names(S2) <- c("orf","id")
S3 <- left_join(S2, S.map, by = c("id"))
S4 <- S3[!is.na(S3$gene),]
write.csv(S4, "SCyc_anno.csv", row.names = F)

V <- read.delim("VFDB_results2.txt",header = F)
V.map <- read.delim("VFDB_setB_anno.txt", sep = "\t", header = F)
names(V.map) <- c("id", "gene", "describe")
V2 <- V[,1:2][which(V$V12>=60),]
names(V2) <- c("orf","id")
V3 <- left_join(V2, V.map, by = c("id"))
V4 <- V3[!is.na(V3$gene),]
write.csv(V4, "VFDB_setB_anno.csv", row.names = F)


ARGs <- read.delim("ARG_results2.txt",header = F)
ARGs.map <- read.csv("ARG_anno.csv")
names(ARGs.map) <- c("gene", "id")
ARGs2 <- ARGs[,1:2][which(ARGs$V12>=60),]
names(ARGs2) <- c("orf","id")
ARGs3 <- left_join(ARGs2, ARGs.map, by = c("id"))
ARGs4 <- ARGs3[!is.na(ARGs3$gene),]
write.csv(ARGs4, "ARGs_anno.csv", row.names = F)

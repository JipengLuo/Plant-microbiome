library(dplyr)
library(stringr)
library(tidyverse)
##Importing the data
s01 <- read.csv("s01_gene_rpkm_coverage.csv", header=T, stringsAsFactors=FALSE)
# 将第一列分裂成两列
s1 <- s01 %>%
  separate(Contig.s01_gene_filter.RPKM, into = c("name", "s01"), sep = "\t")
rownames(s1) <- s1[[1]]
s001 <- as.matrix(s1[,-1])
rownames(s001) <- rownames(s1)
colnames(s001) <- "s01"
write.csv(s001, file = "s01.rpkm.csv")

s02 <- read.csv("s02_gene_rpkm_coverage.csv", header = T)
s2 <- s02 %>%
  separate(Contig.s02_gene_filter.RPKM, into = c("name", "s02"), sep = "\t")
rownames(s2) <- s2[[1]]
s002 <- as.matrix(s2[,-1])
rownames(s002) <- rownames(s2)
colnames(s002) <- "s02"
write.csv(s002, file = "s02.rpkm.csv")

s03 <- read.csv("s03_gene_rpkm_coverage.csv", header = T)
s3 <- s03 %>%
  separate(Contig.s03_gene_filter.RPKM, into = c("name", "s03"), sep = "\t")
rownames(s3) <- s3[[1]]
s003 <- as.matrix(s3[,-1])
rownames(s003) <- rownames(s3)
colnames(s003) <- "s03"
write.csv(s003, file = "s03.rpkm.csv")

s04 <- read.csv("s04_gene_rpkm_coverage.csv", header = T)
s4 <- s04 %>%
  separate(Contig.s04_gene_filter.RPKM, into = c("name", "s04"), sep = "\t")
rownames(s4) <- s4[[1]]
s004 <- as.matrix(s4[,-1])
rownames(s004) <- rownames(s4)
colnames(s004) <- "s04"
write.csv(s004, file = "s04.rpkm.csv")

s05 <- read.csv("s05_gene_rpkm_coverage.csv", header = T)
s5 <- s05 %>%
  separate(Contig.s05_gene_filter.RPKM, into = c("name", "s05"), sep = "\t")
rownames(s5) <- s5[[1]]
s005 <- as.matrix(s5[,-1])
rownames(s005) <- rownames(s5)
colnames(s005) <- "s05"
write.csv(s005, file = "s05.rpkm.csv")

s06 <- read.csv("s06_gene_rpkm_coverage.csv", header = T)
s6 <- s06 %>%
  separate(Contig.s06_gene_filter.RPKM, into = c("name", "s06"), sep = "\t")
rownames(s6) <- s6[[1]]
s006 <- as.matrix(s6[,-1])
rownames(s006) <- rownames(s6)
colnames(s006) <- "s06"
write.csv(s006, file = "s06.rpkm.csv")

s07 <- read.csv("s07_gene_rpkm_coverage.csv", header = T)
s7 <- s07 %>%
  separate(Contig.s07_gene_filter.RPKM, into = c("name", "s07"), sep = "\t")
rownames(s7) <- s7[[1]]
s007 <- as.matrix(s7[,-1])
rownames(s007) <- rownames(s7)
colnames(s007) <- "s07"
write.csv(s007, file = "s07.rpkm.csv")

s08 <- read.csv("s08_gene_rpkm_coverage.csv", header = T)
s8 <- s08 %>%
  separate(Contig.s08_gene_filter.RPKM, into = c("name", "s08"), sep = "\t")
rownames(s8) <- s8[[1]]
s008 <- as.matrix(s8[,-1])
rownames(s008) <- rownames(s8)
colnames(s008) <- "s08"
write.csv(s008, file = "s08.rpkm.csv")

s09 <- read.csv("s09_gene_rpkm_coverage.csv", header = T)
s9 <- s09 %>%
  separate(Contig.s09_gene_filter.RPKM, into = c("name", "s09"), sep = "\t")
rownames(s9) <- s9[[1]]
s009 <- as.matrix(s9[,-1])
rownames(s009) <- rownames(s9)
colnames(s009) <- "s09"
write.csv(s009, file = "s09.rpkm.csv")

s10 <- read.csv("s10_gene_rpkm_coverage.csv", header = T)
s10 <- s10 %>%
  separate(Contig.s10_gene_filter.RPKM, into = c("name", "s10"), sep = "\t")
rownames(s10) <- s10[[1]]
s010 <- as.matrix(s10[,-1])
rownames(s010) <- rownames(s10)
colnames(s010) <- "s10"
write.csv(s010, file = "s10.rpkm.csv")

s11 <- read.csv("s11_gene_rpkm_coverage.csv", header = T)
s11 <- s11 %>%
  separate(Contig.s11_gene_filter.RPKM, into = c("name", "s11"), sep = "\t")
rownames(s11) <- s11[[1]]
s011 <- as.matrix(s11[,-1])
rownames(s011) <- rownames(s11)
colnames(s011) <- "s11"
write.csv(s011, file = "s11.rpkm.csv")

s12 <- read.csv("s12_gene_rpkm_coverage.csv", header = T)
s12 <- s12 %>%
  separate(Contig.s12_gene_filter.RPKM, into = c("name", "s12"), sep = "\t")
rownames(s12) <- s12[[1]]
s012 <- as.matrix(s12[,-1])
rownames(s012) <- rownames(s12)
colnames(s012) <- "s12"
write.csv(s012, file = "s12.rpkm.csv")

s13 <- read.csv("s13_gene_rpkm_coverage.csv", header = T)
s13 <- s13 %>%
  separate(Contig.s13_gene_filter.RPKM, into = c("name", "s13"), sep = "\t")
rownames(s13) <- s13[[1]]
s013 <- as.matrix(s13[,-1])
rownames(s013) <- rownames(s13)
colnames(s013) <- "s13"
write.csv(s013, file = "s13.rpkm.csv")

s14 <- read.csv("s14_gene_rpkm_coverage.csv", header = T)
s14 <- s14 %>%
  separate(Contig.s14_gene_filter.RPKM, into = c("name", "s14"), sep = "\t")
rownames(s14) <- s14[[1]]
s014 <- as.matrix(s14[,-1])
rownames(s014) <- rownames(s14)
colnames(s014) <- "s14"
write.csv(s014, file = "s14.rpkm.csv")

s15 <- read.csv("s15_gene_rpkm_coverage.csv", header = T)
s15 <- s15 %>%
  separate(Contig.s15_gene_filter.RPKM, into = c("name", "s15"), sep = "\t")
rownames(s15) <- s15[[1]]
s015 <- as.matrix(s15[,-1])
rownames(s015) <- rownames(s15)
colnames(s015) <- "s15"
write.csv(s015, file = "s15.rpkm.csv")

s16 <- read.csv("s16_gene_rpkm_coverage.csv", header = T)
s16 <- s16 %>%
  separate(Contig.s16_gene_filter.RPKM, into = c("name", "s16"), sep = "\t")
rownames(s16) <- s16[[1]]
s016 <- as.matrix(s16[,-1])
rownames(s016) <- rownames(s16)
colnames(s016) <- "s16"
write.csv(s016, file = "s16.rpkm.csv")

s17 <- read.csv("s17_gene_rpkm_coverage.csv", header = T)
s17 <- s17 %>%
  separate(Contig.s17_gene_filter.RPKM, into = c("name", "s17"), sep = "\t")
rownames(s17) <- s17[[1]]
s017 <- as.matrix(s17[,-1])
rownames(s017) <- rownames(s17)
colnames(s017) <- "s17"
write.csv(s017, file = "s17.rpkm.csv")

s18 <- read.csv("s18_gene_rpkm_coverage.csv", header = T)
s18 <- s18 %>%
  separate(Contig.s18_gene_filter.RPKM, into = c("name", "s18"), sep = "\t")
rownames(s18) <- s18[[1]]
s018 <- as.matrix(s18[,-1])
rownames(s018) <- rownames(s18)
colnames(s018) <- "s18"
write.csv(s018, file = "s18.rpkm.csv")

s19 <- read.csv("s19_gene_rpkm_coverage.csv", header = T)
s19 <- s19 %>%
  separate(Contig.s19_gene_filter.RPKM, into = c("name", "s19"), sep = "\t")
rownames(s19) <- s19[[1]]
s019 <- as.matrix(s19[,-1])
rownames(s019) <- rownames(s19)
colnames(s019) <- "s19"
write.csv(s019, file = "s19.rpkm.csv")

s20 <- read.csv("s20_gene_rpkm_coverage.csv", header = T)
s20 <- s20 %>%
  separate(Contig.s20_gene_filter.RPKM, into = c("name", "s20"), sep = "\t")
rownames(s20) <- s20[[1]]
s020 <- as.matrix(s20[,-1])
rownames(s020) <- rownames(s20)
colnames(s020) <- "s20"
write.csv(s020, file = "s20.rpkm.csv")

s21 <- read.csv("s21_gene_rpkm_coverage.csv", header = T)
s21 <- s21 %>%
  separate(Contig.s21_gene_filter.RPKM, into = c("name", "s21"), sep = "\t")
rownames(s21) <- s21[[1]]
s021 <- as.matrix(s21[,-1])
rownames(s021) <- rownames(s21)
colnames(s021) <- "s21"
write.csv(s021, file = "s21.rpkm.csv")

s22 <- read.csv("s22_gene_rpkm_coverage.csv", header = T)
s22 <- s22 %>%
  separate(Contig.s21_gene_filter.RPKM, into = c("name", "s22"), sep = "\t")
rownames(s22) <- s22[[1]]
s022 <- as.matrix(s22[,-1])
rownames(s022) <- rownames(s22)
colnames(s022) <- "s22"
write.csv(s022, file = "s22.rpkm.csv")

s23 <- read.csv("s23_gene_rpkm_coverage.csv", header = T)
s23 <- s23 %>%
  separate(Contig.s23_gene_filter.RPKM, into = c("name", "s23"), sep = "\t")
rownames(s23) <- s23[[1]]
s023 <- as.matrix(s23[,-1])
rownames(s023) <- rownames(s23)
colnames(s023) <- "s23"
write.csv(s023, file = "s23.rpkm.csv")

s24 <- read.csv("s24_gene_rpkm_coverage.csv", header = T)
s24 <- s24 %>%
  separate(Contig.s24_gene_filter.RPKM, into = c("name", "s24"), sep = "\t")
rownames(s24) <- s24[[1]]
s024 <- as.matrix(s24[,-1])
rownames(s024) <- rownames(s24)
colnames(s024) <- "s24"
write.csv(s024, file = "s24.rpkm.csv")

s25 <- read.csv("s25_gene_rpkm_coverage.csv", header = T)
s25 <- s25 %>%
  separate(Contig.s25_gene_filter.RPKM, into = c("name", "s25"), sep = "\t")
rownames(s25) <- s25[[1]]
s025 <- as.matrix(s25[,-1])
rownames(s025) <- rownames(s25)
colnames(s025) <- "s25"
write.csv(s025, file = "s25.rpkm.csv")

s26 <- read.csv("s26_gene_rpkm_coverage.csv", header = T)
s26 <- s26 %>%
  separate(Contig.s26_gene_filter.RPKM, into = c("name", "s26"), sep = "\t")
rownames(s26) <- s26[[1]]
s026 <- as.matrix(s26[,-1])
rownames(s026) <- rownames(s26)
colnames(s026) <- "s26"
write.csv(s026, file = "s26.rpkm.csv")

s27 <- read.csv("s27_gene_rpkm_coverage.csv", header = T)
s27 <- s27 %>%
  separate(Contig.s27_gene_filter.RPKM, into = c("name", "s27"), sep = "\t")
rownames(s27) <- s27[[1]]
s027 <- as.matrix(s27[,-1])
rownames(s027) <- rownames(s27)
colnames(s027) <- "s27"
write.csv(s027, file = "s27.rpkm.csv")

s28 <- read.csv("s28_gene_rpkm_coverage.csv", header = T)
s28 <- s28 %>%
  separate(Contig.s28_gene_filter.RPKM, into = c("name", "s28"), sep = "\t")
rownames(s28) <- s28[[1]]
s028 <- as.matrix(s28[,-1])
rownames(s028) <- rownames(s28)
colnames(s028) <- "s28"
write.csv(s028, file = "s28.rpkm.csv")

s29 <- read.csv("s29_gene_rpkm_coverage.csv", header = T)
s29 <- s29 %>%
  separate(Contig.s29_gene_filter.RPKM, into = c("name", "s29"), sep = "\t")
rownames(s29) <- s29[[1]]
s029 <- as.matrix(s29[,-1])
rownames(s029) <- rownames(s29)
colnames(s029) <- "s29"
write.csv(s029, file = "s29.rpkm.csv")

s30 <- read.csv("s30_gene_rpkm_coverage.csv", header = T)
s30 <- s30 %>%
  separate(Contig.s30_gene_filter.RPKM, into = c("name", "30"), sep = "\t")
rownames(s30) <- s30[[1]]
s030 <- as.matrix(s30[,-1])
rownames(s030) <- rownames(s30)
colnames(s030) <- "s30"
write.csv(s030, file = "s30.rpkm.csv")

s31 <- read.csv("s31_gene_rpkm_coverage.csv", header = T)
s31 <- s31 %>%
  separate(Contig.s31_gene_filter.RPKM, into = c("name", "31"), sep = "\t")
rownames(s31) <- s31[[1]]
s031 <- as.matrix(s31[,-1])
rownames(s031) <- rownames(s31)
colnames(s031) <- "s31"
write.csv(s031, file = "s31.rpkm.csv")

s32 <- read.csv("s32_gene_rpkm_coverage.csv", header = T)
s32 <- s32 %>%
  separate(Contig.s32_gene_filter.RPKM, into = c("name", "32"), sep = "\t")
rownames(s32) <- s32[[1]]
s032 <- as.matrix(s32[,-1])
rownames(s032) <- rownames(s32)
colnames(s032) <- "s32"
write.csv(s032, file = "s32.rpkm.csv")

s33 <- read.csv("s33_gene_rpkm_coverage.csv", header = T)
s33 <- s33 %>%
  separate(Contig.s33_gene_filter.RPKM, into = c("name", "33"), sep = "\t")
rownames(s33) <- s33[[1]]
s033 <- as.matrix(s33[,-1])
rownames(s033) <- rownames(s33)
colnames(s033) <- "s33"
write.csv(s033, file = "s33.rpkm.csv")

s34 <- read.csv("s34_gene_rpkm_coverage.csv", header = T)
s34 <- s34 %>%
  separate(Contig.s34_gene_filter.RPKM, into = c("name", "34"), sep = "\t")
rownames(s34) <- s34[[1]]
s034 <- as.matrix(s34[,-1])
rownames(s034) <- rownames(s34)
colnames(s034) <- "s34"
write.csv(s034, file = "s34.rpkm.csv")

s35 <- read.csv("s35_gene_rpkm_coverage.csv", header = T)
s35 <- s35 %>%
  separate(Contig.s35_gene_filter.RPKM, into = c("name", "35"), sep = "\t")
rownames(s35) <- s35[[1]]
s035 <- as.matrix(s35[,-1])
rownames(s035) <- rownames(s35)
colnames(s035) <- "s35"
write.csv(s035, file = "s35.rpkm.csv")

s36 <- read.csv("s36_gene_rpkm_coverage.csv", header = T)
s36 <- s36 %>%
  separate(Contig.s36_gene_filter.RPKM, into = c("name", "36"), sep = "\t")
rownames(s36) <- s36[[1]]
s036 <- as.matrix(s36[,-1])
rownames(s036) <- rownames(s36)
colnames(s036) <- "s36"
write.csv(s036, file = "s36.rpkm.csv")

s37 <- read.csv("s37_gene_rpkm_coverage.csv", header = T)
s37 <- s37 %>%
  separate(Contig.s37_gene_filter.RPKM, into = c("name", "37"), sep = "\t")
rownames(s37) <- s37[[1]]
s037 <- as.matrix(s37[,-1])
rownames(s037) <- rownames(s37)
colnames(s037) <- "s37"
write.csv(s037, file = "s37.rpkm.csv")

s38 <- read.csv("s38_gene_rpkm_coverage.csv", header = T)
s38 <- s38 %>%
  separate(Contig.s38_gene_filter.RPKM, into = c("name", "38"), sep = "\t")
rownames(s38) <- s38[[1]]
s038 <- as.matrix(s38[,-1])
rownames(s038) <- rownames(s38)
colnames(s038) <- "s38"
write.csv(s038, file = "s38.rpkm.csv")

s39 <- read.csv("s39_gene_rpkm_coverage.csv", header = T)
s39 <- s39 %>%
  separate(Contig.s39_gene_filter.RPKM, into = c("name", "39"), sep = "\t")
rownames(s39) <- s39[[1]]
s039 <- as.matrix(s39[,-1])
rownames(s039) <- rownames(s39)
colnames(s039) <- "s39"
write.csv(s039, file = "s39.rpkm.csv")

s40 <- read.csv("s40_gene_rpkm_coverage.csv", header = T)
s40 <- s40 %>%
  separate(Contig.s40_gene_filter.RPKM, into = c("name", "40"), sep = "\t")
rownames(s40) <- s40[[1]]
s040 <- as.matrix(s40[,-1])
rownames(s040) <- rownames(s40)
colnames(s040) <- "s40"
write.csv(s040, file = "s40.rpkm.csv")

s41 <- read.csv("s41_gene_rpkm_coverage.csv", header = T)
s41 <- s41 %>%
  separate(Contig.s41_gene_filter.RPKM, into = c("name", "41"), sep = "\t")
rownames(s41) <- s41[[1]]
s041 <- as.matrix(s41[,-1])
rownames(s041) <- rownames(s41)
colnames(s041) <- "s41"
write.csv(s041, file = "s41.rpkm.csv")

s42 <- read.csv("s42_gene_rpkm_coverage.csv", header = T)
s42 <- s42 %>%
  separate(Contig.s42_gene_filter.RPKM, into = c("name", "42"), sep = "\t")
rownames(s42) <- s42[[1]]
s042 <- as.matrix(s42[,-1])
rownames(s042) <- rownames(s42)
colnames(s042) <- "s42"
write.csv(s042, file = "s42.rpkm.csv")

s43 <- read.csv("s43_gene_rpkm_coverage.csv", header = T)
s43 <- s43 %>%
  separate(Contig.s43_gene_filter.RPKM, into = c("name", "43"), sep = "\t")
rownames(s43) <- s43[[1]]
s043 <- as.matrix(s43[,-1])
rownames(s043) <- rownames(s43)
colnames(s043) <- "s43"
write.csv(s043, file = "s43.rpkm.csv")

s44 <- read.csv("s44_gene_rpkm_coverage.csv", header = T)
s44 <- s44 %>%
  separate(Contig.s44_gene_filter.RPKM, into = c("name", "44"), sep = "\t")
rownames(s44) <- s44[[1]]
s044 <- as.matrix(s44[,-1])
rownames(s044) <- rownames(s44)
colnames(s044) <- "s44"
write.csv(s044, file = "s44.rpkm.csv")

s45 <- read.csv("s45_gene_rpkm_coverage.csv", header = T)
s45 <- s45 %>%
  separate(Contig.s45_gene_filter.RPKM, into = c("name", "45"), sep = "\t")
rownames(s45) <- s45[[1]]
s045 <- as.matrix(s45[,-1])
rownames(s045) <- rownames(s45)
colnames(s045) <- "s45"
write.csv(s045, file = "s45.rpkm.csv")

s46 <- read.csv("s46_gene_rpkm_coverage.csv", header = T)
s46 <- s46 %>%
  separate(Contig.s46_gene_filter.RPKM, into = c("name", "46"), sep = "\t")
rownames(s46) <- s46[[1]]
s046 <- as.matrix(s46[,-1])
rownames(s046) <- rownames(s46)
colnames(s046) <- "s46"
write.csv(s046, file = "s46.rpkm.csv")

s47 <- read.csv("s47_gene_rpkm_coverage.csv", header = T)
s47 <- s47 %>%
  separate(Contig.s47_gene_filter.RPKM, into = c("name", "47"), sep = "\t")
rownames(s47) <- s47[[1]]
s047 <- as.matrix(s47[,-1])
rownames(s047) <- rownames(s47)
colnames(s047) <- "s47"
write.csv(s047, file = "s47.rpkm.csv")

s48 <- read.csv("s48_gene_rpkm_coverage.csv", header = T)
s48 <- s48 %>%
  separate(Contig.s48_gene_filter.RPKM, into = c("name", "48"), sep = "\t")
rownames(s48) <- s48[[1]]
s048 <- as.matrix(s48[,-1])
rownames(s048) <- rownames(s48)
colnames(s048) <- "s48"
write.csv(s048, file = "s48.rpkm.csv")

s49 <- read.csv("s49_gene_rpkm_coverage.csv", header = T)
s49 <- s49 %>%
  separate(Contig.s49_gene_filter.RPKM, into = c("name", "49"), sep = "\t")
rownames(s49) <- s49[[1]]
s049 <- as.matrix(s49[,-1])
rownames(s049) <- rownames(s49)
colnames(s049) <- "s49"
write.csv(s049, file = "s49.rpkm.csv")

s50 <- read.csv("s50_gene_rpkm_coverage.csv", header = T)
s50 <- s50 %>%
  separate(Contig.s50_gene_filter.RPKM, into = c("name", "50"), sep = "\t")
rownames(s50) <- s50[[1]]
s050 <- as.matrix(s50[,-1])
rownames(s050) <- rownames(s50)
colnames(s050) <- "s50"
write.csv(s050, file = "s50.rpkm.csv")

s51 <- read.csv("s51_gene_rpkm_coverage.csv", header = T)
s51 <- s51 %>%
  separate(Contig.s51_gene_filter.RPKM, into = c("name", "51"), sep = "\t")
rownames(s51) <- s51[[1]]
s051 <- as.matrix(s51[,-1])
rownames(s051) <- rownames(s51)
colnames(s051) <- "s51"
write.csv(s051, file = "s51.rpkm.csv")

s52 <- read.csv("s52_gene_rpkm_coverage.csv", header = T)
s52 <- s52 %>%
  separate(Contig.s52_gene_filter.RPKM, into = c("name", "52"), sep = "\t")
rownames(s52) <- s52[[1]]
s052 <- as.matrix(s52[,-1])
rownames(s052) <- rownames(s52)
colnames(s052) <- "s52"
write.csv(s052, file = "s52.rpkm.csv")

s53 <- read.csv("s53_gene_rpkm_coverage.csv", header = T)
s53 <- s53 %>%
  separate(Contig.s53_gene_filter.RPKM, into = c("name", "53"), sep = "\t")
rownames(s53) <- s53[[1]]
s053 <- as.matrix(s53[,-1])
rownames(s053) <- rownames(s53)
colnames(s053) <- "s53"
write.csv(s053, file = "s53.rpkm.csv")

s54 <- read.csv("s54_gene_rpkm_coverage.csv", header = T)
s54 <- s54 %>%
  separate(Contig.s54_gene_filter.RPKM, into = c("name", "54"), sep = "\t")
rownames(s54) <- s54[[1]]
s054 <- as.matrix(s54[,-1])
rownames(s054) <- rownames(s54)
colnames(s054) <- "s54"
write.csv(s054, file = "s54.rpkm.csv")

s55 <- read.csv("s55_gene_rpkm_coverage.csv", header = T)
s55 <- s55 %>%
  separate(Contig.s55_gene_filter.RPKM, into = c("name", "55"), sep = "\t")
rownames(s55) <- s55[[1]]
s055 <- as.matrix(s55[,-1])
rownames(s055) <- rownames(s55)
colnames(s055) <- "s55"
write.csv(s055, file = "s55.rpkm.csv")

s56 <- read.csv("s56_gene_rpkm_coverage.csv", header = T)
s56 <- s56 %>%
  separate(Contig.s56_gene_filter.RPKM, into = c("name", "56"), sep = "\t")
rownames(s56) <- s56[[1]]
s056 <- as.matrix(s56[,-1])
rownames(s056) <- rownames(s56)
colnames(s056) <- "s56"
write.csv(s056, file = "s56.rpkm.csv")

s57 <- read.csv("s57_gene_rpkm_coverage.csv", header = T)
s57 <- s57 %>%
  separate(Contig.s57_gene_filter.RPKM, into = c("name", "57"), sep = "\t")
rownames(s57) <- s57[[1]]
s057 <- as.matrix(s57[,-1])
rownames(s057) <- rownames(s57)
colnames(s057) <- "s57"
write.csv(s057, file = "s57.rpkm.csv")

s58 <- read.csv("s58_gene_rpkm_coverage.csv", header = T)
s58 <- s58 %>%
  separate(Contig.s58_gene_filter.RPKM, into = c("name", "58"), sep = "\t")
rownames(s58) <- s58[[1]]
s058 <- as.matrix(s58[,-1])
rownames(s058) <- rownames(s58)
colnames(s058) <- "s58"
write.csv(s058, file = "s58.rpkm.csv")

s59 <- read.csv("s59_gene_rpkm_coverage.csv", header = T)
s59 <- s59 %>%
  separate(Contig.s59_gene_filter.RPKM, into = c("name", "59"), sep = "\t")
rownames(s59) <- s59[[1]]
s059 <- as.matrix(s59[,-1])
rownames(s059) <- rownames(s59)
colnames(s059) <- "s59"
write.csv(s059, file = "s59.rpkm.csv")

s60 <- read.csv("s60_gene_rpkm_coverage.csv", header = T)
s60 <- s60 %>%
  separate(Contig.s60_gene_filter.RPKM, into = c("name", "60"), sep = "\t")
rownames(s60) <- s60[[1]]
s060 <- as.matrix(s60[,-1])
rownames(s060) <- rownames(s60)
colnames(s060) <- "s60"
write.csv(s060, file = "s60.rpkm.csv")

s61 <- read.csv("s61_gene_rpkm_coverage.csv", header = T)
s61 <- s61 %>%
  separate(Contig.s61_gene_filter.RPKM, into = c("name", "61"), sep = "\t")
rownames(s61) <- s61[[1]]
s061 <- as.matrix(s61[,-1])
rownames(s061) <- rownames(s61)
colnames(s061) <- "s61"
write.csv(s061, file = "s61.rpkm.csv")

s62 <- read.csv("s62_gene_rpkm_coverage.csv", header = T)
s62 <- s62 %>%
  separate(Contig.s62_gene_filter.RPKM, into = c("name", "62"), sep = "\t")
rownames(s62) <- s62[[1]]
s062 <- as.matrix(s62[,-1])
rownames(s062) <- rownames(s62)
colnames(s062) <- "s62"
write.csv(s062, file = "s62.rpkm.csv")

s63 <- read.csv("s63_gene_rpkm_coverage.csv", header = T)
s63 <- s63 %>%
  separate(Contig.s63_gene_filter.RPKM, into = c("name", "63"), sep = "\t")
rownames(s63) <- s63[[1]]
s063 <- as.matrix(s63[,-1])
rownames(s063) <- rownames(s63)
colnames(s063) <- "s63"
write.csv(s063, file = "s63.rpkm.csv")

s64 <- read.csv("s64_gene_rpkm_coverage.csv", header = T)
s64 <- s64 %>%
  separate(Contig.s64_gene_filter.RPKM, into = c("name", "64"), sep = "\t")
rownames(s64) <- s64[[1]]
s064 <- as.matrix(s64[,-1])
rownames(s064) <- rownames(s64)
colnames(s064) <- "s64"
write.csv(s064, file = "s64.rpkm.csv")

s65 <- read.csv("s65_gene_rpkm_coverage.csv", header = T)
s65 <- s65 %>%
  separate(Contig.s65_gene_filter.RPKM, into = c("name", "65"), sep = "\t")
rownames(s65) <- s65[[1]]
s065 <- as.matrix(s65[,-1])
rownames(s065) <- rownames(s65)
colnames(s065) <- "s65"
write.csv(s065, file = "s65.rpkm.csv")

s66 <- read.csv("s66_gene_rpkm_coverage.csv", header = T)
s66 <- s66 %>%
  separate(Contig.s66_gene_filter.RPKM, into = c("name", "66"), sep = "\t")
rownames(s66) <- s66[[1]]
s066 <- as.matrix(s66[,-1])
rownames(s066) <- rownames(s66)
colnames(s066) <- "s66"
write.csv(s066, file = "s66.rpkm.csv")

s67 <- read.csv("s67_gene_rpkm_coverage.csv", header = T)
s67 <- s67 %>%
  separate(Contig.s67_gene_filter.RPKM, into = c("name", "67"), sep = "\t")
rownames(s67) <- s67[[1]]
s067 <- as.matrix(s67[,-1])
rownames(s067) <- rownames(s67)
colnames(s067) <- "s67"
write.csv(s067, file = "s67.rpkm.csv")

s68 <- read.csv("s68_gene_rpkm_coverage.csv", header = T)
s68 <- s68 %>%
  separate(Contig.s68_gene_filter.RPKM, into = c("name", "68"), sep = "\t")
rownames(s68) <- s68[[1]]
s068 <- as.matrix(s68[,-1])
rownames(s068) <- rownames(s68)
colnames(s068) <- "s68"
write.csv(s068, file = "s68.rpkm.csv")

s69 <- read.csv("s69_gene_rpkm_coverage.csv", header = T)
s69 <- s69 %>%
  separate(Contig.s69_gene_filter.RPKM, into = c("name", "69"), sep = "\t")
rownames(s69) <- s69[[1]]
s069 <- as.matrix(s69[,-1])
rownames(s069) <- rownames(s69)
colnames(s069) <- "s69"
write.csv(s069, file = "s69.rpkm.csv")

s70 <- read.csv("s70_gene_rpkm_coverage.csv", header = T)
s70 <- s70 %>%
  separate(Contig.s70_gene_filter.RPKM, into = c("name", "70"), sep = "\t")
rownames(s70) <- s70[[1]]
s070 <- as.matrix(s70[,-1])
rownames(s070) <- rownames(s70)
colnames(s070) <- "s70"
write.csv(s070, file = "s70.rpkm.csv")

s71 <- read.csv("s71_gene_rpkm_coverage.csv", header = T)
s71 <- s71 %>%
  separate(Contig.s71_gene_filter.RPKM, into = c("name", "71"), sep = "\t")
rownames(s71) <- s71[[1]]
s071 <- as.matrix(s71[,-1])
rownames(s071) <- rownames(s71)
colnames(s071) <- "s71"
write.csv(s071, file = "s71.rpkm.csv")

s72 <- read.csv("s72_gene_rpkm_coverage.csv", header = T)
s72 <- s72 %>%
  separate(Contig.s72_gene_filter.RPKM, into = c("name", "72"), sep = "\t")
rownames(s72) <- s72[[1]]
s072 <- as.matrix(s72[,-1])
rownames(s072) <- rownames(s72)
colnames(s072) <- "s72"
write.csv(s072, file = "s72.rpkm.csv")

s73 <- read.csv("s73_gene_rpkm_coverage.csv", header = T)
s73 <- s73 %>%
  separate(Contig.s73_gene_filter.RPKM, into = c("name", "73"), sep = "\t")
rownames(s73) <- s73[[1]]
s073 <- as.matrix(s73[,-1])
rownames(s073) <- rownames(s73)
colnames(s073) <- "s73"
write.csv(s073, file = "s73.rpkm.csv")

s74 <- read.csv("s74_gene_rpkm_coverage.csv", header = T)
s74 <- s74 %>%
  separate(Contig.s74_gene_filter.RPKM, into = c("name", "74"), sep = "\t")
rownames(s74) <- s74[[1]]
s074 <- as.matrix(s74[,-1])
rownames(s074) <- rownames(s74)
colnames(s074) <- "s74"
write.csv(s074, file = "s74.rpkm.csv")

s75 <- read.csv("s75_gene_rpkm_coverage.csv", header = T)
s75 <- s75 %>%
  separate(Contig.s75_gene_filter.RPKM, into = c("name", "75"), sep = "\t")
rownames(s75) <- s75[[1]]
s075 <- as.matrix(s75[,-1])
rownames(s075) <- rownames(s75)
colnames(s075) <- "s75"
write.csv(s075, file = "s75.rpkm.csv")

s76 <- read.csv("s76_gene_rpkm_coverage.csv", header = T)
s76 <- s76 %>%
  separate(Contig.s76_gene_filter.RPKM, into = c("name", "76"), sep = "\t")
rownames(s76) <- s76[[1]]
s076 <- as.matrix(s76[,-1])
rownames(s076) <- rownames(s76)
colnames(s076) <- "s76"
write.csv(s076, file = "s76.rpkm.csv")

s77 <- read.csv("s77_gene_rpkm_coverage.csv", header = T)
s77 <- s77 %>%
  separate(Contig.s77_gene_filter.RPKM, into = c("name", "77"), sep = "\t")
rownames(s77) <- s77[[1]]
s077 <- as.matrix(s77[,-1])
rownames(s077) <- rownames(s77)
colnames(s077) <- "s77"
write.csv(s077, file = "s77.rpkm.csv")

s78 <- read.csv("s78_gene_rpkm_coverage.csv", header = T)
s78 <- s78 %>%
  separate(Contig.s78_gene_filter.RPKM, into = c("name", "78"), sep = "\t")
rownames(s78) <- s78[[1]]
s078 <- as.matrix(s78[,-1])
rownames(s078) <- rownames(s78)
colnames(s078) <- "s78"
write.csv(s078, file = "s78.rpkm.csv")

s87 <- read.csv("s87_gene_rpkm_coverage.csv", header = T)
s87 <- s87 %>%
  separate(Contig.s87_gene_filter.RPKM, into = c("name", "87"), sep = "\t")
rownames(s87) <- s87[[1]]
s087 <- as.matrix(s87[,-1])
rownames(s087) <- rownames(s87)
colnames(s087) <- "s87"
write.csv(s087, file = "s87.rpkm.csv")

s88 <- read.csv("s88_gene_rpkm_coverage.csv", header = T)
s88 <- s88 %>%
  separate(Contig.s88_gene_filter.RPKM, into = c("name", "88"), sep = "\t")
rownames(s88) <- s88[[1]]
s088 <- as.matrix(s88[,-1])
rownames(s088) <- rownames(s88)
colnames(s088) <- "s88"
write.csv(s088, file = "s88.rpkm.csv")

s89 <- read.csv("s89_gene_rpkm_coverage.csv", header = T)
s89 <- s89 %>%
  separate(Contig.s89_gene_filter.RPKM, into = c("name", "89"), sep = "\t")
rownames(s89) <- s89[[1]]
s089 <- as.matrix(s89[,-1])
rownames(s089) <- rownames(s89)
colnames(s089) <- "s89"
write.csv(s089, file = "s89.rpkm.csv")

#################
#s01 <- read.csv("s01.rpkm.csv", header = T, row.names = 1)
#s02 <- read.csv("s02.rpkm.csv", header = T, row.names = 1)
#s03 <- read.csv("s03.rpkm.csv", header = T, row.names = 1)
#s04 <- read.csv("s04.rpkm.csv", header = T, row.names = 1)
#s05 <- read.csv("s05.rpkm.csv", header = T, row.names = 1)
#s06 <- read.csv("s06.rpkm.csv", header = T, row.names = 1)
#s07 <- read.csv("s07.rpkm.csv", header = T, row.names = 1)
#s08 <- read.csv("s08.rpkm.csv", header = T, row.names = 1)
#s09 <- read.csv("s09.rpkm.csv", header = T, row.names = 1)
#s10 <- read.csv("s10.rpkm.csv", header = T, row.names = 1)
#s11 <- read.csv("s11.rpkm.csv", header = T, row.names = 1)
#s12 <- read.csv("s12.rpkm.csv", header = T, row.names = 1)
#s13 <- read.csv("s13.rpkm.csv", header = T, row.names = 1)
#s14 <- read.csv("s14.rpkm.csv", header = T, row.names = 1)
#s15 <- read.csv("s15.rpkm.csv", header = T, row.names = 1)

spe <- cbind(s001,s002,s003,s004,s005,s006,s007,s008,s009,s010,s011,s012,s013,s014,s015,
               s016,s017,s018,s019,s020,s021,s022,s023,s024,s025,s026,s027,s028,s029,s030,
               s031,s032,s033,s034,s035,s036,s037,s038,s039,s040,s041,s042,s043,s044,s045,
               s046,s047,s048,s049,s050,s051,s052,s053,s054,s055,s056,s057,s058,s059,s060,
               s061,s062,s063,s064,s065,s066,s067,s068,s069,s070,s071,s072,s073,s074,s075,
               s076,s077,s078,s087,s088,s087)

library(vegan)
spe.db <- vegdist(t(spe), "bray")
write.csv(as.matrix(spe.db), file="Bray81.csv")








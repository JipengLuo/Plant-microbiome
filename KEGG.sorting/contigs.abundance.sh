#!/bin/bash
#SBATCH -J gene_abund
#SBATCH -p general
#SBATCH -o gene_abund%j.txt
#SBATCH -e gene_abund%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luojip@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --mem=250G
#SBATCH -A r00324

#从ORFs或去冗余的氨基酸序列中提取代表性（核苷酸/氨基酸）序列，生成fasta文件
export PATH=$PATH:/N/project/luojp_CNH2/miniconda3/bin

seqkit seq -n ./dereplication/clusterRes_rep_seq.fasta > ./dereplication/gene_rep.txt
sed -i ' s/ //g' ./dereplication/gene_rep.txt
seqkit grep -f ./dereplication/gene_rep.txt nucl_rename.fnn > ./dereplication/gene_rep.fa

#mapping前代表性核苷酸序列建立索引
export PATH=$PATH:/N/project/luojp_CNH2/applications/bwa
export PATH=$PATH:/N/project/luojp_CNH2/applications/samtools-1.11
export PATH=$PATH:/geode2/home/u050/luojip/Quartz/applications/coverm
#mapping前代表性核苷酸序列建立索引
bwa index -p ./index/gene ./dereplication/gene_rep.fa

#mapping生成bam文件
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s01_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s01_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s01_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s02_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s02_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s02_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s03_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s03_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s03_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s04_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s04_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s04_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s05_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s05_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s05_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s06_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s06_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s06_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s07_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s07_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s07_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s08_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s08_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s08_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s09_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s09_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s09_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s10_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s10_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s10_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s11_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s11_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s11_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s12_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s12_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s12_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s13_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s13_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s13_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s14_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s14_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s14_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s15_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s15_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s15_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s16_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s16_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s16_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s17_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s17_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s17_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s18_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s18_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s18_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s19_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s19_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s19_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s20_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s20_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s20_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s21_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s21_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s21_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s22_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s22_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s22_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s23_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s23_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s23_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s24_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s24_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s24_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s25_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s25_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s25_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s26_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s26_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s26_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s27_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s27_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s27_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s28_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s28_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s28_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s29_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s29_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s29_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s30_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s30_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s30_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s31_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s31_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s31_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s32_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s32_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s32_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s33_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s33_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s33_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s34_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s34_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s34_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s35_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s35_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s35_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s36_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s36_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s36_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s37_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s37_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s37_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s38_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s38_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s38_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s39_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s49_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s39_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s40_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s40_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s40_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s41_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s41_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s41_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s42_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s42_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s42_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s43_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s43_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s43_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s44_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s44_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s44_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s45_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s45_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s45_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s46_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s46_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s46_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s47_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s47_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s47_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s48_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s48_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s48_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s49_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s49_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s49_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s50_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s50_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s50_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s51_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s51_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s51_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s52_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s52_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s52_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s53_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s53_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s53_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s54_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s54_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s54_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s55_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s55_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s55_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s56_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s56_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s56_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s57_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s57_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s57_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s58_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s58_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s58_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s59_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s59_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s59_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s60_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s60_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s60_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s61_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s61_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s61_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s62_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s62_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s62_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s63_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s63_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s63_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s64_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s64_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s64_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s65_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s65_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s65_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s66_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s66_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s66_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s67_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s67_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s67_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s68_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s68_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s68_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s69_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s69_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s69_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s70_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s70_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s70_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s71_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s71_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s71_contigs.bam
bwa mem -t 24 ./index/gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/s72_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/s72_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > ./bam.file/s72_contigs.bam

#选择合适的对齐阈值筛选bam文件,过滤对齐的一致性较低的reads
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s01_contigs.bam -o ./bam.file/s01_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s02_contigs.bam -o ./bam.file/s02_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s03_contigs.bam -o ./bam.file/s03_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s04_contigs.bam -o ./bam.file/s04_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s05_contigs.bam -o ./bam.file/s05_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s06_contigs.bam -o ./bam.file/s06_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s07_contigs.bam -o ./bam.file/s07_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s08_contigs.bam -o ./bam.file/s08_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s09_contigs.bam -o ./bam.file/s09_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/120_contigs.bam -o ./bam.file/s10_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s11_contigs.bam -o ./bam.file/s11_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s12_contigs.bam -o ./bam.file/s12_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s13_contigs.bam -o ./bam.file/s13_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s14_contigs.bam -o ./bam.file/s14_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s15_contigs.bam -o ./bam.file/s15_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s16_contigs.bam -o ./bam.file/s16_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s17_contigs.bam -o ./bam.file/s17_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s18_contigs.bam -o ./bam.file/s18_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s19_contigs.bam -o ./bam.file/s19_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s20_contigs.bam -o ./bam.file/s20_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s21_contigs.bam -o ./bam.file/s21_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s22_contigs.bam -o ./bam.file/s22_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s23_contigs.bam -o ./bam.file/s23_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s24_contigs.bam -o ./bam.file/s24_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s25_contigs.bam -o ./bam.file/s25_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s26_contigs.bam -o ./bam.file/s26_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s27_contigs.bam -o ./bam.file/s27_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s28_contigs.bam -o ./bam.file/s28_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s29_contigs.bam -o ./bam.file/s29_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s30_contigs.bam -o ./bam.file/s30_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s31_contigs.bam -o ./bam.file/s31_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s32_contigs.bam -o ./bam.file/s32_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s33_contigs.bam -o ./bam.file/s33_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s34_contigs.bam -o ./bam.file/s34_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s35_contigs.bam -o ./bam.file/s35_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s36_contigs.bam -o ./bam.file/s36_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s37_contigs.bam -o ./bam.file/s37_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s38_contigs.bam -o ./bam.file/s38_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s39_contigs.bam -o ./bam.file/s39_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s40_contigs.bam -o ./bam.file/s40_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s41_contigs.bam -o ./bam.file/s41_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s42_contigs.bam -o ./bam.file/s42_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s43_contigs.bam -o ./bam.file/s43_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s44_contigs.bam -o ./bam.file/s44_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s45_contigs.bam -o ./bam.file/s45_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s46_contigs.bam -o ./bam.file/s46_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s47_contigs.bam -o ./bam.file/s47_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s48_contigs.bam -o ./bam.file/s48_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s49_contigs.bam -o ./bam.file/s49_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s50_contigs.bam -o ./bam.file/s50_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s51_contigs.bam -o ./bam.file/s51_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s52_contigs.bam -o ./bam.file/s52_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s53_contigs.bam -o ./bam.file/s53_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s54_contigs.bam -o ./bam.file/s54_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s55_contigs.bam -o ./bam.file/s55_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s56_contigs.bam -o ./bam.file/s56_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s57_contigs.bam -o ./bam.file/s57_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s58_contigs.bam -o ./bam.file/s58_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s59_contigs.bam -o ./bam.file/s59_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s60_contigs.bam -o ./bam.file/s60_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s61_contigs.bam -o ./bam.file/s61_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s62_contigs.bam -o ./bam.file/s62_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s63_contigs.bam -o ./bam.file/s63_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s64_contigs.bam -o ./bam.file/s64_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s65_contigs.bam -o ./bam.file/s65_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s66_contigs.bam -o ./bam.file/s66_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s67_contigs.bam -o ./bam.file/s67_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s68_contigs.bam -o ./bam.file/s68_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s69_contigs.bam -o ./bam.file/s69_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s70_contigs.bam -o ./bam.file/s70_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s71_contigs.bam -o ./bam.file/s71_contigs_filter.bam -t 24
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/s72_contigs.bam -o ./bam.file/s72_contigs_filter.bam -t 24

#选择计算覆盖的方法，输出覆盖率结果文件
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s01_contigs_filter.bam > ./abundance/s01_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s02_contigs_filter.bam > ./abundance/s02_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s03_contigs_filter.bam > ./abundance/s03_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s04_contigs_filter.bam > ./abundance/s04_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s05_contigs_filter.bam > ./abundance/s05_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s06_contigs_filter.bam > ./abundance/s06_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s07_contigs_filter.bam > ./abundance/s07_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s08_contigs_filter.bam > ./abundance/s08_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s09_contigs_filter.bam > ./abundance/s09_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s10_contigs_filter.bam > ./abundance/s10_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s11_contigs_filter.bam > ./abundance/s11_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s12_contigs_filter.bam > ./abundance/s12_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s13_contigs_filter.bam > ./abundance/s13_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s14_contigs_filter.bam > ./abundance/s14_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s15_contigs_filter.bam > ./abundance/s15_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s16_contigs_filter.bam > ./abundance/s16_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s17_contigs_filter.bam > ./abundance/s17_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s18_contigs_filter.bam > ./abundance/s18_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s19_contigs_filter.bam > ./abundance/s19_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s20_contigs_filter.bam > ./abundance/s20_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s21_contigs_filter.bam > ./abundance/s21_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s22_contigs_filter.bam > ./abundance/s22_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s23_contigs_filter.bam > ./abundance/s23_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s24_contigs_filter.bam > ./abundance/s24_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s25_contigs_filter.bam > ./abundance/s25_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s26_contigs_filter.bam > ./abundance/s26_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s27_contigs_filter.bam > ./abundance/s27_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s28_contigs_filter.bam > ./abundance/s28_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s29_contigs_filter.bam > ./abundance/s29_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s30_contigs_filter.bam > ./abundance/s30_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s31_contigs_filter.bam > ./abundance/s31_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s32_contigs_filter.bam > ./abundance/s32_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s33_contigs_filter.bam > ./abundance/s33_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s34_contigs_filter.bam > ./abundance/s34_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s35_contigs_filter.bam > ./abundance/s35_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s36_contigs_filter.bam > ./abundance/s36_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s37_contigs_filter.bam > ./abundance/s37_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s38_contigs_filter.bam > ./abundance/s38_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s39_contigs_filter.bam > ./abundance/s39_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s40_contigs_filter.bam > ./abundance/s40_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s41_contigs_filter.bam > ./abundance/s41_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s42_contigs_filter.bam > ./abundance/s42_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s43_contigs_filter.bam > ./abundance/s43_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s44_contigs_filter.bam > ./abundance/s44_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s45_contigs_filter.bam > ./abundance/s45_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s46_contigs_filter.bam > ./abundance/s46_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s47_contigs_filter.bam > ./abundance/s47_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s48_contigs_filter.bam > ./abundance/s48_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s49_contigs_filter.bam > ./abundance/s49_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s50_contigs_filter.bam > ./abundance/s50_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s51_contigs_filter.bam > ./abundance/s51_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s52_contigs_filter.bam > ./abundance/s52_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s53_contigs_filter.bam > ./abundance/s53_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s54_contigs_filter.bam > ./abundance/s54_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s55_contigs_filter.bam > ./abundance/s55_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s56_contigs_filter.bam > ./abundance/s56_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s57_contigs_filter.bam > ./abundance/s57_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s58_contigs_filter.bam > ./abundance/s58_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s59_contigs_filter.bam > ./abundance/s59_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s60_contigs_filter.bam > ./abundance/s60_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s61_contigs_filter.bam > ./abundance/s61_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s62_contigs_filter.bam > ./abundance/s62_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s63_contigs_filter.bam > ./abundance/s63_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s64_contigs_filter.bam > ./abundance/s64_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s65_contigs_filter.bam > ./abundance/s65_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s66_contigs_filter.bam > ./abundance/s66_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s67_contigs_filter.bam > ./abundance/s67_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s68_contigs_filter.bam > ./abundance/s68_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s69_contigs_filter.bam > ./abundance/s69_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s70_contigs_filter.bam > ./abundance/s70_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s71_contigs_filter.bam > ./abundance/s71_contigs_rpkm_coverage.csv
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/s72_contigs_filter.bam > ./abundance/s72_contigs_rpkm_coverage.csv










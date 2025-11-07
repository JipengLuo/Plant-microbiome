#!/bin/bash
#SBATCH -J bwa.contig
#SBATCH -p general
#SBATCH -o bwa.contig%j.txt
#SBATCH -e bwa.contig%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luojip@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=96:00:00
#SBATCH --mem=200G
#SBATCH -A r00324

#mapping前代表性核苷酸序列建立索引
export PATH=$PATH:/N/project/luojp_CNH2/applications/bwa
export PATH=$PATH:/N/project/luojp_CNH2/applications/samtools-1.11
export PATH=$PATH:/geode2/home/u050/luojip/Quartz/applications/coverm
#该步骤会产生五个文件,因为已经生成了该五个文件，因此后续分析不在跑该步了；
#bwa index ./index/s01.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s01/final.contigs.fa
#bwa index ./index/s02.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s02/final.contigs.fa
#bwa index ./index/s03.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s03/final.contigs.fa
#bwa index ./index/s04.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s04/final.contigs.fa
#bwa index ./index/s05.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s05/final.contigs.fa
#bwa index ./index/s06.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s06/final.contigs.fa
#bwa index ./index/s07.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s07/final.contigs.fa
#bwa index ./index/s08.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s08/final.contigs.fa
#bwa index ./index/s09.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s09/final.contigs.fa
#bwa index ./index/s10.gene /N/project/luojp_CNH2/metagenomic02/03.assemble/s10/final.contigs.fa
#将reads文件与contigs库比对
bwa mem -t 24 ./index/all.gene /N/project/luojp_CNH2/metagenomic02/02.cleandata/*_paired_R1.clean.fq.gz /N/project/luojp_CNH2/metagenomic02/02.cleandata/*_paired_R2.clean.fq.gz|samtools sort -O bam -@ 24 -o - > .bam.file/all_contigs.bam

#选择合适的对齐阈值筛选bam文件,过滤对齐的一致性较低的reads
coverm filter --min-read-percent-identity 0.95 --min-read-aligned-percent 0.75 -b ./bam.file/all_contigs.bam -o ./bam.file/all_contigs_filter.bam -t 24

#选择计算覆盖的方法，输出覆盖率结果文件
coverm contig --methods rpkm --trim-max 90 --trim-min 10 --bam-files ./bam.file/all_contigs_filter.bam > ./abundance/all_contigs_rpkm_coverage.csv

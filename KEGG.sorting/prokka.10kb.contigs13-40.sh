#!/bin/bash
#SBATCH -J prokka
#SBATCH -p general
#SBATCH -o prokka%j.txt
#SBATCH -e prokka%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luojip@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=07:00:00
#SBATCH --mem=200G
#SBATCH -A r00324

#使用coverm对比对上的reads进行过滤
export PATH=$PATH:/geode2/home/u050/luojip/Quartz/miniconda3/envs/metagenome/bin
prokka ./metagenomic02/03.assemble/s13/final.contigs.fa --outdir ./metagenomic02/s13 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s14/final.contigs.fa --outdir ./metagenomic02/s14 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s15/final.contigs.fa --outdir ./metagenomic02/s15 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s16/final.contigs.fa --outdir ./metagenomic02/s16 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s17/final.contigs.fa --outdir ./metagenomic02/s17 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s18/final.contigs.fa --outdir ./metagenomic02/s18 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s19/final.contigs.fa --outdir ./metagenomic02/s19 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s20/final.contigs.fa --outdir ./metagenomic02/s20 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s21/final.contigs.fa --outdir ./metagenomic02/s21 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s22/final.contigs.fa --outdir ./metagenomic02/s22 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s23/final.contigs.fa --outdir ./metagenomic02/s23 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s24/final.contigs.fa --outdir ./metagenomic02/s24 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s25/final.contigs.fa --outdir ./metagenomic02/s25 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s26/final.contigs.fa --outdir ./metagenomic02/s26 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s27/final.contigs.fa --outdir ./metagenomic02/s27 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s28/final.contigs.fa --outdir ./metagenomic02/s28 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s29/final.contigs.fa --outdir ./metagenomic02/s29 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s30/final.contigs.fa --outdir ./metagenomic02/s30 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s31/final.contigs.fa --outdir ./metagenomic02/s31 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s32/final.contigs.fa --outdir ./metagenomic02/s32 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s33/final.contigs.fa --outdir ./metagenomic02/s33 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s34/final.contigs.fa --outdir ./metagenomic02/s34 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s35/final.contigs.fa --outdir ./metagenomic02/s35 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s36/final.contigs.fa --outdir ./metagenomic02/s36 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s37/final.contigs.fa --outdir ./metagenomic02/s37 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s38/final.contigs.fa --outdir ./metagenomic02/s38 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s39/final.contigs.fa --outdir ./metagenomic02/s39 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria
prokka ./metagenomic02/03.assemble/s40/final.contigs.fa --outdir ./metagenomic02/s40 --prefix metag --addgenes --addgenes --metagenome --kingdom Bacteria





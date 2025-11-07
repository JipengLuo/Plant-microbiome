#!/bin/bash
#SBATCH -J mmseqs
#SBATCH -p general
#SBATCH -o mmseqs%j.txt
#SBATCH -e mmseqs%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luojip@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --mem=250G
#SBATCH -A r00324

#删除第一行基因名称空格后的字符
#seqkit replace -p " .+" -i nucl.fnn > nucl_rename.fnn
#seqkit replace -p " .+" -i amino.faa > amino_rename.faa

mmseqs easy-linclust amino_rename.faa ./dereplication/clusterRes ./dereplication/tmp --min-seq-id 0.95 -c 0.9 --threads 24 --cluster-mode 2 --cov-mode 1
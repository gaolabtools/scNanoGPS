#! /usr/bin/bash

# https://kb.10xgenomics.com/hc/en-us/articles/4412343032205-Where-can-I-find-the-barcode-whitelist-s-for-Single-Cell-Multiome-ATAC-GEX-product-

#  The ARC Multiome Gene Expression whitelist can be found here:<path_to_cellrangerarc>/cellranger-arc-x.y.z/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz 

cd ~/data/benchmarking/A375/A375_by_blaze

cp ~/tools/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz .

gunzip 737K-arc-v1.txt.gz

python3 ~/tools/BLAZE/bin/blaze.py --full-bc-whitelist 737K-arc-v1.txt --expect-cells 3000 --threads 30 /data/nanopore/2021/NP19_PROM0102_Gao_Kieser_nu_A375_12102021/20211211_0102_1-E7-H7_PAI18615_c649f3e1/fastq_pass/

cp whitelist.csv ~/data/benchmarking/A375_blaze.barcode_list.csv


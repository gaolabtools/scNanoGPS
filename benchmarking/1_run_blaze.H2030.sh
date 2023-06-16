#! /usr/bin/bash

# https://kb.10xgenomics.com/hc/en-us/articles/4412343032205-Where-can-I-find-the-barcode-whitelist-s-for-Single-Cell-Multiome-ATAC-GEX-product-

#  The ARC Multiome Gene Expression whitelist can be found here:<path_to_cellrangerarc>/cellranger-arc-x.y.z/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz 

cd ~/data/benchmarking/H2030/H2030_by_blaze

cp ~/tools/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz .

gunzip 737K-arc-v1.txt.gz

python3 ~/tools/BLAZE/bin/blaze.py --full-bc-whitelist 737K-arc-v1.txt --expect-cells 3000 --threads 30 /data/nanopore/2021/NP18_PROM0102_Gao_Kieser_nu_H2030_12012021/20211202_0228_2-E9-H9_PAI19630_2bbb4022/fastq_pass/

cp whitelist.csv ~/data/benchmarking/H2030_blaze.barcode_list.csv


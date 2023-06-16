#! /usr/bin/bash

snakemake --use-conda --configfile sockeye_config_A375/config.yml --cores 30 -pr all

cp ~/data/benchmarking/A375/A375_by_sockeye/A375/demux/whitelist.tsv ~/data/benchmarking/A375_sockeye.barcode_list.csv


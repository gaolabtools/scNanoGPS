#! /usr/bin/bash

snakemake --use-conda --configfile sockeye_config_H2030/config.yml --cores 30 -pr all

cp ~/data/benchmarking/H2030/H2030_by_sockeye/H2030/demux/whitelist.tsv ~/data/benchmarking/H2030_sockeye.barcode_list.csv


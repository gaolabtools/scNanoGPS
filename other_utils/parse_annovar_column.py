#! /usr/bin/env python3

import time, os, sys, gzip, re, subprocess
import numpy as np
import pandas as pd
from optparse import OptionParser

#===get params===
parser = OptionParser()
parser.add_option("-i", dest = "f_name", nargs = 1, default = "annovar.vcf",
                  help = "Annovar result. "
                         "Default: annovar.vcf")
options, arguments = parser.parse_args()

#===pre-check===
termination = False
if not os.path.isfile(options.f_name):
	print("\nCannot find file at: " + str(options.f_name))
	termination = True

if termination:
	parser.print_help()
	sys.exit(1)

#===start===
print("CHROM\tPOS\tREF\tALT\tGene_name\tRefSeq_ID\tExact_mutation\tFunc_type\tMutation_type\tdbSNP_ID\tCosmic_coding_ID\tCosmic_coding_tumors\tCosmic_noncoding_ID\tCosmic_noncoding_tumors")
df = pd.read_csv(options.f_name, header=None, sep='\t', comment="#")
for idx in range(0, df.shape[0]):
	snp_list = df.iloc[idx, [0, 1, 3, 4]].values
	Gene_name, Refseq_ID, exact_mutation, loc_type, syn_ns_type, snp_id, c_cosmic_id, c_cosmic_tumor, n_cosmic_id, n_cosmic_tumor = "", "", "", "", "", "", "", "", "", ""
	annovar_list = df.iloc[idx, [7]].values[0].split(';')
#	Gene.refGene=AGRN
#	refseq ID:	AAChange.refGene=ISG15:  NM_005101  :exon2:c.A294G:p.V98V
#	exact mutation:	AAChange.refGene=ISG15:NM_005101:exon2:  c.A294G  :p.V98V
#	Func.refGene=intronic
#	ExonicFunc.refGene=nonsynonymous_SNV
#	avsnp150=rs8997
#	cosmic96_coding=ID\x3dCOSV65106735\x3bOCCURENCE\x3d1(kidney);cosmic96_noncoding=ID\x3dCOSV65106735\x3bOCCURENCE\x3d1(kidney)
	for annovar_cont in annovar_list:
		anno_cont = annovar_cont.split('=')
		if anno_cont[0] == "Gene.refGene":
			Gene_name = anno_cont[1]
			Gene_name = Gene_name.replace('\\x3b', ';')
		if anno_cont[0] == "AAChange.refGene":
			tmp_list = anno_cont[1].split(':')
			if len(tmp_list) > 3:
				tmp_list = anno_cont[1].split(':')
				Refseq_ID      = tmp_list[1]
				exact_mutation = tmp_list[3]
			else:
				Refseq_ID = anno_cont[1]
		if anno_cont[0] == "Func.refGene":
			loc_type = anno_cont[1]
		if anno_cont[0] == "ExonicFunc.refGene":
			syn_ns_type = anno_cont[1]
		if anno_cont[0] == "avsnp150":
			snp_id = anno_cont[1]
		if anno_cont[0] == "cosmic96_coding":
			if anno_cont[1] != '.':
				cosmic_list = anno_cont[1].replace('\\x3d', '=').replace('\\x3b', ';').split(';')
				c_cosmic_id = cosmic_list[0].split('=')[1]
				occurrence_list = re.split('\(|\)', cosmic_list[1].split('=')[1])
				c_cosmic_tumor = ';'.join([str(occurrence_list[idx]) for idx in range(1, len(occurrence_list), 2)])
			else:
				c_cosmic_id = '.'
				c_cosmic_tumor = '.'
		if anno_cont[0] == "cosmic96_noncoding":
			if anno_cont[1] != '.':
				cosmic_list = anno_cont[1].replace('\\x3d', '=').replace('\\x3b', ';').split(';')
				n_cosmic_id = cosmic_list[0].split('=')[1]
				occurrence_list = re.split('\(|\)', cosmic_list[1].split('=')[1])
				n_cosmic_tumor = ';'.join([str(occurrence_list[idx]) for idx in range(1, len(occurrence_list), 2)])
			else:
				n_cosmic_id = '.'
				n_cosmic_tumor = '.'

	print("\t".join([str(s) for s in snp_list]), end = "\t")
	print(Gene_name, end = "\t")
	print(Refseq_ID, end = "\t")
	print(exact_mutation, end = "\t")
	print(loc_type, end = "\t")
	print(syn_ns_type, end = "\t")
	print(snp_id, end = "\t")
	print(c_cosmic_id, end = "\t")
	print(c_cosmic_tumor, end = "\t")
	print(n_cosmic_id, end = "\t")
	print(n_cosmic_tumor, end = "\n")

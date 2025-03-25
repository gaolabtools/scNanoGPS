#! /usr/bin/env python3

import time, os, sys, glob, gzip, pysam, subprocess, math
import pandas as pd
import multiprocessing as mp
import scipy.stats as stats

sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../'))

from Bio import SeqIO, pairwise2
from functools import partial
from contextlib import contextmanager
from optparse import OptionParser
from itertools import chain
from gtfparse import read_gtf
from scanner_core import misc
from curator_core import curator_io

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

def add_lv(res_dict, ele):
	if not ele in res_dict.keys():
		res_dict.update({ele: dict()})
	return res_dict

def fusion_isoform_finder(CB, options):
	print("Processing " + CB + "......")
	start_time = time.time()

	blast_res_col = ['qseqid', 'qstart', 'qend', 'qlen', 'sseqid', 'sstart', 'send', 'slen', 'bitscore', 'evalue', 'pident']
	oh = open(os.path.join(options.tmp_dir, CB) + options.blast_o, "wt")
	oh.write("")
	oh.close()

	if os.stat(os.path.join(options.tmp_dir, CB) + options.read_o).st_size > 0:

		cmd = options.blastn + \
		      " -query " + os.path.join(options.tmp_dir, CB) + options.read_o + \
		      " -task blastn -db " + os.path.join(options.tmp_dir, options.db_name) + \
		      " -out " + os.path.join(options.tmp_dir, CB) + options.blast_o + \
		      " -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen bitscore evalue pident'"
		os.system(cmd)

	blast_res_df = pd.read_csv(os.path.join(options.tmp_dir, CB) + options.blast_o, names = blast_res_col, header = None, sep = '\t')

	if blast_res_df.shape[0] > 0:
		#===carve out gid===
		blast_res_df['query_gid_1'] = blast_res_df.apply(lambda row: row['qseqid'].split(';')[2].split(','), axis = 1)
		blast_res_df['query_gid_2'] = blast_res_df.apply(lambda row: row['qseqid'].split(';')[4].split(','), axis = 1)
		blast_res_df['target_gid']  = blast_res_df.apply(lambda row: options.tx_dict[row['sseqid']] if row['sseqid'] in options.tx_dict else "", axis = 1)

		#===data selection===
		blast_res_df = blast_res_df.loc[blast_res_df['target_gid'] != "", :]
		blast_res_df = blast_res_df.loc[blast_res_df.apply(lambda x: (x['target_gid'] in x['query_gid_1']) | (x['target_gid'] in x['query_gid_2']), axis = 1), :]
		blast_res_df['coverage'] = (blast_res_df['qend'] - blast_res_df['qstart']) / blast_res_df['qlen'] * blast_res_df['pident'] / 100

		#===filter out gid < 2===
		blast_res_df_gid_counting = blast_res_df.loc[:, ['qseqid', 'target_gid']].groupby('qseqid').agg({"target_gid": pd.Series.nunique})
		blast_res_df_gid_counting = blast_res_df_gid_counting[blast_res_df_gid_counting['target_gid'] >= 2]
		qseqid_list = blast_res_df_gid_counting.index.unique().tolist()
		blast_res_df = blast_res_df.loc[blast_res_df['qseqid'].isin(qseqid_list), :]

		#===select reads===
		blast_res_df = blast_res_df.sort_values(by=['qseqid', 'target_gid', 'coverage'], ascending = False).drop_duplicates('target_gid')
		res_dict = dict()
		for row in blast_res_df.iterrows():
			for gid1 in row[1]['query_gid_1']:
				res_dict = add_lv(res_dict, gid1)
				for gid2 in row[1]['query_gid_2']:
					res_dict[gid1] = add_lv(res_dict[gid1], gid2)
					if row[1]['target_gid'] == gid1:
						res_dict[gid1][gid2] = add_lv(res_dict[gid1][gid2], "tx1")
						res_dict[gid1][gid2]['tx1'].update({row[1]['sseqid']: row[1]['qseqid'].split(';')[0] + ',' + str(row[1]['sstart']) + ',' + str(row[1]['send']) + ',' + str(round(row[1]['coverage'], 2))})
					if row[1]['target_gid'] == gid2:
						res_dict[gid1][gid2] = add_lv(res_dict[gid1][gid2], "tx2")
						res_dict[gid1][gid2]['tx2'].update({row[1]['sseqid']: row[1]['qseqid'].split(';')[0] + ',' + str(row[1]['sstart']) + ',' + str(row[1]['send']) + ',' + str(round(row[1]['coverage'], 2))})

	res_col = ['gene_id_1', 'gene_id_2', 'gene_name_1', 'gene_name_2', 'tx_id_1', 'tx_id_2', CB]
	res_m = list()
	if blast_res_df.shape[0] > 0:
		for gid1 in res_dict.keys():
			gn1 = options.gn_dict[gid1.strip('\'')]
			for gid2 in res_dict[gid1].keys():
				gn2 = options.gn_dict[gid2.strip('\'')]
				if gid1 == gid2:
					continue
				if (not 'tx1' in res_dict[gid1][gid2].keys()) | (not 'tx2' in res_dict[gid1][gid2].keys()):
					continue
				read_list = []
				for tx1 in res_dict[gid1][gid2]['tx1'].keys():
					read_list.append(tx1 + ',' + res_dict[gid1][gid2]['tx1'][tx1])
				for tx2 in res_dict[gid1][gid2]["tx2"].keys():
					read_list.append(tx2 + ',' + res_dict[gid1][gid2]["tx2"][tx2])
				res_m.append([gid1.strip('\''), gid2.strip('\''), gn1, gn2, tx1.strip('\''), tx2.strip('\''), ';'.join(read_list)])
	res_df = pd.DataFrame(res_m, columns = res_col)

	hours, minutes, seconds = misc.get_time_elapse(start_time)
	print("Fusion isoform finder spent %d : %d : %.2f" % (hours, minutes, seconds))

	return res_df

if __name__ == "__main__":

	#===get params===
	parser = OptionParser()

	parser.add_option("-d",         dest = "o_dir",     default = "scNanoGPS_res",
                          nargs = 1, type = "string",
                          help = "Output directory name. "
                                 "Default: scNanoGPS_res")
	parser.add_option("--tmp_dir",  dest = "tmp_dir",   default = "tmp",
                          nargs = 1, type = "string",
                          help = "Temporary folder name. "
                                 "Default: tmp")
	parser.add_option("-i",         dest = "i_name",    default = "matrix_fusion.tsv.gz",
                          nargs = 1, type = "string",
                          help = "Counting table name. "
                                 "Default: matrix_fusion.tsv.gz")
	parser.add_option("-o",         dest = "o_name",    default = "matrix_fusion_isoform.tsv.gz",
                          nargs = 1, type = "string",
                          help = "Table of major fusion isoform. "
                                 "Default: matrix_fusion_isoform.tsv.gz")
	parser.add_option("-t",         dest = "ncores",    default = 1,
                          nargs = 1, type = "int",
                          help = "Number of cores for program running. "
                                 "Default: 1")
	parser.add_option("--ref_genome", dest = "ref_genome", default = None,
                          nargs = 1, type = "string",
                          help = "* Required ! "
                                 "File for reference genome.")
	parser.add_option("--gtf",      dest = "gtf",       default = None,
                          nargs = 1, type = "string",
                          help = "* Required ! "
                                 "GTF file for expression calling. ")
	parser.add_option("--read_o",   dest = "read_o",    default = ".fusion.fasta",
                          nargs = 1, type = "string",
                          help = "FastA output of fusion reads. "
                                 "Default: .fusion.fasta")
	parser.add_option("--blast_o",   dest = "blast_o",    default = ".blast.tsv",
                          nargs = 1, type = "string",
                          help = "Blast results of fusion reads. "
                                 "Default: .blast.tsv")
	parser.add_option("--tx_sel",   dest = "tx_sel",    default = "select_transcript_ids.txt",
                          nargs = 1, type = "string",
                          help = "IDs list file of potential transcripts. "
                                 "Default: select_transcript_ids.txt")
	parser.add_option("--tx_fa",    dest = "tx_fa",     default = "select_transcript.fa",
                          nargs = 1, type = "string",
                          help = "FastA file of potential transcripts. "
                                 "Default: select_transcript.fa")
	parser.add_option("--db_name",  dest = "db_name",   default = "fusion_isoform_tx",
                          nargs = 1, type = "string",
                          help = "BLAST database name for transcripts. "
                                 "Default: fusion_isoform_tx")
	parser.add_option("--gffread",  dest = "gffread",   default = "gffread",
                          nargs = 1, type = "string",
                          help = "Path to gffread. "
                                 "Default: gffread")
	parser.add_option("--blastn",   dest = "blastn",    default = "blastn",
                          nargs = 1, type = "string",
                          help = "Path to blastn. "
                                 "Default: blastn")
	parser.add_option("--makeblastdb", dest = "makeblastdb",  default = "makeblastdb",
                          nargs = 1, type = "string",
                          help = "Path to makeblastdb. "
                                 "Default: makeblastdb")

	options, arguments = parser.parse_args()

	#===pre-check===
	termination = False
	if not os.path.isdir(options.o_dir):
		print("\nOutput directory is not exist: "    + options.o_dir    + "\n")
		termination = True
	if not options.ref_genome or \
	   not os.path.isfile(options.ref_genome):
		print("\nCannot find reference genome file: " + str(options.ref_genome) + "\n")
		termination = True
	if not options.gtf or \
	   not os.path.isfile(options.gtf):
		print("\nCannot find GTF file: " + str(options.gtf) + "\n")
		termination = True

	cmd = "which " + options.gffread
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	options.gffread = out_msg.decode("utf-8").rstrip()
	cmd = "which " + options.blastn
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	options.blastn = out_msg.decode("utf-8").rstrip()
	cmd = "which " + options.makeblastdb
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	options.makeblastdb = out_msg.decode("utf-8").rstrip()

	if not os.path.isfile(options.gffread):
		print("\nCannot find gffread: " + options.gffread)
		termination = True
	if not os.path.isfile(options.blastn):
		print("\nCannot find blastn: " + options.blastn)
		termination = True
	if not os.path.isfile(options.makeblastdb):
		print("\nCannot find makeblastdb: " + options.makeblastdb)
		termination = True

	if termination:
		parser.print_help()
		sys.exit(1)

	options.tx_sel = os.path.join(options.tmp_dir, options.tx_sel)
	options.tx_fa  = os.path.join(options.tmp_dir, options.tx_fa)

	#=== set env variables
	os.environ["OMP_NUM_THREADS"]        = str(options.ncores)
	os.environ["OPENBLAS_NUM_THREADS"]   = str(options.ncores)
	os.environ["MKL_NUM_THREADS"]        = str(options.ncores)
	os.environ["VECLIB_MAXIMUM_THREADS"] = str(options.ncores)
	os.environ["NUMEXPR_NUM_THREADS"]    = str(options.ncores)

	import numpy as np

	start_time = time.time()
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)

	#===loading annotation===
	print("Loading GTF...", flush = True)
	anno_df = read_gtf(options.gtf)
	tx_df = anno_df.loc[anno_df["feature"] == "transcript", ['gene_id', 'transcript_id']]
	gn_df = anno_df.loc[anno_df["feature"] == "gene", ['gene_id', 'gene_name']]
	gn_df.index = gn_df['gene_id']
	options.gn_dict = gn_df['gene_name'].to_dict()

	#===remove no isoform genes===
	tx_counting_df = tx_df.groupby('gene_id').count()
	tx_counting_df = tx_counting_df[tx_counting_df['transcript_id'] > 1]
	tx_df = tx_df[tx_df['gene_id'].isin(tx_counting_df.index)]

	#===generate isoform gene dictionary===
	options.tx_dict = dict()
	tx_df.apply(lambda row: options.tx_dict.update({row['transcript_id']: row['gene_id']}), axis = 1)
	print("Done", flush = True)

	#===load fusion table===
	compression = None
	if options.i_name.endswith(".gz"):
		compression = 'gzip'
	fus_df   = pd.read_csv(os.path.join(options.o_dir, options.i_name), header = 0, sep = '\t', compression = compression)
	key_list = ['gene1', 'gene2', 'gene_name_1', 'gene_name_2', 'bp1', 'bp2']
	CB_list  = [x for x in fus_df.columns.to_list() if x not in key_list]

	#===fusion gene list===
	fus_gene_list = np.unique(list(chain.from_iterable(fus_df['gene1'].str.split(',').to_list() + fus_df['gene2'].str.split(',').to_list())))
	fus_tx_list   = list(np.unique(tx_df.loc[tx_df['gene_id'].isin(fus_gene_list), 'transcript_id'].to_list()))
	pd.DataFrame(fus_tx_list).to_csv(options.tx_sel, sep = '\t', header = False, index = False)

	#===generate selected transcript db===
	cmd = options.gffread + " -w " + options.tx_fa + " -g " + options.ref_genome + " --ids " + options.tx_sel + " --gtf " + options.gtf
	os.system(cmd)
	cmd = options.makeblastdb + " -dbtype nucl -in " + options.tx_fa + " -input_type fasta -title " + options.db_name + " -out " + os.path.join(options.tmp_dir, options.db_name)
	os.system(cmd)

	#===fusion finder===
	with poolcontext(processes = options.ncores) as pool:
		res_df = pool.map(partial(fusion_isoform_finder, options = options), CB_list)

	compression = None
	if options.o_name.endswith(".gz"):
		compression = 'gzip'
	pd.concat(res_df).reset_index(drop = True).to_csv(os.path.join(options.o_dir, options.o_name), sep='\t', header = True, index = False, compression = compression)

	print("\nBatch fusion finder jobs finished !\n")
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	print("\nBatch fusion finder spends %d : %d : %.2f\n" % (hours, minutes, seconds))


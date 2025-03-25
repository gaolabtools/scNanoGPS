#! /usr/bin/env python3

import time, os, sys, glob, gzip
import numpy as np
import pandas as pd
from optparse import OptionParser
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from scanner_core import misc
from curator_core import curator_io

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

def proc_fc(CB, options):
	import subprocess

	proc_time = time.time()

	print("Generating expression count in " + CB + "...", flush = True)
	cmd = options.featurecounts + " -L" + \
	      " -t gene -g gene_id -f --extraAttributes gene_name" + \
	      " -a " + options.gtf + \
	      " -o " + os.path.join(options.tmp_dir, CB) + ".expr.tsv " + \
	      os.path.join(options.tmp_dir, CB) + ".curated.minimap2.bam"

	proc = subprocess.check_call(cmd.split(' '), stderr = subprocess.DEVNULL)

	hours, minutes, seconds = misc.get_time_elapse(proc_time)
	print(CB + " done with %d : %d : %.2f" % (hours, minutes, seconds), flush = True)

def merge_expr(CB_list, options):
	merging_time = time.time()
	print("Merging expression matrix...", flush = True)
	expr_dict = dict()
	for CB in CB_list:
		expr_dict.update({CB: dict()})
		fh = open(os.path.join(options.tmp_dir, CB) + ".expr.tsv", "rt")
		#===skip cmd line===
		fh.readline()
		#===columns===
		fh.readline()
		#Geneid Chr Start End Strand Length gene_name AAATGCCCAATGGACG
		#     0   1     2   3      4      5         6                7
		while True:
			line = fh.readline()
			if not line:
				break
			line_list = line.rstrip().split("\t")
			expr_dict[CB].update({';'.join(line_list[0:7]): line_list[7]})
		fh.close()
	expr_df_raw = pd.DataFrame.from_dict(expr_dict)
	expr_df_idx = pd.DataFrame([ele.split(';') for ele in expr_df_raw.index.to_list()], columns = options.expr_cols)
	expr_df     = pd.concat([expr_df_idx, expr_df_raw.reset_index(drop = True)], axis = 1)

	hours, minutes, seconds = misc.get_time_elapse(merging_time)
	print("Merging expression matrix spent %d : %d : %.2f" % (hours, minutes, seconds), flush = True)

	return expr_df

if __name__ == "__main__":
	#===get params===
	parser = OptionParser()

	parser.add_option("-d",              dest = "o_dir",         nargs = 1, default = "scNanoGPS_res",
                  help = "Output directory name. "
                         "Default: scNanoGPS_res")
	parser.add_option("--tmp_dir",       dest = "tmp_dir",       nargs = 1, default = "tmp",
                  help = "Temporary folder name. "
                         "Default: tmp")
	parser.add_option("--gtf",           dest = "gtf",           nargs = 1, default = None,
                  help = "* Required ! "
                         "GTF file for expression calling. ")
	parser.add_option("-o",              dest = "o_name",        nargs = 1, default = "gene_exp_matrix.tsv.gz",
                  help = "Counting table name. "
                         "Default: gene_exp_matrix.tsv.gz")
	parser.add_option("--log",           dest = "log_f_name",    nargs = 1, default = 'logs/reporter_expression.log.txt',
                  help = "Log file name."
                         "Default: logs/reporter_expression.log.txt")
	parser.add_option("-t",              dest = "ncores",        nargs = 1, default = 1,
                  help = "Number of cores for program running. "
                         "Default: 1", type = "int")
	parser.add_option("--min_gene_no",   dest = "min_gene_no",   nargs = 1, default = 300,
                  help = "Minimal number of gene per cell. "
                         "Default: 300", type = "int")
	parser.add_option("--sel_bc_o",      dest = "sel_bc_o",      nargs = 1, default = "filtered_barcode_list.txt",
                  help = "Filtered cell barcode list. "
                         "Default: filtered_barcode_list.txt")
	parser.add_option("--keep_meta",     dest = "keep_meta",     nargs = 1, default = None,
                  help = "Keep meta files. Set to 1 to keep meta files."
                         "Default: None")
	parser.add_option("--featurecounts", dest = "featurecounts", nargs = 1, default = "featureCounts",
                  help = "Path to featureCounts."
                         "Default: featureCounts")
	options, arguments = parser.parse_args()

	#===pre-check===
	options.expr_cols  = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'gene_name']
	options.o_compress = ''
	if options.o_name.endswith('.gz'):
		options.o_compress = 'gzip'

	termination = False
	if not options.gtf or \
	   not os.path.isfile(options.gtf):
		print("\nCannot find GTF file: " + str(options.gtf) + "\n")
		termination = True
	if not os.path.isdir(options.o_dir):
		print("\nOutput directory is not exist: "    + options.o_dir   + "\n")
		termination = True
	if not os.path.isdir(options.tmp_dir):
		print("\nTemporary directory is not exist: " + options.tmp_dir + "\n")
		termination = True
	if not curator_io.find_exe(options.featurecounts):
		print("\nCannot find path to featureCounts !\n")
		termination = True

	if termination:
		parser.print_help()
		sys.exit(1)

	logger = open(options.log_f_name, "wt")
	logger.write(" ".join(sys.argv) + "\n\n")
	logger.write("Output directory:    " + options.o_dir + "\n")
	logger.write("Temporary directory: " + options.tmp_dir + "\n")
	logger.close()

	start_time = time.time()
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)

	#===load CB list===
	print("Loading CB list...", flush = True)
	bam_list = glob.glob(os.path.join(options.tmp_dir, "*.curated.minimap2.bam"))
	CB_list = []
	for bam_f in bam_list:
		split_f_name = bam_f.split('/')
		BC = split_f_name[len(split_f_name) - 1].split('.')
		CB_list.append(BC[0])

	print("Generating expression count...", flush = True)
	#===featureCounts===
	with poolcontext(processes = options.ncores) as pool:
		pool.map(partial(proc_fc, options = options), CB_list)

	#===merge df===
	expr_df = merge_expr(CB_list, options)

	#===filtering by gene number===
	counting_list = np.append(np.repeat(options.min_gene_no, 7), expr_df.iloc[:, 7:].astype('int32').apply(lambda x: sum(x > 0), axis = 0))
	res_df = expr_df[expr_df.columns.values[counting_list >= options.min_gene_no]].copy()
	res_df.to_csv(os.path.join(options.o_dir, options.o_name), header = True, index = False, sep = '\t', compression = options.o_compress)

	#===output filtered cell barcode list===
	with open(os.path.join(options.o_dir, options.sel_bc_o), "wt") as cbf:
		cbf.write("\n".join(res_df.columns.values[7:]) + "\n")

	#===remove meta files===
	if not options.keep_meta:
		for cb in CB_list:
			cmd = "rm " + os.path.join(options.tmp_dir, cb) + ".expr.tsv*"
			os.system(cmd)

	print("\nFinished !\n")
	print(str(len(res_df.columns.values[7:])) + " cell barcodes after filtering.")
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)


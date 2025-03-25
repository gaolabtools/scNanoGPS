#! /usr/bin/env python3

import time, os, sys, glob, gzip, subprocess
from multiprocessing.pool import ThreadPool as mp
import numpy as np
import pandas as pd
from functools import partial
from contextlib import contextmanager
from optparse import OptionParser
from scanner_core import misc
from curator_core import curator_io

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp(*args, **kwargs)
	yield pool
	pool.terminate()

def proc_fc(CB, options):
	fc_time = time.time()
	tmp_o = os.path.join(options.tmp_dir, CB) + options.expr_suf
	tmp_l = os.path.join(options.tmp_dir, CB) + options.expr_log
	bam_f = os.path.join(options.tmp_dir, CB) + ".curated.minimap2.bam"

	print("FeatureCounts: Calling expression profile in " + CB + '...', flush = True)
	cmd = options.featurecounts + " -L -t gene -g gene_id -f --extraAttributes gene_name" + \
              " -a " + options.gtf + " -o " + tmp_o + " " + bam_f

	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)

	with open(tmp_l, "wt") as lh:
		lh.write(out_msg.decode("utf-8"))
		lh.write(err_msg.decode("utf-8"))

	hours, minutes, seconds = misc.get_time_elapse(fc_time)
	print("featureCounts spend %d:%d:%.2f on %s" % (hours, minutes, seconds, CB), flush = True)

def combine_expr(CB_list, options):
	expr_columns = ["Geneid", "Chr", "Start", "End", "Strand", "Length", "gene_name"]
	res_df = pd.DataFrame([], columns = expr_columns)
	for CB in CB_list:
		tmp_o = os.path.join(options.tmp_dir, CB) + options.expr_suf
		cb_df = pd.read_csv(tmp_o, header = 0, sep = "\t", comment = '#')
		cb_df = cb_df.rename(columns = {cb_df.columns[7]: CB})

		#===limit gene no===
		if cb_df.loc[cb_df.iloc[:, 7] > 0, :].shape[0] >= options.min_gene_no:
			res_df = res_df.merge(cb_df, left_on = expr_columns, right_on = expr_columns, how = 'outer')

		if not options.keep_meta:
			cmd = "rm " + os.path.join(options.tmp_dir, CB) + options.expr_suf
			os.system(cmd)

	return res_df

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
	parser.add_option("-o",              dest = "o_name",        nargs = 1, default = "matrix.tsv.gz",
                          help = "Counting table name. "
                                 "Default: matrix.tsv.gz")
	parser.add_option("--log",           dest = "log_f_name",    nargs = 1, default = "reporter_expression.log.txt",
                          help = "Log file name."
                                 "Default: reporter_expression.log.txt")
	parser.add_option("-t",              dest = "ncores",        nargs = 1, default = 1,
                          help = "Number of cores for program running. "
                                 "Default: 1", type = "int")
	parser.add_option("--expr_suf",      dest = "expr_suf",      nargs = 1, default = ".expr.tsv",
                          help = "Single cell expression file suffix. "
                                 "Default: .expr.tsv")
	parser.add_option("--expr_log",      dest = "expr_log",      nargs = 1, default = ".expr.log.txt",
                          help = "FeatureCounts log file suffix. "
                                 "Default: .expr.log.txt")
	parser.add_option("--min_gene_no",   dest = "min_gene_no",   nargs = 1, default = 300,
                          help = "Minimal number of gene per cell. "
                                 "Default: 300", type = "int")
	parser.add_option("--keep_meta",     dest = "keep_meta",     nargs = 1, default = None,
                          help = "Set this parameter to 1 to keep meta files. "
                                 "Default: None")
	parser.add_option("--sel_bc_o",      dest = "sel_bc_o",      nargs = 1, default = "filtered_barcode_list.txt",
                          help = "Filtered cell barcode list. "
                                 "Default: filtered_barcode_list.txt")
	parser.add_option("--featurecounts", dest = "featurecounts", nargs = 1, default = "featureCounts",
                          help = "Path to featureCounts."
                                 "Default: featureCounts")
	options, arguments = parser.parse_args()

	#===pre-check===
	options.log_f_name = os.path.join(options.o_dir, options.log_f_name)

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

	#===log===
	logger = open(options.log_f_name, "wt")
	logger.write(" ".join(sys.argv) + "\n\n")
	logger.write("Output directory:    " + options.o_dir + "\n")
	logger.write("Temporary directory: " + options.tmp_dir + "\n")

	start_time = time.time()
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	logger.write("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")
	#===load CB list===
	print("Loading CB list...", flush = True)
	bam_list = glob.glob(os.path.join(options.tmp_dir, "*.curated.minimap2.bam"))
	CB_list = []
	for bam_f in bam_list:
		split_f_name = bam_f.split('/')
		BC = split_f_name[len(split_f_name) - 1].split('.')
		CB_list.append(BC[0])

	#===compute read count===
	print("Generating expression count...", flush = True)
	with poolcontext(processes = options.ncores) as pool:
		pool.map(partial(proc_fc, options = options), CB_list)

	#===combine dtable===
	res_df = combine_expr(CB_list, options)

	#===output result===
	compression = None
	if options.o_name.endswith(".gz"):
		compression = 'gzip'
	res_df.to_csv(os.path.join(options.o_dir, options.o_name), index = False, sep = '\t', compression = compression)

	#===output barcode list===
	with open(os.path.join(options.o_dir, options.sel_bc_o), "wt") as cbf:
		cbf.write("\n".join(res_df.columns[7:].tolist()) + "\n")

	print("\nFinished !\n")
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)
	logger.write("Reporter for expression profile spent: %d : %d : %.2f" % (hours, minutes, seconds) + "\n")
	logger.close()


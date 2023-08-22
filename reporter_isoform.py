#! /usr/bin/env python3

import time, os, sys, glob
import pandas as pd
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from optparse import OptionParser
from scanner_core import misc
from curator_core import curator_io

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

def proc_liqa(CB, options):
	bam_pref = os.path.join(options.tmp_dir, CB)

	if os.path.isfile(os.path.join(bam_pref + options.liqa_o)):
		print("LIQA output exist !!! Skip LIQA on " + CB + "!!!")
	else:
		liqa_time = time.time()
		print("LIQA: Calling isoform in " + CB + '...', flush = True)
		cmd = options.liqa + ' -task quantify -max_distance 20 -f_weight 1' + \
		      ' -refgene ' + options.liqa_ref + \
		      ' -bam ' + bam_pref + '.curated.minimap2.bam' + \
		      ' -out ' + bam_pref + options.liqa_o

		code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
		with open(bam_pref + options.liqa_log, "wt") as fh:
			fh.write(out_msg.decode("utf-8"))
			fh.write(err_msg.decode("utf-8"))
		hours, minutes, seconds = misc.get_time_elapse(liqa_time)
		print("LIQA spend %d:%d:%.2f on %s" % (hours, minutes, seconds, CB), flush = True)

def parse_txid2txname(options):
	tx_mapping = {}
	fh = open(options.gtf, "rt")
	while True:
		line = fh.readline()
		if not line:
			break

		if line.startswith("#"):
			continue

		line_list = line.rstrip().split("\t")
		if line_list[2] != "transcript":
			continue

		attr_list = line_list[8].split('; ')
		txid, txname = "", ""
		for x in attr_list:
			ele_list = x.split(' ')
			if ele_list[0] == "transcript_id":
				txid   = ele_list[1].split('"')[1]
			if ele_list[0] == "transcript_name":
				txname = ele_list[1].split('"')[1]

		tx_mapping.update({txid: txname})
	fh.close()

	return tx_mapping

if __name__ == "__main__":

	#===get params===
	parser = OptionParser()

	parser.add_option("-d",         dest = "o_dir",      nargs = 1, default = "scNanoGPS_res",
	                  help = "Output directory name. "
	                         "Default: scNanoGPS_res")
	parser.add_option("--tmp_dir",  dest = "tmp_dir",    nargs = 1, default = "tmp",
	                  help = "Temporary folder name. "
	                         "Default: tmp")
	parser.add_option("--CB_file",  dest = "CB_file",    nargs = 1, default = "filtered_barcode_list.txt",
                          help = "File name for filtered barcode list. "
                                 "Default: filtered_barcode_list.txt")
	parser.add_option("--gtf",      dest = "gtf",        nargs = 1, default = None,
                          help = "GTF file for obtaining transcript ID. ")
	parser.add_option("--liqa_ref", dest = "liqa_ref",   nargs = 1, default = None,
	                  help = "* Required ! "
	                         "Reference of LIQA. ")
	parser.add_option("-o",         dest = "o_name",     nargs = 1, default = "matrix_isoform.tsv",
	                  help = "Counting table name. "
	                         "Default: matrix_isoform.tsv")
	parser.add_option("--log",      dest = "log_f_name", nargs = 1, default = "reporter_isoform.log.txt",
	                  help = "Log file name. "
	                         "Default: reporter_isoform.log.txt")
	parser.add_option("-t",         dest = "ncores",     nargs = 1, default = 1,
	                  help = "Number of cores for program running. "
	                         "Default: 1", type = "int")
	parser.add_option("--liqa",     dest = "liqa",       nargs = 1, default = "liqa",
	                  help = "Program name of LIQA. "
	                         "Default: liqa")
	parser.add_option("--liqa_log", dest = "liqa_log",   nargs = 1, default = ".liqa.log",
	                  help = "Suffix of LIQA output file. "
	                         "Default: .liqa.log")
	parser.add_option("--liqa_o",   dest = "liqa_o",     nargs = 1, default = ".liqa.tsv",
	                  help = "Suffix of LIQA output file. "
	                         "Default: .liqa.tsv")

	options, arguments = parser.parse_args()

	#===pre-check===
	options.log_f_name = os.path.join(options.o_dir, options.log_f_name)

	termination = False
	if not options.liqa_ref or \
	   not os.path.isfile(options.liqa_ref):
		print("\nCannot find reference of LIQA")
		termination = True
	if not os.path.isdir(options.o_dir):
		print("\nOutput directory is not exist: "    + options.o_dir    + "\n")
		termination = True
	if not os.path.isdir(options.tmp_dir):
		print("\nTemporary directory is not exist: " + options.tmp_dir  + "\n")
		termination = True
	if not curator_io.find_exe(options.liqa):
		print("\nCannot find LIQA !\n")
		termination = True
	if not os.path.isfile(os.path.join(options.o_dir, options.CB_file)):
		print("\nCannot find filter barcode list at: " + os.path.join(options.o_dir, options.CB_file)  + "\n")
		print("All the barcodes under folder " + options.tmp_dir + " will be used insteat.\n")
		options.CB_file = None

	if termination:
		parser.print_help()
		sys.exit(1)

	logger = open(options.log_f_name, "wt")
	logger.write(" ".join(sys.argv) + "\n\n")
	logger.write("Output directory:    " + options.o_dir    + "\n")
	logger.write("Temporary directory: " + options.tmp_dir  + "\n")
	logger.write("Reference of LIQA:   " + options.liqa_ref + "\n")
	logger.close()

	#=== set env variables
	os.environ["OMP_NUM_THREADS"]        = str(options.ncores)
	os.environ["OPENBLAS_NUM_THREADS"]   = str(options.ncores)
	os.environ["MKL_NUM_THREADS"]        = str(options.ncores)
	os.environ["VECLIB_MAXIMUM_THREADS"] = str(options.ncores)
	os.environ["NUMEXPR_NUM_THREADS"]    = str(options.ncores)

	import numpy as np

	start_time = time.time()
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
	#===load CB list===
	print("Loading CB list...", flush = True)
	bam_list, CB_list = [], []
	if not options.CB_file:
		bam_list = glob.glob(os.path.join(options.tmp_dir, "*.curated.minimap2.bam"))
		for bam_f in bam_list:
			split_f_name = bam_f.split('/')
			BC = split_f_name[len(split_f_name) - 1].split('.')
			CB_list.append(BC[0])
	else:
		with open(os.path.join(options.o_dir, options.CB_file), "rt") as CBF:
			while True:
				CB_name = CBF.readline().rstrip()
				if not CB_name:
					break
				CB_list.append(CB_name)
		bam_list = [x + ".curated.minimap2.bam" for x in CB_list]

	#===parse transcript ID===
	if options.gtf:
		print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
		print("Parsing GTF...", flush = True)
		options.tx_mapping = parse_txid2txname(options)

	#===compute isoform===
	with poolcontext(processes = options.ncores) as pool:
		pool.map(partial(proc_liqa, options = options), CB_list)

	print("\nBatch LIQA jobs finished !\n")
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	logger = open(options.log_f_name, "at")
	logger.write("\nBatch LIQA spends %d : %d : %.2f\n" % (hours, minutes, seconds))
	logger.close()

	#===parse & combine results===
	res_tb = dict()
	for CB in CB_list:
		res_tb.update({CB: {}})
		fh = open(os.path.join(options.tmp_dir, CB) + options.liqa_o, "rt")
		#---skip def line---
		fh.readline()
		while True:
			line = fh.readline().rstrip()
			if not line:
				break
			line_list = line.split("\t")
			if options.gtf:
				res_tb[CB].update({options.tx_mapping[line_list[1]] + '_' + line_list[1]: line_list[2]})
			else:
				res_tb[CB].update({line_list[0] + '_' + line_list[1]: line_list[2]})
		fh.close()

	pd.DataFrame(res_tb).to_csv(os.path.join(options.o_dir, options.o_name), sep='\t', na_rep = 0.0, header = True, index = True)

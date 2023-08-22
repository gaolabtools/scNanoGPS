#! /usr/bin/env python3

import time, os, sys, glob, gzip, subprocess
from multiprocessing.pool import ThreadPool as mp
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

def compute_read_count(CB, options):
	comp_time = time.time()
	tmp_o = os.path.join(options.tmp_dir, CB) + options.expr_suf
	bam_f = os.path.join(options.tmp_dir, CB) + ".curated.minimap2.bam"

	print("Computing expression profile in " + CB + '...', flush = True)
	oh = open(tmp_o, "wt")
	oh.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tgene_name\t" + CB + "\n")

	for gtf in options.gtf_data:
		strand_str = ""
		if gtf['strand'] == "+":
			strand_str = "-f 16"
		else:
			strand_str = "-F 16"
		cmd = options.samtools + " view -c " + strand_str + " " +\
		      bam_f + " " + gtf['chr'] + ":" + gtf['start'] + "-" + gtf['end']
		code_msg, out_msg, err_msg = curator_io.sys_run(cmd)

		oh.write(gtf['gene_id'] + "\t" + gtf['chr'] + "\t" + gtf['start'] + "\t" + gtf['end'] + "\t" + gtf['strand'] + "\t" + str(int(gtf['end']) - int(gtf['start']))+ "\t" + gtf['gene_name'] + "\t" + str(out_msg.decode("utf8")).rstrip() + "\n")
	oh.close()
	hours, minutes, seconds = misc.get_time_elapse(comp_time)
	print("Computing read count in %s spends %d:%d:%.2f" % (CB, hours, minutes, seconds), flush = True)

def load_gtf(options):
	gtf_data = list()
	print("Loading GTF ...", flush = True)
	with open(options.gtf, "rt") as fh:
		while True:
			line = fh.readline()
			if not line:
				break
			line_arr = line.split("\t")
			if len(line_arr) < 8:
				continue
			if line_arr[2] != "gene":
				continue
			tmp_dict = dict()
			attr_arr = line_arr[8].split('; ')
			tmp_dict.update({'chr': line_arr[0], 'start': line_arr[3], 'end': line_arr[4], 'strand': line_arr[6]})
			for attr in attr_arr:
				attr_kv = attr.split(' "')
				if attr_kv[0] == 'gene_id':
					tmp_dict.update({'gene_id':   attr_kv[1].split('"')[0]})
				if attr_kv[0] == 'gene_name':
					tmp_dict.update({'gene_name': attr_kv[1].split('"')[0]})
			gtf_data.append(tmp_dict)
	print("Done")
	return(gtf_data)

def combine_expr(CB_list, options):
	print("Combining final matrix...", end = "", flush = True)
	comb_time = time.time()
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
	print("Done")
	hours, minutes, seconds = misc.get_time_elapse(comb_time)
	print("Combining final matrix spends %d:%d:%.2f on %s" % (hours, minutes, seconds, CB), flush = True)
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
	parser.add_option("--samtools",      dest = "samtools",      nargs = 1, default = "samtools",
                          help = "Path to samtools."
                                 "Default: samtools")
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
	if not curator_io.find_exe(options.samtools):
		print("\nCannot find path to samtools !\n")
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

	#===load gtf ===
	options.gtf_data = load_gtf(options)

	#===compute read count===
	print("Generating expression count...", flush = True)
	with poolcontext(processes = options.ncores) as pool:
		pool.map(partial(compute_read_count, options = options), CB_list)

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


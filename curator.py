#! /usr/bin/env python3

import time, os, sys, subprocess, glob
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from scanner_core import misc
from curator_core import curator_io, preprocessing

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

if __name__ == "__main__":

	#===get params===
	parser = preprocessing.getOptions()
	options, arguments = parser.parse_args()

	#===pre-check===
	parser, options = preprocessing.precheck(parser, options)
	logger = open(options.log_f_name, "at")
	logger.write(" ".join(sys.argv) + "\n\n")

	#=== set env variables
	os.environ["OMP_NUM_THREADS"]        = str(options.ncores)
	os.environ["OPENBLAS_NUM_THREADS"]   = str(options.ncores)
	os.environ["MKL_NUM_THREADS"]        = str(options.ncores)
	os.environ["VECLIB_MAXIMUM_THREADS"] = str(options.ncores)
	os.environ["NUMEXPR_NUM_THREADS"]    = str(options.ncores)

	import numpy as np

	#===load CB_counting===
	start_time = time.time()
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
	print("Loading CB_counting: " + options.CB_count + " ...", end = "", flush = True)
	CB_counting_df = curator_io.load_df(options.CB_count)
	#===set index, contributed by Philipp Rentzsch from GitHub===
	CB_counting_df.set_index('idx', inplace = True)
	#===set index, contributed by Philipp Rentzsch from GitHub===
	print("Done", flush = True)

	#===load CB_mrg===
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
	print("Loading CB_mrg: " + options.CB_mrg + " ...", end = "", flush = True)
	CB_dict = dict()
	CB_mrg  = curator_io.open_file(options.CB_mrg, "rt")
	#---skip def line---
	CB_mrg.readline()
	#===contributed by Philipp Rentzsch from GitHub===
	for line in CB_mrg:
		lines  = line.rstrip().split(": ")
		rep_BC = CB_counting_df.loc[int(lines[0]), 'BC']
		all_id = curator_io.extract_id_list(lines[1])
		BC_seq_list = CB_counting_df.loc[all_id, 'BC']
		for ele in BC_seq_list:
			CB_dict[ele] = rep_BC
	#===contributed by Philipp Rentzsch from GitHub===
	CB_mrg.close()
	print("Done", flush = True)
	del CB_counting_df
	print("Time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	#===parse reads===
	time_separation = time.time()
	print("Separation of reads by cell barcodes ...", flush = True)
	BC_list = curator_io.open_file(options.BC_list, "rt")
	# 0   1           2        3  4   5
	# rid orientation BC_start BC UMI mean_BC_quality
	#---skip def line---
	BC_list.readline()

	fastq_f = curator_io.open_file(options.fq_name, "rt")
	while True:
		BC_line = BC_list.readline()
		if not BC_line:
			break

		BC_lines = BC_line.rstrip().split("\t")

		fastq_f.readline()
		seq_line = fastq_f.readline()
		fastq_f.readline()
		qua_line = fastq_f.readline()

		if BC_lines[3] in CB_dict:
			o_name = os.path.join(options.tmp_dir, CB_dict[BC_lines[3]] + ".fastq.gz")
			writer = curator_io.open_file(o_name, "at")
			writer.write("@" + BC_lines[0] + "_" + BC_lines[4] + "\n")
			writer.write(seq_line)
			writer.write("+" + BC_lines[0] + "_" + BC_lines[4] + "\n")
			writer.write(qua_line)
			writer.close()
	BC_list.close()
	fastq_f.close()
	if not options.keep_meta:
		os.system("mv " + options.fq_name  + " " + options.tmp_dir)
		os.system("mv " + options.CB_count + " " + options.tmp_dir)
		os.system("mv " + options.CB_mrg   + " " + options.tmp_dir)

	print("            \rDone\n", flush = True)
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)
	hours, minutes, seconds = misc.get_time_elapse(time_separation)
	logger.write("Separation of reads by cell barcode spent %d : %d : %.2f\n\n" % (hours, minutes, seconds))

	#=== get fastq files list ===
	fq_list = glob.glob(os.path.join(options.tmp_dir, "*.fastq.gz"))
	with poolcontext(processes = options.ncores) as pool:
		logger_res = pool.map(partial(curator_io.curation_master, options = options), fq_list)
	logger.write("\n".join(logger_res) + "\n")

	print(str(len(fq_list)) + " barcodes are done" + " " * 50, flush = True)
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	logger.write("\nCuration process time spent: %d : %d : %.2f\n" % (hours, minutes, seconds))
	logger.close()

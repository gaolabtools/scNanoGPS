#! /usr/bin/env python3

import os, time, gzip
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from scanner_core import scanner_io, preprocessing, searching_core, misc

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

if __name__ == "__main__":

	#===get parameters===
	parser = preprocessing.getOptions()
	options, arguments = parser.parse_args()

	#===precheck===
	parser, options, arguments, counter = preprocessing.precheck(parser, options, arguments)

	#===parameters===
	#---CB list col names---
	df_col = ['rid', 'orientation', 'BC_start', 'BC', 'UMI', 'Seq_end', 'mean_BC_quality']
	with gzip.open(options.bc_f_name, "wt") as BC_output:
		BC_output.write("\t".join(df_col) + "\n")

	#---result counting parameters---
	rid = 0

	#===time benchmarking===
	start_time = time.time()

	print()
	#===read input fastq===
	if options.isFile:
		f_list = list([options.fq_f_name])
		f_idx  = 0
		print("Processing file name: " + str(options.fq_f_name))
		reader = scanner_io.open_file(options.fq_f_name, "rt")

	if options.isDir:
		f_list = preprocessing.get_valid_input_file_list(options.cwd, options.fq_f_name)
		f_idx  = 0
		if not f_list[f_idx].endswith(".fast5"):
			print("Processing file name: " + str(f_list[f_idx]))
		reader = scanner_io.open_file(os.path.join(options.cwd, options.fq_f_name, f_list[f_idx]), "rt")

	writer = scanner_io.open_file(options.fq_o_name, "wt")

	while True:
		batch_data = []
		res_data   = []
		batch_data, eof, rid, reader, f_idx = scanner_io.batch_reading(reader, f_list, f_idx, int(options.batch_no), rid)
		if eof == 1:
			break

		with poolcontext(processes = options.ncores) as pool:
			tmp_data = pool.map(partial(searching_core.ten_nano_workflow, options = options), batch_data)

		#===counting===
		res_data, counter = searching_core.counting_res(res_data, counter, tmp_data)

		#===time benchmarking===
		hours, minutes, seconds = misc.get_time_elapse(start_time)
		print("Processed " + str(rid) + " reads, ", end = "")
		misc.report_time_elapse(hours, minutes, seconds)

		#===output result===
		with gzip.open(options.bc_f_name, "at") as BC_output:
			counter_res_data, counter_batch_data = 0, 0
			for each_row, each_data in zip(res_data, batch_data):
				if each_row[1]:
					BC_output.write("\t".join(str(x) for x in each_row) + "\n")
					scanner_io.output_fastq(writer, each_row, each_data)

	reader.close()
	writer.close()

	print("\nProcessed total " + str(rid) + " reads.")
	print("Found " + str(counter['counter_h_3p_polyT'] + counter['counter_t_3p_polyT']) + " reads have 3'-adaptor + polyT.")
	print("Detecting rate: %.2f" % round((counter['counter_h_3p_polyT'] + counter['counter_t_3p_polyT']) / rid * 100, 2) + '%')

	#===output parameters===
	scanner_io.output_params(options, rid, counter, start_time)

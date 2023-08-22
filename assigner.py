#! /usr/bin/env python3

import os, time, sys
import pandas as pd
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from scanner_core import misc, scanner_io
from assigner_core import preprocessing, wrapping, painting

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

if __name__ == "__main__":

	#===get params===
	parser = preprocessing.getOptions()
	options, arguments = parser.parse_args()

	#===precheck===
	parser, options = preprocessing.precheck(parser, options)

	#=== set env variables
	os.environ["OMP_NUM_THREADS"]        = str(options.ncores)
	os.environ["OPENBLAS_NUM_THREADS"]   = str(options.ncores)
	os.environ["MKL_NUM_THREADS"]        = str(options.ncores)
	os.environ["VECLIB_MAXIMUM_THREADS"] = str(options.ncores)
	os.environ["NUMEXPR_NUM_THREADS"]    = str(options.ncores)

	import numpy as np

	#===loading data===
	start_time = time.time()
	print("Loading file: " + options.input + " ...", flush = True)
	BC_df = pd.read_table(options.input, sep = '\t', header = 0, compression = options.i_compression)
	print("Done\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

	#===counting UMI===
	step_time = time.time()
	print("Counting UMI ...", flush = True)
	UMI_df = wrapping.get_UMI_counting(BC_df, options)
	hours, minutes, seconds = misc.get_time_elapse(step_time)
	print("Counting UMI spend: %d : %d : %.2f" % (hours, minutes, seconds), flush = True)
	BC_df = None

	#===estimate cell number===
	if options.forced_no == 0:
		step_time = time.time()
		print("Estimating cell number ...", flush = True)
		cell_no, cell_no_ext = wrapping.estimate_cell_no(UMI_df, options)
		with open(options.log_f_name, "at") as logger:
			logger.write("\tEstimated cell number(plus extension):" + str(cell_no_ext) + "\n")
			logger.write("\n")
		hours, minutes, seconds = misc.get_time_elapse(step_time)
		print("Estimating cell number spend: %d : %d : %.2f" % (hours, minutes, seconds))
	else:
		print("Warning!\nThe cell number is forcely assigned to " + str(options.forced_no) + " !")

	#===output UMI table===
	step_time = time.time()
	print("\nOutput UMI counting table: " + options.output + "...", end = "", flush = True)
	UMI_df.loc[::, ["idx", "BC", "UMI", "log10_idx", "log10_UMI", "med_log10_slope", "log10_slope"]].to_csv(options.output, sep = '\t', header = True, index = False, compression = options.o_compression)
	print("Done", flush = True)

	hours, minutes, seconds = misc.get_time_elapse(step_time)
	misc.report_time_elapse(hours, minutes, seconds)

	if options.forced_no == 0:
		#===select CB having read no. passing threshold===
		print("Preparing CB table...", flush = True)
		mrg_sel = UMI_df.loc[UMI_df['idx'] <= cell_no_ext, ['idx', 'BC']].sort_values('idx', ascending = True).reset_index(drop = True)
		print("Done", flush = True)

		#===Calc distance===
		print("Computing Levenshtein distance", flush = True)

# mrg_sel:
#          idx                BC
# 0          1  CTACGAAGTGATGAGG
# 1          2  TTGTGCCTCATTGACA

		with poolcontext(processes = options.ncores) as pool:
			pool.map(partial(wrapping.batch_seq_comp, target = mrg_sel, options = options), mrg_sel.iloc[0:(mrg_sel.shape[0] - 1)].values.tolist())

		#===merging CB===
		step_time = time.time()
		print("Merging table...", flush = True)
		res_df = wrapping.merge_cb_new(mrg_sel, options)
		print("Done", flush = True)
		hours, minutes, seconds = misc.get_time_elapse(step_time)
		misc.report_time_elapse(hours, minutes, seconds)

		#===output merged result===
		step_time = time.time()
		print("Writing result to " + options.CB_mrg_o + " ... ", end = "", flush = True)
		res_df = res_df.sort_values(by = ['id2', 'id1'], ascending = True)
		writer = scanner_io.open_file(options.CB_mrg_o, "wt")
		writer.write("Representative cell barcode id:\tincluded cell barcode id")
		cell_no_mrg = 0
		this_id = None
		for idx, row in res_df.iterrows():
			if this_id != row['id2']:
				this_id = row['id2']
				writer.write("\n" + str(row['id2']) + ": " + str(row['id1']))
				cell_no_mrg += 1
			else:
				writer.write("; " + str(row['id1']))
		writer.write("\n")
		writer.close()
		print("Done\n", flush = True)
		hours, minutes, seconds = misc.get_time_elapse(step_time)
		misc.report_time_elapse(hours, minutes, seconds)
	else:
		writer = scanner_io.open_file(options.CB_mrg_o, "wt")
		writer.write("Representative cell barcode id:\tincluded cell barcode id")
		for idx in range(1, (options.forced_no) + 1, 1):
			writer.write(str(idx) + ": " + str(idx) + "\n")
		writer.close()

	#===draw log10 dist plot===
	step_time = time.time()
	print("Plotting read number distribution ...", end = "", flush = True)
	if options.forced_no == 0:
		painting.draw_log_dist_plot(UMI_df, cell_no, cell_no_ext, cell_no_mrg, options)
	else:
		painting.draw_log_dist_plot(UMI_df, options.forced_no, options.forced_no, options.forced_no, options)
	print("Done", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(step_time)
	misc.report_time_elapse(hours, minutes, seconds)

	print("Finish time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), flush = True)

	with open(options.log_f_name, "at") as logger:
		logger.write("Finished time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")

	hours, minutes, seconds = misc.get_time_elapse(start_time)
	print("Assigner time spend: %d : %d : %.2f" % (hours, minutes, seconds))


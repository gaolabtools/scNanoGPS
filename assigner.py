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
	cell_no, cell_no_ext, cell_no_mrg = 1, 1, options.forced_no
	if options.forced_no == 0:
		step_time = time.time()
		print("Estimating cell number ...", flush = True)
		cell_no, cell_no_ext = wrapping.estimate_cell_no(UMI_df, options)
		hours, minutes, seconds = misc.get_time_elapse(step_time)
		print("Estimating cell number spend: %d : %d : %.2f" % (hours, minutes, seconds))
	else:
		cell_no, cell_no_ext = options.forced_no, options.forced_no
		print("Warning!\nThe cell number is forcely assigned to " + str(options.forced_no) + " !")

	#===output UMI table===
	step_time = time.time()
	print("\nOutput UMI counting table: " + options.output + "...", end = "", flush = True)
	UMI_df.loc[::, ["idx", "BC", "UMI", "log10_idx", "log10_UMI", "log10_slope"]].to_csv(options.output, sep = '\t', header = True, index = False, compression = options.o_compression)
	print("Done", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(step_time)
	misc.report_time_elapse(hours, minutes, seconds)

	writer = scanner_io.open_file(options.CB_mrg_o, "wt")
	writer.write("Representative cell barcode id:\tincluded cell barcode id\n")
	#===force cell no===
	if (not options.whitelist) & (options.forced_no != 0):
		for idx in range(1, (options.forced_no + 1) + 1, 1):
			writer.write(str(idx) + ": " + str(idx) + "\n")

	else:
		print("Preparing BC table...\n", flush = True)
		#===having wbc===
		if options.whitelist:
			#===find out exactly matched BC===
			target_df = UMI_df.iloc[0:cell_no_ext, ].loc[UMI_df['BC'].isin(options.wbc['BC']), ['idx', 'BC']].copy()
			print(str(target_df.shape[0]) + " of " + str(cell_no_ext) + " barcodes are exactly matched to whitelist !\n")
			#===mismatched BC===
			query_list = UMI_df.iloc[0:cell_no_ext, ].loc[~UMI_df['BC'].isin(options.wbc['BC']), ['idx', 'BC']].values
			print("Checking the rest " + str(len(query_list)) + " Nanopore barcodes...\n")

		#===no wbc===
		else:
			#===select CB having read no. passing threshold===
			target_df = UMI_df.iloc[0:cell_no, ].loc[:, ['idx', 'BC']].copy()
			print("Estimated " + str(cell_no) + " raw barcodes detected !\n", flush = True)
			query_list = UMI_df.iloc[cell_no:cell_no_ext, ].loc[:, ['idx', 'BC']].values
			print("Checking the extended " + str(len(query_list)) + " Nanopore barcodes...\n")

		#===mapping mismatched BC back to exact matched BC===
		with poolcontext(processes = options.ncores) as pool:
			pool.map(partial(wrapping.batch_seq_comp, target = target_df.rename(columns = {'idx': 'rep_id'}), options = options), query_list)

		step_time = time.time()
		print("Merging distance table...", flush = True)
		CB_mrg_dist_df = wrapping.merge_cb(options)
		print("Done", flush = True)
		hours, minutes, seconds = misc.get_time_elapse(step_time)
		misc.report_time_elapse(hours, minutes, seconds)

# CB_mrg_dist_df:
#    rep_id  Nanopore_id  distance
#0     4127         3837         1

		#===output BC merged list===
		CB_mrg_dist_df_sorted = CB_mrg_dist_df.sort_values(by = ['rep_id', 'Nanopore_id'], ascending = True).reset_index(drop = True)
		this_id = 0

		for idx, row in target_df.iterrows():
			writer.write(str(row['idx']) + ": " + str(row['idx']))
			if CB_mrg_dist_df_sorted.shape[0] > 0:
				while True:
					if this_id > max(CB_mrg_dist_df_sorted.index):
						break
					if CB_mrg_dist_df_sorted.loc[this_id, 'rep_id'] == row['idx']:
						writer.write("; " + str(CB_mrg_dist_df_sorted.loc[this_id, 'Nanopore_id']))
						this_id = this_id + 1
					else:
						break
			cell_no_mrg = cell_no_mrg + 1
			writer.write("\n")
		if not options.whitelist:
			for idx, row in UMI_df.iloc[cell_no:cell_no_ext, ].loc[:, ['idx', 'BC']].iterrows():
				if not row['idx'] in CB_mrg_dist_df_sorted['Nanopore_id'].values:
					writer.write(str(row['idx']) + ": " + str(row['idx']) + "\n")
					cell_no_mrg = cell_no_mrg + 1

	writer.close()
	print("Done\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(step_time)
	misc.report_time_elapse(hours, minutes, seconds)

	#===draw log10 dist plot===
	step_time = time.time()
	print("Plotting read number distribution ...", end = "", flush = True)
	if options.forced_no == 0:
		painting.draw_log_dist_plot(UMI_df, cell_no, cell_no_ext, cell_no_mrg, options)
		with open(options.log_f_name, "at") as logger:
			logger.write("\tEstimated cell number(plus extension):" + str(cell_no_mrg) + "\n")
			logger.write("\n")
	else:
		painting.draw_log_dist_plot(UMI_df, options.forced_no, options.forced_no, options.forced_no, options)
		with open(options.log_f_name, "at") as logger:
			logger.write("\tThe cell number is forcely assigned:  " + str(options.forced_no) + "\n")
			logger.write("\n")
	print("Done", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(step_time)
	misc.report_time_elapse(hours, minutes, seconds)

	print("Finish time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), flush = True)

	with open(options.log_f_name, "at") as logger:
		logger.write("Finished time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")

	hours, minutes, seconds = misc.get_time_elapse(start_time)
	print("Assigner time spend: %d : %d : %.2f" % (hours, minutes, seconds))


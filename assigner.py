#! /usr/bin/env python3

import time
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

	#===loading data===
	print("Time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	print("Loading file: " + options.input + " ...", flush = True)
	start_time = time.time()
	BC_df = pd.read_table(options.input, sep = '\t', header = 0, compression = options.i_compression)
	print("Done\n", flush = True)

	#===counting UMI===
	print("Time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	UMI_df = wrapping.get_UMI_counting(BC_df, options)

	#===estimate cell number===
	print("Time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	cell_no, cell_no_ext = wrapping.estimate_cell_no(UMI_df, options)
	with open(options.log_f_name, "at") as logger:
		logger.write("\tEstimated cell number(plus extension):" + str(cell_no_ext) + "\n")
		logger.write("\n")

	#===output UMI table===
	print("Time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	print("\nOutput UMI counting table: " + options.output + "...", end = "", flush = True)
	UMI_df.loc[::, ["idx", "BC", "UMI", "log10_idx", "log10_UMI", "med_log10_slope", "log10_slope"]].to_csv(options.output, sep = '\t', header = True, index = False, compression = options.o_compression)
	print("Done", flush = True)

	#===select CB having read no. passing threshold===
	print("Time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n")
	mrg_sel = UMI_df.loc[UMI_df['idx'] <= cell_no_ext, ['idx', 'BC']].sort_values('idx', ascending = True).reset_index(drop = True)

	#===calc distance===
	dist_col = ['id1', 'id2', 'distance']
	dist_df  = pd.DataFrame(columns = dist_col)
	for i in range(0, mrg_sel['BC'].count() - 1):
		print("Preparing Cell barcode - batch job: " + str(i + 1) + " of " + str(mrg_sel['BC'].count() - 1) + " ...")
		batch_data = list()

		for j in range(i + 1, mrg_sel['BC'].count()):
			batch_data.append({'id1': mrg_sel.loc[i, 'idx'],
			                   'BC1': mrg_sel.loc[i, 'BC'],
			                   'id2': mrg_sel.loc[j, 'idx'],
			                   'BC2': mrg_sel.loc[j, 'BC']})
		print("Loaded " + str(len(batch_data)) + " cell barcode pairs...")

		print("Calculating Lavenshstein distance...", end = "")
		mrg_batch_time = time.time()
		#===parallel computing===
		with poolcontext(processes = int(options.ncores)) as pool:
			dist_df = pd.concat([dist_df, \
                                     pd.DataFrame(    \
                                        pool.map(partial(wrapping.seq_comp, options = options), batch_data), columns = dist_col)])

		print("Done")
		hours, minutes, seconds = misc.get_time_elapse(mrg_batch_time)
		misc.report_time_elapse(hours, minutes, seconds)
		print()

	dist_df = dist_df.reset_index(drop = True)
	dist_df.to_csv(options.CB_mrg_dist, sep = '\t', header = True, index = False, compression = options.CB_mrg_dist_compression)

	#===merging CB===
	mrg_time = time.time()
	print("Merging CB with Laventhstein distance less than " + str(options.CB_mrg_thr) + " ...", end = "")
	res_df = wrapping.merge_cb(mrg_sel, dist_df, options.CB_mrg_thr)
	print("Done")
	hours, minutes, seconds = misc.get_time_elapse(mrg_time)
	misc.report_time_elapse(hours, minutes, seconds)
	print()

	#===output merged result===
	print("Writing result to " + options.CB_mrg_o + " ... ", end = "")
	res_df = res_df.sort_values(by = ['id2', 'id1'], ascending = True)
	writer = scanner_io.open_file(options.CB_mrg_o, "wt")
	writer.write("Representative cell barcode id:\tincluded cell barcode id")
	cell_no_mrg = 0
	this_id = None
	for idx, row in res_df.iterrows():
		if this_id != row['id2']:
			this_id = row['id2']
			writer.write("\n" + row['id2'] + ": " + row['id1'])
			cell_no_mrg += 1
		else:
			writer.write("; " + row['id1'])
	writer.write("\n")
	writer.close()
	print("Done\n")

	#===draw log10 dist plot===
	print("Plotting read number distribution ...", end = "", flush = True)
	painting.draw_log_dist_plot(UMI_df, cell_no, cell_no_ext, cell_no_mrg, options)
	print("Done", flush = True)

	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)
	print("Finish time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))

	with open(options.log_f_name, "at") as logger:
		logger.write("Finished time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")

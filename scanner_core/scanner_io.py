def open_file(f_name, mode):
	if mode == "rt":
		if f_name.endswith(".fast5"):
			import h5py
			reader = h5py.File(f_name, "r")
		elif f_name.endswith(('.fastq', '.fq')):
			reader = open(f_name, mode)
		elif f_name.endswith(('.fastq.gz', '.fq.gz')):
			import gzip
			reader = gzip.open(f_name, mode)
		return reader

	elif mode == "wt":
		if f_name.endswith('.gz'):
			import gzip
			writer = gzip.open(f_name, mode)
		else:
			writer = open(f_name, mode)
		return writer

	else:
		print("Error occurred during open file: " + str(f_name) + " !\n\n")
		sys.exit(1)

def output_result(o_name, stat_list, df_col):
	import pandas as pd

	stat_df = pd.DataFrame(stat_list, columns = df_col)
	stat_df.to_csv(o_name, header = df_col, index = None, sep = '\t', mode = 'w')

def output_fastq(writer, each_row, each_data):
	if each_row[2]:
		na_seq = each_data['na_seq']
		qu_seq = each_data['qu_seq']
		#---trim off BC & UMI---
		start_pos = int(each_row[2]) + len(each_row[3]) + len(each_row[4]) - 1
		#---no info for timming 5' adaptor---
		end_pos = len(na_seq)
		if each_row[5]:
			end_pos = int(each_row[5])
		if each_row[1] == "T":
			na_seq = each_data['na_seq'].translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))[::-1]
			qu_seq = each_data['qu_seq'][::-1]
			start_pos += 1
		writer.write(each_data['def_line'] + "\n")
		writer.write(na_seq[start_pos:end_pos] + "\n")
		writer.write("+\n")
		writer.write(qu_seq[start_pos:end_pos] + "\n")
	return

def output_params(options, rid, counter, start_time):
	import time
	from scanner_core import misc
	hours, minutes, seconds = misc.get_time_elapse(start_time)

	with open(options.log_f_name, "a+") as logger:
		logger.write("Total " + str(rid) + " reads are processed.\n")
		logger.write("Time elapse: %d : %d : %.2f" % (hours, minutes, seconds) + "\n")
		logger.write("Detecting rate: %.2f" % round((counter['counter_h_3p'] + counter['counter_t_3p']) / rid * 100, 2) + '%' + "\n\n")
		logger.write("Result counting:\n")
		logger.write("\tNumber of 3'-adaptor located on the read head region:           \t" + str(counter['counter_h_3p']) + "\n")
		logger.write("\tNumber of 3'-adaptor + polyT on the read head region:           \t" + str(counter['counter_h_3p_polyT']) + "\n")
		logger.write("\tNumber of 3'-adaptor located on the read tail region:           \t" + str(counter['counter_t_3p']) + "\n")
		logger.write("\tNumber of 3'-adaptor + polyT on the read tail region:           \t" + str(counter['counter_t_3p_polyT']) + "\n\n")

		logger.write("Alignment counting:\n")
		logger.write("\tNumber of 3'-adaptor having no mismatch:                            \t" + str(counter['counter_h_perfect_3p']  + counter['counter_t_perfect_3p'])  + "\n\n")
		logger.write("\tNumber of 3'-adaptor having mismatch at the last one position:      \t" + str(counter['counter_h_last_1_mm']   + counter['counter_t_last_1_mm'])   + "\n")
		logger.write("\tNumber of 3'-adaptor having mismatch at all the last two position:  \t" + str(counter['counter_h_last_12_mm']  + counter['counter_t_last_12_mm'])  + "\n")
		logger.write("\tNumber of 3'-adaptor having mismatch at all the last three position:\t" + str(counter['counter_h_last_123_mm'] + counter['counter_t_last_123_mm']) + "\n\n")
		logger.write("\tNumber of 3'-adaptor having in/del at the last one position:        \t" + str(counter['counter_h_last_1_i']    + counter['counter_t_last_1_i'])    + "\n")
		logger.write("\tNumber of 3'-adaptor having in/del at the last two position:        \t" + str(counter['counter_h_last_2_i']    + counter['counter_t_last_2_i'])    + "\n")
		logger.write("\tNumber of 3'-adaptor having in/del at the last three position:      \t" + str(counter['counter_h_last_3_i']    + counter['counter_t_last_3_i'])    + "\n\n")
		logger.write("\tNumber of rescued truncated 3'-adaptor on the read head region: \t" + str(counter['counter_h_partial_3p']) + "\n")
		logger.write("\tNumber of rescued truncated 3'-adaptor on the read tail region: \t" + str(counter['counter_t_partial_3p']) + "\n\n")

		logger.write("Finish time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")

def printPSAlignment(alignment_obj):
	printAlignment(alignment_obj.seqA[alignment_obj.start:alignment_obj.end], alignment_obj.adaptor_start + alignment_obj.start - alignment_obj.gap_no, "seqA")
	printAlignment(alignment_obj.alignment,                                   alignment_obj.adaptor_start + alignment_obj.start - alignment_obj.gap_no, "aln")
	printAlignment(alignment_obj.seqB[alignment_obj.start:alignment_obj.end], alignment_obj.adaptor_start + alignment_obj.start - alignment_obj.gap_no, "seqB")

def printAlignment(alignment_seq, no_space, label):
	if len(label) > 6:
		print(label[0:6] + ": ", end = "")
	else:
		print(label + ": " + " " * (6 - len(label) - 1), end = " ")
	print(" " * no_space, end = "")
	print(alignment_seq)

def batch_reading(reader, f_list, f_idx, batch_no, rid):
	batch_data = []
	eof = 0

	if f_idx == len(f_list):
		eof = 1
		return batch_data, eof, rid, reader, f_idx
	if f_list[f_idx].endswith(".fast5"):
		# read fast5 at once disregarded to batch_no
		print("Processing file name: " + str(f_list[f_idx]))
		grp_keys = list(reader.keys())
		for my_key in grp_keys:
			rid += 1
			h5_content = list(reader[my_key]['Analyses']['Basecall_1D_000']['BaseCalled_template']['Fastq'][()].decode("utf-8").split("\n"))
			batch_data.append({'rid': rid, 'def_line': h5_content[0], 'na_seq': h5_content[1], 'qu_seq': h5_content[3]})

		f_idx = f_idx + 1
		if f_idx < len(f_list):
			reader.close()
			reader = open_file(f_list[f_idx], "rt")
		return batch_data, eof, rid, reader, f_idx
	else:
		i = 0
		while True:
			i += 1
			if i > batch_no:
				break
			# read fastq
			def_line = reader.readline().rstrip()
			na_seq   = reader.readline().rstrip()
			reader.readline().rstrip()
			qu_seq   = reader.readline().rstrip()

			if not def_line:
				#=== open next file ===
				if f_idx < len(f_list):
					f_idx += 1
					reader.close()

					#=== stop if no further file ===
					if f_idx == len(f_list):
						eof = 1
						return batch_data, eof, rid, reader, f_idx

					print("Processing file name: " + str(f_list[f_idx]))
					reader = open_file(f_list[f_idx], "rt")
					# re-do for loop
					i -= 1
					continue
			else:
				rid += 1
				batch_data.append({'rid': rid, 'def_line': def_line, 'na_seq': na_seq, 'qu_seq': qu_seq})

		return batch_data, eof, rid, reader, f_idx

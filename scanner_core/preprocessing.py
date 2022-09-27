def getOptions():
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("-i", dest = "fq_f_name",         default = None,
                          nargs = 1, type = "string",
                          help = "* Required ! "
                                 "Input FastQ/Fast5 file name, or directory containing multiple input files. "
                                 "Support fastq/fq/fastq.gz/fq.gz/fast5 format.")
	parser.add_option("-o", dest = "fq_o_name",         default = "processed.fastq.gz",
                          nargs = 1, type = "string",
                          help = "Output fastq file name. "
                                 "Could be either in .fastq or .gz format."
                                 "Default: processed.fastq.gz")
	parser.add_option("-d", dest = "o_dir",             default = "scNanoGPS_res",
                          nargs = 1, type = "string",
                          help = "Output directory name. Must give a new directory name to prevent accidental overwriting ! "
                                 "Default: scNanoGPS_res")
	parser.add_option("-b", dest = "bc_f_name",         default = "barcode_list.tsv.gz",
                          nargs = 1, type = "string",
                          help = "Output cell barcode list file. "
                                 "Default: barcode_list.tsv.gz")
	parser.add_option("-t", dest = "ncores",            default = 1,
                          nargs = 1, type = "int",
                          help = "Number of cores for program running. "
                                 "Default: 1")
	parser.add_option("--log", dest = "log_f_name",     default = "scanner.log.txt",
                          nargs = 1, type = "string",
                          help = "Program log file. "
                                 "This file stores program running parameters and counting details. "
                                 "Default: scanner.log.txt")
	parser.add_option("--a5", dest = "adaptor_five_p",  default = "AAGCAGTGGTATCAACGCAGAGTACAT",
                          nargs = 1, type = "string",
                          help = "Sequence of 5'-adaptor. "
                                 "Default: AAGCAGTGGTATCAACGCAGAGTACAT")
	parser.add_option("--a3", dest = "adaptor_three_p", default = "CTACACGACGCTCTTCCGATCT",
                          nargs = 1, type = "string",
                          help = "Sequence of 3'-adaptor. "
                                 "Default: CTACACGACGCTCTTCCGATCT")
	parser.add_option("--pT", dest = "polyT",           default = "TTTTTTTTTTTT",
                          nargs = 1, type = "string",
                          help = "Reverse-complement sequence of polyA. "
                                 "Default: TTTTTTTTTTTT")
	parser.add_option("--lCB", dest = "BC_len",         default = 16,
                          nargs = 1, type = "int",
                          help = "Length of cell barcode. "
                                 "Default: 16")
	parser.add_option("--lUMI", dest = "UMI_len",       default = 12,
                          nargs = 1, type = "int",
                          help = "Length of UMI. "
                                 "Default: 12")

	parser.add_option("--min_read_length", dest = "min_read_length",         default = 200,
                          nargs = 1, type = "int",
	                  help = "Minimal read length. "
	                         "Default: 200")
	parser.add_option("--editing_distance", dest = "allow_editing_distance", default = 2,
                          nargs = 1, type = "int",
                          help = "Editing distance for cell barcode detection. "
                                 "Default: 2")
	parser.add_option("--matching_threshold", dest = "matching_percentage",  default = 0.7,
                          nargs = 1, type = "float",
                          help = "Matching threshold for alignment search. "
                                 "Default: 0.7")
	parser.add_option("--score_threshold", dest = "scoring_threshold",       default = 0.4,
                          nargs = 1, type = "float",
                          help = "Scoring threshold for alignment search. "
                                 "Default: 0.4")
	parser.add_option("--batching_no", dest = "batch_no",                    default = 1000,
                          nargs = 1, type = "int",
                          help = "Number of reads for batch processing. "
                                 "Default: 1000")
	parser.add_option("--debug_mode", dest = "debug_mode",                   default = False,
                          nargs = 1, type = "string",
                          help = "Debug mode switch. "
                                 "Default: False")

	parser.add_option("--penalty_matching",      dest = "dp_ma", default = 2,
                          nargs = 1, type = "int",
                          help = "Dynamic programming matching penalty. "
                                 "Default: 2")
	parser.add_option("--penalty_mismatching",   dest = "dp_mi", default = -3,
                          nargs = 1, type = "int",
                          help = "Dynamic programming mismatching penalty. "
                                 "Default: -3")
	parser.add_option("--penalty_gap_opening",   dest = "dp_go", default = -5,
                          nargs = 1, type = "int",
                          help = "Dynamic programming gap opening penalty. "
                                 "Default: -5")
	parser.add_option("--penalty_gap_extention", dest = "dp_ge", default = -2,
                          nargs = 1, type = "int",
                          help = "Dynamic programming gap extention penalty. "
                                 "Default: -2")

	return parser

def precheck(parser, options, arguments):
	import os, sys, time, math
	from scanner_core.preprocessing import get_valid_input_file_list

	if not options.fq_f_name:
		print("\nInput Fastq file is required !\n")
		parser.print_help()
		sys.exit(1)
	if not options.fq_o_name:
		print("\nOutput Fastq file is required !\n")
		parser.print_help()
		sys.exit(1)
	if not os.path.isdir(options.o_dir):
		os.system("mkdir " + options.o_dir)

	options.fq_o_name  = os.path.join(options.o_dir, options.fq_o_name)
	options.bc_f_name  = os.path.join(options.o_dir, options.bc_f_name)
	options.log_f_name = os.path.join(options.o_dir, options.log_f_name)

	options.isFile = os.path.isfile(options.fq_f_name)
	options.isDir  = os.path.isdir(options.fq_f_name)
	options.cwd    = os.getcwd()
	if options.isDir:
		options.f_no = len(get_valid_input_file_list(options.cwd, options.fq_f_name))

		if options.f_no == 0:
			print("\nCannot find valid input file under directory: " + str(os.path.join(options.cwd, options.fq_f_name)) + "\n")
			parser.print_help()
			sys.exit(1)

	#===output parameters===
	with open(options.log_f_name, "wt") as logger:
		logger.write("Starting time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")
		logger.write("\n")

		logger.write("List of parameters:\n")
		logger.write("\tCurrent working directory:     " + str(options.cwd) + "\n")
		if options.isFile:
			logger.write("\tInput file name:               " + str(options.fq_f_name) + "\n")
		if options.isDir:
			logger.write("\tInput directory:               " + str(options.fq_f_name) + "\n")
			logger.write("\tValid file number:             " + str(options.f_no) + "\n")
		logger.write("\tOutput FastQ file name:        " + str(options.fq_o_name) + "\n")
		logger.write("\tOutput barcode list name:      " + str(options.bc_f_name) + "\n")
		logger.write("\tLog file name:                 " + str(options.log_f_name) + "\n")
		logger.write("\n")

		logger.write("Parameters for pattern search:\n")
		logger.write("\tLength of barcode:             " + str(options.BC_len) + "\n")
		logger.write("\tLength of UMI:                 " + str(options.UMI_len) + "\n")
		logger.write("\t5'-adaptor sequence:           " + str(options.adaptor_five_p) + "\n")
		logger.write("\t3'-adaptor sequence:           " + str(options.adaptor_three_p) + "\n")
		logger.write("\tPolyT sequence:                " + str(options.polyT) + "\n")
		logger.write("\n")

		logger.write("Penalty for dynamic programming:\n")
		logger.write("\tMatching:                      " + str(options.dp_ma) + "\n")
		logger.write("\tMismatching:                   " + str(options.dp_mi) + "\n")
		logger.write("\tGap opening:                   " + str(options.dp_go) + "\n")
		logger.write("\tGap extension:                 " + str(options.dp_ge) + "\n")
		logger.write("\tEditing distance:              " + str(options.allow_editing_distance) + "\n")
		logger.write("\n")

		logger.write("Parameters for computing:\n")
		logger.write("\tNumber of computer cores:      " + str(options.ncores) + "\n")
		logger.write("\tNumber of reads per batch job: " + str(options.batch_no) + "\n")
		logger.write("\tMinimal length of read:        " + str(options.min_read_length) + "\n")
		logger.write("\tMatching threshold:            " + str(options.matching_percentage) + "\n")
		logger.write("\tScoring threshold:             " + str(options.scoring_threshold) + "\n")
		logger.write("\n")

		logger.write("Debug mode switch:             " + str(options.debug_mode) + "\n")
		logger.write("\n")

	#===add 3' half seq of adaptors===
	options.adaptor_five_h  = options.adaptor_five_p[math.floor(len(options.adaptor_five_p) / 2):len(options.adaptor_five_p)]
	options.adaptor_three_h = options.adaptor_three_p[math.floor(len(options.adaptor_three_p) / 2):len(options.adaptor_three_p)]
	options.dp_penalty = [options.dp_ma, options.dp_mi, options.dp_go, options.dp_ge]

	counter = {
		'counter_h_3p':  0,
		'counter_h_3p_polyT': 0,
		'counter_t_3p':  0,
		'counter_t_3p_polyT':  0,
		'counter_ht_3p': 0,
	#---check alignment---
		'counter_h_last_1_mm':   0,
		'counter_h_last_12_mm':  0,
		'counter_h_last_123_mm': 0,
		'counter_t_last_1_mm':   0,
		'counter_t_last_12_mm':  0,
		'counter_t_last_123_mm': 0,
		'counter_h_last_1_i':    0,
		'counter_h_last_2_i':    0,
		'counter_h_last_3_i':    0,
		'counter_t_last_1_i':    0,
		'counter_t_last_2_i':    0,
		'counter_t_last_3_i':    0,
		'counter_h_partial_3p':  0,
		'counter_t_partial_3p':  0,
		'counter_h_perfect_3p':  0,
		'counter_t_perfect_3p':  0}

	return parser, options, arguments, counter

def get_valid_input_file_list(cwd, my_path):
	import os

	f_list = list()
	for f in os.listdir(os.path.join(cwd, my_path)):
		if os.path.isfile(os.path.join(cwd, my_path, f)) and \
		   (f.endswith(".fastq") or f.endswith(".fastq.gz") or
		    f.endswith(".fq") or f.endswith(".fq.gz") or
		    f.endswith(".fast5")):
			f_list.append(os.path.join(cwd, my_path, f))
	return f_list

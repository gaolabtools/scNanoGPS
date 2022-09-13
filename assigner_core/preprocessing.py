def getOptions():
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("-i", dest = "input", nargs = 1, default = "barcode_list.tsv.gz",
                	  help = "Cell barcode list file. Could be either in .tsv or .tsv.gz format. "
                                 "Default: barcode_list.tsv.gz")
	parser.add_option("-o", dest = "output", nargs = 1, default = "CB_counting.tsv.gz",
        	          help = "Cell barcode counting file name. "
	                         "Default: CB_counting.tsv.gz")
	parser.add_option("-d", dest = "o_dir",  nargs = 1, default = "scNanoGPS_res",
                          help = "Output directory name. "
                                 "Default: scNanoGPS_res")
	parser.add_option("-t", dest = "ncores", nargs = 1, default = 1,
                          help = "Number of cores for program running. "
                                 "Default: 1")
	parser.add_option("--log", dest = "log_f_name", nargs = 1, default = "assigner.log.txt",
                          help = "Log file name. "
                                 "Default: assigner.log.txt")
	parser.add_option("--lCB", dest = "BC_len",         default = 16,
                          nargs = 1, type = "int",
                          help = "Length of cell barcode. "
                                 "Default: 16")
	parser.add_option("--CB_no_ext", dest = "CB_no_ext", nargs = 1, default = 0.1,
	                  help = "Increasing CB number on log scale. "
	                         "Default: 0.1")
	parser.add_option("--CB_log10_dist_o", dest = "CB_log10_dist_o", nargs = 1, default = "CB_log10_dist.png",
                          help = "File name used for plotting log10 UMI number distribution. "
                                 "Default: CB_log10_dist.png")
	parser.add_option("--CB_mrg_thr",  dest = "CB_mrg_thr",  nargs = 1, default = 2,
	                  help = "Threshold of distance for merging cell barcodes. "
	                         "Default: 2")
	parser.add_option("--CB_mrg_dist", dest = "CB_mrg_dist", nargs = 1, default = "CB_merged_dist.tsv.gz",
	                  help = "File name for distance matrix of merging cell barcodes. "
	                         "Default: CB_merged_dist.tsv.gz")
	parser.add_option("--CB_mrg_o",    dest = "CB_mrg_o",    nargs = 1, default = "CB_merged_list.tsv.gz",
	                  help = "File name for merged cell barcodes. "
	                         "Default: CB_merged_list.tsv.gz")

	return parser

def precheck(parser, options):
	import os, sys, time

	options.input           = os.path.join(options.o_dir, options.input)
	options.output          = os.path.join(options.o_dir, options.output)
	options.CB_log10_dist_o = os.path.join(options.o_dir, options.CB_log10_dist_o)
	options.CB_mrg_dist     = os.path.join(options.o_dir, options.CB_mrg_dist)
	options.CB_mrg_o        = os.path.join(options.o_dir, options.CB_mrg_o)
	options.log_f_name      = os.path.join(options.o_dir, options.log_f_name)

	if not options.input:
		print("\nCannot find Cell barcode list file: " + options.input + " !\n")
		parser.print_help()
		sys.exit(1)

	if not os.path.isfile(options.input):
		print("\nCannot find Cell barcode list file: " + options.input + " !\n")
		parser.print_help()
		sys.exit(1)

	if not os.path.isdir(options.o_dir):
		print("\nCannot find output dir: " + options.o_dir + " !\n")
		parser.print_help()
		sys.exit(1)

	options.i_compression           = None
	options.o_compression           = None
	options.CB_mrg_dist_compression = None

	if options.input.endswith(".gz"):
		options.i_compression = 'gzip'
	if options.output.endswith(".gz"):
		options.o_compression = 'gzip'
	if options.CB_mrg_dist.endswith(".gz"):
		options.CB_mrg_dist_compression = 'gzip'

	#===output parameters===
	with open(options.log_f_name, "wt") as logger:
		logger.write("Starting time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")
		logger.write("\n")

		logger.write("\tAnalysis result directory:            " + options.o_dir + "\n")
		logger.write("\tCell barcode list file:               " + options.input + "\n")
		logger.write("\tUMI counting file name:               " + options.output + "\n")
		logger.write("\tLog10(UMI) distribution figure name:  " + options.CB_log10_dist_o + "\n")
		logger.write("\tDistance for merging cell barcodes:   " + str(options.CB_mrg_thr) + "\n")
		logger.write("\tMerged cell barcodes distance matrix: " + options.CB_mrg_dist + "\n")
		logger.write("\tMerged cell barcodes file name:       " + options.CB_mrg_o + "\n")
		logger.write("\tLength of barcode:                    " + str(options.BC_len) + "\n")
		logger.write("\tNumber of computer cores:             " + str(options.ncores) + "\n")
		logger.write("\n")

	return parser, options

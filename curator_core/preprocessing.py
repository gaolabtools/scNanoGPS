def getOptions():
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("--fq_name",  dest = "fq_name",  nargs = 1, default = "processed.fastq.gz",
                          help = "Processed fastq file name. "
                                 "Default: processed.fastq.gz")
	parser.add_option("-b",         dest = "BC_list", nargs = 1, default = "barcode_list.tsv.gz",
                          help = "Output cell barcode list file. "
                                 "Default: barcode_list.tsv.gz")
	parser.add_option("--CB_count", dest = "CB_count", nargs = 1, default = "CB_counting.tsv.gz",
                          help = "Cell barcode counting file name. "
                                 "Default: CB_counting.tsv.gz")
	parser.add_option("--CB_list",  dest = "CB_list",  nargs = 1, default = "CB_merged_list.tsv.gz",
                          help = "File name for merged cell barcodes. "
                                 "Default: CB_merged_list.tsv.gz")
	parser.add_option("--ref_genome", dest = "ref_genome", nargs = 1, default = None,
                          help = "* Required ! "
                                 "File for reference genome.")
	parser.add_option("--idx_genome", dest = "idx_genome", nargs = 1, default = None,
                          help = "Path to the Minimap2 genome index. "
                                 "Program will use reference genome if no Minimap2 genome index given. "
                                 "Default: None")
	parser.add_option("-d",         dest = "o_dir",   nargs = 1, default = "scNanoGPS_res",
                          help = "Output directory name. "
                                 "Default: scNanoGPS_res")
	parser.add_option("--tmp_dir",  dest = "tmp_dir", nargs = 1, default = "tmp",
                          help = "Temporary folder name. "
                                 "Default: tmp")
	parser.add_option("-t",         dest = "ncores", nargs = 1, default = 1,
                          help = "Number of cores for computing. "
                                 "Default: 1")
	parser.add_option("--log",      dest = "log_f_name", nargs = 1, default = "curator.log.txt",
                          help = "Log file name. "
                                 "Default: curator.log.txt")
	parser.add_option("--umi_ld",   dest = "umi_ld", nargs = 1, default = 2,
                          help = "Levenshtein distance for merging UMI. "
                                 "Default: 2")
	parser.add_option("--keep_meta", dest = "keep_meta", nargs = 1, default = None,
                          help = "Keep meta data, e.g. bam files, for futher checking. "
                                 "Default: None")
	parser.add_option("--softclipping_thr", dest = "softclipping_thr", nargs = 1, default = 0.8,
                          help = "Threshold for softclipping. "
                                 "Default: 0.8")
	parser.add_option("--minimap2", dest = "minimap2", nargs = 1, default = "minimap2",
                          help = "Path to minimap2. "
                                 "Default: minimap2")
	parser.add_option("--samtools", dest = "samtools", nargs = 1, default = "samtools",
                          help = "Path to samtools. "
                                 "Default: samtools")
	parser.add_option("--spoa",     dest = "spoa", nargs = 1, default = "spoa",
                          help = "Path to spoa. "
                                 "Default: spoa")
	return parser

def precheck(parser, options):
	import sys, os, time
	from curator_core import curator_io

	options.fq_name    = os.path.join(options.o_dir, options.fq_name)
	options.BC_list    = os.path.join(options.o_dir, options.BC_list)
	options.CB_count   = os.path.join(options.o_dir, options.CB_count)
	options.CB_list    = os.path.join(options.o_dir, options.CB_list)
	options.log_f_name = os.path.join(options.o_dir, options.log_f_name)

	termination = False
	if not options.fq_name:
		print("\nProcessed FastQ file is required !\n")
		termination = True
	if not os.path.isfile(options.fq_name):
		print("\nCannot find processed FastQ file: " + options.fq_name + "\n")
		termination = True
	if not os.path.isfile(options.CB_list):
		print("\nCannot find Cell Barcode merged list: " + options.CB_list + "\n")
		termination = True
	if not os.path.isfile(options.CB_count):
		print("\nCannot find Cell Barcode counting file: " + options. CB_count + "\n")
		termination = True
	if not os.path.isfile(options.BC_list):
		print("\nCannot find Barcode list: " + options.BC_list + "\n")
		termination = True
	if os.path.isdir(options.tmp_dir):
		print("\nTemporary directory exist: " + options.tmp_dir)
		termination = True
	if not options.ref_genome:
		print("\nReference genome for minimap2 is required !\n")
		termination = True
	if not curator_io.find_exe(options.minimap2):
		print("\nCannot find path to minimap2 !\n")
		termination = True
	if not curator_io.find_exe(options.samtools):
		print("\nCannot find path to samtools !\n")
		termination = True
	if not curator_io.find_exe(options.spoa):
		print("\nCannot find path to spoa !\n")
		termination = True

	if termination:
		parser.print_help()
		sys.exit(1)

	print("\nCreate new temporary folder: " + options.tmp_dir + "\n")
	os.makedirs(options.tmp_dir)

	#===output parameters===
	with open(options.log_f_name, "wt") as logger:
		logger.write("Starting time stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")
		logger.write("List of parameters:\n")
		logger.write("\tInput fastQ file:           " + str(options.fq_name) + "\n")
		logger.write("\tCell barcode list:          " + str(options.BC_list) + "\n")
		logger.write("\tCell barcode counting file: " + str(options.CB_count) + "\n")
		logger.write("\tMerged cell barcode list:   " + str(options.CB_list) + "\n")
		logger.write("\tReference genome:           " + str(options.ref_genome) + "\n")
		logger.write("\tMinimap2 genome index:      " + str(options.idx_genome) + "\n")
		logger.write("\tOutput directory:           " + str(options.o_dir) + "\n")
		logger.write("\tTemporary directory:        " + str(options.tmp_dir) + "\n")
		logger.write("\tNumber of computer cores:   " + str(options.ncores) + "\n")
		logger.write("\tLog file name:              " + str(options.log_f_name) + "\n")
		logger.write("\tThreshold for softclipping: " + str(options.softclipping_thr) + "\n")
		logger.write("\tPath of minimap2:           " + str(options.minimap2) + "\n")
		logger.write("\tPath of samtools:           " + str(options.samtools) + "\n")
		logger.write("\tPath spoa:                  " + str(options.spoa) + "\n")
		logger.write("\n")

	return parser, options

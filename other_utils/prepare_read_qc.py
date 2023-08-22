#! /usr/bin/env python3

import os, sys, gzip
from optparse import OptionParser
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../'))
from scanner_core import preprocessing, scanner_io

def getTailna(na_seq, options):
	return na_seq[-options.read_len:len(na_seq)].translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))[::-1]

def getTailqu(qu_seq, options):
	return qu_seq[-options.read_len:len(qu_seq)][::-1]

#===get params===
parser = OptionParser()

parser.add_option("-i",      dest = "fq_f_name", nargs = 1, default = None,
                  help = "File or folder name of reads. ")
parser.add_option("-d",      dest = "o_dir",     nargs = 1, default = "scNanoGPS_res",
                  help = "Output directory name. "
                         "Default: scNanoGPS_res")
parser.add_option("-l",      dest = "read_len",  nargs = 1, default = 100, type = "int",
                  help = "Length of the extracted first and last read nucleotides. "
                         "Default: 100")
parser.add_option("--o1",    dest = "o1",        nargs = 1, default = "first_tail.fastq.gz",
                  help = "First given length of reads. "
                         "Default: first_tail.fastq.gz")
parser.add_option("--o2",    dest = "o2",        nargs = 1, default = "last_tail.fastq.gz",
                  help = "Last given length of reads. "
                         "Default: last_tail.fastq.gz")

options, arguments = parser.parse_args()

#===pre-check===
termination = False
if not os.path.isdir(options.o_dir):
	os.system("mkdir " + options.o_dir)
if not options.fq_f_name:
	print("\nPlease give valid file or folder name: " + str(options.fq_f_name) + "\n")
	termination = True
if termination:
	parser.print_help()
	sys.exit(1)

#===precheck===
options.isFile   = os.path.isfile(options.fq_f_name)
options.isDir    = os.path.isdir(options.fq_f_name)
options.cwd      = os.getcwd()
options.batch_no = 4000
if options.isDir:
	options.f_no = len(preprocessing.get_valid_input_file_list(options.cwd, options.fq_f_name))
	if options.f_no == 0:
		print("\nCannot find valid input file under directory: " + str(os.path.join(options.cwd, options.fq_f_name)) + "\n")
		parser.print_help()
		sys.exit(1)

#===count read length===
if options.isFile:
	f_list = list([options.fq_f_name])
	f_idx  = 0
	print("Parsing fastQ file: " + str(options.fq_f_name) + " ...")
	reader = scanner_io.open_file(options.fq_f_name, "rt")

if options.isDir:
	f_list = preprocessing.get_valid_input_file_list(options.cwd, options.fq_f_name)
	f_idx  = 0
	if not f_list[f_idx].endswith(".fast5"):
		print("Processing file name: " + str(f_list[f_idx]))
	reader = scanner_io.open_file(os.path.join(options.cwd, options.fq_f_name, f_list[f_idx]), "rt")

rid = 0
while True:
	batch_data = []
	batch_data, eof, rid, reader, f_idx = scanner_io.batch_reading(reader, f_list, f_idx, int(options.batch_no), rid)
	if eof == 1:
		break

	o1h = gzip.open(os.path.join(options.o_dir, options.o1), "at")
	o2h = gzip.open(os.path.join(options.o_dir, options.o2), "at")
	for read in batch_data:
		o1h.write("@" + read['def_line'] + "\n")
		o1h.write(read['na_seq'][0:options.read_len] + "\n")
		o1h.write("+" + read['def_line'] + "\n")
		o1h.write(read['qu_seq'][0:options.read_len] + "\n")

		o2h.write("@" + read['def_line'] + "\n")
		o2h.write(getTailna(read['na_seq'], options) + "\n")
		o2h.write("+" + read['def_line'] + "\n")
		o2h.write(getTailqu(read['qu_seq'], options) + "\n")

	batch_data = []
	o1h.close()
	o2h.close()
print("\nDone.")


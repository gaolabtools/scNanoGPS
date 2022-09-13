#! /usr/bin/env python3

import os, sys, gzip
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../'))
from scanner_core import preprocessing, scanner_io

#===get params===
parser = OptionParser()

parser.add_option("-i",      dest = "fq_f_name", nargs = 1, default = None,
                  help = "File or folder name of reads. ")
parser.add_option("-d",           dest = "o_dir",      nargs = 1, default = "scNanoGPS_res",
                  help = "Output directory name. "
                         "Default: scNanoGPS_res")
parser.add_option("-f",      dest = "fig_name",  nargs = 1, default = "read_length.png",
                  help = "Read length histogram file name. "
                         "Default: read_length.png")
parser.add_option("-o",      dest = "o_name",    nargs = 1, default = "read_length.tsv.gz",
                  help = "Read length file name. "
                         "Default: read_length.tsv.gz")
parser.add_option("--fig_w", dest = "fig_w",     nargs = 1, default = 12,
                  help = "Width of figure (inch). "
                         "Default: 12")
parser.add_option("--fig_h", dest = "fig_h",     nargs = 1, default = 7,
                  help = "Height of figure (inch). "
                         "Default: 7")

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
seq_len_list = []
while True:
	batch_data = []
	batch_data, eof, rid, reader, f_idx = scanner_io.batch_reading(reader, f_list, f_idx, int(options.batch_no), rid)
	if eof == 1:
		break

	for read in batch_data:
		seq_len_list.append(len(read['na_seq']))

	batch_data = []

print("\nDone.")

sns.set_theme(style = "white")
plt.rcParams["figure.figsize"]    = [options.fig_w, options.fig_h]
plt.rcParams['figure.dpi']        = 300
plt.rcParams["figure.autolayout"] = True
p = sns.histplot(seq_len_list)
p.set(xlim = (0, 5000))
p.set_xlabel("Read length")
plt.show()
plt.savefig(os.path.join(options.o_dir, options.fig_name))
plt.figure().clear()

oh = gzip.open(os.path.join(options.o_dir, options.o_name), "wt")
for rlen in seq_len_list:
	oh.write(str(rlen) + "\n")
oh.close()

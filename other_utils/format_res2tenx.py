#! /usr/bin/env python3

import os, sys
import gzip
from optparse import OptionParser

#===get params===
parser = OptionParser()

parser.add_option("-i", dest = "i_name", nargs = 1, default = None,
                  help = "scNanoGPS matrix table file name.")
parser.add_option("-o", dest = "o_dir",  nargs = 1, default = None,
                  help = "Output folder name.")

parser.add_option("-b", dest = "b_name", nargs = 1, default = "barcodes.tsv.gz",
                  help = "Cell barcode file name.")
parser.add_option("-f", dest = "f_name", nargs = 1, default = "features.tsv.gz",
                  help = "Feature file name.")
parser.add_option("-m", dest = "m_name", nargs = 1, default = "matrix.mtx.gz",
                  help = "Matrix file name.")

options, arguments = parser.parse_args()

#===pre-check===
termination = False
if not options.i_name or \
   not os.path.isfile(options.i_name):
	print("\nCannot find scNanoGPS matrix file: " + str(options.i_name) + "\n")
	termination = True
if os.path.isdir(options.o_dir):
	print("\nOutput directory exist: " + options.o_dir   + "\n")
	print("\nPlease assign new folder name to prevent accidental over-writing.")
	termination = True

if termination:
	parser.print_help()
	sys.exit(1)

print("Create output folder: " + options.o_dir)
os.system("mkdir " + options.o_dir)

options.b_name = os.path.join(options.o_dir, options.b_name)
options.f_name = os.path.join(options.o_dir, options.f_name)
options.m_name = os.path.join(options.o_dir, options.m_name)

print("\nConverting...", end = "")
#===parse data===
ih = open(options.i_name, "rt")

#---skip command---
ih.readline()

#---parse header---
header_list = ih.readline().rstrip().split("\t")

bh = gzip.open(options.b_name, "wt")
for bc_idx, bc in enumerate(header_list[7::]):
	bc_list = bc.split('/')
	bc_str = bc_list[len(bc_list) - 1].split('.')[0]
	bh.write(bc_str + "\n")
bh.close()

tmp_fname = os.path.join(options.o_dir, "tmp")
th = open(tmp_fname, "wt")
fh = gzip.open(options.f_name, "wt")

gene_idx = 0
count_idx = 0
bc_idx = 0
while True:
	line = ih.readline().rstrip()
	if not line:
		break
	gene_idx += 1
	line_list = line.split("\t")

	fh.write(line_list[0] + "\t" + line_list[6] + "\tGene Expression\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\n")

	for bc_idx, count in enumerate(line_list[7::]):
		count_idx += 1
		th.write(str(gene_idx) + "\t" + str(bc_idx + 1) + "\t" + str(count) + "\n")
th.close()
fh.close()
ih.close()

mh = gzip.open(options.m_name, "wt")
mh.write("%%MatrixMarket matrix coordinate integer general\n")
mh.write("%metadata_json: {\"software_version\": \"scNanoGPS\", \"format_version\": 2}\n")
mh.write(str(gene_idx) + "\t" + str(bc_idx + 1) + "\t" + str(count_idx) + "\n")
th = open(tmp_fname, "rt")
while True:
	line = th.readline()
	if not line:
		break
	mh.write(line)
th.close()
mh.close()

os.system("rm " + tmp_fname)

print("Done.\n")

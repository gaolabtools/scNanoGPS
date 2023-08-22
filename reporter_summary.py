#! /usr/bin/env python3

import time, os, sys, gzip, re, subprocess, glob
import numpy as np
import pandas as pd
from optparse import OptionParser
from curator_core import curator_io

#===get params===
parser = OptionParser()

parser.add_option("-d",            dest = "o_dir",       nargs = 1, default = "scNanoGPS_res",
                  help = "Output directory name. "
                         "Default: scNanoGPS_res")
parser.add_option("--tmp_dir",     dest = "tmp_dir",     nargs = 1, default = "tmp",
                  help = "Temporary folder name. "
                         "Default: tmp")
parser.add_option("--scanner_log", dest = "scanner_log", nargs = 1, default = "scanner.log.txt",
                  help = "Scanner log file name. "
                         "Default: scanner.log.txt")
parser.add_option("--bc_f",        dest = "bc_f",        nargs = 1, default = "barcode_list.tsv.gz",
                  help = "Cell barcode list file. "
                         "Default: barcode_list.tsv.gz")
parser.add_option("--read_len_f",  dest = "read_len_f",  nargs = 1, default = "read_length.tsv.gz",
                  help = "Scanner log file name. "
                         "Default: read_length.tsv.gz")
parser.add_option("--CB_file",     dest = "CB_file",     nargs = 1, default = "filtered_barcode_list.txt",
                  help = "File name for filtered barcode list. "
                         "Default: filtered_barcode_list.txt")
parser.add_option("--exp_tb",      dest = "exp_tb",      nargs = 1, default = "matrix.tsv",
                  help = "Counting table name. "
                         "Default: matrix.tsv")
parser.add_option("--ref_genome",  dest = "ref_genome",  nargs = 1, default = None,
                  help = "* Required ! "
                         "File for reference genome.")
parser.add_option("--gtf",         dest = "gtf",         nargs = 1, default = None,
                  help = "* Required ! "
                         "Genome annotation file GTF.")
parser.add_option("--log",         dest = "log_f_name",  nargs = 1, default = "summary.txt",
                  help = "Log file name. "
                         "Default: summary.txt")
parser.add_option("--samtools",    dest = "samtools",    nargs = 1, default = "samtools",
                  help = "Path to samtools. "
                         "Default: samtools")
parser.add_option("--qualimap",    dest = "qualimap",    nargs = 1, default = "qualimap",
                  help = "Path to qualimap. "
                         "Default: qualimap")
options, arguments = parser.parse_args()

#===pre-check===
options.exp_tb = os.path.join(options.o_dir, options.exp_tb)
options.compression = None
if options.exp_tb.endswith('.gz'):
	options.compression = 'gzip'
options.log_f_name = os.path.join(options.o_dir, options.log_f_name)

termination = False
if not options.ref_genome or \
   not os.path.isfile(options.ref_genome):
	print("\nCannot find reference genome file: " + str(options.ref_genome) + "\n")
	termination = True
if not options.gtf or \
   not os.path.isfile(options.gtf):
	print("\nCannot find genome annotation file: " + str(options.gtf) + "\n")
	termination = True
if not os.path.isdir(options.o_dir):
	print("\nOutput directory is not exist: "    + options.o_dir   + "\n")
	termination = True
if not os.path.isdir(options.tmp_dir):
	print("\nTemporary directory is not exist: " + options.tmp_dir + "\n")
	print(options.tmp_dir)
	termination = True
if not os.path.isfile(os.path.join(options.o_dir, options.scanner_log)):
	print("\nCannot find scanner log file at: " + str(os.path.join(options.o_dir, options.scanner_log)))
	termination = True
if not os.path.isfile(os.path.join(options.o_dir, options.bc_f)):
	print("\nCannot find barcode list file at: " + str(os.path.join(options.o_dir, options.bc_f)))
	termination = True
if not os.path.isfile(os.path.join(options.o_dir, options.read_len_f)):
	print("\nCannot find read length file at: " + str(os.path.join(options.o_dir, options.read_len_f)))
	print("Please run read_length_profiler.py first\n")
	termination = True
if not os.path.isfile(os.path.join(options.o_dir, options.CB_file)):
	options.CB_file = None
else:
	options.CB_file = os.path.join(options.o_dir, options.CB_file)

cmd = "which " + options.samtools
code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
options.samtools = out_msg.decode("utf-8").rstrip()
if not os.path.isfile(options.samtools):
	print("\nCannot find samtools: " + options.samtools)
	termination = True

cmd = "which " + options.qualimap
code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
options.qualimap = out_msg.decode("utf-8").rstrip()
if not os.path.isfile(options.qualimap):
	print("\nCannot find qualimap: " + options.qualimap)
	termination = True

if termination:
	parser.print_help()
	sys.exit(1)

#===parse scanner log===
print("\n\nParsing scanner log ...", flush = True)
fh = open(os.path.join(options.o_dir, options.scanner_log), "rt")
total_read_no, pass_read_no, detecting_rate = 0, 0, 0
while True:
	line = fh.readline()
	if not line:
		break

	match = re.match("Total (\d+) reads are processed.", line)
	if match:
		total_read_no = int(match[1])

	match = re.match("Detecting rate: ([\d\.]+)\%", line)
	if match:
		detecting_rate = float(match[1])

	match = re.match("\tNumber of 3\'-adaptor located on the read head region:\s+(\d+)", line)
	if match:
		pass_read_no  = int(match[1])

	match = re.match("\tNumber of 3\'-adaptor located on the read tail region:\s+(\d+)", line)
	if match:
		pass_read_no += int(match[1])
fh.close()
print("Done.\n", flush = True)

#===calc mean read length===
print("Computing read length ...", flush = True)
df = pd.read_csv(os.path.join(options.o_dir, options.read_len_f), header = None, sep = '\t', compression = 'gzip')

len_median = df.iloc[:, 0].median()
len_mean   = round(df.iloc[:, 0].mean(), 2)
len_max    = df.iloc[:, 0].max()
print("Done.\n", flush = True)

#===calc quality===
print("Calculating quality score ...", flush = True)
df = pd.read_csv(os.path.join(options.o_dir, options.bc_f), header = 0, sep = '\t', compression = 'gzip')
qua_median = round(df['mean_BC_quality'].median(), 2)
qua_mean   = round(df['mean_BC_quality'].mean(), 2)
print("Done.\n", flush = True)

#===loading CB list===
bam_list, CB_list = [], []
if not options.CB_file:
	bam_list = glob.glob(os.path.join(options.tmp_dir, "*.curated.minimap2.bam"))
	for bam_f in bam_list:
		split_f_name = bam_f.split('/')
		BC = split_f_name[len(split_f_name) - 1].split('.')
		CB_list.append(BC[0])
else:
	with open(options.CB_file, "rt") as CBF:
		while True:
			CB_name = CBF.readline().rstrip()
			if not CB_name:
				break
			CB_list.append(CB_name)
	bam_list = [x + ".curated.minimap2.bam" for x in CB_list]

#===merge bam files===
print("\nMerging all bam file for qualimap...\n", flush = True)
bam_counter, umi_per_cell = 0, list()
for bam_f in bam_list:
	bam_counter += 1
	print(str(bam_counter) + " of " + str(len(bam_list)) + " files...", end = "\r", flush = True)
	cmd = options.samtools + " view " + os.path.join(options.tmp_dir, bam_f) + \
	      " >> " + os.path.join(options.tmp_dir, "master.sam")
	os.system(cmd)

	#===calc median UMI no per cell==
	cmd = options.samtools + " view " + os.path.join(options.tmp_dir, bam_f) + \
	      ' | cut -f1 | cut -d \'_\' -f1 | sort -u | wc -l'
	proc = subprocess.Popen(cmd, shell = True,
	       stdout = subprocess.PIPE,
	       stderr = subprocess.PIPE)
	out_msg, err_msg = proc.communicate()
	umi_per_cell.append(int(out_msg.decode("utf-8").rstrip()))

print()
cmd = options.samtools + " view -Sb " + os.path.join(options.tmp_dir, "master.sam") + \
      " -T " + options.ref_genome + \
      " -o " + os.path.join(options.tmp_dir, "master.unsorted.bam")
os.system(cmd)

cmd = "rm " + os.path.join(options.tmp_dir, "master.sam")
os.system(cmd)

cmd = options.samtools + " sort " + os.path.join(options.tmp_dir, "master.unsorted.bam") + \
      " -o " + os.path.join(options.tmp_dir, "master.bam")
os.system(cmd)

cmd = "rm " + os.path.join(options.tmp_dir, "master.unsorted.bam")
os.system(cmd)

cmd = options.samtools + " index " + os.path.join(options.tmp_dir, "master.bam")
os.system(cmd)
print("Done.\n", flush = True)

#===count confidently mapped reads===
cmd = options.samtools + " view " + os.path.join(options.tmp_dir, "master.bam") + \
      ' | cut -f1 | cut -d \'_\' -f1 | uniq | sort -u | wc -l'
proc = subprocess.Popen(cmd, shell = True,
       stdout = subprocess.PIPE,
       stderr = subprocess.PIPE)
out_msg, err_msg = proc.communicate()

UMI_no = int(out_msg.decode("utf-8").rstrip())

#===count mean & median expressed gene number===
expr_df = pd.read_csv(options.exp_tb, header = 0, sep = '\t', skiprows =1, compression = options.compression)
expr_list = list()
for i in range(7, expr_df.shape[1]):
	expr_list.append(sum(expr_df.iloc[:, i] > 0))
expr_mean = np.mean(expr_list)
expr_med  = np.median(expr_list)

#===count median UMI number per cell===
#cmd = options.samtools + " view " + os.path.join(options.tmp_dir, "master.bam") + \
#      ' | cut -f1 | cut -d \'_\' -f2'
#proc = subprocess.Popen(cmd, shell = True,
#       stdout = subprocess.PIPE,
#       stderr = subprocess.PIPE)
#out_msg, err_msg = proc.communicate()

#UMI_no_list = out_msg.decode("utf-8").rstrip().split("\n")
#unique, counts = np.unique(UMI_no_list, return_counts = True)
median_umi_no = round(np.median(umi_per_cell), 2)

#===qualimap===
print("Calling Qualimap...", flush = True)
cmd = options.qualimap + " rnaseq -bam " + \
      os.path.join(options.tmp_dir, "master.bam") + " -gtf " + options.gtf
os.system(cmd)

exonic_ratio, intronic_ratio, intergenic_ratio, overlapping_exon_ratio = 0, 0, 0, 0
fh = open(os.path.join(options.tmp_dir, "master_rnaseq_qc", "rnaseq_qc_results.txt"), "rt")
while True:
	line = fh.readline()
	if not line:
		break

	match = re.match("\s+exonic \=\s+[\d,]+\s\(([\d\.]+)\%\)", line)
	if match:
		exonic_ratio = float(match[1])

	match = re.match("\s+intronic \=\s+[\d,]+\s\(([\d\.]+)\%\)", line)
	if match:
		intronic_ratio = float(match[1])

	match = re.match("\s+intergenic \=\s+[\d,]+\s\(([\d\.]+)\%\)", line)
	if match:
		intergenic_ratio = float(match[1])

#	match = re.match("\s+overlapping exon \=\s+[\d,]+\s\(([\d\.]+)\%\)", line)
#	if match:
#		overlapping_exon_ratio = float(match[1])
fh.close()
print("Done.\n", flush = True)

#===Init. file handler, start report===
print("Generating summary...", flush = True)
logger = open(options.log_f_name, "wt")

logger.write("Read yield:                  " + str(total_read_no) + "\n")
logger.write("Valid read number:           " + str(pass_read_no) + "\n")
logger.write("Detecting rate:              " + str(detecting_rate) + '%' + "\n")
logger.write("\n")
logger.write("Median read length:          " + str(len_median) + "\n")
logger.write("Mean read length:            " + str(len_mean) + "\n")
logger.write("Maximal read length:         " + str(len_max) + "\n")
logger.write("Median cell barcode quality: " + str(qua_median) + "\n")
logger.write("Mean cell barcode quality:   " + str(qua_mean) + "\n")
logger.write("\n")
logger.write("Cell number:                 " + str(len(CB_list)) + "\n")
logger.write("Raw reads per cell:          " + str(round(total_read_no / len(CB_list), 2)) + "\n")
logger.write("UMI counts:                  " + str(UMI_no) + "\n")
logger.write("Median UMI counts per cell:  " + str(median_umi_no) + "\n")
logger.write("Mean UMI counts per cell:    " + str(round(UMI_no / len(CB_list), 2)) + "\n")
logger.write("Median gene number:          " + str(expr_med) + "\n")
logger.write("Mean gene number:            " + str(round(expr_mean, 2)) + "\n")
logger.write("\n")
logger.write("Exonic:                      " + str(exonic_ratio) + '%' + "\n")
logger.write("Intronic:                    " + str(intronic_ratio) + '%' + "\n")
logger.write("Intergenic:                  " + str(intergenic_ratio) + '%' + "\n")

logger.close()
print("Done.\n", flush = True)

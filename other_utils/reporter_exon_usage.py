#! /usr/bin/env python3

import time, os, sys, glob
from optparse import OptionParser
from scanner_core import misc
from curator_core import curator_io

#===get params===
parser = OptionParser()

parser.add_option("-d",              dest = "o_dir",         nargs = 1, default = "scNanoGPS_res",
                  help = "Output directory name. "
                         "Default: scNanoGPS_res")
parser.add_option("--tmp_dir",       dest = "tmp_dir",       nargs = 1, default = "tmp",
                  help = "Temporary folder name. "
                         "Default: tmp")
parser.add_option("--gtf",           dest = "gtf",           nargs = 1, default = None,
                  help = "* Required ! "
                         "GTF file for expression calling. ")
parser.add_option("-o",              dest = "o_name",        nargs = 1, default = "matrix_exon_usage.tsv",
                  help = "Counting table name. "
                         "Default: matrix_exon_usage.tsv")
parser.add_option("--log",           dest = "log_f_name",    nargs = 1, default = "reporter_exon_usage.log.txt",
                  help = "Log file name. "
                         "Default: reporter_exon_usage.log.txt")
parser.add_option("-t",              dest = "ncores",        nargs = 1, default = 1,
                  help = "Number of cores for program running. "
                         "Default: 1")
parser.add_option("--featurecounts", dest = "featurecounts", nargs = 1, default = "featureCounts",
                  help = "Path to featureCounts. "
                         "Default: featureCounts")
parser.add_option("--strandness",    dest = "strandness",    nargs = 1, default = "2",
                  help = "Strand of reads. Setting for featureCounts. "
                         "Default is reversely stranded because scNanoGPS scans from 3'-adaptor. "
                         "0 (unstranded), 1 (stranded) and 2 (reversely stranded). "
                         "Default: 2")

options, arguments = parser.parse_args()

#===pre-check===
options.log_f_name = os.path.join(options.o_dir, options.log_f_name)

termination = False
if not options.gtf or \
   not os.path.isfile(options.gtf):
	print("\nCannot find GTF file: " + str(options.gtf) + "\n")
	termination = True
if not os.path.isdir(options.o_dir):
	print("\nOutput directory is not exist: "    + options.o_dir   + "\n")
	termination = True
if not os.path.isdir(options.tmp_dir):
	print("\nTemporary directory is not exist: " + options.tmp_dir + "\n")
	termination = True
if not curator_io.find_exe(options.featurecounts):
	print("\nCannot find path to featureCounts !\n")
	termination = True

if termination:
	parser.print_help()
	sys.exit(1)

logger = open(options.log_f_name, "wt")
logger.write(" ".join(sys.argv) + "\n\n")
logger.write("Output directory:    " + options.o_dir + "\n")
logger.write("Temporary directory: " + options.tmp_dir + "\n")
logger.close()

start_time = time.time()
print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
#===load CB list===
print("Loading CB list...", flush = True)
bam_list = glob.glob(os.path.join(options.tmp_dir, "*.curated.minimap2.bam"))
CB_list = []
for bam_f in bam_list:
	split_f_name = bam_f.split('/')
	BC = split_f_name[len(split_f_name) - 1].split('.')
	CB_list.append(BC[0])

print("Generating exon usage matrix...", flush = True)

cmd = options.featurecounts + " -L -O -T " + str(options.ncores) + \
      " -t exon -g exon_id -f --extraAttributes gene_name,transcript_id,exon_number" + \
      " -s " + options.strandness + \
      " -a " + options.gtf + \
      " -o " + os.path.join(options.o_dir, options.o_name)

for CB in CB_list:
	cmd = cmd + " " + os.path.join(options.tmp_dir, CB) + ".curated.minimap2.bam"

code_msg, stdout, stderr = curator_io.sys_run(cmd)

with open(os.path.join(options.o_dir, options.o_name), "at") as oh:
	oh.write(stdout.decode("utf-8"))

with open(options.log_f_name, "at") as logger:
	logger.write(stderr.decode("utf-8"))

#===combine exon number to TxID===
fh = open(os.path.join(options.o_dir, options.o_name),          "rt")
oh = open(os.path.join(options.o_dir, options.o_name) + ".tmp", "wt")
oh.write(fh.readline())

line_list = fh.readline().rstrip().split("\t")
#---extract BC---
BC_list = list()
for i in line_list[9:]:
	line_list_list = i.split('.')
	pre_BC = line_list_list[len(line_list_list) - 4].split('/')
	BC_list.append(pre_BC[len(pre_BC) - 1])
#---extract BC---
oh.write("\t".join(line_list[0:7]) + "_" + line_list[7] + "_" + line_list[8] + "\t" + "\t".join(BC_list) + "\n")

exon_dict = dict()
while True:
	line = fh.readline()
	if not line:
		break
	line_list = line.split("\t")
	if line_list[0] in exon_dict:
		continue
	else:
		oh.write("\t".join(line_list[0:7]) + "_" + line_list[7] + "_exon_" + line_list[8] + "\t" + "\t".join(line_list[9:]))
		exon_dict[line_list[0]] = 1
oh.close()
fh.close()

os.system("mv " + os.path.join(options.o_dir, options.o_name) + ".tmp " + os.path.join(options.o_dir, options.o_name))

print("\nFinished !\n")
print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
hours, minutes, seconds = misc.get_time_elapse(start_time)
misc.report_time_elapse(hours, minutes, seconds)


#! /usr/bin/env python3

import time, os, sys, glob
import pandas as pd
import multiprocessing as mp
from functools import partial
from contextlib import contextmanager
from optparse import OptionParser
from scanner_core import misc
from curator_core import curator_io

@contextmanager
def poolcontext(*args, **kwargs):
	pool = mp.Pool(*args, **kwargs)
	yield pool
	pool.terminate()

def load_CB(options):
	#===load CB list===
	print("Loading CB list...", flush = True)
	bam_list, CB_list = [], []
	if not options.CB_file:
		bam_list = glob.glob(os.path.join(options.tmp_dir, "*" + options.bam_suf))
		for bam_f in bam_list:
			split_f_name = bam_f.split('/')
			BC = split_f_name[len(split_f_name) - 1].split('.')
			CB_list.append(BC[0])
	else:
		with open(os.path.join(options.o_dir, options.CB_file), "rt") as CBF:
			while True:
				CB_name = CBF.readline().rstrip()
				if not CB_name:
					break
				CB_list.append(CB_name)
		bam_list = [x + options.bam_suf for x in CB_list]
	return bam_list, CB_list

def add_RG(CB, options):
	scTime = time.time()
	print("Adding RG on " + CB + "...", flush = True)
	#===Add RG into bam===
	cmd = options.samtools + " addreplacerg -r ID:" + CB + \
              " -o " + os.path.join(options.isoquant_d, CB) + '.RG.sam' + \
              " " + os.path.join(options.tmp_dir, CB) + options.bam_suf
	os.system(cmd)
	hours, minutes, seconds = misc.get_time_elapse(scTime)
	print("Adding RG on " + CB + " spent %d : %d : %.2f" % (hours, minutes, seconds), flush = True)

def pre_sort(i, start_no, end_no, CB_list, options):
	mrg_time = time.time()
	tmp_bam_pref = os.path.join(options.isoquant_d, "merged_tmp_") + str(i)
	for CB in CB_list[start_no[i]:end_no[i]]:
		cmd = 'grep -v "^@" ' + os.path.join(options.isoquant_d, CB) + ".RG.sam >> " + tmp_bam_pref + ".unsorted.sam"
		os.system(cmd)
		cmd = 'rm ' + os.path.join(options.isoquant_d, CB) + ".RG.sam"
		os.system(cmd)
	hours, minutes, seconds = misc.get_time_elapse(mrg_time)
	print("Merging batch %d spent %d : %d : %.2f\n" % (i, hours, minutes, seconds), flush = True)

def merge_bam(CB_list, options):
	mrg_oa_time = time.time()
	mrg_pref = os.path.join(options.isoquant_d, "merged")
	start_no = list(range(0, len(CB_list), options.batch_no))
	end_no   = [x + options.batch_no for x in start_no]
	if end_no[-1] > len(CB_list):
		end_no[-1] = len(CB_list)

	with poolcontext(processes = options.ncores) as pool:
		pool.map(partial(pre_sort, start_no = start_no, end_no = end_no, CB_list = CB_list, options = options), list(range(0, len(start_no))))

	mrg_time = time.time()
	cmd = "cat " + options.hdr_f + " > " + mrg_pref + ".unsorted.sam"
	os.system(cmd)
	for i in list(range(0, len(start_no))):
		tmp_sam_f = os.path.join(options.isoquant_d, "merged_tmp_") + str(i) + ".unsorted.sam"
		cmd = 'cat ' + tmp_sam_f + " >> " + mrg_pref + ".unsorted.sam"
		os.system(cmd)
		cmd = 'rm ' + tmp_sam_f
		os.system(cmd)

	hours, minutes, seconds = misc.get_time_elapse(mrg_time)
	print("Sam files merging spent %d : %d : %.2f\n" % (hours, minutes, seconds), flush = True)

	mrg_time = time.time()
	cmd = options.samtools + ' sort -@ ' + str(options.ncores) + ' -m 2G -o ' + mrg_pref + ".curated.minimap2.bam " + mrg_pref + ".unsorted.sam"
	os.system(cmd)
	hours, minutes, seconds = misc.get_time_elapse(mrg_time)
	print("Samtools sort spent %d : %d : %.2f\n" % (hours, minutes, seconds), flush = True)

	cmd = 'rm ' + mrg_pref + ".unsorted.sam"
	os.system(cmd)

	mrg_time = time.time()
	cmd = options.samtools + " index " + mrg_pref + ".curated.minimap2.bam"
	os.system(cmd)

	hours, minutes, seconds = misc.get_time_elapse(mrg_oa_time)
	print("Merging bam files overall spent %d : %d : %.2f\n" % (hours, minutes, seconds), flush = True)

def parse_txid2txname(options):
	tx_mapping = {}
	fh = open(options.gtf, "rt")
	while True:
		line = fh.readline()
		if not line:
			break
		if line.startswith("#"):
			continue
		line_list = line.rstrip().split("\t")
		if line_list[2] != "transcript":
			continue
		attr_list = line_list[8].split('; ')
		txid, txname = "", ""
		for x in attr_list:
			ele_list = x.split(' ')
			if ele_list[0] == "transcript_id":
				txid   = ele_list[1].split('"')[1]
			if ele_list[0] == "transcript_name":
				txname = ele_list[1].split('"')[1]
		tx_mapping.update({txid: txname})
	fh.close()
	return tx_mapping

if __name__ == "__main__":

	#===get params===
	parser = OptionParser()

	parser.add_option("-d",           dest = "o_dir",      nargs = 1, default = "scNanoGPS_res",
	                  help = "Output directory name. "
	                         "Default: scNanoGPS_res")
	parser.add_option("--tmp_dir",    dest = "tmp_dir",    nargs = 1, default = "tmp",
	                  help = "Temporary folder name. "
	                         "Default: tmp")
	parser.add_option("--CB_file",    dest = "CB_file",    nargs = 1, default = "filtered_barcode_list.txt",
	                  help = "File name for filtered barcode list. "
	                         "Default: filtered_barcode_list.txt")
	parser.add_option("--ref_genome", dest = "ref_genome", nargs = 1, default = None,
	                  help = "* Required ! "
	                         "File for reference genome.")
	parser.add_option("--gtf",        dest = "gtf",        nargs = 1, default = None,
	                  help = "* Required ! "
	                         "File for genome annotation.")
	parser.add_option("--bam_suf",    dest = "bam_suf",    nargs = 1, default = ".curated.minimap2.bam",
	                  help = "Suffix of the bam files. "
	                         "Default: .curated.minimap2.bam")
	parser.add_option("-o",           dest = "o_name",     nargs = 1, default = "isoform_exp_matrix.tsv.gz",
	                  help = "Counting table name. "
	                         "Default: isoform_exp_matrix.tsv.gz")
	parser.add_option("--log",        dest = "log_f_name", nargs = 1, default = 'logs/reporter_isoform.log.txt',
	                  help = "Log file name. "
	                         "Default: logs/reporter_isoform.log.txt")
	parser.add_option("-t",           dest = "ncores",     nargs = 1, default = 1,
	                  help = "Number of cores for program running. "
	                         "Default: 1", type = "int")
	parser.add_option("--batch_no",   dest = "batch_no",   nargs = 1, default = 500,
	                  help = "Batch number for merging bam files. "
	                         "Default: 500", type = "int")
	parser.add_option("--samtools",   dest = "samtools",   nargs = 1, default = "samtools",
	                  help = "Path to samtools. "
	                         "Default: samtools")
	parser.add_option("--isoquant",   dest = "isoquant",   nargs = 1, default = "isoquant.py",
	                  help = "Program name of IsoQuant. "
	                         "Default: isoquant.py")
	parser.add_option("--isoquant_d", dest = "isoquant_d", nargs = 1, default = "scNanoGPS_res/IsoQuant_res",
	                  help = "IsoQuant output directory. "
	                         "Default: scNanoGPS_res/IsoQuant_res")
	parser.add_option("--isoquant_o", dest = "isoquant_o",   nargs = 1, default = 'OUT/OUT.transcript_model_grouped_counts.tsv',
	                  help = "Suffix of IsoQuant output file. "
	                         "Default: OUT/OUT.transcript_model_grouped_counts.tsv")

	options, arguments = parser.parse_args()

	#===pre-check===
	termination = False
	if not options.ref_genome and not os.path.isdir(options.ref_genome):
		print("\nCannot find reference genome\n")
		termination = True
	if not options.gtf and not os.path.isdir(options.gtf):
		print("\nCannot find genome annotation\n")
		termination = True
	if not os.path.isdir(options.o_dir):
		print("\nOutput directory is not exist: "    + options.o_dir    + "\n")
		termination = True
	if not os.path.isdir(options.tmp_dir):
		print("\nTemporary directory is not exist: " + options.tmp_dir  + "\n")
		termination = True
	if not os.path.isdir(options.isoquant_d):
		print("\nCreate IsoQuant output directory: " + options.isoquant_d + "\n")
		os.system("mkdir " + options.isoquant_d)
	else:
		print("\nFind previous IsoQuant output directory !!! \n")
	if not os.path.isfile(os.path.join(options.o_dir, options.CB_file)):
		print("\nCannot find filter barcode list at: " + os.path.join(options.o_dir, options.CB_file)  + "\n")
		print("All the barcodes under folder " + options.tmp_dir + " will be used insteat.\n")
		options.CB_file = None
	if not curator_io.find_exe(options.samtools):
		print("\nCannot find path to samtools !\n")
		termination = True
	if not options.isoquant and not curator_io.find_exe(options.isoquant):
		print("\nCannot find IsoQuant !\n")
		termination = True
	options.hdr_f = os.path.join(options.isoquant_d, "merged.hdr")

	if termination:
		parser.print_help()
		sys.exit(1)

	logger = open(options.log_f_name, "wt")
	logger.write(" ".join(sys.argv) + "\n\n")
	logger.write("Output directory:    " + options.o_dir    + "\n")
	logger.write("Temporary directory: " + options.tmp_dir  + "\n")
	logger.write("Reference genome:    " + options.ref_genome + "\n")
	logger.close()

	#=== set env variables
	os.environ["OMP_NUM_THREADS"]        = str(options.ncores)
	os.environ["OPENBLAS_NUM_THREADS"]   = str(options.ncores)
	os.environ["MKL_NUM_THREADS"]        = str(options.ncores)
	os.environ["VECLIB_MAXIMUM_THREADS"] = str(options.ncores)
	os.environ["NUMEXPR_NUM_THREADS"]    = str(options.ncores)

	import numpy as np

	start_time = time.time()
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)

	#===load CB list===
	bam_list, CB_list = load_CB(options)

	#===Extracting header===
		#===extract header from first bam===
	cmd = options.samtools + " view -H " + os.path.join(options.tmp_dir, CB_list[0]) + options.bam_suf
	hdr, sm_str = "", os.getcwd().split(os.sep)[-1] 
	code_msg, out_msg, err_msg = curator_io.sys_run(cmd)
	if err_msg.decode("utf-8") == "":
		hdr = out_msg.decode("utf-8")
	else:
		print("Something wrong doing \"" + cmd + "\"\n")
		sys.exit()

	with open(options.hdr_f, "wt") as oh:
		for line in hdr.splitlines():
			#===remove @PG===
			if line.startswith('@PG'):
				continue
			oh.write(line + "\n")
		#===Add RG to header===
		oh.write("\n".join(["@RG\tID:" + CB + "\tSM:" + sm_str for CB in CB_list]) + "\n")

	#===add RG, extract bam, generate to master bam file===
	#===Add RG into bam===
	with poolcontext(processes = options.ncores) as pool:
		pool.map(partial(add_RG, options = options), CB_list)

	hours, minutes, seconds = misc.get_time_elapse(start_time)
	print("Adding RG takes %d : %d : %.2f\n" % (hours, minutes, seconds), flush = True)

	merge_bam(CB_list, options)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	print("Merging BAM files (include adding RG) takes %d : %d : %.2f\n" % (hours, minutes, seconds), flush = True)

	#===Run IsoQuant===
	isoquant_time = time.time()
	cmd = "python3 " + options.isoquant + " -d nanopore --complete_genedb" + \
	      " -r " + options.ref_genome + \
	      " -g " + options.gtf + \
	      " -t " + str(options.ncores) + \
	      " --read_group tag:RG" + \
	      " --bam " + os.path.join(options.isoquant_d, "merged.curated.minimap2.bam") + \
	      " -o " + options.isoquant_d
	os.system(cmd)
	hours, minutes, seconds = misc.get_time_elapse(isoquant_time)
	print("IsoQuant spent %d : %d : %.2f\n" % (hours, minutes, seconds), flush = True)

	hours, minutes, seconds = misc.get_time_elapse(start_time)
	print("Total time spent %d : %d : %.2f\n" % (hours, minutes, seconds), flush = True)

	#===parse transcript ID===
	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), "\n", flush = True)
	print("Parsing GTF...", flush = True)
	options.tx_mapping = parse_txid2txname(options)

	#===generate isoform matrix===
	res_df = pd.read_csv(os.path.join(options.isoquant_d, options.isoquant_o), header = 0, sep = '\t')
	res_df = res_df.rename(columns={'#feature_id': 'transcript_id'})
	res_df = res_df.loc[res_df['transcript_id'].isin(options.tx_mapping.keys()), :].copy()
	res_df['transcript_name'] = [options.tx_mapping[tid] for tid in res_df['transcript_id'].values]

	res_df = res_df[['transcript_id', 'transcript_name'] + res_df.columns[1:len(res_df.columns)-1].to_list()].copy()

	if res_df.shape[0] > 0:
		compression = None
		if options.o_name.endswith(".gz"):
			compression = 'gzip'
		res_df.to_csv(os.path.join(options.o_dir, options.o_name), sep='\t', header = True, index = False, compression = compression)
	else:
		if options.o_name.endswith(".gz"):
			oh = gzip.open(os.path.join(options.o_dir, options.o_name), "wt")
		else:
			oh = open(os.path.join(options.o_dir, options.o_name), "wt")
		oh.write("\t".join(key_list) + "\n")
		oh.close()

	print("\nTime stamp: " + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n", flush = True)
	hours, minutes, seconds = misc.get_time_elapse(start_time)
	misc.report_time_elapse(hours, minutes, seconds)

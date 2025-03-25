def open_file(f_name, mode):
	if f_name.endswith(".gz"):
		import gzip
		reader = gzip.open(f_name, mode)
	else:
		reader = open(f_name, mode)
	return reader

def load_df(f_name):
	import pandas as pd

	compression = None
	if f_name.endswith(".gz"):
		compression = 'gzip'
	res_df = pd.read_table(f_name, sep = '\t', header = 0, compression = compression)
	return res_df

def extract_id_list(line):
	return [int(element) for element in line.split("; ")]

def find_exe(p_name):
	import os

	for path in os.environ['PATH'].split(os.pathsep):
		if os.path.isfile(os.path.join(path, p_name)):
			return(os.path.join(path, p_name))

	return False

def sys_run(cmd):
	import subprocess

	proc = subprocess.Popen(cmd.split(' '),
		stdout = subprocess.PIPE,
		stderr = subprocess.PIPE)

	stdout, stderr = proc.communicate()
	proc.wait()

	return proc.returncode, stdout, stderr

def filter_softclipping(softclipping_thr, input_bam, output_bam, drop_bam):
	import os, pysam

	bam_i = pysam.AlignmentFile(input_bam,  "rb")
	bam_o = pysam.AlignmentFile(output_bam, "wb", template = bam_i)

	clipped_rid_dict = dict()

	for read in bam_i.fetch():
		no_m, no_s = 0, 0

		for cigar in read.cigartuples:
			if cigar[0] == 0:
				no_m += cigar[1]
			if cigar[0] == 4 or cigar[0] == 5:
				no_s += cigar[1]

#       M       BAM_CMATCH      0
#       I       BAM_CINS        1
#       D       BAM_CDEL        2
#       N       BAM_CREF_SKIP   3
#       S       BAM_CSOFT_CLIP  4
#       H       BAM_CHARD_CLIP  5
#       P       BAM_CPAD        6
#       =       BAM_CEQUAL      7
#       X       BAM_CDIFF       8
#       B       BAM_CBACK       9

		if (no_m / (no_m + no_s)) >= softclipping_thr:
			bam_o.write(read)
		else:
			clipped_rid_dict.update({read.query_name: 1})

	bam_i.close()
	bam_o.close()

	bam_i = pysam.AlignmentFile(input_bam,  "rb")
	bam_d = pysam.AlignmentFile(drop_bam,   "wb", template = bam_i)

	for read in bam_i.fetch():
		if read.query_name in clipped_rid_dict:
			bam_d.write(read)

	bam_i.close()
	bam_d.close()

def split_reads_by_pos(bam_i, seq_dict):
	i = 0
	buf_list, compiled_list = [], []
	pre_chr   = ""
	pre_start = 0
	for read in bam_i.fetch():
		i += 1
		if (pre_chr != read.reference_name or \
		    pre_start + 5 < read.reference_start) and i > 1:
			compiled_list.append(buf_list)
			buf_list  = []
			pre_chr   = read.reference_name
			pre_start = read.reference_start

		buf_list.append({'query_name': read.query_name, \
		                 'chr_name':   read.reference_name, \
		                 'start':      read.reference_start, \
		                 'sequence':   seq_dict[read.query_name]})
		if i == 1:
			pre_chr   = read.reference_name
			pre_start = read.reference_start

	
	compiled_list.append(buf_list)
	buf_list = []

	return compiled_list

def rec_umi_lookup(dist_df, this_umi, umi_array):
	if this_umi in umi_array:
		return this_umi
	else:
		for ele in dist_df:
			if ele[1] == this_umi:
				return rec_umi_lookup(dist_df, ele[0], umi_array)

def rec_swap_dist_df(test_df):
	#===check bc1===
	if len(test_df[test_df['bc1'].duplicated(keep = 'last')].values) > 0:
		for idx in test_df[test_df['bc1'].duplicated(keep = 'last')].index.values:
			tmp_val = test_df.iloc[idx]['bc1']
			test_df.iloc[idx]['bc1'] = test_df.iloc[idx]['bc2']
			test_df.iloc[idx]['bc2'] = tmp_val
			rec_swap_dist_df(test_df)
	#===check bc2===
	if len(test_df[test_df['bc2'].duplicated(keep = 'last')].values) > 0:
		for idx in test_df[test_df['bc2'].duplicated(keep = 'last')].index.values:
			tmp_val = test_df.iloc[idx]['bc1']
			test_df.iloc[idx]['bc1'] = test_df.iloc[idx]['bc2']
			test_df.iloc[idx]['bc2'] = tmp_val
			rec_swap_dist_df(test_df)
	return test_df

def collapse_UMI(BC, reads, options):
	import os, distance
	import pandas as pd

	consensus, duplicates, singletons = {}, [], []
	umi_list = dict()

	#===format reads to umi_list
	for read in reads:
		umi_split_list = read['query_name'].split("_")
		if umi_split_list[len(umi_split_list) - 1] in umi_list:
			umi_list[umi_split_list[len(umi_split_list) - 1]].append(read)
		else:
			umi_list[umi_split_list[len(umi_split_list) - 1]] = [read]

	#===use LD to merge UMI===
	if int(options.umi_ld) > 0:
		umi_seq = list(umi_list.keys())
		dist_df = list()
		if len(umi_seq) > 1:
			for ui in range(0, len(umi_seq) - 1):
				for uj in range(ui + 1, len(umi_seq)):
					ld = distance.levenshtein(umi_seq[ui], umi_seq[uj])
					if ld <= int(options.umi_ld):
						dist_df.append([umi_seq[ui], umi_seq[uj]])
			if len(dist_df) > 0:
				for dist_idx in range(0, len(dist_df)):
					if (dist_df[dist_idx][0] in umi_list) & (dist_df[dist_idx][1] in umi_list):
						umi_list[dist_df[dist_idx][0]].extend(umi_list[dist_df[dist_idx][1]])
						del(umi_list[dist_df[dist_idx][1]])

	#===need parallel computing===
	for umi in umi_list:
		if len(umi_list[umi]) > 1:
			#===build consensus===
			sorted_list = sorted(umi_list[umi], key = lambda key: key['query_name'])
			tmp_fname = os.path.join(options.tmp_dir, "tmp_") + BC + "_" + sorted_list[0]['chr_name'] + "_" + str(sorted_list[0]['start']) + "_" + sorted_list[0]['query_name'] + "_consensus.fa"
			con_t = open(tmp_fname, "wt")
			counter = 0
			for read in sorted_list:
				counter = counter + 1
				duplicates.append(read['query_name'] + ";" + \
				                  read['chr_name']   + ";" + \
				              str(read['start']))
				con_t.write(">" + read['query_name'] + "\n")
				con_t.write(read['sequence'] + "\n")
				if counter >= options.max_umi_duplicates:
					break
			con_t.close()

			cmd = options.spoa + " " + tmp_fname
			code_msg, out_msg, err_msg = sys_run(cmd)
			cmd = "rm " + tmp_fname
			os.system(cmd)

			consensus[sorted_list[0]['query_name']] = out_msg.decode("utf-8").split("\n")[1]
		else:
			#===output===
			singletons.append(umi_list[umi][0]['query_name'] + ";" + \
			                  umi_list[umi][0]['chr_name']   + ";" + \
			              str(umi_list[umi][0]['start']))

	return [consensus, duplicates, singletons]

def build_read_seq_dict(sam_content):
	seq_dict = dict()
	for read in sam_content.split("\n"):
		read_list = read.split("\t")
		if len(read_list) < 11:
			continue
		if read_list[0] not in seq_dict and \
		   read_list[9] != '*':
			seq_dict[read_list[0]] = read_list[9]
	return seq_dict

def check_mapping_db(options):
	if options.idx_genome:
		return options.idx_genome
	else:
		return options.ref_genome

def filter_bam_inc(fq_pref, options):
	import os

	if options.inc_bed:
		cmd = options.samtools + " view -b -L " + options.inc_bed + " " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam -o " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.selected.bam"
		os.system(cmd)

		cmd = ""
		if not options.keep_meta:
			cmd += "cp "
		else:
			cmd += "mv "
		cmd += os.path.join(options.tmp_dir, fq_pref) + ".minimap2.selected.bam " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam"
		os.system(cmd)

		cmd = options.samtools + " index " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam"
		os.system(cmd)

def filter_bam_exc(fq_pref, options):
	import os

	if options.exc_bed:
		cmd = options.samtools + " view -b -L " + options.exc_bed + " " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam -o " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.selected.bam -U " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.unselected.bam"
		os.system(cmd)

		cmd = ""
		if not options.keep_meta:
			cmd += "cp "
		else:
			cmd += "mv "
		cmd += os.path.join(options.tmp_dir, fq_pref) + ".minimap2.unselected.bam " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam"
		os.system(cmd)

		if not options.keep_meta:
			cmd = "rm " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.selected.bam"
			os.system(cmd)

		cmd = options.samtools + " index " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam"
		os.system(cmd)

def curation_master(fq, options):
	import time, os, pysam, subprocess
	from scanner_core import misc

	time_curation = time.time()

	fq_pref = fq.split(".fastq.gz")[0].split(options.tmp_dir)[1].split('/')[1]
	stat_msg = "Processing " + fq_pref + " ... "

	#---mapping---
	cmd = options.minimap2 + " -ax splice " + \
	      check_mapping_db(options) + " " + fq
	code_msg, out_msg, err_msg = sys_run(cmd)

	seq_dict = dict()
	with open(os.path.join(options.tmp_dir, fq_pref) + ".sam", "wt") as SAM:
		#---extract sequence before filtering---
		seq_dict = build_read_seq_dict(out_msg.decode("utf-8"))
		if options.inc_contig:
			SAM.write(out_msg.decode("utf-8"))
		#---filter out non-autosome---
		else:   
			lines = out_msg.decode("utf-8").split("\n")
			for line in lines:
				if line == "":
					break
				if line.split("\t")[2].startswith("chr"):
					SAM.write(line + "\n")

	with open(os.path.join(options.tmp_dir, fq_pref) + ".minimap2.log.txt", "wt") as MINIMAP2:
		MINIMAP2.write(err_msg.decode("utf-8"))

	#---samtools conversion---
	   #--- with header ---
	cmd = options.samtools + " view -Sb " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".sam " + \
	      "-T " + options.ref_genome + " " + \
	      "-o " + os.path.join(options.tmp_dir, fq_pref) + ".unsorted.bam"
	code_msg, out_msg, err_msg = sys_run(cmd)

	if not options.keep_meta:
		cmd = "rm " + os.path.join(options.tmp_dir, fq_pref) + ".sam"
		os.system(cmd)

	cmd = options.samtools + " sort " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".unsorted.bam " + \
	      "-o " + os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam"
	code_msg, out_msg, err_msg = sys_run(cmd)

	if not options.keep_meta:
		cmd = "rm " + os.path.join(options.tmp_dir, fq_pref) + ".unsorted.bam"
		os.system(cmd)

	cmd = options.samtools + " index " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam"
	code_msg, out_msg, err_msg = sys_run(cmd)

	#===filter in/out regions, ex: rRNA/tRNA===
	filter_bam_inc(fq_pref, options)
	filter_bam_exc(fq_pref, options)

	#===filter softclipping===
	filter_softclipping(options.softclipping_thr, \
	      os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam", \
	      os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam", \
	      os.path.join(options.tmp_dir, fq_pref) + ".high_softclipping.bam")

	cmd = options.samtools + " view " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam | " + \
	      "wc -l"
	no_raw = subprocess.check_output(cmd, shell = True).decode("utf-8").strip()

	cmd = options.samtools + " view " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam | " + \
	      "wc -l"
	no_filt = subprocess.check_output(cmd, shell = True).decode("utf-8").strip()

	cmd = options.samtools + " view " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".high_softclipping.bam | " + \
	      "wc -l"
	no_h_soft = subprocess.check_output(cmd, shell = True).decode("utf-8").strip()

	cmd = options.samtools + " index " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam"
	code_msg, out_msg, err_msg = sys_run(cmd)

	if not options.skip_curation:
		#===UMI collapse===
		bam_i = pysam.AlignmentFile(os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam",  "rb")
		compiled_list = split_reads_by_pos(bam_i, seq_dict)
		bam_i.close()

		uniq_consensus = dict()
		uniq_duplicates = dict()
		singleton_list = []

		for pos_list in compiled_list:
			consensus, duplicates, singletons = collapse_UMI(fq_pref, pos_list, options)
			for key in consensus:
				uniq_consensus[key] = consensus[key]
			for key in duplicates:
				uniq_duplicates[key] = 1
			singleton_list.extend(singletons)

		con_o = open(os.path.join(options.tmp_dir, fq_pref) + ".consensus.fasta", "wt")
		for readID in uniq_consensus:
			con_o.write(">" + readID + "\n")
			con_o.write(uniq_consensus[readID])
		con_o.close()

		bam_i = pysam.AlignmentFile(os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam",  "rb")
		bam_o = pysam.AlignmentFile(os.path.join(options.tmp_dir, fq_pref) + ".singleton.bam", "wb", template = bam_i)
		for read in bam_i.fetch():
			RID_pos_str = read.query_name     + ";" + \
			              read.reference_name + ";" + \
			          str(read.reference_start)
			if RID_pos_str in uniq_duplicates:
				continue
			read_str = read.query_name     + ";" + \
			           read.reference_name + ";" + \
			           str(read.reference_start)
			if read_str in singleton_list:
				bam_o.write(read)
		bam_i.close()
		bam_o.close()

		if not options.keep_meta:
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam"
			os.system(cmd)
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".minimap2.bam.bai"
			os.system(cmd)

		cmd = "grep '>' " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".consensus.fasta | " \
		      "wc -l "
		no_dup = subprocess.check_output(cmd, shell = True).decode("utf-8").strip()

		cmd = options.samtools + " view " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".singleton.bam | " + \
		      "wc -l"
		no_singleton = subprocess.check_output(cmd, shell = True).decode("utf-8").strip()

		if not options.keep_meta:
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam"
			os.system(cmd)
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam.bai"
			os.system(cmd)

		#---mapping consensus reads---
		cmd = options.minimap2 + " -ax splice " + \
		      check_mapping_db(options) + " " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".consensus.fasta"
		code_msg, out_msg, err_msg = sys_run(cmd)

		with open(os.path.join(options.tmp_dir, fq_pref) + ".consensus.sam", "wt") as SAM:
			SAM.write(out_msg.decode("utf-8"))

		with open(os.path.join(options.tmp_dir, fq_pref) + ".consensus.minimap2.log.txt", "wt") as MINIMAP2:
			MINIMAP2.write(err_msg.decode("utf-8"))

		if not options.keep_meta:
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".consensus.fasta"
			os.system(cmd)

		#---merge collaped reads---
			#---with header---
		cmd = "cat " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".consensus.sam >> " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".curated.sam"
		os.system(cmd)

		if not options.keep_meta:
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".consensus.?am"
			os.system(cmd)
			#---without header---
		cmd = options.samtools + " view " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".singleton.bam " + \
		      "-o " + os.path.join(options.tmp_dir, fq_pref) + ".singleton.sam"
		os.system(cmd)
		cmd = "cat " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".singleton.sam >> " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".curated.sam"
		os.system(cmd)

		if not options.keep_meta:
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".singleton.?am"
			os.system(cmd)

			#--- with partial header ---
		cmd = options.samtools + " view -Sb " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".curated.sam " + \
		      "-T " + options.ref_genome + " " + \
		      "-o " + os.path.join(options.tmp_dir, fq_pref) + ".curated.unsorted.bam"
		code_msg, out_msg, err_msg = sys_run(cmd)

		if not options.keep_meta:
			cmd = "rm " + \
			      os.path.join(options.tmp_dir, fq_pref) + ".curated.sam"
			os.system(cmd)

		cmd = options.samtools + " sort " + \
		      os.path.join(options.tmp_dir, fq_pref) + ".curated.unsorted.bam " + \
		      "-o " + os.path.join(options.tmp_dir, fq_pref) + ".curated.minimap2.bam"
		code_msg, out_msg, err_msg = sys_run(cmd)

		if not options.keep_meta:
			cmd = "rm " + os.path.join(options.tmp_dir, fq_pref) + ".curated.unsorted.bam"
			os.system(cmd)
	else:
		no_singleton, no_dup = "NA", "NA"
		cmd = "mv " + os.path.join(options.tmp_dir, fq_pref) + ".filtered.bam " + os.path.join(options.tmp_dir, fq_pref) + ".curated.minimap2.bam"
		os.system(cmd)

	cmd = options.samtools + " index " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".curated.minimap2.bam"
	code_msg, out_msg, err_msg = sys_run(cmd)

	cmd = options.samtools + " view " + \
	      os.path.join(options.tmp_dir, fq_pref) + ".curated.minimap2.bam | " + \
	      "wc -l"
	no_curated = subprocess.check_output(cmd, shell = True).decode("utf-8").strip()

	hours, minutes, seconds = misc.get_time_elapse(time_curation)
	time_spent = "%d : %d : %.2f" % (hours, minutes, seconds)

	print(stat_msg + " Done. Spent " + time_spent, flush = True)
	return("\t".join([fq_pref, no_raw, no_filt, no_h_soft, no_singleton, no_dup, no_curated, time_spent]))


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

def collapse_UMI(reads, options):
	import os, copy, distance

	con_res, con_reads, bam_res = {}, [], []
	umi_list = dict()

	for read in reads:
		umi_split_list = read['query_name'].split("_")

		if umi_split_list[len(umi_split_list) - 1] in umi_list:
			umi_list[umi_split_list[len(umi_split_list) - 1]].append(read)
		else:
			umi_list[umi_split_list[len(umi_split_list) - 1]] = [read]

	#===use 2 LD to merge UMI===
	if int(options.umi_ld) > 0:
		umi_seq      = list(umi_list.keys())
		umi_mrg_dict = dict()
		dist_df      = list()
		for seq_idx in range(0, len(umi_seq)):
			umi_mrg_dict.update({seq_idx: [seq_idx]})
		if len(umi_seq) > 1:
			for ui in range(0, len(umi_seq) - 1):
				for uj in range(ui + 1, len(umi_seq)):
					ld = distance.levenshtein(umi_seq[ui], umi_seq[uj])
					if ld <= int(options.umi_ld):
						dist_df.append([umi_seq[ui], umi_seq[uj]])
			if len(dist_df) > 0:
				for dist_idx in range(len(dist_df) - 1, -1, -1):
					umi_list[dist_df[dist_idx][0]].extend(umi_list[dist_df[dist_idx][1]])
				for dist_idx in range(0, len(dist_df)):
					if dist_df[dist_idx][1] in umi_list:
						del(umi_list[dist_df[dist_idx][1]])

	#===need parallel computing===
	for umi in umi_list:
		if len(umi_list[umi]) > 1:
			#===build consensus===
			sorted_list = sorted(umi_list[umi], key = lambda key: key['query_name'])
			fname_str = sorted_list[0]['chr_name'] + "_" + str(sorted_list[0]['start']) + "_" + sorted_list[0]['query_name'] + "_consensus.fa"
			tmp_fname = os.path.join(options.tmp_dir, "tmp_") + fname_str
			con_t = open(tmp_fname, "wt")
			for read in sorted_list:
				con_reads.append(read['query_name'] + ";" + \
				                 read['chr_name']   + ";" + \
				             str(read['start']))
				con_t.write(">" + read['query_name'] + "\n")
				con_t.write(read['sequence'] + "\n")
			con_t.close()

			fname = os.path.join(options.tmp_dir, fname_str)
			cmd = options.spoa + " " + tmp_fname + " > " + fname
			os.system(cmd)
			cmd = "rm " + tmp_fname
			os.system(cmd)

			con_i = open(fname, "rt")
			con_head = con_i.readline()
			con_seq  = con_i.readline()
			con_i.close()

			cmd = "rm " + fname
			os.system(cmd)

			con_res[sorted_list[0]['query_name']] = con_seq

		else:
			#===output===
			bam_res.append(umi_list[umi][0]['query_name'] + ";" + \
			               umi_list[umi][0]['chr_name']   + ";" + \
			           str(umi_list[umi][0]['start']))

	return con_res, con_reads, bam_res

def build_read_seq_dict(bam_name):
	import pysam

	bam_f = pysam.AlignmentFile(bam_name, "rb")

	seq_dict = dict()
	for read in bam_f.fetch():
		if read.query_name not in seq_dict and \
		   read.query_sequence is not None:
			seq_dict[read.query_name] = read.query_sequence

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


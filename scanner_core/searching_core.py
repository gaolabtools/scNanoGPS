def getHeader(na_seq, options):
	return na_seq[0: options.scan_region]

def getTail(na_seq, options):
	return na_seq[-options.scan_region:len(na_seq)].translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))[::-1]

def gen_alignment(seq1, seq2, options):
	from Bio import Align

	aligner                  = Align.PairwiseAligner()
	aligner.mode             = 'local'
	aligner.match_score      = options.dp_penalty[0]
	aligner.mismatch_score   = options.dp_penalty[1]
	aligner.open_gap_score   = options.dp_penalty[2]
	aligner.extend_gap_score = options.dp_penalty[3]

	return aligner.align(seq1, seq2)

def adaptor_search(na_seq, adaptor, options):
	alignment_res = gen_alignment(na_seq, adaptor, options)

	if len(alignment_res) == 0:
	        return None

	aln_res = alignment_res[0]
	for aln in alignment_res:
		if aln.score >= aln_res.score:
			aln_res = aln

	if aln_res.score < len(adaptor) * options.dp_penalty[0] * options.scoring_threshold:
		return None

	start = int(aln_res.indices[0][0])
	end   = int(aln_res.indices[0][-1]) + 1
	seqA, seqB = aln_res
	alignment, gap_no = [], 0
	for a, b in zip(seqA, seqB):
		if a == '-' or b == '-':
			alignment.append('-')
			gap_no += 1
		elif a == b:
			alignment.append('|')
		else:
			alignment.append(' ')

	return mappingRes(seqA, seqB, start, end, aln_res.score, "".join(alignment), gap_no)

def polyT_search(na_seq, options):
	pos_list = list()
	alignment_res = gen_alignment(na_seq, options.polyT, options)

	for aln in alignment_res:
		if aln.score >= len(options.polyT) * options.dp_penalty[0] * options.scoring_threshold:
			pos_list.append([aln.indices[0][0], aln.indices[0][1]])

	return pos_list

def ten_nano_workflow(read_data, options):
	from scanner_core import scanner_io

	#===check read length===
	if len(read_data['na_seq']) < options.min_read_length:
		return {'rid': read_data['def_line'].split(' ')[0].split('@')[1], 'orientation': None, 'BC_start': None,
			'BC_seq': None, 'UMI_seq': None, 'Seq_end': None, 'mean_quality': None}

	#===init result dict===
	res_data = {}

	#===print serial number under debug_mode===
	if options.debug_mode:
		print(str(read_data['rid']))

	#===get head/tail region of na_seq===
	na_seq_header = getHeader(read_data['na_seq'], options)
	na_seq_tail   = getTail(read_data['na_seq'],   options)
	qu_seq_header = getHeader(read_data['qu_seq'], options)
	qu_seq_tail   = getTail(read_data['qu_seq'],   options)

	#===Step 1: Brute force median search===
	ht_res = polyT_search(na_seq_header, options)
	tt_res = polyT_search(na_seq_tail,   options)

	#===Step 2: Precisely search===
	h5_ps_res, h3_ps_res, t5_ps_res, t3_ps_res = None, None, None, None
	if len(ht_res) > 0:
		h3_ps_res = adaptor_search(na_seq_header, options.adaptor_three_p, options)
		t5_ps_res = adaptor_search(na_seq_tail,   options.adaptor_five_p,  options)
	if len(tt_res) > 0:
		t3_ps_res = adaptor_search(na_seq_tail,   options.adaptor_three_p, options)
		h5_ps_res = adaptor_search(na_seq_header, options.adaptor_five_p,  options)

	#===Step 3: tie breaker for co-existence of h3 and t3===
	#===use BC+UMI+polyT to break tie===
	if h3_ps_res and t3_ps_res:
		h_in_order = check_adaptor_BC_UMI_polyT_in_distance(h3_ps_res, ht_res, options)
		t_in_order = check_adaptor_BC_UMI_polyT_in_distance(t3_ps_res, tt_res, options)

		if h_in_order and t_in_order:
			h3_ps_res = None
			t3_ps_res = None
		elif h_in_order:
			t3_ps_res = None
		elif t_in_order:
			h3_ps_res = None
		else:
			h3_ps_res = None
			t3_ps_res = None

	#===Step 4: check reads status===
	if h3_ps_res:
		res_data['counter_h_3p'] = 1
		#===check boundaries===
		if h3_ps_res.end < len(options.adaptor_three_p):
			res_data['counter_h_partial_3p']  = 1
		if ((h3_ps_res.end - h3_ps_res.start) == len(options.adaptor_three_p)) and \
		   h3_ps_res.gap_no == 0 and \
		   len([s for s in h3_ps_res.alignment if s == '|']) == len(options.adaptor_three_p):
			res_data['counter_h_perfect_3p']  = 1
		if h3_ps_res.alignment[-1] == " ":
			res_data['counter_h_last_1_mm']   = 1
		if h3_ps_res.alignment[-1] == " " and h3_ps_res.alignment[-2] == " ":
			res_data['counter_h_last_12_mm']  = 1
		if h3_ps_res.alignment[-1] == " " and h3_ps_res.alignment[-2] == " " and h3_ps_res.alignment[-3] == " ":
			res_data['counter_h_last_123_mm'] = 1
		#===check boundaries===

	if t3_ps_res:
		res_data['counter_t_3p'] = 1
		#===check boundaries===
		if t3_ps_res.end < len(options.adaptor_three_p):
			res_data['counter_t_partial_3p']  = 1
		if ((t3_ps_res.end - t3_ps_res.start) == len(options.adaptor_three_p)) and \
		   t3_ps_res.gap_no == 0 and \
		   len([s for s in t3_ps_res.alignment if s == '|']) == len(options.adaptor_three_p):
			res_data['counter_t_perfect_3p']  = 1
		if t3_ps_res.alignment[-1] == " ":
			res_data['counter_t_last_1_mm']   = 1
		if t3_ps_res.alignment[-1] == " " and t3_ps_res.alignment[-2] == " ":
			res_data['counter_t_last_12_mm']  = 1
		if t3_ps_res.alignment[-1] == " " and t3_ps_res.alignment[-2] == " " and t3_ps_res.alignment[-3] == " ":
			res_data['counter_t_last_123_mm'] = 1
		#===check boundaries===

	#===Step 5: Extract BC, UMI===
	orientation, BC_start, UMI_start, BC_seq, UMI_seq, Seq_end, mean_quality = None, None, None, None, None, len(read_data['na_seq']), None
	if h3_ps_res and check_adaptor_BC_UMI_polyT_in_distance(h3_ps_res, ht_res, options):
		res_data['counter_h_3p_polyT'] = 1
		BC_start     = h3_ps_res.end
		rc_seq       = read_data['na_seq']
		BC_seq       = rc_seq[BC_start:BC_start + options.BC_len]
		UMI_start    = BC_start + options.BC_len
		UMI_seq      = rc_seq[UMI_start:UMI_start + options.UMI_len]

		if t5_ps_res:
			Seq_end = len(read_data['na_seq']) - t5_ps_res.end

		BC_quality   = read_data['qu_seq'][BC_start:BC_start + options.BC_len]
		sum_quality  = 0
		for i in range(0, len(BC_quality)):
			sum_quality += (ord(BC_quality[i]) - 33)
		mean_quality = round(sum_quality / len(BC_quality), 2)
		orientation  = "H"

	if t3_ps_res and check_adaptor_BC_UMI_polyT_in_distance(t3_ps_res, tt_res, options):
		res_data['counter_t_3p_polyT'] = 1
		BC_start     = t3_ps_res.end
		rc_seq       = read_data['na_seq'].translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))[::-1]
		BC_seq       = rc_seq[BC_start:BC_start + options.BC_len]
		UMI_start    = BC_start + options.BC_len
		UMI_seq      = rc_seq[UMI_start:UMI_start + options.UMI_len]

		if h5_ps_res:
			Seq_end = len(read_data['na_seq']) - h5_ps_res.end

		BC_quality   = read_data['qu_seq'][::-1][BC_start:BC_start + options.BC_len]
		sum_quality  = 0
		for i in range(0, len(BC_quality)):
			sum_quality += (ord(BC_quality[i]) - 33)
		mean_quality = round(sum_quality / len(BC_quality), 2)
		orientation  = "T"

	res_data['rid']          = read_data['def_line'].split(' ')[0].split('@')[1]
	res_data['orientation']  = orientation
	res_data['BC_start']     = str(BC_start)
	res_data['BC_seq']       = BC_seq
	res_data['UMI_seq']      = UMI_seq
	res_data['Seq_end']      = Seq_end
	res_data['mean_quality'] = mean_quality

	#===debug===
	if options.debug_mode:
		print("header: " + na_seq_header)
		if h3_ps_res:
			scanner_io.printPSAlignment(h3_ps_res)
			BC_start = h3_ps_res.end
			print("h3_ps_res.start: " + str(h3_ps_res.start + 1))
			print("h3_ps_res.end:   " + str(h3_ps_res.end))
			print("h3_ps_res.gap_n: " + str(h3_ps_res.gap_no))
			print("BC_start:        " + str(BC_start + 1))
			scanner_io.printAlignment(na_seq_header[BC_start:                  BC_start + options.BC_len],
			                                        BC_start,                  "BC")
			scanner_io.printAlignment(na_seq_header[BC_start + options.BC_len: BC_start + options.BC_len + options.UMI_len],
			                                        BC_start + options.BC_len, "UMI")
		if len(ht_res) > 0:
			scanner_io.printPT(ht_res)
		if h5_ps_res:
			scanner_io.printPSAlignment(h5_ps_res)
			sequence_end = h5_ps_res.end
			print("h5_ps_res.start: " + str(h5_ps_res.start + 1))
			print("h5_ps_res.end:   " + str(h5_ps_res.end))
			print("h5_ps_res.gap_n: " + str(h5_ps_res.gap_no))
			print("sequence_end:    " + str(sequence_end))
		print()

		print("tail:   " + na_seq_tail)
		if t3_ps_res:
			scanner_io.printPSAlignment(t3_ps_res)
			BC_start = t3_ps_res.end
			print("t3_ps_res.start: " + str(t3_ps_res.start + 1))
			print("t3_ps_res.end:   " + str(t3_ps_res.end))
			print("t3_ps_res.gap_n: " + str(t3_ps_res.gap_no))
			print("BC_start:        " + str(BC_start + 1))
			scanner_io.printAlignment(na_seq_tail[BC_start:                  BC_start + options.BC_len],
			                                      BC_start,                  "BC")
			scanner_io.printAlignment(na_seq_tail[BC_start + options.BC_len: BC_start + options.BC_len + options.UMI_len],
			                                      BC_start + options.BC_len, "UMI")
		if len(tt_res) > 0:
			scanner_io.printPT(tt_res)
		if t5_ps_res:
			scanner_io.printPSAlignment(t5_ps_res)
			sequence_end = t5_ps_res.end
			print("t5_ps_res.start: " + str(t5_ps_res.start + 1))
			print("t5_ps_res.end:   " + str(t5_ps_res.end))
			print("t5_ps_res.gap_n: " + str(t5_ps_res.gap_no))
			print("sequence_end:    " + str(sequence_end))
		print()

		if res_data['orientation']:
			if h3_ps_res:
				print("Count by header!!!")
			if t3_ps_res:
				print("Count by tail!!!")

			if res_data['BC_seq']:
				print("BC_seq:      " + res_data['BC_seq'])
				print("UMI_seq:     " + res_data['UMI_seq'])
			print()
	#===debug===

	return res_data

def counting_res(res_data, counter, tmp_data):
	for each_row in tmp_data:
		res_data.append([each_row['rid'], each_row['orientation'], each_row['BC_start'],
		                 each_row['BC_seq'], each_row['UMI_seq'], each_row['Seq_end'], each_row['mean_quality']])

		for row_key in each_row:
			if row_key.startswith('counter_'):
				if each_row[row_key]:
					counter[row_key] += each_row[row_key]
	return res_data, counter

def check_adaptor_BC_UMI_polyT_in_distance(ps_res, t_res, options):
	BC_end = ps_res.end
	if len(t_res) > 0:
		for pt_list in t_res:
			if (pt_list[0] - BC_end) >= (options.BC_len + options.UMI_len * (1 - 0.75)) and \
			   (pt_list[0] - BC_end) <= (options.BC_len + options.UMI_len * (1 + 0.75)):
				return True
	return False

class mappingRes:
	def __init__(self, seqA, seqB, start, end, score, alignment, gap_no):
		self.seqA      = seqA
		self.seqB      = seqB
		self.start     = start
		self.end       = end
		self.score     = score
		self.alignment = alignment
		self.gap_no    = gap_no


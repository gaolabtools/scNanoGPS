def getHeader(na_seq):
	return na_seq[0: 100]

def getTail(na_seq):
	return na_seq[-100:len(na_seq)].translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))[::-1]

def mySort(posList, rev = False):
	posList = list(filter(None, posList))
	if len(posList) > 0:
		return sorted(posList, reverse = rev)[0]
	else:
		return None

def median_search(na_seq, adaptor, m_threshold):
	min_mm = len(adaptor) + 1
	best_pos  = -1
	for i in range(len(adaptor), len(na_seq)):
		matching_p = ''
		mm_no = 0
		for j in range(0, len(adaptor)):
			na_pos  = i - j
			ada_pos = len(adaptor) - 1 - j

			if na_seq[na_pos] != adaptor[ada_pos]:
				mm_no = mm_no + 1
				matching_p = matching_p + ' '
			else:
				matching_p = matching_p + '|'

			if mm_no > min_mm:
				matching_p = ''
				break

		if matching_p != '':
			if mm_no < min_mm:
				min_mm = mm_no
				best_pos = i

	if min_mm <= (len(adaptor) * (1 - m_threshold)):
		return best_pos - len(adaptor) + 1
	else:
		return None

def precise_search(full_na_seq, adaptor, start_pos, end_pos, scoring_threshold, dp_penalty):
	from Bio import pairwise2

	if start_pos < 0:
		start_pos = 0
	if end_pos >= len(full_na_seq):
		end_pos = len(full_na_seq) - 1
	na_seq = full_na_seq[start_pos:end_pos]

	alignment_res = pairwise2.align.localms(na_seq, adaptor, dp_penalty[0], dp_penalty[1], dp_penalty[2], dp_penalty[3], one_alignment_only = True)

	if len(alignment_res) == 0:
		return None

	seqA  = str(alignment_res[0].seqA)
	seqB  = str(alignment_res[0].seqB)
	start = 0
	end   = len(alignment_res[0].seqA) - 1
	while seqB[start] == '-':
		start += 1
	while seqB[end] == '-':
		end -= 1

	if alignment_res[0].score < len(adaptor) * dp_penalty[0] * scoring_threshold and (start_pos != 0 or start != 0):
		return False

	alignment, gap_no = "", 0
	for i in range(start, end + 1):
		if seqA[i] == seqB[i]:
			alignment += "|"
		else:
			alignment += " "
		if seqA[i] == '-':
			gap_no += 1

	return mappingRes(seqA, seqB, start, end + 1, alignment_res[0].score, alignment, gap_no, start_pos)

def chech_alignment(na_seq, adaptor, start_pos, end_pos, a):
	print("na_seq: ", end = "")
	if (start_pos + a.start - a.seqA[a.start:a.end].count('-')) > 0:
		print(" " * (start_pos + a.start - a.seqA[a.start:a.end].count('-')), end = "")
	print(a.seqA[a.start:a.end])

	print("        ", end = "")
	print(" " * (start_pos + a.start - a.seqA[a.start:a.end].count('-')), end = "")
	for i in range(a.start, a.end):
		if a.seqA[i] == a.seqB[i]:
			print("|", end = "")
		else:
			print(" ", end = "")
	print()

	print("adaptor:", end = "")
	print(" " * (a.start_pos + start - seqA[start:end].count('-')), end = "")
	print(seqB[start:end])

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

	#===get head/tail 100 of na_seq===
	na_seq_header = getHeader(read_data['na_seq'])
	na_seq_tail   = getTail(read_data['na_seq'])
	qu_seq_header = getHeader(read_data['qu_seq'])
	qu_seq_tail   = getTail(read_data['qu_seq'])

	#===Step 1: Brute force median search===
	ht_res = median_search(na_seq_header, options.polyT, options.matching_percentage)
	tt_res = median_search(na_seq_tail,   options.polyT, options.matching_percentage)

	#===Step 2: Precisely search===
	h5_ps_res, h3_ps_res, t5_ps_res, t3_ps_res = None, None, None, None
	if ht_res:
		h3_ps_res = precise_search(na_seq_header, options.adaptor_three_p, 0, 100, options.scoring_threshold, options.dp_penalty)
		t5_ps_res = precise_search(na_seq_tail,   options.adaptor_five_p,  0, 100, options.scoring_threshold, options.dp_penalty)
	if tt_res:
		t3_ps_res = precise_search(na_seq_tail,   options.adaptor_three_p, 0, 100, options.scoring_threshold, options.dp_penalty)
		h5_ps_res = precise_search(na_seq_header, options.adaptor_five_p,  0, 100, options.scoring_threshold, options.dp_penalty)

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
		if h3_ps_res.adaptor_start + h3_ps_res.end - h3_ps_res.gap_no < len(options.adaptor_three_p):
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

		if h3_ps_res.seqA[h3_ps_res.end - 1] == "-" or h3_ps_res.seqB[h3_ps_res.end - 1] == "-":
			res_data['counter_h_last_1_i'] = 1
		if h3_ps_res.seqA[h3_ps_res.end - 2] == "-" or h3_ps_res.seqB[h3_ps_res.end - 2] == "-":
			res_data['counter_h_last_2_i'] = 1
		if h3_ps_res.seqA[h3_ps_res.end - 3] == "-" or h3_ps_res.seqB[h3_ps_res.end - 3] == "-":
			res_data['counter_h_last_3_i'] = 1
		#===check boundaries===

	if t3_ps_res:
		res_data['counter_t_3p'] = 1
		#===check boundaries===
		if t3_ps_res.adaptor_start + t3_ps_res.end - t3_ps_res.gap_no < len(options.adaptor_three_p):
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

		if t3_ps_res.seqA[t3_ps_res.end - 1] == "-" or t3_ps_res.seqB[t3_ps_res.end - 1] == "-":
			res_data['counter_t_last_1_i'] = 1
		if t3_ps_res.seqA[t3_ps_res.end - 2] == "-" or t3_ps_res.seqB[t3_ps_res.end - 2] == "-":
			res_data['counter_t_last_2_i'] = 1
		if t3_ps_res.seqA[t3_ps_res.end - 3] == "-" or t3_ps_res.seqB[t3_ps_res.end - 3] == "-":
			res_data['counter_t_last_3_i'] = 1
		#===check boundaries===

	#===Step 5: Extract BC, UMI===
	orientation, BC_start, UMI_start, BC_seq, UMI_seq, Seq_end, mean_quality = None, None, None, None, None, len(read_data['na_seq']), None
	if h3_ps_res and ht_res and check_adaptor_BC_UMI_polyT_in_distance(h3_ps_res, ht_res, options):
		res_data['counter_h_3p_polyT'] = 1
		BC_start    = h3_ps_res.adaptor_start + h3_ps_res.end - h3_ps_res.gap_no
		rc_seq      = na_seq_header
		if BC_start > 50:
			rc_seq = read_data['na_seq']
		BC_seq      = rc_seq[BC_start:BC_start + options.BC_len]

		UMI_start   = BC_start + options.BC_len
		UMI_seq     = rc_seq[UMI_start:UMI_start + options.UMI_len]

		if t5_ps_res:
			Seq_end = len(read_data['na_seq']) - (t5_ps_res.adaptor_start + t5_ps_res.end - t5_ps_res.gap_no)

		BC_quality  = read_data['qu_seq'][BC_start:BC_start + options.BC_len]
		sum_quality = 0
		for i in range(0, len(BC_quality)):
			sum_quality += (ord(BC_quality[i]) - 33)
		mean_quality = round(sum_quality / len(BC_quality), 2)

		orientation = "H"

	if t3_ps_res and tt_res and check_adaptor_BC_UMI_polyT_in_distance(t3_ps_res, tt_res, options):
		res_data['counter_t_3p_polyT'] = 1
		BC_start    = t3_ps_res.adaptor_start + t3_ps_res.end - t3_ps_res.gap_no
		rc_seq      = na_seq_tail
		if BC_start > 50:
			rc_seq     = read_data['na_seq'].translate(str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}))[::-1]
		BC_seq      = rc_seq[BC_start:BC_start + options.BC_len]

		UMI_start   = BC_start + options.BC_len
		UMI_seq     = rc_seq[UMI_start:UMI_start + options.UMI_len]

		if h5_ps_res:
			Seq_end = len(read_data['na_seq']) - (h5_ps_res.adaptor_start + h5_ps_res.end - h5_ps_res.gap_no)

		BC_quality  = read_data['qu_seq'][::-1][BC_start:BC_start + options.BC_len]
		sum_quality = 0
		for i in range(0, len(BC_quality)):
			sum_quality += (ord(BC_quality[i]) - 33)
		mean_quality = round(sum_quality / len(BC_quality), 2)

		orientation = "T"

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
			BC_start = h3_ps_res.adaptor_start + h3_ps_res.end - h3_ps_res.gap_no
			print("h3_ps_res.adaptor_start: " + str(h3_ps_res.adaptor_start))
			print("h3_ps_res.end:           " + str(h3_ps_res.end))
			print("h3_ps_res.gap_n:         " + str(h3_ps_res.gap_no))
			print("BC_start:                " + str(BC_start))
			scanner_io.printAlignment(na_seq_header[BC_start:                  BC_start + options.BC_len],
			                                        BC_start,                  "BC")
			scanner_io.printAlignment(na_seq_header[BC_start + options.BC_len: BC_start + options.BC_len + options.UMI_len],
			                                        BC_start + options.BC_len, "UMI")
		if ht_res:
			scanner_io.printAlignment(options.polyT, ht_res, "polyT")
		if h5_ps_res:
			scanner_io.printPSAlignment(h5_ps_res)
			sequence_end = h5_ps_res.adaptor_start + h5_ps_res.end - h5_ps_res.gap_no
			print("h5_ps_res.adaptor_start: " + str(h5_ps_res.adaptor_start))
			print("h5_ps_res.end:           " + str(h5_ps_res.end))
			print("h5_ps_res.gap_n:         " + str(h5_ps_res.gap_no))
			print("sequence_end:            " + str(sequence_end))
		print()

		print("tail:   " + na_seq_tail)
		if t3_ps_res:
			scanner_io.printPSAlignment(t3_ps_res)
			BC_start = t3_ps_res.adaptor_start + t3_ps_res.end - t3_ps_res.gap_no
			print("t3_ps_res.adaptor_start: " + str(t3_ps_res.adaptor_start))
			print("t3_ps_res.end:           " + str(t3_ps_res.end))
			print("t3_ps_res.gap_n:         " + str(t3_ps_res.gap_no))
			print("BC_start:                " + str(BC_start))
			scanner_io.printAlignment(na_seq_tail[BC_start:                  BC_start + options.BC_len],
			                                      BC_start,                  "BC")
			scanner_io.printAlignment(na_seq_tail[BC_start + options.BC_len: BC_start + options.BC_len + options.UMI_len],
			                                      BC_start + options.BC_len, "UMI")
		if tt_res:
			scanner_io.printAlignment(options.polyT, tt_res, "polyT")
		if t5_ps_res:
			scanner_io.printPSAlignment(t5_ps_res)
			sequence_end = t5_ps_res.adaptor_start + t5_ps_res.end - t5_ps_res.gap_no
			print("t5_ps_res.adaptor_start: " + str(t5_ps_res.adaptor_start))
			print("t5_ps_res.end:           " + str(t5_ps_res.end))
			print("t5_ps_res.gap_n:         " + str(t5_ps_res.gap_no))
			print("sequence_end:            " + str(sequence_end))
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
	BC_end = ps_res.adaptor_start + ps_res.end - ps_res.gap_no
	if (t_res - BC_end) >= (options.BC_len + options.UMI_len - options.UMI_len * 0.75) and \
	   (t_res - BC_end) <= (options.BC_len + options.UMI_len + options.UMI_len * 0.75):
		return True
	else:
		return False

class mappingRes:
	def __init__(self, seqA, seqB, start, end, score, alignment, gap_no, adaptor_start):
		self.seqA  = seqA
		self.seqB  = seqB
		self.start = start
		self.end   = end
		self.score = score
		self.alignment = alignment
		self.gap_no    = gap_no
		self.adaptor_start = adaptor_start


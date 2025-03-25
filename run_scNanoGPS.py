#! /usr/bin/env python3

import glob, os, sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i",              dest = "fq_f_name",     default = None,
                  nargs = 1, type = "string",
                  help = "* Required ! "
                         "Input FastQ/Fast5 file name, or directory containing multiple input files. "
                         "Support fastq/fq/fastq.gz/fq.gz/fast5 format.")
parser.add_option("-d",              dest = "o_dir",         default = "scNanoGPS_res",
                  nargs = 1, type = "string",
                  help = "Output directory name. "
                         "Default: scNanoGPS_res")
parser.add_option("--tmp_dir",       dest = "tmp_dir",       default = "tmp",
                  nargs = 1, type = "string",
                  help = "Temporary folder name. "
                         "Default: tmp")
parser.add_option("-p",              dest = "protocol",      default = "3p",
                  nargs = 1, type = "string",
                  help = "10x barcoding protocol. (3p / 5p / spatial) "
                         "Default: 3p")
parser.add_option("-t",              dest = "ncores",        default = 1,
                  nargs = 1, type = "string",
                  help = "Number of cores for program running. "
                         "Default: 1")
parser.add_option("--gtf",           dest = "gtf",           default = None,
                  nargs = 1, type = "string",
                  help = "* Required ! "
                         "GTF file for expression calling. ")
parser.add_option("--ref_genome",    dest = "ref_genome",    default = None,
                  nargs = 1, type = "string",
                  help = "* Required ! "
                         "File for reference genome.")
parser.add_option("--idx_genome",    dest = "idx_genome",    default = None,
                  nargs = 1, type = "string",
                  help = "Path to the Minimap2 genome index. "
                         "Program will use reference genome if no Minimap2 genome index given. "
                         "Default: None")
parser.add_option("--whitelist",     dest = "whitelist",     default = None,
                  nargs = 1, type = "string",
                  help = "Path to the cell barcode whitelist. "
                         "Default: None")
parser.add_option("--exc_bed",       dest = "exc_bed",       default = None,
                  nargs = 1, type = "string",
                  help = "Exclude specific regions (BED) in file. "
                         "Default: None")
parser.add_option("--isoquant",      dest = "isoquant",      default = None,
                  nargs = 1, type = "string",
                  help = "Provide path to IsoQuant to conduct isoform calling. "
                         "Default: None")
parser.add_option("--annovar",       dest = "annovar",       default = None,
                  nargs = 1, type = "string",
                  help = "Provide directory path to ANNOVAR to conduct SNP calling. "
                         "Default: None")
parser.add_option("--annovardb",     dest = "annovardb",     default = "hg38db",
                  nargs = 1, type = "string",
                  help = "Name of ANNOVAR database. "
                         "Default: hg38db")
parser.add_option("--annovargv",     dest = "annovargv",     default = "hg38",
                  nargs = 1, type = "string",
                  help = "Version of ANNOVAR genome version. "
                         "Default: hg38")
parser.add_option("--annovarprot",   dest = "annovarprot",   default = "refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c,cosmic96_coding,cosmic96_noncoding",
                  nargs = 1, type = "string",
                  help = "Analysis protocol of ANNOVAR. "
                         "Default: refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c,cosmic96_coding,cosmic96_noncoding")
parser.add_option("--annovarop",     dest = "annovarop",     default = "gx,r,f,f,f,f,f",
                  nargs = 1, type = "string",
                  help = "Analysis operation of ANNOVAR. "
                         "Default: gx,r,f,f,f,f,f")
parser.add_option("--annovar_xref",  dest = "annovar_xref",  default = "hg38db/omim/gene_xref.txt",
                  nargs = 1, type = "string",
                  help = "Path to cross-reference genome of ANNOVAR. "
                         "Default: hg38db/omim/gene_xref.txt")
options, arguments = parser.parse_args()

#===pre-check===
polyT_seq = "TTTTTTTTTTTT"
#===10x TSO for 5 prime protocol, https://kb.10xgenomics.com/hc/en-us/articles/360001493051-What-is-a-template-switch-oligo-TSO===
TSO_seq   = "TTTCTTATATGGG"
Termination = False

if not options.fq_f_name:
	print("FastQ file(s) is required!")
	Termination = True
else:
	options.fq_f_name = os.path.abspath(os.path.expanduser(options.fq_f_name))
	if os.path.isdir(options.fq_f_name):
		fq_list = glob.glob(os.path.join(options.fq_f_name, "*fq*"))
		if len(fq_list) > 0:
			print("Found " + str(len(fq_list)) + " FastQ files under " + options.fq_f_name + "\n")
	elif os.path.isfile(options.fq_f_name):
		print("Found FastQ file: " + str(len(options.fq_f_name)) + "\n")
	else:
		print("Cannot find given FastQ file: " + str(len(options.fq_f_name)) + "\n")
		Termination = True

if not options.ref_genome:
	print("Reference genome is required!")
	Termination = True
else:
	options.ref_genome = os.path.abspath(os.path.expanduser(options.ref_genome))
	if not os.path.isfile(options.ref_genome):
		print("Cannot find given reference genome file: " + str(len(options.ref_genome)) + "\n")
		Termination = True

if not options.idx_genome:
	print("Reference genome index file is required!")
	Termination = True
else:
	options.idx_genome = os.path.abspath(os.path.expanduser(options.idx_genome))
	if not os.path.isfile(options.idx_genome):
		print("Cannot find given reference genome index file: " + str(len(options.idx_genome)) + "\n")
		Termination = True

if options.exc_bed:
	options.exc_bed = os.path.abspath(os.path.expanduser(options.exc_bed))
	if not os.path.isfile(options.exc_bed):
		print("Cannot find given exclusive BED file: " + str(len(options.exc_bed)) + "\n")
		Termination = True

if options.isoquant:
	options.isoquant = os.path.abspath(os.path.expanduser(options.isoquant))
	if not os.path.isfile(options.isoquant):
		print("Cannot find IsoQuant at: " + options.isoquant)
		Termination = True

if options.annovar:
	options.annovar = os.path.abspath(os.path.expanduser(options.annovar))
	if not os.path.isdir(options.annovar):
		print("Cannot find ANNOVAR directory at: " + options.annovar)
		Termination = True

if str.lower(options.protocol) == "3p":
	options.pT = polyT_seq
elif str.lower(options.protocol) == "5p":
	options.pT = TSO_seq
elif str.lower(options.protocol) == "spatial":
	if not options.whitelist:
		print("Spatial-seq whitelist is required.\n")
		Termination = True
		if not os.path.isfile(options.whitelist):
			print("Cannot find whitelist file at: " + options.whitelist)
			Termination = True
else:
	print("Valid protocol (5p / 3p / spatial) is required.\n")
	Termination = True

if Termination:
	parser.print_help()
	sys.exit(1)

if not os.path.isdir('logs'):
	os.mkdir('logs')
wbc_cmd, exc_cmd = '', ''

oh = open("run_scNanoGPS.sh", "wt")
oh.write('#! /bin/bash' + "\n\n")
oh.write('P_DIR="' + os.path.dirname(os.path.abspath(__file__)) + '"' + "\n")
oh.write('FASTQ="' + os.path.abspath(os.path.expanduser(options.fq_f_name)) + '"' + "\n")
oh.write('REF_GENOME="' + os.path.abspath(os.path.expanduser(options.ref_genome)) + '"' + "\n")
oh.write('IND_GENOME="' + os.path.abspath(os.path.expanduser(options.idx_genome)) + '"' + "\n")
oh.write('GENOME_ANNOTATION="' + os.path.abspath(os.path.expanduser(options.gtf)) + '"' + "\n")
if options.whitelist:
	wbc_cmd = '--whitelist $WBC'
	oh.write('WBC="' + os.path.abspath(os.path.expanduser(options.whitelist)) + '"' + "\n")
if options.exc_bed:
	exc_cmd = '--exc_bed $EXC_BED'
	oh.write('EXC_BED="' + os.path.abspath(os.path.expanduser(options.exc_bed)) + '"' + "\n")
oh.write('ncores=' + options.ncores + "\n")
if options.isoquant:
	oh.write('ISOQUANT="' + os.path.abspath(os.path.expanduser(options.isoquant)) + '"' + "\n")
if options.annovar:
	oh.write('ANNOVAR="' + options.annovar + '"' + "\n")
	oh.write('ANNOVAR_DB="' + os.path.join(options.annovar, options.annovardb) + '"' + "\n")
	oh.write('ANNOVAR_GV="' + options.annovargv + '"' + "\n")
	oh.write('ANNOVAR_PROTOCOL="' + options.annovarprot + '"' + "\n")
	oh.write('ANNOVAR_OP="' + options.annovarop + '"' + "\n")
	oh.write('ANNOVAR_XREF="' + os.path.join(options.annovar, options.annovar_xref) + '"' + "\n")
oh.write('PT_SEQ="' + options.pT + '"' + "\n")
oh.write("\n")
oh.write('python3 $P_DIR/other_utils/read_length_profiler.py -i $FASTQ &> logs/run_read_length_profiler.log.txt &' + "\n")
oh.write('python3 $P_DIR/scanner.py -t $ncores -i $FASTQ --pT $PT_SEQ &> logs/run_scanner.log.txt' + "\n")
oh.write('python3 $P_DIR/assigner.py -t $ncores ' + wbc_cmd + ' &> logs/run_assigner.log.txt' + "\n")
oh.write('python3 $P_DIR/curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME ' + exc_cmd + ' &> logs/run_curator.log.txt' + "\n")
oh.write('python3 $P_DIR/reporter_expression.py -t $ncores --gtf $GENOME_ANNOTATION &> logs/run_reporter_expression.log.txt' + "\n")
if options.isoquant:
	oh.write('python3 $P_DIR/reporter_isoform.py -t $ncores --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION --isoquant $ISOQUANT &> logs/run_reporter_isoform.log.txt' + "\n")
if options.annovar:
	oh.write('python3 $P_DIR/reporter_SNV.py -t $ncores --ref_genome $REF_GENOME --annovar $ANNOVAR --annovar_db $ANNOVAR_DB --annovar_gv $ANNOVAR_GV --annovar_protocol $ANNOVAR_PROTOCOL --annovar_operation $ANNOVAR_OP --annovar_xref $ANNOVAR_XREF &> logs/run_reporter_SNV.log.txt' + "\n")
	oh.write('python3 $P_DIR/other_utils/parse_annovar_column.py -i scNanoGPS_res/annovar.hg38_multianno.vcf > scNanoGPS_res/annovar.hg38_multianno.tsv' + "\n")
mrg_bam_param = ""
if options.isoquant:
	mrg_bam_param = "--mrg_bam scNanoGPS_res/IsoQuant_res/merged.curated.minimap2.bam"
oh.write('python3 $P_DIR/reporter_summary.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION ' + mrg_bam_param + ' --qualimap_param "--java-mem-size=300G" &> logs/run_reporter_summary.log.txt' + "\n")
oh.close()

print("\n *** Please review or edit the script file 'run_scNanoGPS.sh' and then run it ! *** \n")


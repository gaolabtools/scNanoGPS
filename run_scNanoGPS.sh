#! /bin/bash

P_DIR="~/data/scNanoGPS_v1.1"
FASTQ="~/data/scNanoGPS_v1.1/example/fastq/"
REF_GENOME="~/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
IND_GENOME="~/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.mmi"
GENOME_ANNOTATION="~/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"
REF_LIQA="~/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/liqa/genes.refgene"
ncores=20
ANNOVAR="~/tools/annovar"
ANNOVAR_DB="~/tools/annovar/hg38db/"
ANNOVAR_GV="hg38"
ANNOVAR_PROTOCOL="refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c,cosmic96_coding,cosmic96_noncoding"
ANNOVAR_OP="gx,r,f,f,f,f,f"
ANNOVAR_XREF="~/tools/annovar/hg38db/omim/gene_xref.txt"

python3 $P_DIR/scanner.py -i $FASTQ -t $ncores
python3 $P_DIR/other_utils/read_length_profiler.py -i $FASTQ
python3 $P_DIR/assigner.py -t $ncores
python3 $P_DIR/curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME
python3 $P_DIR/reporter_expression.py -t $ncores --gtf $GENOME_ANNOTATION
python3 $P_DIR/reporter_isoform.py -t $ncores --liqa_ref $REF_LIQA
python3 $P_DIR/reporter_SNV.py -t $ncores --ref_genome $REF_GENOME --annovar $ANNOVAR --annovar_db $ANNOVAR_DB --annovar_gv $ANNOVAR_GV --annovar_protocol $ANNOVAR_PROTOCOL --annovar_operation $ANNOVAR_OP --annovar_xref $ANNOVAR_XREF
#python3 $P_DIR/other_utils/parse_annovar_column.py -i scNanoGPS_res/annovar.hg38_multianno.vcf > scNanoGPS_res/annovar.hg38_multianno.tsv
python3 $P_DIR/reporter_summary.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION


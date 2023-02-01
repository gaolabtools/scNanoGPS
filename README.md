# scNanoGPS
scNanoGPS is a computational toolkit for analyzing high throughput single cell nanopore sequencing data to detect Genotypes and Phenotype Simultaneously from same cells.  scNanoGPS includes 5 major steps: 1) **NanoQC** to perform quality control of the raw seqeucning data; 2) **Scanner** to scan and filter out reads that do not have expected adapater sequence patterns, i.e., TrueSeq Read 1 adapter sequence, TSO adaper sequence, poly (A/T)n block sequence,  Cell Barcodes (CB) and unique molecule identifier (UMI) sequence blocks; 3) **Assigner** to detect the list of true cell barcodes, merge cell barcodes with sequencing errors and assign raw reads into single cells; 4) **Curator** to detect reads with true UMIs and collapse them to make consensus sequences of individual molecules to curate sequencing errors on gene bodies; 5) **Reporter** to detect single cell transcriptomes, single cell gene isoforms and single cell mutations from consensus single cell long reads data.

# Index
- [Installation](#installation)
- [Step 1: NanoQC](#step-1-nanoqc)
- [Step 2: Scanner](#step-2-scanner)
- [Step 3: Assigner](#step-3-assigner)
- [Step 4: Curator](#step-4-curator)
- [Step 5: Reporter](#step-5-reporter)

# Installation
The scNanoGPS pipeline is built with python3. We recommend users to use anaconda/miniconda virtual environment to install it. Refer to [Anaconda turorial](https://docs.anaconda.com/anaconda/user-guide/tasks/switch-environment/) for environment building. <br />

### Build python3 virtual environment
- Example codes for creating and activating python3 environment on Linux-based OS:
  ```
  conda create -n scNanoGPS python=3 numpy scipy
  source activate scNanoGPS
  ```

### Install scNanoGPS and dependencies
- The scNanoGPS requires the following dependencies to work:
  ```
  - biopython 1.79
  - distance 0.1.3
  - matplotlib 3.5.2
  - pandas 1.4.2
  - pysam 0.19.0
  - seaborn 0.11.2
  ```

- Example codes for obtaining scNanoGPS from GitHub and installation of dependencies:
  ```
  git clone https://github.com/gaolabtools/scNanoGPS/
  cd scNanoGPS
  pip3 install -r requirements.txt
  ```

### Install other essential tools
scNanoGPS uses the following third party tools for mapping again genome reference, collapsing reads with same UMIs, and sumamrizing single cell gene expression, isoform, and SNV profiles.

- Example codes for installation of third party tools
  - [minimap2](https://lh3.github.io/minimap2/) ([GitHub](https://github.com/lh3/minimap2#installation), [Anaconda](https://anaconda.org/bioconda/minimap2))
    ```
    conda install -c bioconda minimap2
    ```
  - [Samtools](https://www.htslib.org/) ([GitHub](https://github.com/samtools/samtools), [Anaconda](https://anaconda.org/bioconda/samtools))
    ```
    conda install -c bioconda samtools
    ```
  - tabix ([Anaconda](https://anaconda.org/bioconda/tabix))
    ```
    conda install -c bioconda tabix
    ```
  - SPOA ([GitHub](https://github.com/rvaser/spoa), [Anaconda](https://anaconda.org/bioconda/spoa))
    ```
    conda install -c bioconda spoa
    ```
  - SubRead featureCounts ([SourceForge](http://subread.sourceforge.net/), [Anaconda](https://anaconda.org/bioconda/subread))
    ```
    # You can download and unzip pre-compiled binary file from https://sourceforge.net/projects/subread/files/subread-2.0.3/
    
    tar -xzf subread-2.0.3-<platform>.tar.gz
    
    # or install subread via anaconda
    
    conda install -c bioconda subread
    ```
  - LIQA ([GitHub](https://github.com/WGLab/LIQA))
    ```
    pip3 install liqa
    ```
  - Longshot ([GitHub](https://github.com/pjedge/longshot), [Anaconda](https://anaconda.org/bioconda/longshot))
    ```
    conda install -c bioconda longshot
    ```
  - [BCFtools](https://samtools.github.io/bcftools/howtos/index.html) ([GitHub](https://github.com/samtools/bcftools), [Anaconda](https://anaconda.org/bioconda/bcftools))
    ```
    conda install -c bioconda bcftools
    ```
  - [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/)
    ```
    # You can download ANNOVAR from https://www.openbioinformatics.org/annovar/annovar_download_form.php

    tar -xvf annovar.latest.tar.gz
    ```
  - (optional) [gffread](https://ccb.jhu.edu/software/stringtie/gff.shtml) ([GitHub](https://github.com/gpertea/gffread))
    ```
    git clone https://github.com/gpertea/gffread
    cd gffread
    make release
    ```
  - [Qualimap](http://qualimap.conesalab.org/) ([Anaconda](https://anaconda.org/bioconda/qualimap))
    ```
    # You can download Qualimap from http://qualimap.conesalab.org/
    
    unzip qualimap_v2.2.1.zip
    
    # or install qualimap via anaconda
    
    conda install -c bioconda qualimap
    ```

### Prepare reference genome and annotations
  - Reference genome<br />
    Users can obtain reference genome from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/), [Ensembl](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/), or any other autorities
    ```
    wget https://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    ```
  - Gene annotation (GTF/GFF)<br />
    Users can obtain reference gene annotations from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/), [Ensembl](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/), or any other autorities
    ```
    wget https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
    ```
    - Note: Many tools cannot use compressed GTF file. Please try to gunzip compress GTF.gz beforehand.
  
  - Index reference genome for minimap2<br />
    Prepare indexed genome for minimap2 to boost mapping.  Refer to the [Minimap2 instruction](https://github.com/lh3/minimap2#getting-started).<br />
    - Example code:
      ```
      minimap2 -x map-ont -d example/GRCh38_chr22.mmi example/GRCh38_chr22.fa.gz
      ```
  - Index reference genome for LIQA<br />
    This version of scNanoGPS uses LIQA as a default tool to calculate single cell isoforms. Refer to the [LIQA's manual](https://github.com/WGLab/LIQA/blob/master/doc/Usage.md#step-1-transforming-isoforms-to-compatible-matrix-based-on-reference-annotation-file) for complete instructions<br />
    - Example code:
      ```
      liqa -task refgene -ref example/GRCh38_chr22.gtf -format gtf -out example/GRCh38_chr22.liqa.refgene
      ```
  - Annotation tables for ANNOVAR<br />
    This version of scNanoGPS uses ANNOVAR to annotate single cell SNVs results, please refer to [ANNOVAR's webpage](https://annovar.openbioinformatics.org/en/latest/) for more information.
     - Example codes:
    ```
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene hg38db/
    perl annotate_variation.pl -buildver hg38 -downdb cytoBand hg38db/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome hg38db/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 hg38db/
    perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c hg38db/
    ```

# Step 1: NanoQC
### Read length distribution <br />
scNanoGPS contains a script named “read_length_profiler.py” to compute the raw read lengths of all reads.  The script can either read through individual FastQ/Fast5 files or all FastQ/Fast5 files under a given folder. The raw read length histogram is drawn accordingly.

  - Manual of read_length_profiler.py
    ```
    python3 other_utils/read_length_profiler.py -h
    Usage: read_length_profiler.py [options]

    Options:
      -h, --help     show this help message and exit
      -i FQ_F_NAME   File or folder name of reads.
      -d O_DIR       Output directory name. Default: scNanoGPS_res
      -f FIG_NAME    Read length histogram file name. Default: read_length.png
      -o O_NAME      Read length file name. Default: read_length.tsv.gz
      --fig_w=FIG_W  Width of figure (inch). Default: 12
      --fig_h=FIG_H  Height of figure (inch). Default: 7
    ```

  - Example code:
    ```
    python3 other_utils/read_length_profiler.py -i example/fastq/
    ```

### FastQC (optional) <br />
Per the experimental design of read architecture, the TruSeq Read1, cell barcode (CB), unique molecular identifier (UMI), and polyA tail are expected to locate either in the first or the last 100 nucleotide range of each Nanopore read.  Users can use FastQC to check the qualities of the first and last 100 nucleotides of individual Nanopore reads to draw the per-base quality score boxplot.

# Step 2: Scanner
This step of scNanoGPS pipeline is executed by a python script called “scanner.py”. This script scans for both TruSeq Read1 and polyA tail of the reads. Scanning for other sequence modules are optional.
To boost the scanning speed, we scan the first and last 100 nucleotides of reads to recognize TruSeq Read 1 and PolyA. Following by recognition, Scanner extracts CBs and UMIs which are neighbored by TruSeq Read 1 and Poly(A/T)n sequence blocks. Then the Scanner outputs two different files. One is a processed FastQ file holding the insert sequences without TruSeq Read 1, CB, and UMI sequences. The CBs of individual reads are moved to the read names as tags.  The other file is table named, “barcode_list.tsv”, storing reads information including read names, CBs, UMIs, and others.

- Manual of scanner.py
  ```
  python3 scanner.py -h
  Usage: scanner.py [options]

  Options:
    -h, --help            show this help message and exit
    -i FQ_F_NAME          * Required ! Input FastQ/Fast5 file name, or directory
                          containing multiple input files. Support
                          fastq/fq/fastq.gz/fq.gz/fast5 format.
    -o FQ_O_NAME          Output fastq file name. Could be either in .fastq or
                          .gz format.Default: processed.fastq.gz
    -d O_DIR              Output directory name. Must give a new directory name
                          to prevent accidental overwriting ! Default:
                          scNanoGPS_res
    -b BC_F_NAME          Output cell barcode list file. Default:
                          barcode_list.tsv.gz
    -t NCORES             Number of cores for program running. Default: 1
    --log=LOG_F_NAME      Program log file. This file stores program running
                          parameters and counting details. Default:
                          scanner.log.txt
    --a5=ADAPTOR_FIVE_P   Sequence of 5'-adaptor. Default:
                          AAGCAGTGGTATCAACGCAGAGTACAT
    --a3=ADAPTOR_THREE_P  Sequence of 3'-adaptor. Default:
                          CTACACGACGCTCTTCCGATCT
    --pT=POLYT            Reverse-complement sequence of polyA. Default:
                          TTTTTTTTTTTT
    --lCB=BC_LEN          Length of cell barcode. Default: 16
    --lUMI=UMI_LEN        Length of UMI. Default: 12
    --min_read_length=MIN_READ_LENGTH
                          Minimal read length. Default: 200
    --editing_distance=ALLOW_EDITING_DISTANCE
                          Editing distance for cell barcode detection. Default:
                          2
    --matching_threshold=MATCHING_PERCENTAGE
                          Matching threshold for alignment search. Default: 0.7
    --score_threshold=SCORING_THRESHOLD
                          Scoring threshold for alignment search. Default: 0.4
    --batching_no=BATCH_NO
                          Number of reads for batch processing. Default: 1000
    --scanning_region=SCAN_REGION
                          Region length for adaptor scanning. Default: 100
    --debug_mode=DEBUG_MODE
                          Debug mode switch. Default: False
    --penalty_matching=DP_MA
                          Dynamic programming matching penalty. Default: 2
    --penalty_mismatching=DP_MI
                          Dynamic programming mismatching penalty. Default: -3
    --penalty_gap_opening=DP_GO
                          Dynamic programming gap opening penalty. Default: -5
    --penalty_gap_extention=DP_GE
                          Dynamic programming gap extention penalty. Default: -2
  ```

- Example code:
  ```
  python3 scanner.py -i example/fastq/ -t 2
  ```

# Step 3: Assigner
This step of scNanoGPS pipeline is executed by a python script called “assigner.py”. This script is designed for CB collapsing and estimation of the optimal CB number without guidance of 10X short-read sequencing data or any CB whitelist.
To estimate the number of optimal CB, we use edge detection strategy to find out the point where has dramatical signal dropping (Fig. 1b). The detailed method is that the assigner first calculates the supporting UMI number to every CB, and sorts the CB list by UMI number in decreasing order. Following by computing the partial derivatives (slopes) per CB in log10 scale, the medium number of slope changes in log10 scale are computed per 0.001 log10 tick. Then the maximal medium log10 slope change is selected, and where is the crude anchoring for following processes. To fathom the fully signal dropping point and include more useful CBs, we allow 10% more signal in log10 scale.
Next, the script collapses CBs which have similar sequences. Previous study shows that the most accurate criterial for CB and UMI collapsing in Illumina samples are three and two Levenshtein Distance (LD), respectively. Here we use two LD to merge similar CB as Refinery Local Optimization. Then a list of representative CB having sufficient supporting read is generated.

- Manual of assigner.py
  ```
  python3 assigner.py -h
  Usage: assigner.py [options]

  Options:
    -h, --help            show this help message and exit
    -i INPUT              Cell barcode list file. Could be either in .tsv or
                          .tsv.gz format. Default: barcode_list.tsv.gz
    -o OUTPUT             Cell barcode counting file name. Default:
                          CB_counting.tsv.gz
    -d O_DIR              Output directory name. Default: scNanoGPS_res
    -t NCORES             Number of cores for program running. Default: 1
    --log=LOG_F_NAME      Log file name. Default: assigner.log.txt
    --lCB=BC_LEN          Length of cell barcode. Default: 16
    --CB_no_ext=CB_NO_EXT
                          Increasing CB number on log scale. Default: 0.1
    --CB_log10_dist_o=CB_LOG10_DIST_O
                          File name used for plotting log10 UMI number
                          distribution. Default: CB_log10_dist.png
    --CB_mrg_thr=CB_MRG_THR
                          Threshold of distance for merging cell barcodes.
                          Default: 2
    --CB_mrg_dist=CB_MRG_DIST
                          File name for distance matrix of merging cell
                          barcodes. Default: CB_merged_dist.tsv.gz
    --CB_mrg_o=CB_MRG_O   File name for merged cell barcodes. Default:
                          CB_merged_list.tsv.gz
  ```

- Example code:
  ```
  python3 assigner.py -t 2
  ```

# Step 4: Curator
This step of scNanoGPS pipeline is executed by a python script called “curator.py”. This script is used for demultiplexing, filtering, reference genome mapping and re-mapping, and UMI collapsing.

The master FastQ file of all cells is demultiplexed according to the true CB list determined by Assigner into single cell FastQ files each representing one cell. Curator then maps individual FastQ files onto given reference genome by Mimimap2 under splice mode. Chimeric reads from different chromosomes are filtered out in this step (fusion gene detection functionality is under development). Next, Curator scans UMIs through their full length reads by their mapped genomic orders.  UMIs that within 2 LD and mapping to same genomic coordinates are considered as same UMI barcodes. To further accommodate possible small indels (<5bp) that causes minor drifting of mapping coordinates, Curator allows 5bp differences to buffer these sequencing errors. To perform parallel computing, the reads scanning is placed into batches based on their genomic coordinates.  The reads that share same UMIs are collapsed to generate consensus sequences of individual molecules using software, SPOA.  Finally, we re-map the consensus sequences of single cells onto reference genome by using Minimap2 under splice mode. There is a portion of reads that are singletons having only one UMI, which is mapped previsouly.  We merged both singleton BAMs with consensus BAMs to formal a final BAMs as curated data.

- Manual of curator.py
  ```
  python3 curator.py -h 
  Usage: curator.py [options]

  Options:
    -h, --help            show this help message and exit
    --fq_name=FQ_NAME     Processed fastq file name. Default: processed.fastq.gz
    -b BC_LIST            Output cell barcode list file. Default:
                          barcode_list.tsv.gz
    --CB_count=CB_COUNT   Cell barcode counting file name. Default:
                          CB_counting.tsv.gz
    --CB_list=CB_LIST     File name for merged cell barcodes. Default:
                          CB_merged_list.tsv.gz
    --ref_genome=REF_GENOME
                          * Required ! File for reference genome.
    --idx_genome=IDX_GENOME
                          Path to the Minimap2 genome index. Program will use
                          reference genome if no Minimap2 genome index given.
                          Default: None
    -d O_DIR              Output directory name. Default: scNanoGPS_res
    --tmp_dir=TMP_DIR     Temporary folder name. Default: tmp
    -t NCORES             Number of cores for computing. Default: 1
    --log=LOG_F_NAME      Log file name. Default: curator.log.txt
    --umi_ld=UMI_LD       Levenshtein distance for merging UMI. Default: 2
    --keep_meta=KEEP_META
                          Keep meta data, e.g. bam files, for futher checking.
                          Default: None
    --softclipping_thr=SOFTCLIPPING_THR
                          Threshold for softclipping. Default: 0.8
    --minimap2=MINIMAP2   Path to minimap2. Default: minimap2
    --samtools=SAMTOOLS   Path to samtools. Default: samtools
    --spoa=SPOA           Path to spoa. Default: spoa
  ```

- Example code:
  ```
  python3 curator.py -t 2 --ref_genome example/GRCh38_chr22.fa.gz --idx_genome example/GRCh38_chr22.mmi
  ```

# Step 5: Reporter
Lastly, scNanoGPS contains a set of reporter scripts for generating multi-omics profiles from same single cells with Nanopore long-read sequencing data.
This version of scNanoGPS detects the gene expression, isoform, and single nucleotide variations (SNVs) profiles by using FeatureCounts, LIQA, and longshot, respectively.

### 5.1 Single cell gene expression profile
- Manual of reporter_expression.py
  ```
  python3 reporter_expression.py -h
  Usage: reporter_expression.py [options]

  Options:
    -h, --help            show this help message and exit
    -d O_DIR              Output directory name. Default: scNanoGPS_res
    --tmp_dir=TMP_DIR     Temporary folder name. Default: tmp
    --gtf=GTF             * Required ! GTF file for expression calling.
    -o O_NAME             Counting table name. Default: matrix.tsv
    --log=LOG_F_NAME      Log file name.Default: reporter_expression.log.txt
    -t NCORES             Number of cores for program running. Default: 1
    --min_gene_no=MIN_GENE_NO
                          Minimal number of gene per cell. Default: 300
    --sel_bc_o=SEL_BC_O   Filtered cell barcode list. Default:
                          filtered_barcode_list.txt
    --featurecounts=FEATURECOUNTS
                          Path to featureCounts.Default: featureCounts
  ```

- Example code:
  ```
  # Please add featureCounts into your path to use which command
  
  python3 reporter_expression.py -t 2 --gtf example/GRCh38_chr22.gtf --featurecounts $(which featureCounts)
  
  # or
  
  python3 reporter_expression.py -t 2 --gtf example/GRCh38_chr22.gtf --featurecounts /path/to/subread_folder/bin/featureCounts
  ```

### 5.2 Single cell isoform profile
- Manual of reporter_isoform.py
  ```
  python3 reporter_isoform.py -h
  Usage: reporter_isoform.py [options]

  Options:
    -h, --help           show this help message and exit
    -d O_DIR             Output directory name. Default: scNanoGPS_res
    --tmp_dir=TMP_DIR    Temporary folder name. Default: tmp
    --CB_file=CB_FILE    File name for filtered barcode list. Default:
                         filtered_barcode_list.txt
    --liqa_ref=LIQA_REF  * Required ! Reference of LIQA.
    -o O_NAME            Counting table name. Default: matrix_isoform.tsv
    --log=LOG_F_NAME     Log file name. Default: reporter_isoform.log.txt
    -t NCORES            Number of cores for program running. Default: 1
    --liqa=LIQA          Program name of LIQA. Default: liqa
    --liqa_log=LIQA_LOG  Suffix of LIQA output file. Default: .liqa.log
    --liqa_o=LIQA_O      Suffix of LIQA output file. Default: .liqa.tsv
  ```

- Example code:
  ```
  # Please add liqa into your path to use which command
  
  python3 reporter_isoform.py -t 2 --liqa_ref example/GRCh38_chr22.liqa.refgene
  
  # or
  
  python3 reporter_isoform.py -t 2 --liqa $(which liqa) --liqa_ref example/GRCh38_chr22.liqa.refgene
  
  # or
  
  python3 reporter_isoform.py -t 2 --liqa /path/to/liqa --liqa_ref example/GRCh38_chr22.liqa.refgene
  ```

### 5.3 single cell SNV profile
- Manual of reporter_SNV.py
  ```
  python3 reporter_SNV.py -h
  Usage: reporter_SNV.py [options]

  Options:
    -h, --help            show this help message and exit
    -d O_DIR              Output directory name. Default: scNanoGPS_res
    --tmp_dir=TMP_DIR     Temporary folder name. Default: tmp
    --CB_file=CB_FILE     File name for filtered barcode list. Default:
                          filtered_barcode_list.txt
    --ref_genome=REF_GENOME
                          * Required ! File for reference genome.
    --longshot_min_cov=LONGSHOT_MIN_COV
                          Minimal coverage for longshot. Default: 2
    --longshot_min_alt_count=LONGSHOT_MIN_ALT_COUNT
                          Minimal alternative count for longshot. Default: 2
    --longshot_o=LONGSHOT_O
                          Prefix of LongShot output VCF file. Default:
                          longshot.output
    -o O_NAME             Result SNV matrix file name. Must be ended with
                          .vcf.gz. Default: matrix_SNV.vcf.gz
    --log=LOG_F_NAME      Log file name. Default: reporter_SNV.log.txt
    --o_snv_l=O_SNV_L     Filtered SNVs position list. Default:
                          filtered_SNV_position_list.tsv.gz
    --o_snv_dp=O_SNV_DP   SNVs depth matrix. Default: matrix_SNV_dp.tsv.gz
    -t NCORES             Number of cores for program running. Default: 1
    --samtools=SAMTOOLS   Path to samtools. Default: samtools
    --bcftools=BCFTOOLS   Path to bcftools. Default: bcftools
    --tabix=TABIX         Path to tabix. Default: tabix
    --longshot=LONGSHOT   Path to longshot. Default: None
    --prevalence=PREVALENCE
                          SNV prevalence rate. Default: 0.01
    --min_read_quality=MIN_READ_QUALITY
                          Minimal read quality. Default: 0
    --annovar=ANNOVAR     Path to Annovar. Default: None
    --annovar_db=ANNOVAR_DB
                          Path to Annovar database. Default: None
    --annovar_gv=ANNOVAR_GV
                          Annovar database genome version. Default: hg38
    --annovar_protocol=ANNOVAR_PROTOCOL
                          Protocols of Annovar database. Default:
                          refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c
    --annovar_operation=ANNOVAR_OPERATION
                          Operation of Annovar database. Default: gx,r,f,f,f
    --annovar_xref=ANNOVAR_XREF
                          Path to Omim xref. Default: None
  ```

- Example code for variant detection with LongShot and annotation with ANNOVAR (please make sure to include ANNOVAR in your path):
  ```
  # Please add annovar into your path to use which command
  
  python3 reporter_SNV.py -t 2 --ref_genome example/GRCh38_chr22.fa.gz --annovar $(which table_annovar.pl) --annovar_db path/to/annovar_db/
  
  # or
  
  python3 reporter_SNV.py -t 2 --ref_genome example/GRCh38_chr22.fa.gz --annovar /path/to/annovar/table_annovar.pl --annovar_db path/to/annovar_db/
  ```

- Example code for variant calling without ANNOVAR annotations (optional, not suggested):
  ```
  python3 reporter_SNV.py -t 2 --ref_genome example/GRCh38_chr22.fa.gz
  ```

### 5.4 Generate final summary table
- Manual of reporter_summary.py
  ```
  python3 reporter_summary.py -h
  Usage: reporter_summary.py [options]

  Options:
    -h, --help            show this help message and exit
    -d O_DIR              Output directory name. Default: scNanoGPS_res
    --tmp_dir=TMP_DIR     Temporary folder name. Default: tmp
    --scanner_log=SCANNER_LOG
                          Scanner log file name. Default: scanner.log.txt
    --bc_f=BC_F           Cell barcode list file. Default: barcode_list.tsv.gz
    --read_len_f=READ_LEN_F
                          Scanner log file name. Default: read_length.tsv.gz
    --CB_file=CB_FILE     File name for filtered barcode list. Default:
                          filtered_barcode_list.txt
    --ref_genome=REF_GENOME
                          * Required ! File for reference genome.
    --gtf=GTF             * Required ! Genome annotation file GTF.
    --log=LOG_F_NAME      Log file name. Default: summary.txt
    --samtools=SAMTOOLS   Path to samtools. Default: samtools
    --qualimap=QUALIMAP   Path to qualimap. Default: qualimap
  ```

- Example code for generating final summary
  ```
  python3 reporter_summary.py --ref_genome example/GRCh38_chr22.fa.gz --gtf example/GRCh38_chr22.gtf --qualimap path/to/qualimap_folder/qualimap
  ```

- Example summary table
  ```
  Read yield:                  98349656
  Valid read number:           76630302
  Detecting rate:              77.92%

  Median read length:          921.0
  Mean read length:            1146.76
  Maximal read length:         190885
  Median read quality:         21.12
  Mean read quality:           20.47

  Cell number:                 3470
  Raw reads per cell:          28342.84
  UMI counts:                  27292859
  Mean UMI counts per cell:    7865.38
  Median UMI counts per cell:  5212.5

  Exonic:                      13.91%
  Intronic:                    70.2%
  Intergenic:                  15.89%
  ```

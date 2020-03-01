Custom analysis software for ChIP-seq of DNA repair proteins after Cas9-mediated DSBs
====

### Software requirements
- python 2.7 (https://www.anaconda.com/distribution/)
    - Anaconda distribution comes with the required numpy and scipy libraries
- pysam (https://pysam.readthedocs.io/en/latest/installation.html)
- bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- samtools (http://www.htslib.org/download/)


### Generate processed BAM files (if starting from raw FASTQ files)
1. Download sequencing reads in FASTQ format from SRA
2. Download either the human (hg38) or mouse (mm10) prebuilt bowtie2 indices
    - hg38 (ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)
    - mm10 (ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip)
    a. move to the corresponding folders named `hg38_bowtie2/` or `mm10_bowtie2/`
3. Download either human (hg38) or mouse (mm10) genome assembly in FASTA format
    - hg38.fa (https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
    - mm10.fa (http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)
    a. extract from gz files
    b. move to the corresponding folders named `hg38_bowtie2/` or `mm10_bowtie2/`
4. Generate FASTA file indices
    - `samtools faidx hg38_bowtie2/hg38.fa`
    - `samtools faidx mm10_bowtie2/mm10.fa`
5. Align reads to either hg38 or mm10 using bowtie2
    - will only use hg38 as demonstration from now on, simply replace with mm10 as desired.
   ```
   bowtie2 -q --local -X 1000 hg38_bowtie2/hg38 \
   -1 /path/to/read1.fastq \
   -2 /path/to/read2.fastq \ 
   -S /path/to/output.sam
   ```
6. Convert from SAM to BAM, filter for mapping quality >=25 and singleton reads
    ```
   samtools view -h -S -b -F0x08 -q25 \
   /path/to/output.sam \
   /path/to/output_unsorted.bam
   ```
7. Add mate score tags to ensure that paired-end reads contain correct information about the mate reads
   ```
   samtools fixmate -m \
   /path/to/output_unsorted.bam \
   /path/to/output_fixmate.bam
   ```
8. Sort BAM file entries by genomic position
   ```
   samtools sort \
   /path/to/output_fixmate.bam \
   /path/to/output_sorted.bam
   ```
9. Remove potential PCR duplicates
   ```
   samtools markdup -r \
   /path/to/output_sorted.bam \
   /path/to/output_rmdup.bam
   ```
10. (optional) Subset BAM files to normalize by # of reads in a set  

11. Index sorted BAM file
    ``` samtools index /path/to/output_final.bam ```

### Start from pre-processed BAM files
We have included pre-processed sequencing reads in BAM format, which is the output of the previous section.
1. Download pre-processed paired-end reads in BAM format from SRA.
2. Move the downloaded data for MRE11 and gH2AX ChIP-seq to the desired folder.
3. Open the following files:
    a. `chipseq_mre11.py`, which runs the analysis code for MRE11 ChIP-seq at ACTB.
    b. `chipseq_h2ax.py`, which runs the analysis code for gH2AX ChIP-seq at ACTB.
4. For each file, change the `base` variable to be the path to the directory that holds the BAM files.
5. Create a new folder to hold the output of the analysis, set its path to `base_a`.
6. Ensure that the file paths in `chipseq_mre11.py` and `chipseq_h2ax.py` are correct before running both scripts.

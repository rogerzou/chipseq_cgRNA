Custom analysis software for ChIP-seq of DNA repair proteins after CRISPR/Cas9-mediated DSBs
====

## Software requirements
- [python 2.7](https://www.anaconda.com/distribution/) (Anaconda's python distribution comes with the required numpy and scipy libraries)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- Ensure that both `samtools` and `bowtie2` are added to path and can be called directly from bash

## Data requirements
- The data for MRE11 and γH2AX ChIP-seq before/after Cas9 activation with light can be downloaded from [Here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA609749).
- The data from DISCOVER-seq ([Wienert & Wyman et al, Science, 2019](https://www.ncbi.nlm.nih.gov/pubmed/31000663)) for comparison can be downloaded from [Here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA509652).

## Generate processed BAM files (if starting from raw FASTQ files)
1. Download sequencing reads in FASTQ format from SRA
2. Download either the human or mouse prebuilt bowtie2 indices
    - [Human hg38](http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)
    - Mouse mm10 (ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip)
    - move to the corresponding folders named `hg38_bowtie2/` or `mm10_bowtie2/`
3. Download either human (hg38) or mouse (mm10) genome assembly in FASTA format
    - [hg38.fa](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
    - [mm10.fa](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)
    - extract from gz files
    - move to the corresponding folders named `hg38_bowtie2/` or `mm10_bowtie2/`
4. Generate FASTA file indices
    - `samtools faidx hg38_bowtie2/hg38.fa`
    - `samtools faidx mm10_bowtie2/mm10.fa`
5. Modify the first lines of bash script `process_reads.sh`
    - Fill the list variable `filelist` with paths to each **sample name** for processing. Note that the actual file paths include the **sample name** followed by "\_1.fastq" or "\_2.fastq", denoting read1 or read2 of paired-end reads, respectively.
    - Fill in the path to the indexed genome denoted by variable `genomepath`.
6. Run the following code snippet, where `-p` denotes the number of samples to process in parallel; modify accordingly. This performs genome alignment, filtering, sorting, removal of PCR duplicates, indexing, and sample statistics output. This step takes less than one day to complete on our Intel i7-8700K, 32GB RAM desktop, though speed appears to be bottlenecked by read/writes to disk.
    ```
    bash process_reads.sh -p 6
    ```
7. *(optional)* For fair comparison between different time points in a timeseries, we subset reads from all relevant samples to the sample with the fewest reads of the set. Run the following code snippet, where `-s` inputs the number of mapped reads to subset for each sample; modify accordingly. This step takes less than 10 minutes.
    ```
    bash process_reads.sh -p 6 -s 24400000
    ```
    - for MRE11, subsetted to 24,400,000 reads for both replicates.
    - for γH2AX, subsetted to 43,275,829 reads for replicate 1 and 62,271,079 reads for replicate 2.

## Start from pre-processed BAM files
In addition to raw paired-end reads in FASTQ format, we have also uploaded pre-processed sequencing reads in BAM format to SRA. These are the output of the previous section. It is highly recommended to start from these BAM files.
1. Download pre-processed paired-end reads in BAM format from SRA.
2. Move the downloaded data for MRE11 and γH2AX ChIP-seq to the desired folder.
3. If not already done, index the downloaded BAM files:
    ```
    samtools index /path/to/output.bam
    ```

## chipseq_mre11.py
#### This script runs the code for analyzing MRE11 ChIP-seq data after activation of Cas9/cgRNA targeting ACTB
1. Open `chipseq_mre11.py`, set `base` variable to be the path to the directory that holds the BAM files.
2. Create a new folder to hold the output of the analysis, set its path to `base_a`.
3. Ensure that all file names are correct (if the BAM files were directly downloaded from SRA, they should be), then run script.

## chipseq_h2ax.py
#### This script runs the code for analyzing γH2AX ChIP-seq data after activation of Cas9/cgRNA targeting ACTB
1. Open `chipseq_h2ax.py`, set `base` variable to be the path to the directory that holds the BAM files.
2. Create a new folder to hold the output of the analysis, set its path to `base_a`.
3. Ensure that all file names are correct (if the BAM files were directly downloaded from SRA, they should be), then run script.

## discoverseq_mre11.py
#### This scripts runs the code for analyzing MRE11 ChIP-seq data from [Wienert & Wyman et al (Science, 2019)](https://www.ncbi.nlm.nih.gov/pubmed/31000663)
1. FASTQ reads with the following SRA run accession codes (SRR) were downloaded from [Here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA509652).
    ```
    SRR8550692, SRR8550673, SRR8550703, SRR8550680, SRR8550681, SRR8550704, SRR8550684, SRR8550705, SRR8550693, SRR8550695, SRR8553800, SRR8553810, SRR8553804, SRR8553806
    ```
2. Generate BAM files from raw FASTQ reads following instructions from the previous sections. These files will be used in this script.
3. Ensure that the file names are correctly referenced in script, then run script.

## discoverseq_others.py
#### This scripts runs the code for analyzing ChIP-seq against multiple repair factors from [Wienert & Wyman et al (Science, 2019)](https://www.ncbi.nlm.nih.gov/pubmed/31000663)
1. FASTQ reads with the following SRA run accession codes (SRR) were downloaded from [Here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA509652).
    ```
    SRR8550677, SRR8550696, SRR8550679, SRR8550694, SRR8550682, SRR8550699, SRR8550678, SRR8550697, SRR8550698, SRR8550690
    ```
2. Generate BAM files from raw FASTQ reads following instructions from the previous sections. These files will be used in this script.
3. Ensure that the file names are correctly referenced in script, then run script.

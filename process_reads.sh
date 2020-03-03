#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
"/mnt/d/cgRNA_SRA/mre11_actb_00m_rep1" \
"/mnt/d/cgRNA_SRA/mre11_actb_02m_rep1" \
"/mnt/d/cgRNA_SRA/mre11_actb_05m_rep1" \
"/mnt/d/cgRNA_SRA/mre11_actb_15m_rep1" \
"/mnt/d/cgRNA_SRA/mre11_actb_30m_rep1" \
"/mnt/d/cgRNA_SRA/mre11_actb_60m_rep1" \
"/mnt/d/cgRNA_SRA/mre11_actb_00m_rep2" \
"/mnt/d/cgRNA_SRA/mre11_actb_02m_rep2" \
"/mnt/d/cgRNA_SRA/mre11_actb_05m_rep2" \
"/mnt/d/cgRNA_SRA/mre11_actb_15m_rep2" \
"/mnt/d/cgRNA_SRA/mre11_actb_30m_rep2" \
"/mnt/d/cgRNA_SRA/mre11_actb_60m_rep2" \
"/mnt/d/cgRNA_SRA/h2ax_actb_00m_rep1" \
"/mnt/d/cgRNA_SRA/h2ax_actb_02m_rep1" \
"/mnt/d/cgRNA_SRA/h2ax_actb_05m_rep1" \
"/mnt/d/cgRNA_SRA/h2ax_actb_15m_rep1" \
"/mnt/d/cgRNA_SRA/h2ax_actb_30m_rep1" \
"/mnt/d/cgRNA_SRA/h2ax_actb_60m_rep1" \
"/mnt/d/cgRNA_SRA/h2ax_actb_00m_rep2" \
"/mnt/d/cgRNA_SRA/h2ax_actb_02m_rep2" \
"/mnt/d/cgRNA_SRA/h2ax_actb_05m_rep2" \
"/mnt/d/cgRNA_SRA/h2ax_actb_15m_rep2" \
"/mnt/d/cgRNA_SRA/h2ax_actb_30m_rep2" \
"/mnt/d/cgRNA_SRA/h2ax_actb_60m_rep2" \
)
# Enter path to indexed genome
genomepath="/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38"
##########################################


# processing arguments, proceed with bioinformatic pipeline
main() {
  numthreads=1                      # default number of parallel processes
  while getopts 'p:s:' opt; do      # pass number of reads to subset via -s option
    case "$opt" in
      p) numthreads="$OPTARG"
        ;;
      s) numsubset="$OPTARG"
        ;;
      \?) echo "Usage: $(basename $0) [-p number of threads] [-s number of reads to subset]"
      exit 1
        ;;
    esac
  done
  echo "Number of parallel processes: $numthreads"
  for file in "${filelist[@]}"; do  # process each file
    ((i=i%numthreads)); ((i++==0)) && wait
    process ${genomepath} ${file} ${numsubset} &
  done
}


# subset BAM file only if [-s] argument is included, otherwise,perform main bioinformatics pipeline
process() {
  if [ ! -z ${3+x} ] ; then
    subset $3 $2
  else
    align2bam $1 $2
  fi
}


# main bioinformatics pipeline (alignment to indexing and read statistics)
align2bam() {

  # Align reads to either hg38 or mm10 using bowtie2
  bowtie2 -p 6 -q --local -X 1000 -x $1 \
  -1 "$2_1.fastq" -2 "$2_2.fastq" -S "$2.sam" ;

  # Convert from SAM to BAM, filter for mapping quality >=25 and singleton reads
  samtools view -h -S -b -F 0x08 -q 25 \
  "$2.sam" > "$2_unsorted.bam" ;

  # Add mate score tags to ensure that paired-end reads contain correct information about the mate reads
  samtools fixmate -m \
  "$2_unsorted.bam" "$2_fixmate.bam" ;

  # Sort BAM file entries by genomic position
  samtools sort \
  "$2_fixmate.bam" > "$2_sorted.bam" ;

  # Remove potential PCR duplicates
  samtools markdup -r \
  "$2_sorted.bam" "$2_rmdup.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$2_rmdup.bam" ;

  # Retrieve read count statistics
  samtools flagstat \
  "$2_rmdup.bam" > "$2_flagstats.txt"

}


# normalize by mapped read counts, only if script has been run before without [-s] flag so that "*_rmdup.bam" files have been generated in align2bam()
subset() {

  if test -f "$2_rmdup.bam"; then       # only run program if *_rmdup.bam exists

    # Subset BAM file by number of mapped reads
    total=$( samtools view -F 0x04 -c "$2_rmdup.bam" )
    frac=$(echo "$1/$total" | bc -l)
    if (( $(echo "$frac >= 1" | bc -l) )) ; then
      cp "$2_rmdup.bam" "$2_subset.bam"
    else
      samtools view -bs ${frac} \
      "$2_rmdup.bam" > "$2_subset.bam"
    fi

    # Index the subsetted BAM file
    samtools index \
    "$2_subset.bam" ;

    # Delete index output for rmdup
    rm "$2_rmdup.bam.bai"

    # Retrieve read count statistics
    samtools flagstat \
    "$2_subset.bam" > "$2_flagstats.txt"

  else
    echo "ERROR: *_rmdup.bam does not exist. First run process_reads.sh without the -s flag"
  fi
}


main "$@"; exit

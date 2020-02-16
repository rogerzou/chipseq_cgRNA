# -*- coding: utf-8 -*-
""" Useful functions for ChIP-seq analysis after Cas9 cleavage
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

from collections import defaultdict
import pysam
import os
import re
import numpy as np
import csv
from scipy import stats


def status_statement(current, final, count, chr=None):
    """ Print progress statements for long processes

    :param current: current iteration number
    :param final: total number of iterations
    :param count: number of print statements total
    :param chr: print chromosome information of current progress

    """
    if current % int(final/count) == 0:
        if chr is None:
            print("Processed %i out of %i" % (current, final))
        else:
            print("Processed %i out of %i in %s" % (current, final, chr))


def hg38_generator():
    """ Generate the number of base pairs for each chromosome

    :return generator that outputs the number of base pairs for each chromosome
             in the format: [chr7, 159345973]

    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/hg38.sizes", 'r') as f:   # parse sizes for each chromosome of hg38
        for row in csv.reader(f, dialect='excel', delimiter='\t'):  # iterate over each chromosome
            yield row


def read_pair_generator(bam, region_string=None):
    """ Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.

    :param bam: pysam AlignmentFile loaded with BAM file that contains paired-end reads
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :return generator of read pairs with the following format: [read1, read2]

    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_pair_align(read1, read2):
    """ Extract read pair locations as a fragment oriented in increasing chromosome coordinates

    :param read1: read #1 of pair in pysam AlignedSegment format
    :param read2: read #2 of pair in pysam AlignedSegment format
    :return 4-item array in the following format: [fragA-start, fragA-end, fragB-start, fragB-end]
            with monotonically increasing chromosome coordinates

    """
    r1pos = [x+1 for x in read1.positions]
    r2pos = [x+1 for x in read2.positions]
    if read1.mate_is_reverse and r1pos[0] < r2pos[0]:  # read1 is earlier
        read = [r1pos[0], r1pos[-1], r2pos[0], r2pos[-1]]
        # print("%s--%s>  <%s--%s" % tuple(read))
    elif read2.mate_is_reverse and r2pos[0] < r1pos[0]:  # read2 is earlier
        read = [r2pos[0], r2pos[-1], r1pos[0], r1pos[-1]]
        # print("%s--%s>  <%s--%s" % tuple(read))
    else:
        read = []
        print("FAIL")
    return read


def get_read_subsets(filein, fileout, region_string, target):
    """ Categorize fragments (paired reads) by its relationship to the target (cut) site

    For a target (cut) site, subset the fragments within region_range to
    - fragments that span the cut site
    - fragments that do not -> fragments that start from the cut site and span left vs right

    :param filein: BAM file containing paired-end reads
    :param fileout: base output file name with extension (.bam) omitted - function generates
                    multiple files where each file only contains reads from a specific category
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :param target: coordinate within region_string that denotes the 4th nucleotide from PAM

    """
    bam = pysam.AlignmentFile(filein, 'rb')
    countM = 0      # fragments that span
    countN = 0      # fragments that don't span
    countL = 0      # fragments that don't span, on left
    countR = 0      # fragments that don't span, on right
    countS1 = 0     # fragments where read1 is sequenced
    countS2 = 0     # fragments where read2 is sequenced
    countI = 0      # reads with indels
    counter = 0     # count total number of reads
    fileM = fileout + "_M.bam"
    fileN = fileout + "_N.bam"
    fileL = fileout + "_L.bam"
    fileR = fileout + "_R.bam"
    fileS1 = fileout + "_S1.bam"
    fileS2 = fileout + "_S2.bam"
    fileI = fileout + "_I.bam"
    readsM = pysam.AlignmentFile(fileM, "wb", template=bam)
    readsN = pysam.AlignmentFile(fileN, "wb", template=bam)
    readsL = pysam.AlignmentFile(fileL, "wb", template=bam)
    readsR = pysam.AlignmentFile(fileR, "wb", template=bam)
    readsS1 = pysam.AlignmentFile(fileS1, "wb", template=bam)
    readsS2 = pysam.AlignmentFile(fileS2, "wb", template=bam)
    readsI = pysam.AlignmentFile(fileI, "wb", template=bam)
    for read1, read2 in read_pair_generator(bam, region_string):
        counter += 1
        read = read_pair_align(read1, read2)
        if read[0] < target < read[3]:          # fragments that span
            countM += 1
            readsM.write(read1)
            readsM.write(read2)
            sequencedend = False
            if read[0] < target < read[1]:      # fragments where read1 is sequenced
                countS1 += 1
                readsS1.write(read1)
                readsS1.write(read2)
                sequencedend = True
            if read[2] < target < read[3]:      # fragments where read2 is sequenced
                countS2 += 1
                readsS2.write(read1)
                readsS2.write(read2)
                sequencedend = True
            # get indel-containing fragments
            if sequencedend and ('I' in read1.cigarstring or 'I' in read2.cigarstring):
                countI += 1
                readsI.write(read1)
                readsI.write(read2)
        else:                                       # fragments that don't span
            countN += 1
            readsN.write(read1)
            readsN.write(read2)
            # separate reads by bias left vs right
            if read[0] >= target:
                countR += 1
                readsR.write(read1)
                readsR.write(read2)
            elif read[-1] <= target:
                countL += 1
                readsL.write(read1)
                readsL.write(read2)
    readsM.close()
    readsN.close()
    readsL.close()
    readsR.close()
    readsS1.close()
    readsS2.close()
    readsI.close()
    bam.close()
    pysam.sort("-o", fileM, fileM)
    pysam.sort("-o", fileN, fileN)
    pysam.sort("-o", fileL, fileL)
    pysam.sort("-o", fileR, fileR)
    pysam.sort("-o", fileS1, fileS1)
    pysam.sort("-o", fileS2, fileS2)
    pysam.sort("-o", fileI, fileI)
    os.system("samtools index " + fileM)
    os.system("samtools index " + fileN)
    os.system("samtools index " + fileL)
    os.system("samtools index " + fileR)
    os.system("samtools index " + fileS1)
    os.system("samtools index " + fileS2)
    os.system("samtools index " + fileI)
    print("%i span | %i read1seq | %i read2seq | %i indels | %i !span | %i left | %i right | %i total"
          % (countM, countS1, countS2, countI, countN, countL, countR, counter))


def to_wiggle_pairs(filein, fileout, region_string):
    """ Constructs fragment pile-ups in wiggle format using paired-end information

    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.wig) omitted
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160

    """
    [chr, sta, end] = re.split('[:-]', region_string)
    sta = int(sta)
    end = int(end)
    wlist = [0] * (end-sta+1)
    wig = open(fileout + ".wig", "w")
    wig.write("variableStep\tchrom=%s\n" % chr)
    bam = pysam.AlignmentFile(filein, 'rb')
    for read1, read2 in read_pair_generator(bam, region_string):
        read = read_pair_align(read1, read2)
        wlist = [x+1 if read[0]-sta <= i <= read[-1]-sta else x for i, x in enumerate(wlist)]
        # print("%i\t%i" % (read[0]-sta, read[-1]-sta))
    for i, x in enumerate(wlist):
        wig.write("%i\t%i\n" % (sta + i, x))
    wig.close()
    bam.close()


def to_wiggle_windows(filein, fileout, window, chr=None):
    """ Outputs wiggle file that counts number of reads in each window

    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.wig) omitted
    :param window: size of window (i.e. if 500, genome is divided into 500bp windows, function
                outputs the number of reads in each window)
    :param chr: array of chromosome strings to limit analysis to particular chromosomes
                (i.e. ['chr7', 'chr8'])

    """
    bam = pysam.AlignmentFile(filein, 'rb')
    wig = open(fileout + ".wig", "w")
    for row in hg38_generator():  # iterate over each chromosome
        if chr is None or (chr is not None and row[0] in chr):      # checks for chr #
            chr_i = row[0]
            wig.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
            numbins = int(int(row[1]) / window)  # number of bins
            for i in range(numbins):              # iterate over each bin (total - 1)
                start1 = i * window
                finish1 = (i + 1) * window - 1
                wig.write("%i\n" % bam.count(chr_i, start1, finish1))
                status_statement(i, numbins, 50, chr_i)
    wig.close()
    bam.close()


def to_bins(filein, fileout, window, numbins, chr=None):
    """ Outputs CSV file that counts number of reads in each bin (column) for each window (row)

    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.csv) omitted; each row corresponds
                    to each window, each column corresponds to the number of reads for each bin
    :param window: size of window (i.e. if 500, genome is divided into 500bp windows)
    :param numbins: number of bins for each window
    :param chr: array of chromosome strings to limit analysis to particular chromosomes
                (i.e. ['chr7', 'chr8'])

    """
    bam = pysam.AlignmentFile(filein, 'rb')
    for row in hg38_generator():                                    # iterate over each chromosome
        if chr is None or (chr is not None and row[0] in chr):      # checks for chr #
            count = int(int(row[1]) / window)                       # number of windows
            res = int(window / numbins)
            cm = np.zeros((count, 3 + numbins), dtype=object)       # array to hold bin counts info
            chr_i = row[0]
            for i in range(count):                                  # iterate over each window
                win_start = i * window
                win_finish = (i + 1) * window - 1
                for j in range(numbins):                            # iterate over each bin
                    cm[i, 0] = chr_i
                    cm[i, 1] = win_start
                    cm[i, 2] = win_finish
                    bin_start = win_start + j * res
                    bin_finish = win_start + (j + 1) * res - 1
                    cm[i, j + 3] = bam.count(chr_i, bin_start, bin_finish)
                status_statement(i, numbins, 50, chr_i)
            np.savetxt(fileout + ".csv", cm, fmt='%s', delimiter=',')
    bam.close()


def ttest_two(samp_file, ctrl_file, fileout, p=0.01):
    """ From to_bins() outputs, perform significance testing between the bins of each window,
    between experimental and control conditions

    - T-test is passed if bin values of experimental condition windows is significantly higher than
    those of the control condition (one-sided with Bonferroni correction)
    - Outputs CSV file of binning results with t-test result
    - Outputs WIG file indicating windows with significant t-test results

    :param samp_file: output of to_bins() as CSV file for experimental condition
    :param ctrl_file: output of to_bins() as CSV file for control as comparison
    :param fileout: base output file name with extension (.csv) omitted; each row corresponds
                    to each window, each column corresponds to the number of reads for each bin
    :param p: p value cut-off

    """
    np_samp = np.loadtxt(samp_file, delimiter=',', dtype=object)
    np_ctrl = np.loadtxt(ctrl_file, delimiter=',', dtype=object)
    count = np_samp.shape[0]
    if np_samp.shape[0] != np_ctrl.shape[0]:
        raise ValueError("Number of windows between sample and control files are different!")
    ot = np.zeros((count, 9), dtype=object)         # array to hold info for each bin
    wig = open(fileout + ".wig", "w")
    chr_i = None
    prev = 0
    chr_ctrl = np_ctrl[:, 0]                        #
    chr_samp = np_samp[:, 0]
    ran_ctrl = np_ctrl[:, 1:3].astype(int)
    ran_samp = np_samp[:, 1:3].astype(int)
    bin_ctrl = np_ctrl[:, 3:].astype(int)
    bin_samp = np_samp[:, 3:].astype(int)
    for i in range(count):
        if chr_samp[i] != chr_ctrl[i] or ran_samp[i, 0] != ran_ctrl[i, 0] \
                or ran_samp[i, 1] != ran_ctrl[i, 1]:
            raise ValueError("Between sample and control files, the chr, start, or end coordinates "
                             "are different!")
        else:
            if chr_samp[i] != chr_i:
                chr_i = chr_samp[i]
                wig.write("variableStep\tchrom=%s\n" % chr_i)
            ot[i, 0:3] = np_samp[i, 0:3]
            ot[i, 3] = np.sum(bin_ctrl[i, :])
            ot[i, 4] = np.sum(bin_samp[i, :])
            ttest = stats.ttest_rel(bin_samp[i, :], bin_ctrl[i, :])
            ot[i, 5] = ttest[0]  # t test
            ot[i, 6] = ttest[1]  # t test
            start_i = ran_samp[i, 0]
            if ot[i, 5] > 0 and ot[i, 6] / 2 < p / count:
                ot[i, 7] = 1
                ot[i, 8] = start_i - prev
                prev = start_i
                wig.write("%i\t%i\n" % (start_i, 1))
            else:
                ot[i, 7] = np.nan
                ot[i, 8] = np.nan
                wig.write("%i\t%i\n" % (start_i, 0))
            if i % 100 == 0:
                print("Processed %i of %i total windows in %s" % (i, count, chr_i))
    np.savetxt(fileout + "_ttest.csv", ot, fmt='%s', delimiter=',')

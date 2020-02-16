# -*- coding: utf-8 -*-
""" gH2AX ChIP-seq analysis for ACTB-targeting Cas9/cgRNA
"""

import src.chipseq as c

# home directory of BAM files and 'analysis' output directory
base = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/191017_chipseq/"
base_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/191017_chipseq/analysis/"

win = 50000         # window span in base pairs
numbins = 50            # number of bins per window for significance testing

# Covert BAM file to wiggle file that counts the number of reads in each window span
c.to_wiggle_windows(base+"h2ax_actb_nolight_rep1_final.bam", base_a+"00m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_2m_rep1_final.bam", base_a+"02m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_5m_rep1_final.bam", base_a+"05m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_15m_rep1_final.bam", base_a+"15m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_30m_rep1_final.bam", base_a+"30m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_1h_rep1_final.bam", base_a+"60m_rep1", win, ['chr7'])

c.to_wiggle_windows(base+"h2ax_actb_nolight_rep2_final.bam", base_a+"00m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_2m_rep2_final.bam", base_a+"02m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_5m_rep2_final.bam", base_a+"05m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_15m_rep2_final.bam", base_a+"15m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_30m_rep2_final.bam", base_a+"30m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"h2ax_actb_1h_rep2_final.bam", base_a+"60m_rep2", win, ['chr7'])

# For each window span, count number of reads in each bin
c.to_bins(base+"h2ax_actb_nolight_rep1_final.bam", base_a+"00m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_2m_rep1_final.bam", base_a+"02m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_5m_rep1_final.bam", base_a+"05m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_15m_rep1_final.bam", base_a+"15m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_30m_rep1_final.bam", base_a+"30m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_1h_rep1_final.bam", base_a+"60m_rep1", win, numbins, ['chr7'])

c.to_bins(base+"h2ax_actb_nolight_rep2_final.bam", base_a+"00m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_2m_rep2_final.bam", base_a+"02m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_5m_rep2_final.bam", base_a+"05m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_15m_rep2_final.bam", base_a+"15m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_30m_rep2_final.bam", base_a+"30m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"h2ax_actb_1h_rep2_final.bam", base_a+"60m_rep2", win, numbins, ['chr7'])

# Perform T-test on bins by comparing each time point to no-light samples
c.ttest_two(base_a+"00m_rep1.csv", base_a+"00m_rep1.csv", base_a+"ttest-00m_rep1", p=0.05)
c.ttest_two(base_a+"02m_rep1.csv", base_a+"02m_rep1.csv", base_a+"ttest-02m_rep1", p=0.05)
c.ttest_two(base_a+"05m_rep1.csv", base_a+"05m_rep1.csv", base_a+"ttest-05m_rep1", p=0.05)
c.ttest_two(base_a+"15m_rep1.csv", base_a+"15m_rep1.csv", base_a+"ttest-15m_rep1", p=0.05)
c.ttest_two(base_a+"30m_rep1.csv", base_a+"30m_rep1.csv", base_a+"ttest-30m_rep1", p=0.05)
c.ttest_two(base_a+"60m_rep1.csv", base_a+"60m_rep1.csv", base_a+"ttest-60m_rep1", p=0.05)

c.ttest_two(base_a+"00m_rep2.csv", base_a+"00m_rep2.csv", base_a+"ttest-00m_rep2", p=0.05)
c.ttest_two(base_a+"02m_rep2.csv", base_a+"02m_rep2.csv", base_a+"ttest-02m_rep2", p=0.05)
c.ttest_two(base_a+"05m_rep2.csv", base_a+"05m_rep2.csv", base_a+"ttest-05m_rep2", p=0.05)
c.ttest_two(base_a+"15m_rep2.csv", base_a+"15m_rep2.csv", base_a+"ttest-15m_rep2", p=0.05)
c.ttest_two(base_a+"30m_rep2.csv", base_a+"30m_rep2.csv", base_a+"ttest-30m_rep2", p=0.05)
c.ttest_two(base_a+"60m_rep2.csv", base_a+"60m_rep2.csv", base_a+"ttest-60m_rep2", p=0.05)

# Converts BAM to WIG format in 40kb window around cut site to visualize sub kilobase-scale features
tr = 5529660                    # ACTB cleavage site
rr = "chr7:5509660-5549660"     # 40kb window centered at cut site
c.to_wiggle_pairs(base+"h2ax_actb_nolight_rep1_final.bam", base_a+"40kb_00m_rep1", rr)
c.to_wiggle_pairs(base+"h2ax_actb_2m_rep1_final.bam", base_a+"40kb_02m_rep1", rr)
c.to_wiggle_pairs(base+"h2ax_actb_5m_rep1_final.bam", base_a+"40kb_05m_rep1", rr)
c.to_wiggle_pairs(base+"h2ax_actb_15m_rep1_final.bam", base_a+"40kb_15m_rep1", rr)
c.to_wiggle_pairs(base+"h2ax_actb_30m_rep1_final.bam", base_a+"40kb_30m_rep1", rr)
c.to_wiggle_pairs(base+"h2ax_actb_1h_rep1_final.bam", base_a+"40kb_60m_rep1", rr)

c.to_wiggle_pairs(base+"h2ax_actb_nolight_rep1_final.bam", base_a+"40kb_00m_rep2", rr)
c.to_wiggle_pairs(base+"h2ax_actb_2m_rep2_final.bam", base_a+"40kb_02m_rep2", rr)
c.to_wiggle_pairs(base+"h2ax_actb_5m_rep2_final.bam", base_a+"40kb_05m_rep2", rr)
c.to_wiggle_pairs(base+"h2ax_actb_15m_rep2_final.bam", base_a+"40kb_15m_rep2", rr)
c.to_wiggle_pairs(base+"h2ax_actb_30m_rep2_final.bam", base_a+"40kb_30m_rep2", rr)
c.to_wiggle_pairs(base+"h2ax_actb_1h_rep2_final.bam", base_a+"40kb_60m_rep2", rr)

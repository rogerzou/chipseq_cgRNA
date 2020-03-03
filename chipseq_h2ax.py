# -*- coding: utf-8 -*-
""" gamma H2AX ChIP-seq analysis for ACTB-targeting Cas9/cgRNA
"""

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
base = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/cgRNA_SRA/"
base_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/cgRNA_SRA/analysis/"

win = 50000             # window span in base pairs
numbins = 50            # number of bins per window for significance testing


""" Convert BAM file to WIG file that counts the number of reads in each window span. """
c.to_wiggle_windows(base+"gh2ax_actb_00m_rep1.bam", base_a+"gh2ax_00m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_02m_rep1.bam", base_a+"gh2ax_02m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_05m_rep1.bam", base_a+"gh2ax_05m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_15m_rep1.bam", base_a+"gh2ax_15m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_30m_rep1.bam", base_a+"gh2ax_30m_rep1", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_60m_rep1.bam", base_a+"gh2ax_60m_rep1", win, ['chr7'])

c.to_wiggle_windows(base+"gh2ax_actb_00m_rep2.bam", base_a+"gh2ax_00m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_02m_rep2.bam", base_a+"gh2ax_02m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_05m_rep2.bam", base_a+"gh2ax_05m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_15m_rep2.bam", base_a+"gh2ax_15m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_30m_rep2.bam", base_a+"gh2ax_30m_rep2", win, ['chr7'])
c.to_wiggle_windows(base+"gh2ax_actb_60m_rep2.bam", base_a+"gh2ax_60m_rep2", win, ['chr7'])


""" For each window span, count number of reads in each bin. """
c.to_bins(base+"gh2ax_actb_00m_rep1.bam", base_a+"gh2ax_00m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_02m_rep1.bam", base_a+"gh2ax_02m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_05m_rep1.bam", base_a+"gh2ax_05m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_15m_rep1.bam", base_a+"gh2ax_15m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_30m_rep1.bam", base_a+"gh2ax_30m_rep1", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_60m_rep1.bam", base_a+"gh2ax_60m_rep1", win, numbins, ['chr7'])

c.to_bins(base+"gh2ax_actb_00m_rep2.bam", base_a+"gh2ax_00m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_02m_rep2.bam", base_a+"gh2ax_02m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_05m_rep2.bam", base_a+"gh2ax_05m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_15m_rep2.bam", base_a+"gh2ax_15m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_30m_rep2.bam", base_a+"gh2ax_30m_rep2", win, numbins, ['chr7'])
c.to_bins(base+"gh2ax_actb_60m_rep2.bam", base_a+"gh2ax_60m_rep2", win, numbins, ['chr7'])


""" Perform T-test on bins by comparing each time point to no-light samples. """
c.ttest_two(base_a+"gh2ax_00m_rep1.csv", base_a+"gh2ax_00m_rep1.csv", base_a+"ttest-00m_rep1", p=0.05)
c.ttest_two(base_a+"gh2ax_02m_rep1.csv", base_a+"gh2ax_00m_rep1.csv", base_a+"ttest-02m_rep1", p=0.05)
c.ttest_two(base_a+"gh2ax_05m_rep1.csv", base_a+"gh2ax_00m_rep1.csv", base_a+"ttest-05m_rep1", p=0.05)
c.ttest_two(base_a+"gh2ax_15m_rep1.csv", base_a+"gh2ax_00m_rep1.csv", base_a+"ttest-15m_rep1", p=0.05)
c.ttest_two(base_a+"gh2ax_30m_rep1.csv", base_a+"gh2ax_00m_rep1.csv", base_a+"ttest-30m_rep1", p=0.05)
c.ttest_two(base_a+"gh2ax_60m_rep1.csv", base_a+"gh2ax_00m_rep1.csv", base_a+"ttest-60m_rep1", p=0.05)

c.ttest_two(base_a+"gh2ax_00m_rep2.csv", base_a+"gh2ax_00m_rep2.csv", base_a+"ttest-00m_rep2", p=0.05)
c.ttest_two(base_a+"gh2ax_02m_rep2.csv", base_a+"gh2ax_00m_rep2.csv", base_a+"ttest-02m_rep2", p=0.05)
c.ttest_two(base_a+"gh2ax_05m_rep2.csv", base_a+"gh2ax_00m_rep2.csv", base_a+"ttest-05m_rep2", p=0.05)
c.ttest_two(base_a+"gh2ax_15m_rep2.csv", base_a+"gh2ax_00m_rep2.csv", base_a+"ttest-15m_rep2", p=0.05)
c.ttest_two(base_a+"gh2ax_30m_rep2.csv", base_a+"gh2ax_00m_rep2.csv", base_a+"ttest-30m_rep2", p=0.05)
c.ttest_two(base_a+"gh2ax_60m_rep2.csv", base_a+"gh2ax_00m_rep2.csv", base_a+"ttest-60m_rep2", p=0.05)

# Converts BAM to WIG format in 40kb window around cut site to visualize sub kilobase-scale features
tr = 5529660                    # ACTB cleavage site
rr = "chr7:5509660-5549660"     # 40kb window centered at cut site
c.to_wiggle_pairs(base+"gh2ax_actb_00m_rep1_final.bam", base_a+"40kb_00m_rep1", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_02m_rep1_final.bam", base_a+"40kb_02m_rep1", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_05m_rep1_final.bam", base_a+"40kb_05m_rep1", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_15m_rep1_final.bam", base_a+"40kb_15m_rep1", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_30m_rep1_final.bam", base_a+"40kb_30m_rep1", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_60m_rep1_final.bam", base_a+"40kb_60m_rep1", rr)

c.to_wiggle_pairs(base+"gh2ax_actb_00m_rep2_final.bam", base_a+"40kb_00m_rep2", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_02m_rep2_final.bam", base_a+"40kb_02m_rep2", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_05m_rep2_final.bam", base_a+"40kb_05m_rep2", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_15m_rep2_final.bam", base_a+"40kb_15m_rep2", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_30m_rep2_final.bam", base_a+"40kb_30m_rep2", rr)
c.to_wiggle_pairs(base+"gh2ax_actb_60m_rep2_final.bam", base_a+"40kb_60m_rep2", rr)

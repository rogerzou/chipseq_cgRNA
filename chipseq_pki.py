# -*- coding: utf-8 -*-
""" MRE11 ChIP-seq analysis for ACTB-targeting Cas9/cgRNA
"""

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
base = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/cgRNA_SRA/"
base_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/cgRNA_SRA/analysis/"

tr = 5529660                    # ACTB cleavage site
rr = "chr7:5527160-5532160"     # 5kb window centered at cut site


""" Subset BAM files (in 5kb window centered at cut site) by characteristics of each paired-end read
    in relation to the cut site. """
c.get_read_subsets(base+"mre11_actb_noD_rep1.bam", base_a+"mre11_noD_rep1", rr, tr)
c.get_read_subsets(base+"mre11_actb_noD_rep2.bam", base_a+"mre11_noD_rep2", rr, tr)
c.get_read_subsets(base+"mre11_actb_PKi_rep1.bam", base_a+"mre11_PKi_rep1", rr, tr)
c.get_read_subsets(base+"mre11_actb_PKi_rep2.bam", base_a+"mre11_PKi_rep2", rr, tr)


""" Convert all fragments (in 5kb window centered at cut site) to wiggle format. """
c.to_wiggle_pairs(base+"mre11_actb_noD_rep1.bam", base_a+"mre11_noD_rep1", rr)
c.to_wiggle_pairs(base+"mre11_actb_noD_rep2.bam", base_a+"mre11_noD_rep2", rr)
c.to_wiggle_pairs(base+"mre11_actb_PKi_rep1.bam", base_a+"mre11_PKi_rep1", rr)
c.to_wiggle_pairs(base+"mre11_actb_PKi_rep2.bam", base_a+"mre11_PKi_rep2", rr)


""" Convert spanning fragments to wiggle format. """
c.to_wiggle_pairs(base_a+"mre11_noD_rep1_M.bam", base_a+"mre11_noD_rep1_M", rr)
c.to_wiggle_pairs(base_a+"mre11_noD_rep2_M.bam", base_a+"mre11_noD_rep2_M", rr)
c.to_wiggle_pairs(base_a+"mre11_PKi_rep1_M.bam", base_a+"mre11_PKi_rep1_M", rr)
c.to_wiggle_pairs(base_a+"mre11_PKi_rep2_M.bam", base_a+"mre11_PKi_rep2_M", rr)


""" Convert fragments that start/end 5bp away from cut site to wiggle format. """
c.to_wiggle_pairs(base_a+"mre11_noD_rep1_N.bam", base_a+"mre11_noD_rep1_N", rr)
c.to_wiggle_pairs(base_a+"mre11_noD_rep2_N.bam", base_a+"mre11_noD_rep2_N", rr)
c.to_wiggle_pairs(base_a+"mre11_PKi_rep1_N.bam", base_a+"mre11_PKi_rep1_N", rr)
c.to_wiggle_pairs(base_a+"mre11_PKi_rep2_N.bam", base_a+"mre11_PKi_rep2_N", rr)


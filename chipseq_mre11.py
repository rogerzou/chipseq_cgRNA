# -*- coding: utf-8 -*-
""" MRE11 ChIP-seq analysis for ACTB-targeting Cas9/cgRNA
"""

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
base = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/191010_chipseq/"
base_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/191010_chipseq/analysis/"

tr = 5529660                    # ACTB cleavage site
rr = "chr7:5527160-5532160"     # 5kb window centered at cut site


""" Subset BAM files (in 5kb window centered at cut site) by characteristics of each paired-end read
    in relation to the cut site. """
c.get_read_subsets(base+"actb_nolight_rep1_final.bam", base_a+"actb_00m_rep1_final", rr, tr)
c.get_read_subsets(base+"actb_2m_rep1_final.bam", base_a+"actb_02m_rep1_final", rr, tr)
c.get_read_subsets(base+"actb_5m_rep1_final.bam", base_a+"actb_05m_rep1_final", rr, tr)
c.get_read_subsets(base+"actb_15m_rep1_final.bam", base_a+"actb_15m_rep1_final", rr, tr)
c.get_read_subsets(base+"actb_30m_rep1_final.bam", base_a+"actb_30m_rep1_final", rr, tr)
c.get_read_subsets(base+"actb_1h_rep1_final.bam", base_a+"actb_60m_rep1_final", rr, tr)

c.get_read_subsets(base+"actb_nolight_rep2_final.bam", base_a+"actb_00m_rep2_final", rr, tr)
c.get_read_subsets(base+"actb_2m_rep2_final.bam", base_a+"actb_02m_rep2_final", rr, tr)
c.get_read_subsets(base+"actb_5m_rep2_final.bam", base_a+"actb_05m_rep2_final", rr, tr)
c.get_read_subsets(base+"actb_15m_rep2_final.bam", base_a+"actb_15m_rep2_final", rr, tr)
c.get_read_subsets(base+"actb_30m_rep2_final.bam", base_a+"actb_30m_rep2_final", rr, tr)
c.get_read_subsets(base+"actb_1h_rep2_final.bam", base_a+"actb_60m_rep2_final", rr, tr)


""" Convert all fragments (in 5kb window centered at cut site) to wiggle format. """
c.to_wiggle_pairs(base+"actb_nolight_rep1_final.bam", base_a+"actb_00m_rep1_final", rr)
c.to_wiggle_pairs(base+"actb_2m_rep1_final.bam", base_a+"actb_02m_rep1_final", rr)
c.to_wiggle_pairs(base+"actb_5m_rep1_final.bam", base_a+"actb_05m_rep1_final", rr)
c.to_wiggle_pairs(base+"actb_15m_rep1_final.bam", base_a+"actb_15m_rep1_final", rr)
c.to_wiggle_pairs(base+"actb_30m_rep1_final.bam", base_a+"actb_30m_rep1_final", rr)
c.to_wiggle_pairs(base+"actb_1h_rep1_final.bam", base_a+"actb_60m_rep1_final", rr)

c.to_wiggle_pairs(base+"actb_nolight_rep2_final.bam", base_a+"actb_00m_rep2_final", rr)
c.to_wiggle_pairs(base+"actb_2m_rep2_final.bam", base_a+"actb_02m_rep2_final", rr)
c.to_wiggle_pairs(base+"actb_5m_rep2_final.bam", base_a+"actb_05m_rep2_final", rr)
c.to_wiggle_pairs(base+"actb_15m_rep2_final.bam", base_a+"actb_15m_rep2_final", rr)
c.to_wiggle_pairs(base+"actb_30m_rep2_final.bam", base_a+"actb_30m_rep2_final", rr)
c.to_wiggle_pairs(base+"actb_1h_rep2_final.bam", base_a+"actb_60m_rep2_final", rr)


""" Convert spanning fragments to wiggle format. """
c.to_wiggle_pairs(base_a+"actb_00m_rep1_final_M.bam", base_a+"actb_00m_rep1_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_02m_rep1_final_M.bam", base_a+"actb_02m_rep1_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_05m_rep1_final_M.bam", base_a+"actb_05m_rep1_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_15m_rep1_final_M.bam", base_a+"actb_15m_rep1_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_30m_rep1_final_M.bam", base_a+"actb_30m_rep1_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_60m_rep1_final_M.bam", base_a+"actb_60m_rep1_final_M", rr)

c.to_wiggle_pairs(base_a+"actb_00m_rep2_final_M.bam", base_a+"actb_00m_rep2_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_02m_rep2_final_M.bam", base_a+"actb_02m_rep2_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_05m_rep2_final_M.bam", base_a+"actb_05m_rep2_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_15m_rep2_final_M.bam", base_a+"actb_15m_rep2_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_30m_rep2_final_M.bam", base_a+"actb_30m_rep2_final_M", rr)
c.to_wiggle_pairs(base_a+"actb_60m_rep2_final_M.bam", base_a+"actb_60m_rep2_final_M", rr)


""" Convert fragments that start/end 5bp away from cut site to wiggle format. """
c.to_wiggle_pairs(base_a+"actb_00m_rep1_final_N.bam", base_a+"actb_00m_rep1_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_02m_rep1_final_N.bam", base_a+"actb_02m_rep1_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_05m_rep1_final_N.bam", base_a+"actb_05m_rep1_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_15m_rep1_final_N.bam", base_a+"actb_15m_rep1_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_30m_rep1_final_N.bam", base_a+"actb_30m_rep1_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_60m_rep1_final_N.bam", base_a+"actb_60m_rep1_final_N", rr)

c.to_wiggle_pairs(base_a+"actb_00m_rep2_final_N.bam", base_a+"actb_00m_rep2_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_02m_rep2_final_N.bam", base_a+"actb_02m_rep2_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_05m_rep2_final_N.bam", base_a+"actb_05m_rep2_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_15m_rep2_final_N.bam", base_a+"actb_15m_rep2_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_30m_rep2_final_N.bam", base_a+"actb_30m_rep2_final_N", rr)
c.to_wiggle_pairs(base_a+"actb_60m_rep2_final_N.bam", base_a+"actb_60m_rep2_final_N", rr)

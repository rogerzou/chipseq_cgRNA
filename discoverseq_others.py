# -*- coding: utf-8 -*-
""" ChIP-seq analysis for multiple repair factors after Cas9 targeting VEGFA.
    Data was downloaded from SRA (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA509652) and was part
    of the DISCOVER-seq paper (Wienert & Wyman et al., Science, 2019).
"""

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
base = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/data_discoverseq/"
base_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/data_discoverseq/analysis/"

# [cut site coordinate, chromosome and 5kb window flanking cut site]
vegfa = [43770825, "chr6:43768325-43773325"]

ku70_r = "SRR8550677_rmdup.bam"         # Ku70 ChIP-seq, experimental sample
ku70_c = "SRR8550696_rmdup.bam"         # Ku70 ChIP-seq, negative control sample
nbs1_r = "SRR8550679_rmdup.bam"         # NBS1 ChIP-seq, experimental sample
nbs1_c = "SRR8550694_rmdup.bam"         # NBS1 ChIP-seq, negative control sample
rad50_r = "SRR8550682_rmdup.bam"        # RAD50 ChIP-seq, experimental sample
rad50_c = "SRR8550699_rmdup.bam"        # RAD50 ChIP-seq, negative control sample
gh2ax_r = "SRR8550678_rmdup.bam"        # gamma H2AX ChIP-seq, experimental sample
gh2ax_c = "SRR8550697_rmdup.bam"        # gamma H2AX ChIP-seq, negative control sample
fancd2_r = "SRR8550698_rmdup.bam"       # FANCD2 ChIP-seq, experimental sample
fancd2_c = "SRR8550690_rmdup.bam"       # FANCD2 ChIP-seq, negative control sample


""" Subset BAM files (in 5kb window centered at cut site) by characteristics of each paired-end read
    in relation to the cut site. """
c.get_read_subsets(base + ku70_r, base_a + "VEGFA_ku70_r_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + ku70_c, base_a + "VEGFA_ku70_c_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + nbs1_r, base_a + "VEGFA_nbs1_r_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + nbs1_c, base_a + "VEGFA_nbs1_c_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + rad50_r, base_a + "VEGFA_rad50_r_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + rad50_c, base_a + "VEGFA_rad50_c_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + gh2ax_r, base_a + "VEGFA_gh2ax_r_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + gh2ax_c, base_a + "VEGFA_gh2ax_c_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + fancd2_r, base_a + "VEGFA_fancd2_r_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + fancd2_c, base_a + "VEGFA_fancd2_c_sub", vegfa[1], vegfa[0])


""" Convert all fragments (in 5kb window centered at cut site) to wiggle format. """
c.to_wiggle_pairs(base + ku70_r, base_a + "VEGFA_ku70_r_pile", vegfa[1])
c.to_wiggle_pairs(base + ku70_c, base_a + "VEGFA_ku70_c_pile", vegfa[1])
c.to_wiggle_pairs(base + nbs1_r, base_a + "VEGFA_nbs1_r_pile", vegfa[1])
c.to_wiggle_pairs(base + nbs1_c, base_a + "VEGFA_nbs1_c_pile", vegfa[1])
c.to_wiggle_pairs(base + rad50_r, base_a + "VEGFA_rad50_r_pile", vegfa[1])
c.to_wiggle_pairs(base + rad50_c, base_a + "VEGFA_rad50_c_pile", vegfa[1])
c.to_wiggle_pairs(base + gh2ax_r, base_a + "VEGFA_gh2ax_r_pile", vegfa[1])
c.to_wiggle_pairs(base + gh2ax_c, base_a + "VEGFA_gh2ax_c_pile", vegfa[1])
c.to_wiggle_pairs(base + fancd2_r, base_a + "VEGFA_fancd2_r_pile", vegfa[1])
c.to_wiggle_pairs(base + fancd2_c, base_a + "VEGFA_fancd2_c_pile", vegfa[1])


""" Convert spanning fragments to wiggle format. """
c.to_wiggle_pairs(base_a + "VEGFA_ku70_r_sub_M.bam", base_a + "VEGFA_ku70_r_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_ku70_c_sub_M.bam", base_a + "VEGFA_ku70_c_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_nbs1_r_sub_M.bam", base_a + "VEGFA_nbs1_r_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_nbs1_c_sub_M.bam", base_a + "VEGFA_nbs1_c_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_rad50_r_sub_M.bam", base_a + "VEGFA_rad50_r_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_rad50_c_sub_M.bam", base_a + "VEGFA_rad50_c_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_gh2ax_r_sub_M.bam", base_a + "VEGFA_gh2ax_r_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_gh2ax_c_sub_M.bam", base_a + "VEGFA_gh2ax_c_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_fancd2_r_sub_M.bam", base_a + "VEGFA_fancd2_r_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_fancd2_c_sub_M.bam", base_a + "VEGFA_fancd2_c_pile_M", vegfa[1])


""" Convert fragments that start/end 5bp away from cut site to wiggle format. """
c.to_wiggle_pairs(base_a + "VEGFA_ku70_r_sub_N.bam", base_a + "VEGFA_ku70_r_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_ku70_c_sub_N.bam", base_a + "VEGFA_ku70_c_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_nbs1_r_sub_N.bam", base_a + "VEGFA_nbs1_r_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_nbs1_c_sub_N.bam", base_a + "VEGFA_nbs1_c_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_rad50_r_sub_N.bam", base_a + "VEGFA_rad50_r_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_rad50_c_sub_N.bam", base_a + "VEGFA_rad50_c_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_gh2ax_r_sub_N.bam", base_a + "VEGFA_gh2ax_r_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_gh2ax_c_sub_N.bam", base_a + "VEGFA_gh2ax_c_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_fancd2_r_sub_N.bam", base_a + "VEGFA_fancd2_r_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_fancd2_c_sub_N.bam", base_a + "VEGFA_fancd2_c_pile_N", vegfa[1])

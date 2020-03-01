# -*- coding: utf-8 -*-
""" MRE11 ChIP-seq analysis for ACTB-targeting Cas9/cgRNA
"""

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
base = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/data_discoverseq/"
base_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/data_discoverseq/analysis/"

hbb = [5226984, "chr11:5224484-5229484"]
vegfa = [43770825, "chr6:43768325-43773325"]
emx1 = [72933869, "chr2:72931369-72936369"]
rnf2 = [185087641, "chr1:185085141-185090141"]
pcsk9 = [106463862, "chr4:106461362-106466362"]
hbb_1 = "SRR8550692_rmdup.bam"
hbb_2 = "SRR8550673_rmdup.bam"
vegfa_1 = "SRR8550703_rmdup.bam"
vegfa_2 = "SRR8550680_rmdup.bam"
emx1_1 = "SRR8550681_rmdup.bam"
emx1_2 = "SRR8550704_rmdup.bam"
rnf2_1 = "SRR8550684_rmdup.bam"
rnf2_2 = "SRR8550705_rmdup.bam"
neghg_1 = "SRR8550693_rmdup.bam"
neghg_2 = "SRR8550695_rmdup.bam"
pcsk9_1 = "SRR8553800_rmdup.bam"
pcsk9_2 = "SRR8553810_rmdup.bam"
negmm_1 = "SRR8553804_rmdup.bam"
negmm_2 = "SRR8553806_rmdup.bam"


""" Subset BAM files (in 5kb window centered at cut site) by characteristics of each paired-end read
    in relation to the cut site. """
c.get_read_subsets(base + hbb_1, base_a + "HBB_MRE11_r1_sub", hbb[1], hbb[0])
c.get_read_subsets(base + hbb_2, base_a + "HBB_MRE11_r2_sub", hbb[1], hbb[0])
c.get_read_subsets(base + vegfa_1, base_a + "VEGFA_MRE11_r1_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + vegfa_2, base_a + "VEGFA_MRE11_r2_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + emx1_1, base_a + "EMX1_MRE11_r1_sub", emx1[1], emx1[0])
c.get_read_subsets(base + emx1_2, base_a + "EMX1_MRE11_r2_sub", emx1[1], emx1[0])
c.get_read_subsets(base + rnf2_1, base_a + "RNF2_MRE11_r1_sub", rnf2[1], rnf2[0])
c.get_read_subsets(base + rnf2_2, base_a + "RNF2_MRE11_r2_sub", rnf2[1], rnf2[0])
c.get_read_subsets(base + pcsk9_1, base_a + "PCSK9_MRE11_r1_sub", pcsk9[1], pcsk9[0])
c.get_read_subsets(base + pcsk9_2, base_a + "PCSK9_MRE11_r2_sub", pcsk9[1], pcsk9[0])

c.get_read_subsets(base + neghg_1, base_a + "HBB_MRE11_c1_sub", hbb[1], hbb[0])
c.get_read_subsets(base + neghg_2, base_a + "HBB_MRE11_c2_sub", hbb[1], hbb[0])
c.get_read_subsets(base + neghg_1, base_a + "VEGFA_MRE11_c1_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + neghg_2, base_a + "VEGFA_MRE11_c2_sub", vegfa[1], vegfa[0])
c.get_read_subsets(base + neghg_1, base_a + "EMX1_MRE11_c1_sub", emx1[1], emx1[0])
c.get_read_subsets(base + neghg_2, base_a + "EMX1_MRE11_c2_sub", emx1[1], emx1[0])
c.get_read_subsets(base + neghg_1, base_a + "RNF2_MRE11_c1_sub", rnf2[1], rnf2[0])
c.get_read_subsets(base + neghg_2, base_a + "RNF2_MRE11_c2_sub", rnf2[1], rnf2[0])
c.get_read_subsets(base + negmm_1, base_a + "PCSK9_MRE11_c1_sub", pcsk9[1], pcsk9[0])
c.get_read_subsets(base + negmm_2, base_a + "PCSK9_MRE11_c2_sub", pcsk9[1], pcsk9[0])


""" Convert all fragments (in 5kb window centered at cut site) to wiggle format. """
c.to_wiggle_pairs(base + hbb_1, base_a + "HBB_MRE11_r1_pile", hbb[1], True)
c.to_wiggle_pairs(base + hbb_2, base_a + "HBB_MRE11_r2_pile", hbb[1], True)
c.to_wiggle_pairs(base + vegfa_1, base_a + "VEGFA_MRE11_r1_pile", vegfa[1], True)
c.to_wiggle_pairs(base + vegfa_2, base_a + "VEGFA_MRE11_r2_pile", vegfa[1], True)
c.to_wiggle_pairs(base + emx1_1, base_a + "EMX1_MRE11_r1_pile", emx1[1], True)
c.to_wiggle_pairs(base + emx1_2, base_a + "EMX1_MRE11_r2_pile", emx1[1], True)
c.to_wiggle_pairs(base + rnf2_1, base_a + "RNF2_MRE11_r1_pile", rnf2[1], True)
c.to_wiggle_pairs(base + rnf2_2, base_a + "RNF2_MRE11_r2_pile", rnf2[1], True)
c.to_wiggle_pairs(base + pcsk9_1, base_a + "PCSK9_MRE11_r1_pile", pcsk9[1], True)
c.to_wiggle_pairs(base + pcsk9_2, base_a + "PCSK9_MRE11_r2_pile", pcsk9[1], True)

c.to_wiggle_pairs(base + neghg_1, base_a + "HBB_MRE11_c1_pile", hbb[1], True)
c.to_wiggle_pairs(base + neghg_2, base_a + "HBB_MRE11_c2_pile", hbb[1], True)
c.to_wiggle_pairs(base + neghg_1, base_a + "VEGFA_MRE11_c1_pile", vegfa[1], True)
c.to_wiggle_pairs(base + neghg_2, base_a + "VEGFA_MRE11_c2_pile", vegfa[1], True)
c.to_wiggle_pairs(base + neghg_1, base_a + "EMX1_MRE11_c1_pile", emx1[1], True)
c.to_wiggle_pairs(base + neghg_2, base_a + "EMX1_MRE11_c2_pile", emx1[1], True)
c.to_wiggle_pairs(base + neghg_1, base_a + "RNF2_MRE11_c1_pile", rnf2[1], True)
c.to_wiggle_pairs(base + neghg_2, base_a + "RNF2_MRE11_c2_pile", rnf2[1], True)
c.to_wiggle_pairs(base + negmm_1, base_a + "PCSK9_MRE11_c1_pile", pcsk9[1], True)
c.to_wiggle_pairs(base + negmm_2, base_a + "PCSK9_MRE11_c2_pile", pcsk9[1], True)


""" Convert spanning fragments to wiggle format. """
c.to_wiggle_pairs(base_a + "HBB_MRE11_r1_sub_M.bam", base_a + "HBB_MRE11_r1_pile_M", hbb[1])
c.to_wiggle_pairs(base_a + "HBB_MRE11_r2_sub_M.bam", base_a + "HBB_MRE11_r2_pile_M", hbb[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_r1_sub_M.bam", base_a + "VEGFA_MRE11_r1_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_r2_sub_M.bam", base_a + "VEGFA_MRE11_r2_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_r1_sub_M.bam", base_a + "EMX1_MRE11_r1_pile_M", emx1[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_r2_sub_M.bam", base_a + "EMX1_MRE11_r2_pile_M", emx1[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_r1_sub_M.bam", base_a + "RNF2_MRE11_r1_pile_M", rnf2[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_r2_sub_M.bam", base_a + "RNF2_MRE11_r2_pile_M", rnf2[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_r1_sub_M.bam", base_a + "PCSK9_MRE11_r1_pile_M", pcsk9[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_r2_sub_M.bam", base_a + "PCSK9_MRE11_r2_pile_M", pcsk9[1])

c.to_wiggle_pairs(base_a + "HBB_MRE11_c1_sub_M.bam", base_a + "HBB_MRE11_c1_pile_M", hbb[1])
c.to_wiggle_pairs(base_a + "HBB_MRE11_c2_sub_M.bam", base_a + "HBB_MRE11_c2_pile_M", hbb[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_c1_sub_M.bam", base_a + "VEGFA_MRE11_c1_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_c2_sub_M.bam", base_a + "VEGFA_MRE11_c2_pile_M", vegfa[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_c1_sub_M.bam", base_a + "EMX1_MRE11_c1_pile_M", emx1[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_c2_sub_M.bam", base_a + "EMX1_MRE11_c2_pile_M", emx1[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_c1_sub_M.bam", base_a + "RNF2_MRE11_c1_pile_M", rnf2[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_c2_sub_M.bam", base_a + "RNF2_MRE11_c2_pile_M", rnf2[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_c1_sub_M.bam", base_a + "PCSK9_MRE11_c1_pile_M", pcsk9[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_c2_sub_M.bam", base_a + "PCSK9_MRE11_c2_pile_M", pcsk9[1])


""" Convert fragments that start/end 5bp away from cut site to wiggle format. """
c.to_wiggle_pairs(base_a + "HBB_MRE11_r1_sub_N.bam", base_a + "HBB_MRE11_r1_pile_N", hbb[1])
c.to_wiggle_pairs(base_a + "HBB_MRE11_r2_sub_N.bam", base_a + "HBB_MRE11_r2_pile_N", hbb[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_r1_sub_N.bam", base_a + "VEGFA_MRE11_r1_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_r2_sub_N.bam", base_a + "VEGFA_MRE11_r2_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_r1_sub_N.bam", base_a + "EMX1_MRE11_r1_pile_N", emx1[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_r2_sub_N.bam", base_a + "EMX1_MRE11_r2_pile_N", emx1[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_r1_sub_N.bam", base_a + "RNF2_MRE11_r1_pile_N", rnf2[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_r2_sub_N.bam", base_a + "RNF2_MRE11_r2_pile_N", rnf2[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_r1_sub_N.bam", base_a + "PCSK9_MRE11_r1_pile_N", pcsk9[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_r2_sub_N.bam", base_a + "PCSK9_MRE11_r2_pile_N", pcsk9[1])

c.to_wiggle_pairs(base_a + "HBB_MRE11_c1_sub_N.bam", base_a + "HBB_MRE11_c1_pile_N", hbb[1])
c.to_wiggle_pairs(base_a + "HBB_MRE11_c2_sub_N.bam", base_a + "HBB_MRE11_c2_pile_N", hbb[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_c1_sub_N.bam", base_a + "VEGFA_MRE11_c1_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "VEGFA_MRE11_c2_sub_N.bam", base_a + "VEGFA_MRE11_c2_pile_N", vegfa[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_c1_sub_N.bam", base_a + "EMX1_MRE11_c1_pile_N", emx1[1])
c.to_wiggle_pairs(base_a + "EMX1_MRE11_c2_sub_N.bam", base_a + "EMX1_MRE11_c2_pile_N", emx1[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_c1_sub_N.bam", base_a + "RNF2_MRE11_c1_pile_N", rnf2[1])
c.to_wiggle_pairs(base_a + "RNF2_MRE11_c2_sub_N.bam", base_a + "RNF2_MRE11_c2_pile_N", rnf2[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_c1_sub_N.bam", base_a + "PCSK9_MRE11_c1_pile_N", pcsk9[1])
c.to_wiggle_pairs(base_a + "PCSK9_MRE11_c2_sub_N.bam", base_a + "PCSK9_MRE11_c2_pile_N", pcsk9[1])
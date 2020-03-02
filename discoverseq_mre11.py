# -*- coding: utf-8 -*-
""" MRE11 ChIP-seq analysis for Cas9 targeting multiple targets. Data was downloaded from SRA
    (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA509652) and was part of the DISCOVER-seq paper
    (Wienert & Wyman et al., Science, 2019).
"""

import src.chipseq as c

""" Home directory of BAM files and 'analysis' output directory; MODIFY AS APPROPRIATE. """
base = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/data_discoverseq/"
base_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/data_discoverseq/analysis/"

# [cut site coordinate, chromosome and 5kb window flanking cut site]
hbb = [5226984, "chr11:5224484-5229484"]            # coordinates for HBB cut site in hg38
vegfa = [43770825, "chr6:43768325-43773325"]        # coordinates for VEGFA cut site in hg38
emx1 = [72933869, "chr2:72931369-72936369"]         # coordinates for EMX1 cut site in hg38
rnf2 = [185087641, "chr1:185085141-185090141"]      # coordinates for RNF2 cut site in hg38
pcsk9 = [106463862, "chr4:106461362-106466362"]     # coordinates for PCSK9 cut site in mm10

hbb_1 = "SRR8550692_rmdup.bam"                      # Cas9 targeting HBB, replicate 1
hbb_2 = "SRR8550673_rmdup.bam"                      # Cas9 targeting HBB, replicate 2
vegfa_1 = "SRR8550703_rmdup.bam"                    # Cas9 targeting VEGFA, replicate 1
vegfa_2 = "SRR8550680_rmdup.bam"                    # Cas9 targeting VEGFA, replicate 2
emx1_1 = "SRR8550681_rmdup.bam"                     # Cas9 targeting EMX1, replicate 1
emx1_2 = "SRR8550704_rmdup.bam"                     # Cas9 targeting EMX1, replicate 2
rnf2_1 = "SRR8550684_rmdup.bam"                     # Cas9 targeting RNF2, replicate 1
rnf2_2 = "SRR8550705_rmdup.bam"                     # Cas9 targeting RNF2, replicate 2
neghg_1 = "SRR8550693_rmdup.bam"                    # non-targeting ctrl in K562, replicate 1
neghg_2 = "SRR8550695_rmdup.bam"                    # non-targeting ctrl in K562, replicate 2
pcsk9_1 = "SRR8553800_rmdup.bam"                    # Cas9 targeting PCSK9, replicate 1
pcsk9_2 = "SRR8553810_rmdup.bam"                    # Cas9 targeting PCSK9, replicate 2
negmm_1 = "SRR8553804_rmdup.bam"                    # non-targeting ctrl in mouse liver, replicate 1
negmm_2 = "SRR8553806_rmdup.bam"                    # non-targeting ctrl in mouse liver, replicate 2


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
c.to_wiggle_pairs(base + hbb_1, base_a + "HBB_MRE11_r1_pile", hbb[1])
c.to_wiggle_pairs(base + hbb_2, base_a + "HBB_MRE11_r2_pile", hbb[1])
c.to_wiggle_pairs(base + vegfa_1, base_a + "VEGFA_MRE11_r1_pile", vegfa[1])
c.to_wiggle_pairs(base + vegfa_2, base_a + "VEGFA_MRE11_r2_pile", vegfa[1])
c.to_wiggle_pairs(base + emx1_1, base_a + "EMX1_MRE11_r1_pile", emx1[1])
c.to_wiggle_pairs(base + emx1_2, base_a + "EMX1_MRE11_r2_pile", emx1[1])
c.to_wiggle_pairs(base + rnf2_1, base_a + "RNF2_MRE11_r1_pile", rnf2[1])
c.to_wiggle_pairs(base + rnf2_2, base_a + "RNF2_MRE11_r2_pile", rnf2[1])
c.to_wiggle_pairs(base + pcsk9_1, base_a + "PCSK9_MRE11_r1_pile", pcsk9[1])
c.to_wiggle_pairs(base + pcsk9_2, base_a + "PCSK9_MRE11_r2_pile", pcsk9[1])

c.to_wiggle_pairs(base + neghg_1, base_a + "HBB_MRE11_c1_pile", hbb[1])
c.to_wiggle_pairs(base + neghg_2, base_a + "HBB_MRE11_c2_pile", hbb[1])
c.to_wiggle_pairs(base + neghg_1, base_a + "VEGFA_MRE11_c1_pile", vegfa[1])
c.to_wiggle_pairs(base + neghg_2, base_a + "VEGFA_MRE11_c2_pile", vegfa[1])
c.to_wiggle_pairs(base + neghg_1, base_a + "EMX1_MRE11_c1_pile", emx1[1])
c.to_wiggle_pairs(base + neghg_2, base_a + "EMX1_MRE11_c2_pile", emx1[1])
c.to_wiggle_pairs(base + neghg_1, base_a + "RNF2_MRE11_c1_pile", rnf2[1])
c.to_wiggle_pairs(base + neghg_2, base_a + "RNF2_MRE11_c2_pile", rnf2[1])
c.to_wiggle_pairs(base + negmm_1, base_a + "PCSK9_MRE11_c1_pile", pcsk9[1])
c.to_wiggle_pairs(base + negmm_2, base_a + "PCSK9_MRE11_c2_pile", pcsk9[1])


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

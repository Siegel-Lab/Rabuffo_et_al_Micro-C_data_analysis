### This script was written by Markus Schmidt (LMU)
### It takes a .pairs file generated after 3C-like methods as input and corrects it for ploidy.
### A matrix file has at least the following 5 columns:  chr1    pos1    chr2    pos2   value
### The script duplicates each line (i.e. each interaction pair) and split the interactions coming from the cores in
### interactions coming from either coreA or coreB.
### It makes 2 assumptions:
### - core interactions occur only between A-A or B-B
### - core A interacts always with A arms (3 or 5), core B interacts always with B arms
### Both the assumptions remain true also with interchromosomal interactions.
### The best way to display heatmaps from the ploidy-corrected matrix file is to separate the 2 haplotipes.

import sys


DIPLOID_CONTIGS = {
    "Chr1_core_Tb427v12",
    "Chr2_core_Tb427v12",
    "Chr3_core_Tb427v12",
    "Chr3_3A_Tb427v12",
    "Chr4_core_Tb427v12",
    "Chr5_core_Tb427v12",
    "Chr6_core_Tb427v12",
    "Chr7_core_Tb427v12",
    "Chr7_5A_Tb427v12",
    "Chr8_core_Tb427v12",
    "Chr9_core_Tb427v12",
    "Chr10_core_Tb427v12",
    "Chr11_core_Tb427v12",
    "BES15_Tb427v12",
}

DIPLOID_A_CONTIG_NAMES = {
    "Chr1_core_Tb427v12": "Chr1_coreA_Tb427v12",
    "Chr2_core_Tb427v12": "Chr2_coreA_Tb427v12",
    "Chr3_core_Tb427v12": "Chr3_coreA_Tb427v12",
    "Chr3_3A_Tb427v12": "Chr3_3A.I_Tb427v12",
    "Chr4_core_Tb427v12": "Chr4_coreA_Tb427v12",
    "Chr5_core_Tb427v12": "Chr5_coreA_Tb427v12",
    "Chr6_core_Tb427v12": "Chr6_coreA_Tb427v12",
    "Chr7_core_Tb427v12": "Chr7_coreA_Tb427v12",
    "Chr7_5A_Tb427v12": "Chr7_5A.I_Tb427v12",
    "Chr8_core_Tb427v12": "Chr8_coreA_Tb427v12",
    "Chr9_core_Tb427v12": "Chr9_coreA_Tb427v12",
    "Chr10_core_Tb427v12": "Chr10_coreA_Tb427v12",
    "Chr11_core_Tb427v12": "Chr11_coreA_Tb427v12",
    "BES15_Tb427v12":     "BES15.I_Tb427v12",
}

DIPLOID_B_CONTIG_NAMES = {
    "Chr1_core_Tb427v12": "Chr1_coreB_Tb427v12",
    "Chr2_core_Tb427v12": "Chr2_coreB_Tb427v12",
    "Chr3_core_Tb427v12": "Chr3_coreB_Tb427v12",
    "Chr3_3A_Tb427v12": "Chr3_3A.II_Tb427v12",
    "Chr4_core_Tb427v12": "Chr4_coreB_Tb427v12",
    "Chr5_core_Tb427v12": "Chr5_coreB_Tb427v12",
    "Chr6_core_Tb427v12": "Chr6_coreB_Tb427v12",
    "Chr7_core_Tb427v12": "Chr7_coreB_Tb427v12",
    "Chr7_5A_Tb427v12": "Chr7_5A.II_Tb427v12",
    "Chr8_core_Tb427v12": "Chr8_coreB_Tb427v12",
    "Chr9_core_Tb427v12": "Chr9_coreB_Tb427v12",
    "Chr10_core_Tb427v12": "Chr10_coreB_Tb427v12",
    "Chr11_core_Tb427v12": "Chr11_coreB_Tb427v12",
    "BES15_Tb427v12":     "BES15.II_Tb427v12",

}

HAPLOID_A_CONTIGS = {
    "Chr1_5B_Tb427v12",
    "Chr1_3A_Tb427v12",
    "Chr2_5A_Tb427v12",
    "BES12_Tb427v12",
    "Chr3_5A_Tb427v12",
    "Chr4_5A_Tb427v12",
    "Chr4_3A_Tb427v12",
    "BES5_Tb427v12",
    "Chr5_3B_Tb427v12",
    "BES7_Tb427v12",
    "Chr6_3A_Tb427v12",
    "Chr8_5B_Tb427v12",
    "Chr8_3B_Tb427v12",
    "Chr9_5A_Tb427v12",
    "Chr9_3A_Tb427v12",
    "Chr10_5B_Tb427v12",
    "Chr10_3B_Tb427v12",
    "Chr11_5B_Tb427v12",
    "Chr11_3B_Tb427v12",
    "BES2_Tb427v12",
    "BES4_Tb427v12",
    "BES10_Tb427v12",
    "BES11_Tb427v12",
    "BES13_Tb427v12",
    "BES14_Tb427v12",
    "BES17_Tb427v12"
}

HAPLOID_B_CONTIGS = { # unused - just for completeness sake
    "Chr1_5A_Tb427v12",
    "Chr1_3B_Tb427v12",
    "Chr2_coreB_Tb427v12",
    "Chr3_5B_Tb427v12",
    "Chr3_coreB_Tb427v12",
    "BES3_Tb427v12",
    "Chr4_5B_Tb427v12",
    "Chr4_coreB_Tb427v12",
    "Chr4_3B_Tb427v12",
    "Chr5_coreB_Tb427v12",
    "Chr5_3A_Tb427v12",
    "BES1_Tb427v12",
    "Chr6_3B_Tb427v12",
    "Chr8_5A_Tb427v12",
    "Chr8_3A_Tb427v12",
    "Chr9_5B_Tb427v12",
    "Chr9_3B_Tb427v12",
    "Chr10_5A_Tb427v12",
    "Chr10_3A_Tb427v12",
    "Chr11_5A_Tb427v12",
    "Chr11_3A_Tb427v12"
}

HAPLOID_UNKNOWN_CONTIGS = {
}

with open(sys.argv[1]) as in_file:
    # iterate through all lines
    for line in in_file.readlines():
        # split line into columns -> we only need the first 5, 
        # the rest we store in extra and print them the same way we got them
        chr_1, pos_1, chr_2, pos_2, value = line[:-1].split()
        # both contigs diploid: -> split into A and B
        if chr_1 in DIPLOID_CONTIGS and chr_2 in DIPLOID_CONTIGS:
            # print line for A
            print(DIPLOID_A_CONTIG_NAMES[chr_1], pos_1, 
                  DIPLOID_A_CONTIG_NAMES[chr_2], pos_2, value, sep="\t")
            # print line for B
            print(DIPLOID_B_CONTIG_NAMES[chr_1], pos_1, 
                  DIPLOID_B_CONTIG_NAMES[chr_2], pos_2, value, sep="\t")
        
        # neither diploid: -> print twice
        elif chr_1 not in DIPLOID_CONTIGS and chr_2 not in DIPLOID_CONTIGS:
            print(chr_1, pos_1, chr_2, pos_2, (int(value)*2), sep="\t")
        
        # only one diploid: -> figure out proper contig to put in & print twice
        else:
            fst_hapl = chr_1 not in DIPLOID_CONTIGS
            # get diploid contig name
            chr_dipl = chr_2 if fst_hapl else chr_1
            
            # get haploid contig name
            chr_hapl = chr_1 if fst_hapl else chr_2
            # check if we should turn diploid contig name into A or B variant
            if chr_hapl in HAPLOID_UNKNOWN_CONTIGS:
                # adapt diploid contig name to A
                chr_dipl_a = DIPLOID_A_CONTIG_NAMES[chr_dipl]
                # save back diploid contig name to the correct "position" in the line
                if fst_hapl:
                    chr_2 = chr_dipl_a
                else:
                    chr_1 = chr_dipl_a
                # print line
                print(chr_1, pos_1, chr_2, pos_2, value, sep="\t")
                chr_dipl_b = DIPLOID_B_CONTIG_NAMES[chr_dipl]
                # save back diploid contig name to the correct "position" in the line
                if fst_hapl:
                    chr_2 = chr_dipl_b
                else:
                    chr_1 = chr_dipl_b
                # print line
                print(chr_1, pos_1, chr_2, pos_2, value, sep="\t")
            else:
                if chr_hapl in HAPLOID_A_CONTIGS:
                    # adapt diploid contig name to A
                    chr_dipl = DIPLOID_A_CONTIG_NAMES[chr_dipl]
                else:
                    # adapt diploid contig name to B
                    chr_dipl = DIPLOID_B_CONTIG_NAMES[chr_dipl]
                
                # save back diploid contig name to the correct "position" in the line
                if fst_hapl:
                    chr_2 = chr_dipl
                else:
                    chr_1 = chr_dipl
                
                # print line twice
                print(chr_1, pos_1, chr_2, pos_2, (int(value)*2), sep="\t")

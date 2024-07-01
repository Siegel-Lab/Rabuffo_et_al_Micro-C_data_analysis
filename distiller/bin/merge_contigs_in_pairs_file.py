# created by Markus Schmidt and Claudia Rabuffo @ LMU 2023
# 
# merges contigs in a pairs file (i.e. combines contigs by placing their interactions consecutively)
# as input give the pairs file and a text file containing the sizes of the contigs
# the text file (.ctsizes) should be in the following format:
# ### new_contig_name
# old_contig_name_a contig_size_a
# old_contig_name_b contig_size_b


import fileinput
import sys

def merge_contigs(ctsizes_file, pairs_file):
    translate_dict = {}
    with open(ctsizes_file, "r") as ctsizes_file:
        curr_new_contig_name = None
        curr_size_add = 0
        for line in ctsizes_file:
            if line[:4] == "### ":
                curr_new_contig_name = line[4:-1]
                curr_size_add = 0
                continue
            else:
                contig_name, pos = line[:-1].split()
                translate_dict[contig_name] = (curr_new_contig_name, curr_size_add)
                curr_size_add += int(pos)
    for line in fileinput.input(pairs_file):
        if len(line) == 0 or line[0] == "#":
            print(line[:-1])
            continue
        read_name, contig_1, pos_1, contig_2, pos_2, *extra = line[:-1].split("\t")

        pos_1 = int(pos_1)
        pos_2 = int(pos_2)

        if contig_1 in translate_dict:
            new_c, add_pos = translate_dict[contig_1]
            contig_1 = new_c
            pos_1 += add_pos
        if contig_2 in translate_dict:
            new_c, add_pos = translate_dict[contig_2]
            contig_2 = new_c
            pos_2 += add_pos

        print(read_name, contig_1, pos_1, contig_2, pos_2, *extra, sep="\t")

if __name__ == "__main__":
    merge_contigs(sys.argv[1], sys.argv[2])


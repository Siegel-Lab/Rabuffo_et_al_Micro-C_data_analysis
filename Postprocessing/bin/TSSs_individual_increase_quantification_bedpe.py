import cooler
import sys
import fileinput
import math
from scipy.stats import ks_2samp
import random


# requirements: 
#   cooler
#
# usage: 
#   python3 individual_increase_quantification.py <cooler_path> <bedpe_path> <ctrl_offset> > <outfile>
#
# cooler_path is a .cool file (single resolution, must be balanced 'cooler balance')
# gff_path is a .gff file
# ctrl_offset is the offset to the control region. default: 20000
#
#
# This script quantifies the increase of interactions around a set of annotations
# for this it averages the balanced number of interactions at every intersection of the annotations
# additionally for each intersection it calculates a control by shifting the region downstream by ctrl_offset

def anno_from_line(line):
    chr_x, start_x, end_x, chr_y, start_y, end_y = line.strip().split()
    return chr_x, int(start_x), int(end_x), chr_y, int(start_y), int(end_y)

def load_annotations(bedpe_path, bedpe_ctrl):
    annotations = []
    for line, ctrl_line in zip(fileinput.input(bedpe_path), fileinput.input(bedpe_ctrl)):
        annotations.append((anno_from_line(line), anno_from_line(ctrl_line)))
    return annotations

def get_idx(bin_chr, bin_pos, clr, bin_size):
    return clr.offset(bin_chr) + bin_pos // bin_size

def get_chr_size(clr, chr_name):
    return clr.chromsizes[chr_name]

def extract_value_cool(clr, bin_size, chr_x, start_x, end_x, chr_y, start_y, end_y):
    start_idx_x = get_idx(chr_x, start_x, clr, bin_size)
    end_idx_x = get_idx(chr_x, end_x, clr, bin_size)
    start_idx_y = get_idx(chr_y, start_y, clr, bin_size)
    end_idx_y = get_idx(chr_y, end_y, clr, bin_size)
    val = []
    vals = clr.matrix(balance=True)[start_idx_x : end_idx_x + 1, start_idx_y : end_idx_y + 1]
    for v in vals.flatten():
        if not math.isnan(v):
            val.append(v)
    return val

def mean(xs):
    return sum(xs) / len(xs)

def individual_increase_quantification(cooler_path, bedpe_path, bedpe_ctrl, low_interactions_out):

    clr = cooler.Cooler(cooler_path)
    bin_size = clr.binsize

    annotations = load_annotations(bedpe_path, bedpe_ctrl)
    # print(annotations, file=sys.stderr)

    print("#individual increase quantification")
    print("#parameter", "cooler_path:", cooler_path)
    print("#parameter", "bedpe_path:", bedpe_path)
    print("#columns:", "chr_x", "start_x", "end_x", "chr_y", "start_y", "end_y", "val", "val_ctrl", sep="\t")
    count_interactions_intra_chr = 0
    count_interactions_inter_chr = 0
    count_enriched_intra_chr = 0
    count_enriched_inter_chr = 0
    count_NaN = 0
    with open(low_interactions_out, "w") as f_out:
        print("#columns:", "chr_x", "start_x", "end_x", "chr_y", "start_y", "end_y", sep="\t", file=f_out)
        for idx, (anno_pos, ctrl_pos) in enumerate(annotations):
            print("", idx, "/", len(annotations), end="\r", file=sys.stderr, flush=True)
            local_data = extract_value_cool(clr, bin_size, *anno_pos)

            global_data = extract_value_cool(clr, bin_size, *ctrl_pos)

            if len(local_data) == 0 or len(global_data) == 0:
                count_NaN += 1
                continue

            if mean(global_data) == 0:
                count_NaN += 1
                continue

            enrichment = mean(local_data) / mean(global_data)

            chr_x, start_x, end_x, chr_y, start_y, end_y = anno_pos
            if enrichment > 1.25:
                if chr_x == chr_y:
                    count_enriched_intra_chr += 1
                else:
                    count_enriched_inter_chr += 1
            if chr_x == chr_y:
                count_interactions_intra_chr += 1
            else:
                count_interactions_inter_chr += 1

            print(chr_x, start_x, end_x, chr_y, start_y, end_y, enrichment, sep="\t")
            if enrichment < 1:
                print(chr_x, start_x, end_x, chr_y, start_y, end_y, sep="\t", file=f_out)
    print("#without data or mean of control == 0:", count_NaN, file=sys.stderr)
    print("#enriched_intra_chr:", count_enriched_intra_chr, file=sys.stderr)
    print("#count_interactions_intra_chr:", count_interactions_intra_chr, file=sys.stderr)
    print(file=sys.stderr)
    print("#enriched_inter_chr:", count_enriched_inter_chr, file=sys.stderr)
    print("#count_interactions_inter_chr:", count_interactions_inter_chr, file=sys.stderr)


if __name__ == "__main__":
    individual_increase_quantification(*sys.argv[1:])
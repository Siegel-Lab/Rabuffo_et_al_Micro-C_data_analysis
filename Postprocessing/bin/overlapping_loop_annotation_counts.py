#### written by Markus Schmidt (LMU)
# usage: python3 overlapping_peak_loop_counts.py <loop_input_file> <require_peaks_on_both_axes> <annotation_input_file> <annotation_type>

import sys
from intervaltree import Interval, IntervalTree

annotation_trees = {}


total_annotations = 0
with open(sys.argv[3], "r") as in_file:
    for line in in_file:
        if len(line) == 0 or line[0] == "#":
            continue
        else:
            chrom, _, anno, start, end, *extra = line[:-1].split()
            if anno in sys.argv[4].split(","):
                if chrom not in annotation_trees:
                    annotation_trees[chrom] = IntervalTree()
                start = int(start)
                end = int(end)
                if end >= start:
                    annotation_trees[chrom].add(Interval(start, end))
                    total_annotations += 1


on_both_axes = sys.argv[2].lower() == "true"
print("total number of annotations of type", sys.argv[4], total_annotations)
total_overlaps = 0
total_loops = 0
with open(sys.argv[1], "r") as in_file:
    for line in in_file:
        if len(line) == 0 or line[0] == "#":
            continue
        else:
            chrom_1, start_1, end_1, chrom_2, start_2, end_2, *extra = line[:-1].split()
            start_1 = int(start_1)
            end_1 = int(end_1)
            start_2 = int(start_2)
            end_2 = int(end_2)
            num_overalaps = 0
            a = []
            b = []
            if chrom_1 in annotation_trees:
                if end_1 >= start_1:
                    a = annotation_trees[chrom_1].overlap(start_1, end_1)
                    num_overalaps += min(1, len(a))
            if chrom_2 in annotation_trees:
                if end_2 >= start_2:
                    b = annotation_trees[chrom_2].overlap(start_2, end_2)
                    num_overalaps += min(1, len(b))
            # print(line[:-1], num_overalaps, a, b)
            if on_both_axes and num_overalaps == 2:
                print(line[:-1])
                total_overlaps += 1
            elif not on_both_axes and num_overalaps >= 1:
                print(line[:-1])
                total_overlaps += 1
            total_loops += 1

print("total number of loops :", total_loops)

print("total number of annotations of type", sys.argv[4], "overlapping at least one loop", "on both axes:" if on_both_axes else "on either axis:", total_overlaps)

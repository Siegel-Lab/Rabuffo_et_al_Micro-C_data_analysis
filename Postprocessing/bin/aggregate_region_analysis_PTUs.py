import fileinput
import sys

# Usage:
# python3 aggregate_region_analysis.py in.pairs regions > out.pairs
#
# Use ge-zipped files as such:
# zcat in.pairs.gz | python3 aggregate_region_analysis.py - regions | gzip > out.pairs.gz
#
#
# the in_pairs file should be a upper-triangle pairs file.
#
# the regions file should be a tab-delimited file containing the regions.
# Each region is basically a stretch on the genome that gets translated into a stretch on the aggregated contigs.
# the file has the following columns:
# 1. contig name on the genome
# 1.5 the index of the region
# 2. start position of this region in the genome
# 3. end position of this region in the genome
# 4. start position this region in the aggregated contig
# 5. end position of this region in the aggregated contig
# 6. the orientation of the region (+/-)
#
#
# Version 2
#



def aggregate_region_analysis(in_pairs, regions, out_inter_region, out_inter_chrom, out_intra_chrom):
    region_dict = {}
    with fileinput.input(regions) as region_in:
        for line in region_in:
            if len(line) <= 1 or line[0] == '#':
                continue
            contig, idx, start_from, end_from, start_to, end_to, strand = line.strip().split()
            if contig not in region_dict:
                region_dict[contig] = []

            region_dict[contig].append([int(idx), int(start_from), int(end_from), int(start_to), int(end_to), strand])
    
    for _, value in region_dict.items():
        value.sort()
        for a, b in zip(value[:-1], value[1:]):
            if a > b:
                print("cannot have overlapping regions. Aborting", file=sys.stderr)
                exit()


    have_warned = False
    def fix_pos(ctg, pos, chrom_size):
        nonlocal have_warned
        for idx, start_from, end_from, start_to, end_to, strand in region_dict[ctg]:
            if start_from <= pos and pos < end_from:
                p_ret = int((end_to - start_to) * (pos - start_from) / (end_from - start_from) + start_to)
                assert p_ret >= start_to
                assert p_ret < end_to
                if strand == "-":
                    p_ret = chrom_size - p_ret - 1
                return idx, p_ret
        if not have_warned:
            print("Found interaction that is outside of any defined region. Ignoring this interaction", file=sys.stderr)
            print(line[:-1], file=sys.stderr)
            have_warned = True
        return 0, None

    with open(out_inter_region, "w") as out_inter_region_file, \
         open(out_inter_chrom, "w") as our_inter_chrom_file, \
         open(out_intra_chrom, "w") as out_intra_chrom_file:
        chrom_size = max(end_to for v in region_dict.values() for _, _, _, _, end_to, _ in v)
        for f in [out_inter_region_file, our_inter_chrom_file, out_intra_chrom_file]:
            print("## pairs format v1.0.0", file=f)
            print("#shape: upper-triangle", file=f)
            print("#chromsize: chromosome", chrom_size, file=f)
        with fileinput.input(in_pairs) as reads_file:
            for line in reads_file:
                if len(line) == 0 or line[0] == '#':
                    if line.startswith("#columns:") or line.startswith("#samheader:"):
                        for f in [out_inter_region_file, our_inter_chrom_file, out_intra_chrom_file]:
                            print(line[:-1], file=f)
                    continue
                readname, ctg1, pos1, ctg2, pos2, *extra = line[:-1].strip().split()
                
                while len(extra) < 12:
                    extra.append("")

                if ctg1 in region_dict and ctg2 in region_dict:
                    pos1 = int(pos1)
                    pos2 = int(pos2)
                    idx1, pos1 = fix_pos(ctg1, pos1, chrom_size)
                    idx2, pos2 = fix_pos(ctg2, pos2, chrom_size)
                    if pos1 is None or pos2 is None:
                        continue
                    if pos1 > pos2:
                        pos1, pos2 = pos2, pos1

                    if ctg1 == ctg2:
                        if idx1 == idx2:
                            print(readname, "chromosome", pos1, "chromosome", pos2, *extra, sep="\t", 
                                file=out_inter_region_file)
                        else:
                            print(readname, "chromosome", pos1, "chromosome", pos2, *extra, sep="\t", 
                                  file=our_inter_chrom_file)
                    else:
                        print(readname, "chromosome", pos1, "chromosome", pos2, *extra, sep="\t", 
                             file=out_intra_chrom_file)



if __name__ == '__main__':
    aggregate_region_analysis(*sys.argv[1:])
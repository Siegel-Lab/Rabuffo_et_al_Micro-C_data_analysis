import fileinput
import sys
import math
from scipy.stats import ks_2samp

# REQUIREMENTS:
# scipy needs to be installed -> pip install scipy
#
# USAGE:
# python3 compute_p_values.py <in_file> <local_size> <pov_chromosome> <enrichment_threshold>
#
# in_file: input file with the following columns "chromosome start end value"
# local_sizes: size of the local region to consider for the p-value computation (can be a comma separated list for multiple sizes)
# pov_chromosome: chromosome to consider as point of view (pov), data from this chromosome will be exluded. can be given as a comma separated list for multiple chromosomes.
# enrichment_threshold: minimal enrichment to test for


def compute_p_values(in_file, local_sizes="50000,100000,500000,1000000,5000000", 
                     pov_chromosome="chromosome_2B", anno_file="/dev/null", enrichment_threshold=2):
    local_sizes=[int(x) for x in local_sizes.split(',')]
    enrichment_threshold=float(enrichment_threshold)
    
    data = []
    for idx, line in enumerate(fileinput.input(in_file)):
        # if len(line) == 0:
        #     continue
        if line[0] == '#':
            print(line.strip())
            continue
        cols = line.strip().split()
        if cols[0] in pov_chromosome.split(","):
            cols[3] = float('NaN')
        data.append([cols[0], int(cols[1]), int(cols[2]), float(cols[3])])

    with open(anno_file, "w") as f_out:
        for chr_name, start_pos, end_pos, val in data:
            print(chr_name, start_pos, end_pos, sep='\t', end='')
            significant = False
            very_significant = False
            for local_size in local_sizes:
                local_data = [row[3] / enrichment_threshold for row in data if row[0] == chr_name and row[1] >= start_pos - local_size and row[2] <= end_pos + local_size and not math.isnan(row[3])]
                global_data = [row[3] for row in data if row[0] == chr_name and (row[1] < start_pos - local_size or row[2] > end_pos + local_size) and not math.isnan(row[3])]

                if len(local_data) == 0 or len(global_data) == 0:
                    print("\tNaN", end='')
                    continue

                p_val = ks_2samp(local_data, global_data, alternative="less").pvalue
                print("\t" + str(p_val), end='')
                significant = significant or p_val < 0.05
                very_significant = very_significant or p_val < 0.01
            print()
            if very_significant:
                print(chr_name, start_pos, end_pos, "**", file=f_out)
            elif significant:
                print(chr_name, start_pos, end_pos, "*", file=f_out)



if __name__ == '__main__':
    compute_p_values(*sys.argv[1:])

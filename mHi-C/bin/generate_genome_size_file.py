import sys
from Bio import SeqIO
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    print("{}\t{}".format(seq_record.id, len(seq_record)))

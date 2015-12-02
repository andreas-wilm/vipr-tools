import sys
import random
from Bio import SeqIO

fasta_in = sys.argv[1]
outfh = sys.stdout

seqrecs = list(SeqIO.parse(fasta_in, "fasta"))
               
num_contigs = 10
contig_no = 0
for sr in seqrecs:
    for i in range(num_contigs):
        start_pos = random.randint(0, len(sr.seq)-1)
        end_pos = random.randint(start_pos+1, len(sr.seq))
        seq = sr[start_pos:end_pos]
        if random.choice([0, 1]):
            seq = seq.reverse_complement()
            ori = "-"
        else:
            ori="+"
        contig_no += 1
        seqname = "contig-{}-{}-{}-ori{}".format(contig_no, start_pos, end_pos, ori)
        outfh.write(">{}\n{}\n".format(seqname, seq.seq))

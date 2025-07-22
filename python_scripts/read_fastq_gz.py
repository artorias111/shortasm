# read gzip files and print as a fasta

import gzip
from Bio import SeqIO

k = 31  # Set your desired k-mer length here - might have to be an input later


def reverse_complement(sequence):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca') # this is cool
    return sequence.translate(complement)[::-1]


fastq_file = '../tests/fastq_files/6_Swamp_S1_18S_2019_minq7.fastq.gz'
with gzip.open(fastq_file, 'rt') as f:
    for record in SeqIO.parse(f, "fastq"):
        seq = str(record.seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            revcomp = reverse_complement(kmer)
            canonical = min(kmer, revcomp)
            print(f"kmer: {kmer}\trevcomp: {revcomp}\tcanonical: {canonical}")
# read gzip files and print as a fasta

import gzip
from Bio import SeqIO

k = 6  # Set your desired k-mer length here - might have to be an input later


def reverse_complement(sequence):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca') # this is cool
    return sequence.translate(complement)[::-1]


fastq_file = '../tests/fastq_files/6_Swamp_S1_18S_2019_minq7.fastq.gz'
kmer_counts = {}
with gzip.open(fastq_file, 'rt') as f:
    for record in SeqIO.parse(f, "fastq"):
        seq = str(record.seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            revcomp = reverse_complement(kmer)
            canonical = min(kmer, revcomp)
            # print(f"kmer: {kmer}\trevcomp: {revcomp}\tcanonical: {canonical}")
            if canonical in kmer_counts:
                kmer_counts[canonical] += 1
            else:
                kmer_counts[canonical] = 1

# Print a sample of the kmer_counts dictionary
# note that these are only canonical k-mers, not all k-mers
for kmer, count in list(kmer_counts.items())[:10]:
    print(f"{kmer}: {count}")

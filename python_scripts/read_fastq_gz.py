# read gzip files and print as a fasta

import gzip
from Bio import SeqIO

fastq_file = '../tests/fastq_files/6_Swamp_S1_18S_2019_minq7.fastq.gz'
with gzip.open(fastq_file, 'rt') as f:
    for record in SeqIO.parse(f, "fastq"):
        print(">"+record.id)
        print(record.seq)
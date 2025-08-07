import gzip
from Bio import SeqIO
import argparse


class FastqProcessor:
    def __init__(self, reads_file, kmer_length=6):
        self.reads_file = reads_file
        self.kmer_length = kmer_length
        self.kmer_counts = {}
    
    def reverse_complement(self, sequence):
        complement = str.maketrans('ACGTacgt', 'TGCAtgca')
        return sequence.translate(complement)[::-1]
    
    def process_fastq(self):
        with gzip.open(self.reads_file, 'rt') as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                for i in range(len(seq) - self.kmer_length + 1):
                    kmer = seq[i:i+self.kmer_length]
                    revcomp = self.reverse_complement(kmer)
                    canonical = min(kmer, revcomp)
                    
                    if canonical in self.kmer_counts:
                        self.kmer_counts[canonical] += 1
                    else:
                        self.kmer_counts[canonical] = 1
    
    def print_sample_kmers(self, sample_size=10):
        print("=== K-mer Counts Sample ===")
        for kmer, count in list(self.kmer_counts.items())[:sample_size]:
            print(f"{kmer}: {count}")
    
    def get_kmer_counts(self):
        return self.kmer_counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', type=str, required=True)
    parser.add_argument('--kmer', type=int, default=6)
    args = parser.parse_args()
    
    processor = FastqProcessor(args.reads, args.kmer)
    
    print(f"Processing FASTQ file: {args.reads}")
    print(f"K-mer length: {args.kmer}")
    processor.process_fastq()
    
    processor.print_sample_kmers()


if __name__ == "__main__":
    main()

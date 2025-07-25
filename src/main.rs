use clap::Parser;
use needletail::parse_fastx_file;
use needletail::Sequence;
use std::str;

/// Produce De Bruijn graphs from short reads
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]

struct Args {
    /// Path to short reads file
    #[arg(short, long)]
    reads: String,

    /// k-mer size
    #[arg(short, long, default_value_t = 7)]
    kmer: u8,
}

fn main() {
    let args = Args::parse();

    let kmer_size = args.kmer;
    // Path to the test fastq file :
    // /Users/sbhat/Documents/projects/shortasm/tests/fastq_files/6_Swamp_S1_18S_2019_minq7.fastq.gz

    println!("Your short reads file is {}, and your k-mer size is {}", args.reads, args.kmer);
    let mut reader = parse_fastx_file(&args.reads).expect("the fastq path is not valid");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let seqid = seqrec.id();
        let seqvec = seqrec.seq();
        let id_str = str::from_utf8(seqid).unwrap();
        let seq_str = str::from_utf8(&seqvec).unwrap();

        // k-merize the sequence

        println!("k-mers for {id_str}");

        for k in seqrec.kmers(kmer_size) {
            let kmer_str = str::from_utf8(k).unwrap();
            println!("{kmer_str}");
        }

        // println!("{seq_str}");

    }
}

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
    // /Users/sbhat/Documents/projects/shortasm/tests/fastq_files/6_Swamp_S1_18S_2019_minq7.fastq.gz - path to the test fastq file

    println!("Your short reads file is {}, and your k-mer size is {}", args.reads, args.kmer);
    let mut reader = parse_fastx_file(&args.reads).expect("the fastq path is not valid");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let seqid = seqrec.id();
        let seqvec = seqrec.seq();
        let id_str = str::from_utf8(seqid).unwrap();
        let _seq_str = str::from_utf8(&seqvec).unwrap(); // unused variable

        // normalize the sequence to ensure all nucleotides are uppercase, and all special
        // character from the fastx file is removed


        let seq_norm = seqrec.normalize(false);
        // get the reverse complement
        let rc = seq_norm.reverse_complement();

        for (_, kmer_obj, rc_flag) in seq_norm.canonical_kmers(kmer_size, &rc) { // I'm sure this
                                                                                 // part can be a
                                                                                 // lot more
                                                                                 // idiomatic, but
                                                                                 // hey, I'm still
                                                                                 // learning Rust
                                                                                 // :D
            if rc_flag == true {
                let og_seq = kmer_obj;
                let rc_kmer = kmer_obj.reverse_complement();
                let og_seq_str = str::from_utf8(&og_seq).unwrap();
                let rc_kmer_str = str::from_utf8(&rc_kmer).unwrap();
                let kmer_str = str::from_utf8(kmer_obj).unwrap();
                println!("orignal sequence: {}\tReverse complement: {}\tCanonical k-mer:{}", og_seq_str, rc_kmer_str, kmer_str);

            } else {
                let og_seq = kmer_obj.reverse_complement();
                let rc_kmer = kmer_obj;
                let og_seq_str = str::from_utf8(&og_seq).unwrap();
                let rc_kmer_str = str::from_utf8(&rc_kmer).unwrap();
                let kmer_str = str::from_utf8(kmer_obj).unwrap();
                println!("orignal sequence: {}\tReverse complement: {}\tCanonical k-mer:{}", og_seq_str, rc_kmer_str, kmer_str);
            }

        }

        /*
        // k-merize the sequence

        println!("k-mers for {id_str}");

        for k in seqrec.kmers(kmer_size) {
            let kmer_str = str::from_utf8(k).unwrap();
            println!("{kmer_str}");
        }

        // println!("{seq_str}");
        */
    }
}

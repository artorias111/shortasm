use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use std::str;
use std::collections::HashMap;

/// Produce De Bruijn graphs from short reads
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]

struct Args {
    /// Path to short reads file
    #[arg(short, long)]
    reads: String,

    /// k-mer size
    #[arg(short, long, default_value_t = 6)]
    kmer: u8,
}

fn main() {
    let args = Args::parse();
    let kmer_size = args.kmer;

    println!("Your short reads file is {}, and your k-mer size is {}", args.reads, args.kmer);
    let mut reader = parse_fastx_file(&args.reads).expect("the fastq path is not valid");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let seqid = seqrec.id();
        let seqvec = seqrec.seq();
        let _id_str = str::from_utf8(seqid).unwrap();
        let _seq_str = str::from_utf8(&seqvec).unwrap(); 

        // normalize the sequence to ensure all nucleotides are uppercase, and all special
        // character from the fastx file is removed


        let seq_norm = seqrec.normalize(false);
        // get the reverse complement
        let rc = seq_norm.reverse_complement();

        // store canonical kmers
        let mut kmer_count:HashMap<&[u8], usize> = HashMap::new();

        for (_, kmer_obj, rc_flag) in seq_norm.canonical_kmers(kmer_size, &rc) { // I'm sure this
                                                                                 // part can be a
                                                                                 // lot more
                                                                                 // idiomatic, but
                                                                                 // hey, I'm still
                                                                                 // learning Rust
                                                                                 // :D
            kmer_count.entry(kmer_obj)
                .and_modify(|e| { *e += 1 })
                .or_insert(1);


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

        // get your k-mer count hashmap, prints once for each record though, need to fix that

        for (k, v) in &kmer_count {
            let kmer_str = str::from_utf8(k).unwrap();
            let count = v;
            println!("{kmer_str}: {count}");
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

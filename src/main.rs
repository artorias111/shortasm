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

fn count_kmers(file: &str, kmer_size: u8) ->HashMap<String, usize> {
    let mut reader = parse_fastx_file(&file).expect("the fastq path is not valid");
    let mut kmer_count:HashMap<String, usize> = HashMap::new();

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

        for (_, kmer_obj, rc_flag) in seq_norm.canonical_kmers(kmer_size, &rc) { // I'm sure this
                                                                                 // block can be a
                                                                                 // lot more
                                                                                 // idiomatic, but
                                                                                 // hey, I'm still
                                                                                 // learning Rust:D
            kmer_count.entry(str::from_utf8(kmer_obj).unwrap().to_string())
                .and_modify(|e| { *e += 1 })
                .or_insert(1);

            if rc_flag {
                let og_seq = kmer_obj;
                let rc_kmer = kmer_obj.reverse_complement();
                let _og_seq_str = str::from_utf8(&og_seq).unwrap();
                let _rc_kmer_str = str::from_utf8(&rc_kmer).unwrap();
                let _kmer_str = str::from_utf8(kmer_obj).unwrap();
                // println!("orignal sequence: {}\tReverse complement: {}\tCanonical k-mer:{}", og_seq_str, rc_kmer_str, kmer_str);

            } else {
                let og_seq = kmer_obj.reverse_complement();
                let rc_kmer = kmer_obj;
                let _og_seq_str = str::from_utf8(&og_seq).unwrap();
                let _rc_kmer_str = str::from_utf8(&rc_kmer).unwrap();
                let _kmer_str = str::from_utf8(kmer_obj).unwrap();
                // println!("orignal sequence: {}\tReverse complement: {}\tCanonical k-mer:{}", og_seq_str, rc_kmer_str, kmer_str);
            }
        }
    }

    kmer_count
}

// new main function, much cleaner :D
fn main() {
    let args = Args::parse();
    // let mut kmer_count:HashMap<String, usize> = HashMap::new();

    println!("Your short reads file is {}, and your k-mer size is {}", args.reads, args.kmer);
    
    let kmer_count = count_kmers(&args.reads, args.kmer);
    let mut counts = 0;
    for (k, v) in &kmer_count {
        println!("{k}: {v}");
        counts +=1;
        if counts >10 {
            break;
        }
    }
}



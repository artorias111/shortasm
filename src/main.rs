mod cli;
use clap::Parser;
use cli::Args;
use needletail::{Sequence, parse_fastx_file};
use std::collections::HashMap;
use std::str;

fn count_kmers(file: &str, kmer_size: u8) -> anyhow::Result<HashMap<String, usize>> {
    let mut reader = parse_fastx_file(file)?;
    let mut kmer_count = HashMap::new();

    while let Some(record) = reader.next() {
        let seqrec = record?;
        let seq_norm = seqrec.normalize(false);
        let rc = seq_norm.reverse_complement();

        // store canonical k-mers
        for (_, kmer, _) in seq_norm.canonical_kmers(kmer_size, &rc) {
            let kmer_str = str::from_utf8(kmer)?.to_string();
            *kmer_count.entry(kmer_str).or_insert(0) += 1;
        }
    }

    Ok(kmer_count)
}

// new main function, much cleaner :D
fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    // let mut kmer_count:HashMap<String, usize> = HashMap::new();

    if args.kmer % 2 == 0 || args.kmer > 31 {
        anyhow::bail!("k must be odd and at most 31, got {}", args.kmer);
    }

    println!(
        "Your short reads file is {}, and your k-mer size is {}",
        args.reads, args.kmer
    );

    let kmer_count = count_kmers(&args.reads, args.kmer)?;

    for (kmer, count) in kmer_count.iter().take(10) {
        println!("{kmer} : {count}");
    }

    Ok(())
}

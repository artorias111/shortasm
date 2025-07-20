use clap::Parser;

/// PRoduce De Bruijn graphs from short reads
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]

struct Args {
    /// Directory containing all your short reads
    #[arg(short, long)]
    reads: String,

    /// k-mer size
    #[arg(short, long, default_value_t = 7)]
    kmer: u8,
}

fn main() {
    let args = Args::parse();

    println!("Your short reads are in {}, and your k-mer size is {}", args.reads, args.kmer);
}

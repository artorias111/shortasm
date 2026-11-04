use clap::Parser;

/// Produce De Bruijn graphs from short reads
/// Nodes: k-mers, edges: overlapping k-mers (k-1)
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Path to short reads file
    #[arg(short, long)]
    pub reads: String, // maybe i can make it a glob of fastq(gz) in the future

    /// k-mer size
    #[arg(short, long, default_value_t = 31)]
    pub kmer: u8,
}

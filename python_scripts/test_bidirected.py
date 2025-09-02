#!/usr/bin/env python3
"""
Test script to demonstrate bidirected DNA functionality in de Bruijn graphs.
"""

from de_bruijn_graph_builder import DeBruijnGraphBuilder

def test_reverse_complement():
    """Test reverse complement functionality with simple examples."""
    print("=== Testing Reverse Complement Functionality ===")
    
    # Create a builder instance (we won't actually build a graph for this test)
    builder = DeBruijnGraphBuilder("dummy.fastq", kmer_length=4)
    
    test_kmers = ["ATGC", "GCTA", "AAAA", "ATAT", "GCGC"]
    
    for kmer in test_kmers:
        rev_comp = builder.get_reverse_complement(kmer)
        canonical = builder.get_canonical_kmer(kmer)
        
        print(f"K-mer: {kmer}")
        print(f"  Reverse complement: {rev_comp}")
        print(f"  Canonical form: {canonical}")
        print(f"  Self-complementary: {kmer == rev_comp}")
        print()

def test_bidirected_edge_logic():
    """Test the logic for creating bidirected edges."""
    print("=== Testing Bidirected Edge Logic ===")
    
    builder = DeBruijnGraphBuilder("dummy.fastq", kmer_length=4)
    
    # Test with a 4-mer
    kmer = "ATGC"
    print(f"Original k-mer: {kmer}")
    
    # Get prefix and suffix
    prefix = kmer[:-1]  # "ATG"
    suffix = kmer[1:]   # "TGC"
    print(f"Prefix: {prefix}, Suffix: {suffix}")
    
    # Forward edge: ATG -> TGC
    print(f"Forward edge: {prefix} -> {suffix}")
    
    # Reverse complement edge
    rev_comp_prefix = builder.get_reverse_complement(prefix)  # "CAT"
    rev_comp_suffix = builder.get_reverse_complement(suffix)  # "GCA"
    print(f"Reverse complement edge: {rev_comp_suffix} -> {rev_comp_prefix}")
    
    print(f"Note: {rev_comp_suffix} is reverse complement of {suffix}")
    print(f"      {rev_comp_prefix} is reverse complement of {prefix}")
    print()

def demonstrate_dna_bidirectionality():
    """Demonstrate the bidirected nature of DNA with examples."""
    print("=== DNA Bidirectionality Demonstration ===")
    
    print("""
DNA is double-stranded and antiparallel:
- Strand 1: 5'-ATGC-3'
- Strand 2: 3'-TACG-5'

The reverse complement of ATGC is TACG:
- Reverse: ATGC -> CGTA
- Complement: CGTA -> TACG

In a de Bruijn graph:
- When we see k-mer ATGC, we create edge ATG -> TGC
- We also need to consider the reverse complement TACG
- This creates edge GCA -> CAT (reverse complement of TGC -> ATG)

This ensures our graph captures both DNA strands!
    """)

def test_canonical_kmers():
    """Test canonical k-mer selection."""
    print("=== Testing Canonical K-mer Selection ===")
    
    builder = DeBruijnGraphBuilder("dummy.fastq", kmer_length=4)
    
    test_pairs = [
        ("ATGC", "GCAT"),  # ATGC < GCAT lexicographically
        ("GCAT", "ATGC"),  # Same as above, should both return ATGC
        ("AAAA", "TTTT"),  # AAAA < TTTT
        ("ATAT", "ATAT"),  # Self-complementary
    ]
    
    for kmer1, kmer2 in test_pairs:
        canonical1 = builder.get_canonical_kmer(kmer1)
        canonical2 = builder.get_canonical_kmer(kmer2)
        
        print(f"K-mer 1: {kmer1} -> Canonical: {canonical1}")
        print(f"K-mer 2: {kmer2} -> Canonical: {canonical2}")
        print(f"Same canonical: {canonical1 == canonical2}")
        print()

if __name__ == "__main__":
    test_reverse_complement()
    test_bidirected_edge_logic()
    demonstrate_dna_bidirectionality()
    test_canonical_kmers()
    
    print("=== Summary ===")
    print("""
The bidirected de Bruijn graph implementation:

1. **Handles DNA bidirectionality**: Each k-mer is processed along with its reverse complement
2. **Maintains edge symmetry**: If edge A->B exists, then rev_comp(B)->rev_comp(A) also exists
3. **Uses canonical k-mers**: Lexicographically smaller of k-mer and its reverse complement
4. **Captures both DNA strands**: Essential for accurate genome assembly
5. **Validates structure**: Can check if the graph properly represents bidirected DNA

This implementation ensures that the de Bruijn graph correctly represents the 
double-stranded, antiparallel nature of DNA molecules.
    """)

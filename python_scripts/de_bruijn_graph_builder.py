import networkx as nx
from read_fastq_gz import FastqProcessor


class DeBruijnGraphBuilder:
    def __init__(self, reads_file, kmer_length=6):
        self.reads_file = reads_file
        self.kmer_length = kmer_length
        self.graph = nx.DiGraph()
        self.kmer_counts = {}
        
        # DNA complement mapping
        self.complement_map = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
        }
    
    def get_reverse_complement(self, sequence):
        """
        Compute the reverse complement of a DNA sequence.
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            str: Reverse complement sequence
        """
        # First reverse the sequence, then complement each base
        reversed_seq = sequence[::-1]
        complement = ''.join(self.complement_map.get(base, base) for base in reversed_seq)
        return complement
    
    def get_canonical_kmer(self, kmer):
        """
        Get the canonical form of a k-mer (lexicographically smaller of k-mer and its reverse complement).
        
        Args:
            kmer (str): DNA k-mer
            
        Returns:
            str: Canonical k-mer
        """
        reverse_complement = self.get_reverse_complement(kmer)
        return min(kmer, reverse_complement)
    
    def build_graph_from_kmers(self):
        processor = FastqProcessor(self.reads_file, self.kmer_length)
        processor.process_fastq()
        self.kmer_counts = processor.get_kmer_counts()

        kmer_count_len = len(self.kmer_counts)
        
        print(f"Building bidirected de Bruijn graph from {kmer_count_len} k-mers...")
        print("Including both canonical and reverse complement edges for DNA bidirectionality...")
        
        # Process each k-mer and its reverse complement
        for kmer, count in self.kmer_counts.items():
            # Get the canonical form of the k-mer
            canonical_kmer = self.get_canonical_kmer(kmer)
            
            # Get reverse complement of the k-mer
            reverse_complement_kmer = self.get_reverse_complement(kmer)
            
            # Process canonical k-mer edges
            self._add_kmer_edges(canonical_kmer, count)
            
            # Process reverse complement k-mer edges
            self._add_kmer_edges(reverse_complement_kmer, count)
    
    def _add_kmer_edges(self, kmer, count):
        """
        Add edges for a k-mer to the graph.
        
        Args:
            kmer (str): DNA k-mer
            count (int): Count of this k-mer
        """
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        # Add the forward edge (prefix -> suffix)
        if self.graph.has_edge(prefix, suffix):
            current_count = self.graph[prefix][suffix]['count']
            self.graph[prefix][suffix]['count'] = current_count + count
        else:
            self.graph.add_edge(prefix, suffix, count=count)
        
        # For bidirected DNA, we also need to consider the reverse complement edge
        # The reverse complement of prefix -> suffix is (reverse_complement of suffix) -> (reverse_complement of prefix)
        rev_comp_suffix = self.get_reverse_complement(suffix)
        rev_comp_prefix = self.get_reverse_complement(prefix)
        
        # Add the reverse complement edge
        if self.graph.has_edge(rev_comp_suffix, rev_comp_prefix):
            current_count = self.graph[rev_comp_suffix][rev_comp_prefix]['count']
            self.graph[rev_comp_suffix][rev_comp_prefix]['count'] = current_count + count
        else:
            self.graph.add_edge(rev_comp_suffix, rev_comp_prefix, count=count)
    
    def print_graph_stats(self):
        print("\n=== Bidirected De Bruijn Graph Statistics ===")
        print(f"Number of nodes: {self.graph.number_of_nodes()}")
        print(f"Number of edges: {self.graph.number_of_edges()}")
        
        print("\nSample edges with counts:")
        edge_count = 0
        for u, v, data in self.graph.edges(data=True):
            if edge_count < 10:
                print(f"  {u} -> {v} (count: {data['count']})")
                edge_count += 1
            else:
                break
        
        in_degrees = dict(self.graph.in_degree())
        out_degrees = dict(self.graph.out_degree())
        
        print(f"\nNode degree analysis:")
        print(f"  Average in-degree: {sum(in_degrees.values()) / len(in_degrees):.2f}")
        print(f"  Average out-degree: {sum(out_degrees.values()) / len(out_degrees):.2f}")
        
        max_in_degree = max(in_degrees.values()) if in_degrees else 0
        max_out_degree = max(out_degrees.values()) if out_degrees else 0
        
        print(f"  Max in-degree: {max_in_degree}")
        print(f"  Max out-degree: {max_out_degree}")
        
        # Show bidirected nature
        print(f"\nBidirected DNA Analysis:")
        bidirectional_edges = 0
        for u, v, data in self.graph.edges(data=True):
            if self.graph.has_edge(v, u):
                bidirectional_edges += 1
        
        print(f"  Bidirectional edge pairs: {bidirectional_edges // 2}")
        print(f"  Total edges: {self.graph.number_of_edges()}")
        print(f"  Bidirectionality ratio: {bidirectional_edges / self.graph.number_of_edges():.2f}")
    
    def get_high_weight_edges(self, threshold=10):
        high_weight_edges = []
        for u, v, data in self.graph.edges(data=True):
            if data['count'] >= threshold:
                high_weight_edges.append((u, v, data['count']))
        
        return sorted(high_weight_edges, key=lambda x: x[2], reverse=True)
    
    def print_high_weight_edges(self, threshold=10, max_show=20):
        high_weight_edges = self.get_high_weight_edges(threshold)
        
        print(f"\n=== High Weight Edges (count >= {threshold}) ===")
        for i, (u, v, count) in enumerate(high_weight_edges[:max_show]):
            print(f"  {i+1}. {u} -> {v} (count: {count})")
        
        if len(high_weight_edges) > max_show:
            print(f"  ... and {len(high_weight_edges) - max_show} more edges")
    
    def demonstrate_bidirected_nature(self, max_examples=5):
        """
        Demonstrate the bidirected nature of the DNA graph by showing
        reverse complement relationships between edges.
        """
        print(f"\n=== Bidirected DNA Demonstration ===")
        print("Showing reverse complement edge pairs:")
        
        examples_shown = 0
        for u, v, data in self.graph.edges(data=True):
            if examples_shown >= max_examples:
                break
                
            # Get reverse complement of the edge
            rev_comp_u = self.get_reverse_complement(u)
            rev_comp_v = self.get_reverse_complement(v)
            
            # Check if the reverse complement edge exists
            if self.graph.has_edge(rev_comp_v, rev_comp_u):
                rev_data = self.graph[rev_comp_v][rev_comp_u]
                print(f"  Forward:  {u} -> {v} (count: {data['count']})")
                print(f"  Reverse:  {rev_comp_v} -> {rev_comp_u} (count: {rev_data['count']})")
                print(f"  Note: {rev_comp_v} is reverse complement of {u}")
                print(f"        {rev_comp_u} is reverse complement of {v}")
                print()
                examples_shown += 1
    
    def validate_bidirected_structure(self):
        """
        Validate that the graph properly represents bidirected DNA structure.
        Returns True if the graph is properly bidirected.
        """
        print(f"\n=== Validating Bidirected Structure ===")
        
        total_edges = self.graph.number_of_edges()
        bidirectional_pairs = 0
        missing_complements = 0
        
        for u, v, data in self.graph.edges(data=True):
            # Get reverse complement of the edge
            rev_comp_u = self.get_reverse_complement(u)
            rev_comp_v = self.get_reverse_complement(v)
            
            if self.graph.has_edge(rev_comp_v, rev_comp_u):
                bidirectional_pairs += 1
            else:
                missing_complements += 1
                print(f"  Missing reverse complement: {rev_comp_v} -> {rev_comp_u}")
        
        print(f"  Total edges: {total_edges}")
        print(f"  Bidirectional pairs: {bidirectional_pairs // 2}")
        print(f"  Missing complements: {missing_complements}")
        
        # Each edge should have a reverse complement (except self-complementary edges)
        is_valid = missing_complements == 0
        print(f"  Graph is properly bidirected: {is_valid}")
        
        return is_valid
    
    def get_reverse_complement_analysis(self, kmer):
        """
        Analyze a k-mer and show its reverse complement relationships.
        
        Args:
            kmer (str): DNA k-mer to analyze
        """
        print(f"\n=== Reverse Complement Analysis for '{kmer}' ===")
        
        rev_comp = self.get_reverse_complement(kmer)
        canonical = self.get_canonical_kmer(kmer)
        
        print(f"  Original k-mer: {kmer}")
        print(f"  Reverse complement: {rev_comp}")
        print(f"  Canonical form: {canonical}")
        
        if kmer == rev_comp:
            print(f"  Self-complementary: Yes")
        else:
            print(f"  Self-complementary: No")
        
        # Show the edges that would be created
        prefix = kmer[:-1]
        suffix = kmer[1:]
        rev_prefix = self.get_reverse_complement(prefix)
        rev_suffix = self.get_reverse_complement(suffix)
        
        print(f"  Forward edge: {prefix} -> {suffix}")
        print(f"  Reverse complement edge: {rev_suffix} -> {rev_prefix}")
        
        # Check if these edges exist in the graph
        if self.graph.has_edge(prefix, suffix):
            print(f"  Forward edge exists in graph: Yes (count: {self.graph[prefix][suffix]['count']})")
        else:
            print(f"  Forward edge exists in graph: No")
            
        if self.graph.has_edge(rev_suffix, rev_prefix):
            print(f"  Reverse edge exists in graph: Yes (count: {self.graph[rev_suffix][rev_prefix]['count']})")
        else:
            print(f"  Reverse edge exists in graph: No")
    
    def get_graph(self):
        return self.graph
    
    def get_kmer_counts(self):
        return self.kmer_counts


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Build bidirected de Bruijn graph from FASTQ reads")
    parser.add_argument('--reads', type=str, required=True, help="Input FASTQ file")
    parser.add_argument('--kmer', type=int, default=6, help="K-mer length")
    parser.add_argument('--threshold', type=int, default=10, help="Threshold for high-weight edges")
    parser.add_argument('--analyze-kmer', type=str, help="Analyze a specific k-mer and its reverse complement")
    parser.add_argument('--validate', action='store_true', help="Validate bidirected structure")
    parser.add_argument('--demonstrate', action='store_true', help="Demonstrate bidirected nature")
    args = parser.parse_args()
    
    builder = DeBruijnGraphBuilder(args.reads, args.kmer)
    
    print(f"Processing {args.reads} with k-mer length {args.kmer}")
    print("Building bidirected de Bruijn graph with reverse complement handling...")
    builder.build_graph_from_kmers()
    
    builder.print_graph_stats()
    
    # Demonstrate bidirected nature if requested
    if args.demonstrate:
        builder.demonstrate_bidirected_nature()
    
    # Validate bidirected structure if requested
    if args.validate:
        builder.validate_bidirected_structure()
    
    # Analyze specific k-mer if provided
    if args.analyze_kmer:
        if len(args.analyze_kmer) == args.kmer:
            builder.get_reverse_complement_analysis(args.analyze_kmer)
        else:
            print(f"Error: K-mer length must be {args.kmer}, but '{args.analyze_kmer}' has length {len(args.analyze_kmer)}")
    
    builder.print_high_weight_edges(args.threshold)
    
    print(f"\n=== Journal Entry: Bidirected DNA Implementation ===")
    print("""
**Learning: Understanding Bidirected Graphs for DNA**

1. **DNA Bidirectionality**: DNA is double-stranded and antiparallel. Each strand has a reverse complement.
   - A pairs with T, G pairs with C
   - Strands run in opposite directions (5' to 3' vs 3' to 5')

2. **Reverse Complement Logic**:
   - For a k-mer like 'ATG', the reverse complement is 'CAT'
   - Process: reverse the sequence ('ATG' -> 'GTA'), then complement each base ('GTA' -> 'CAT')

3. **Bidirected Graph Representation**:
   - When adding edge ATG -> TGC, also add reverse complement edge GCA -> CAT
   - This ensures the graph captures both DNA strands
   - Each edge should have a corresponding reverse complement edge

4. **Implementation Details**:
   - Canonical k-mers: Use lexicographically smaller of k-mer and its reverse complement
   - Edge symmetry: If A->B exists, then rev_comp(B)->rev_comp(A) should also exist
   - Self-complementary k-mers: Palindromic sequences like 'ATAT' are their own reverse complement

5. **Benefits**:
   - More accurate representation of DNA structure
   - Captures both strands of double-stranded DNA
   - Essential for proper genome assembly
   - Handles the natural bidirectionality of DNA replication and transcription
    """)


if __name__ == "__main__":
    main() 
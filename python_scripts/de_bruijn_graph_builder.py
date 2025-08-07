import networkx as nx
from read_fastq_gz import FastqProcessor


class DeBruijnGraphBuilder:
    def __init__(self, reads_file, kmer_length=6):
        self.reads_file = reads_file
        self.kmer_length = kmer_length
        self.graph = nx.DiGraph()
        self.kmer_counts = {}
    
    def build_graph_from_kmers(self):
        processor = FastqProcessor(self.reads_file, self.kmer_length)
        processor.process_fastq()
        self.kmer_counts = processor.get_kmer_counts()

        kmer_count_len = len(self.kmer_counts)
        
        print(f"Building de Bruijn graph from {kmer_count_len} k-mers...")
        
        for kmer, count in self.kmer_counts.items():
            prefix = kmer[:-1]
            suffix = kmer[1:]
            
            if self.graph.has_edge(prefix, suffix):
                current_count = self.graph[prefix][suffix]['count']
                self.graph[prefix][suffix]['count'] = current_count + count
            else:
                self.graph.add_edge(prefix, suffix, count=count)
    
    def print_graph_stats(self):
        print("\n=== De Bruijn Graph Statistics ===")
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
    
    def get_graph(self):
        return self.graph
    
    def get_kmer_counts(self):
        return self.kmer_counts


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Build de Bruijn graph from FASTQ reads")
    parser.add_argument('--reads', type=str, required=True, help="Input FASTQ file")
    parser.add_argument('--kmer', type=int, default=6, help="K-mer length")
    parser.add_argument('--threshold', type=int, default=10, help="Threshold for high-weight edges")
    args = parser.parse_args()
    
    builder = DeBruijnGraphBuilder(args.reads, args.kmer)
    
    print(f"Processing {args.reads} with k-mer length {args.kmer}")
    builder.build_graph_from_kmers()
    
    builder.print_graph_stats()
    
    builder.print_high_weight_edges(args.threshold)


if __name__ == "__main__":
    main() 
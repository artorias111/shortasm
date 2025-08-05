# NetworkX directed graph experiments
import networkx as nx


class GraphExperimenter:
    def __init__(self):
        self.graph = nx.DiGraph()
    
    def create_basic_graph(self):
        """Create a basic directed graph with nodes and edges."""
        # Add nodes
        self.graph.add_node("A")
        self.graph.add_node("B")
        self.graph.add_node("C")
        self.graph.add_node("D")
        
        # Add edges
        self.graph.add_edge("A", "B")
        self.graph.add_edge("B", "C")
        self.graph.add_edge("C", "D")
        self.graph.add_edge("D", "A")  # Creates a cycle
        self.graph.add_edge("A", "C")  # Creates a shortcut
    
    def print_graph_info(self):
        """Print basic information about the graph."""
        print("=== NetworkX Directed Graph Experiment ===")
        print(f"Number of nodes: {self.graph.number_of_nodes()}")
        print(f"Number of edges: {self.graph.number_of_edges()}")
        print(f"Nodes: {list(self.graph.nodes())}")
        print(f"Edges: {list(self.graph.edges())}")
    
    def analyze_graph(self):
        """Perform basic graph analysis."""
        # Check if graph has cycles
        has_cycles = nx.is_directed_acyclic_graph(self.graph)
        print(f"Is acyclic: {has_cycles}")
        
        # Basic graph analysis
        print(f"Node degrees: {dict(self.graph.degree())}")
        print(f"In-degrees: {dict(self.graph.in_degree())}")
        print(f"Out-degrees: {dict(self.graph.out_degree())}")
    
    def create_custom_graph(self, nodes, edges):
        """Create a custom graph with specified nodes and edges."""
        self.graph.clear()
        
        # Add nodes
        for node in nodes:
            self.graph.add_node(node)
        
        # Add edges
        for edge in edges:
            self.graph.add_edge(edge[0], edge[1])
    
    def get_graph(self):
        """Return the graph object."""
        return self.graph


def main():
    """Main function to demonstrate NetworkX functionality."""
    experimenter = GraphExperimenter()
    
    # Create and analyze basic graph
    experimenter.create_basic_graph()
    experimenter.print_graph_info()
    experimenter.analyze_graph()
    
    print("\n=== Custom Graph Example ===")
    # Create a custom graph
    custom_nodes = ["X", "Y", "Z"]
    custom_edges = [("X", "Y"), ("Y", "Z"), ("Z", "X")]
    
    experimenter.create_custom_graph(custom_nodes, custom_edges)
    experimenter.print_graph_info()
    experimenter.analyze_graph()


if __name__ == "__main__":
    main() 
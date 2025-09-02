import networkx as nx


class GraphExperimenter:
    def __init__(self):
        self.graph = nx.DiGraph()
    
    def create_basic_graph(self):
        self.graph.add_node("A")
        self.graph.add_node("B")
        self.graph.add_node("C")
        self.graph.add_node("D")
        
        self.graph.add_edge("A", "B")
        self.graph.add_edge("B", "C")
        self.graph.add_edge("C", "D")
        self.graph.add_edge("D", "A")
        self.graph.add_edge("A", "C")
    
    def print_graph_info(self):
        print("=== NetworkX Directed Graph Experiment ===")
        print(f"Number of nodes: {self.graph.number_of_nodes()}")
        print(f"Number of edges: {self.graph.number_of_edges()}")
        print(f"Nodes: {list(self.graph.nodes())}")
        print(f"Edges: {list(self.graph.edges())}")
    
    def analyze_graph(self):
        has_cycles = nx.is_directed_acyclic_graph(self.graph)
        print(f"Is acyclic: {has_cycles}")
        
        print(f"Node degrees: {dict(self.graph.degree())}")
        print(f"In-degrees: {dict(self.graph.in_degree())}")
        print(f"Out-degrees: {dict(self.graph.out_degree())}")
    
    def create_custom_graph(self, nodes, edges):
        self.graph.clear()
        
        for node in nodes:
            self.graph.add_node(node)
        
        for edge in edges:
            self.graph.add_edge(edge[0], edge[1])
    
    def get_graph(self):
        return self.graph


def main():
    experimenter = GraphExperimenter()
    
    experimenter.create_basic_graph()
    experimenter.print_graph_info()
    experimenter.analyze_graph()
    
    print("\n=== Custom Graph Example ===")
    custom_nodes = ["X", "Y", "Z"]
    custom_edges = [("X", "Y"), ("Y", "Z"), ("Z", "X")] # list of tuples for edges, are they directed? 
                                                        # I don't think so, but tuples are ordered, so I'm not sure, maybe they ARE directed
    
    experimenter.create_custom_graph(custom_nodes, custom_edges)
    experimenter.print_graph_info()
    experimenter.analyze_graph()


if __name__ == "__main__":
    main() 

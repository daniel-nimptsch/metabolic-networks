import queue
import networkx as nx
from networkx.algorithms import isomorphism
import parse_smiles
import traversal
import constants
# import util
# from pyvis.network import Network
from queue import SimpleQueue

# from bio_pathway_comparison import find_synthesis_paths, reaction_paths_unique_test
# organism_name = "ecoli_cimIV"
queue = SimpleQueue()

# Function to print if amino acids or cofactors are in the graph
def print_AS_cofactor_info(G):
    list_of_all_products = traversal.breadth_first_search(G, "D-glucose")
    for amino_acid in constants.amino_acid_list:
        if amino_acid in list_of_all_products:
            print(amino_acid + " is there")
        else:
            print(amino_acid + " is NOT there")
    for cofactor in constants.all_cofactors:
        if cofactor in list_of_all_products:
            print(cofactor + " is there")
        else:
            print(cofactor + " is NOT there")

# Functio to print reactions associated with a metabolite
def print_metabolite_reactions(G, metabolite):
    for pre in G.predecessors(metabolite):
        print(G.nodes[pre]['description'])
    for neigh in G.neighbors(metabolite):
        print(G.nodes[neigh]['description'])

# Main analysis for WP1
# TODO:
# Make alternative paths finding faster!
# Save results as files!
def analysis():
    for organism_name in constants.path_dict:
        G = parse_smiles.parse_smiles(organism_name, add_leucine=True)
        print(f'{organism_name} full:', len(list(G.nodes)))
        H = nx.DiGraph()
        P = traversal.get_subgraph_glucose_S(G)
        for amino_acid in constants.amino_acid_list:
            H = nx.compose(H, traversal.breadth_first_search_reverse_amino_acid(G, amino_acid))
        R = P.copy()
        R.remove_nodes_from(n for n in P if n not in H)
        R.remove_edges_from(e for e in P.edges if e not in H.edges)
        print(f'{organism_name} AS-bio-pw:', len(list(R.nodes)))
        print("Metabolite info full:")
        print_AS_cofactor_info(G)
        print("Metabolite info AS-bio-pw:")
        print_AS_cofactor_info(R)

        for amino_acid in constants.amino_acid_list:
            print(amino_acid)    
            # T, paths = find_synthesis_paths(R, amino_acid)
        # Test
        # T, paths = find_synthesis_paths(R, 'L-leucine')

        # Visualization
        # nt = Network('1000px', '1000px', select_menu=True, filter_menu=True)
        # nt.from_nx(paths[0])
        # nt.show('output/nx.html')

def reaction_paths_unique_test(path_graph, g, synthesis_paths, i):
    graph_found = False
    if len(list(path_graph.nodes)) < len(list(g.nodes)):
        if isomorphism.DiGraphMatcher(G1=g, G2=path_graph, node_match=None, edge_match=None).subgraph_is_isomorphic():
            graph_found = True
            return [synthesis_paths, graph_found]
    elif len(list(path_graph.nodes)) > len(list(g.nodes)):
        if isomorphism.DiGraphMatcher(G1=path_graph, G2=g, node_match=None, edge_match=None).subgraph_is_isomorphic():
            graph_found = True
            synthesis_paths[i] = path_graph
            return [synthesis_paths, graph_found]
    else:
        if isomorphism.DiGraphMatcher(G1=path_graph, G2=g, node_match=None, edge_match=None).is_isomorphic():
            graph_found = True
            return [synthesis_paths, graph_found]
    return[synthesis_paths, graph_found]

def find_synthesis_paths(G: nx.DiGraph, amino_acid: str):
    queue.__init__()
    # Initialize a list to store all paths that synthesize the amino acid
    synthesis_paths = []
    # Mark all nodes as not available
    for node in G.nodes:
        G.nodes[node]['available'] = False
    # Add AS
    G.nodes[amino_acid]['available'] = True
    # Add cofactors
    for cofactor in constants.cofactors:
        try:
            G.nodes[cofactor]['available'] = True
        except KeyError:
            print(f'{cofactor} is not in network')
    # First queue input. Empty digraph and AS
    queue.put([amino_acid, nx.DiGraph()])

    while not queue.empty():
        current_node, current_path = queue.get()
        G.nodes[current_node]['available'] = True
        current_path.add_node(current_node)
        # Check if all required preds are available for the current node
        preds_available = all(G.nodes[pred]['available'] for pred in G.predecessors(current_node))
        # If all preds are available, add the current path to the synthesis_paths list
        if preds_available:
            graph_found = False
            for pred in G.predecessors(current_node):
                current_path.add_node(pred)
                current_path.add_edge(pred, current_node)
            path_graph = G.subgraph(current_path).copy()
            for i, g in enumerate(synthesis_paths):
                synthesis_paths, graph_found = reaction_paths_unique_test(path_graph, g, synthesis_paths, i)
                if graph_found:
                    break
            # If the path_graph is not already in the synthesis_paths list, add it
            if not graph_found:
                synthesis_paths.append(path_graph)
        else: # If not all preds are available, continue the search by enqueueing the preds of the current node
            for pred in G.predecessors(current_node):
                new_path = current_path.copy()
                new_path.add_node(pred)
                new_path.add_edge(pred, current_node)
                if G.nodes[pred]['available']: # or pred in current_path 
                    continue
                queue.put([pred, new_path])
    
    # Return the list of all paths that synthesize the amino acid
    ret_list  = []
    for ret in G.nodes:
        if G.nodes[ret]['available'] == True:
            ret_list.append(ret)
    return [G.subgraph(ret_list).copy(), synthesis_paths]


import queue
import networkx as nx
# from networkx.algorithms import isomorphism
import parse_smiles
import traversal
import constants
# import util
from pyvis.network import Network
from queue import SimpleQueue

organism_name = "ecoli_cimIV"

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
        # print("Metabolite info full:")
        # print_AS_cofactor_info(G)
        print("Metabolite info AS-bio-pw:")
        print_AS_cofactor_info(R)
        # for amino_acid in constants.amino_acid_list:
        #     print(amino_acid)    
        #     T, paths = find_synthesis_paths(R, amino_acid)

        # Test
        # paths = find_synthesis_paths(R, 'L-leucine')
        paths = find_synthesis_paths(R, 'L-leucine', max_length=1000)
        # len(paths)
        #
        # to_remove = []
        # for i, path in enumerate(paths):
        #     if not proof_path(R, paths[i], 'L-leucine'):
        #         to_remove.append(i)
        # len(to_remove)
        # 
        # for i in range(0, 10):
        #     print(set(list(paths[i].nodes)))
        #
        # # Visualization
        # nt = Network('1000px', '1000px', select_menu=True, filter_menu=True)
        # nt.from_nx(paths[10])
        # nt.show('test.html')


# For a found path proof if the path is legal by
# walking again through the nodes of the main path G
def proof_path(G:nx.DiGraph, path:nx.DiGraph, amino_acid:str):
    queue = SimpleQueue()
    queue.__init__()
    for node in path.nodes():
        path.nodes[node]['preds_avail'] = False
        path.nodes[node]['visited'] = False
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
    queue.put(amino_acid)
    while not queue.empty():
        current_node = queue.get()
        # print(current_node)
        # print(queue.qsize())
        # print(queue.empty())
        path.nodes[current_node]['visited'] = True
        # Check if all required preds are available for the current node
        preds_available = all(G.nodes[pred]['available'] for pred in G.predecessors(current_node))
        # print(f'preds avail: {preds_available}')
        # If all preds are available, add the current path to the synthesis_paths list
        if preds_available:
            path.nodes[current_node]['preds_avail'] = True
        else:
            # if there are no preds
            no_preds = len(list(path.predecessors(current_node))) == 0 
            # print(f'no preds: {no_preds}')
            if not no_preds:
                # print(list(path.predecessors(current_node)))
                for pred in path.predecessors(current_node):
                    visited = path.nodes[pred]['visited']
                    # print(f'visited: {visited}')
                    if not visited:
                        # print(f'put pred in queue: {pred}')
                        queue.put(pred)
    proof_list = []
    for node in path.nodes():
        proof_list.append(path.nodes[node]['preds_avail'])
    if all(proof_list):
        # print('exit: True')
        return True
    else:
        # print('exit: False')
        return False   

# def find_synthesis_paths(G, amino_acid, max_length=None):
#     """
#     Find all possible paths from which an amino acid can be synthesized
#     given a directed networkx graph with reactions and metabolites as nodes.
#
#     Args:
#     - G: networkx.DiGraph, the directed graph
#     - amino_acid: str, the amino acid to be synthesized
#     - max_length: int, optional, the maximum length of a path
#
#     Returns:
#     - synthesis_paths: list of lists, each inner list contains the nodes along a
#     path that synthesizes the amino acid
#     """
#     # Mark all nodes as not available
#     for node in G.nodes:
#         G.nodes[node]['available'] = False
#     # Add D-glucose
#     G.nodes['D-glucose']['available'] = True
#     # Add cofactors
#     for cofactor in constants.cofactors:
#         try:
#             G.nodes[cofactor]['available'] = True
#         except KeyError:
#             print(f'{cofactor} is not in network')
#
#     # initialize synthesis_paths and seen_nodes
#     synthesis_paths = []
#     seen_nodes = set()
#
#     # define helper function to recursively search for synthesis paths
#     def search_path(node, current_path):
#         # add node to current path
#         current_path.append(node)
#         print(f'node : {node}')
#
#         # if the current node is not the amino_acid and not available, return
#         # if node != amino_acid and not G.nodes[node]['available']:
#         #     return
#
#         # if the current node has already been seen, return to avoid infinite loops
#         if node in seen_nodes:
#             return
#
#         # mark the current node as seen
#         seen_nodes.add(node)
#
#         # if the current node is the amino acid, add the current path to synthesis_paths
#         if node == amino_acid:
#             synthesis_paths.append(current_path[:])
#             return
#
#         # recursively search for synthesis paths through the predecessors of the current node
#         for predecessor in G.predecessors(node):
#             # check if the maximum path length has been exceeded
#             if max_length is not None and len(current_path) >= max_length:
#                 continue
#             # # if the predecessor is a reversible reaction, check both directions
#             # if G.nodes[predecessor].get('reversible', False):
#             #     search_path(predecessor, current_path)
#             #     reverse_reaction = predecessor + '_reverse'
#             #     if reverse_reaction in G:
#             #         search_path(reverse_reaction, current_path)
#             # # if the predecessor is not a reversible reaction, search in the forward direction
#             # else:
#             search_path(predecessor, current_path)
#
#         # remove the current node from the current path before returning
#         current_path.pop()
#
#     # find all possible synthesis paths by searching from the amino acid node
#     for predecessor in G.predecessors(amino_acid):
#         search_path(predecessor, [])
#
#     return synthesis_paths


def find_synthesis_paths(G: nx.DiGraph, amino_acid: str):
    queue = SimpleQueue()
    queue.__init__()
    # Initialize a list to store all paths that synthesize the amino acid
    synthesis_paths = []
    # Mark all nodes as not available
    for node in G.nodes:
        G.nodes[node]['available'] = False
    # Add AS
    # G.nodes[amino_acid]['available'] = True
    # Add D-glucose
    G.nodes['D-glucose']['available'] = True
    # Add cofactors
    for cofactor in constants.cofactors:
        try:
            G.nodes[cofactor]['available'] = True
        except KeyError:
            print(f'{cofactor} is not in network')
    # First queue input. Empty digraph and AS
    queue.put([amino_acid, [], G.copy()])

    while not queue.empty():
        current_node, current_path, main_graph = queue.get()
        # print(f'current_node: {current_node}')
        main_graph.nodes[current_node]['available'] = True
        current_path.append(current_node)
        # Check if all required preds are available for the current node
        preds_available = all(G.nodes[pred]['available'] for pred in G.predecessors(current_node))
        # If all preds are available, add the current path to the synthesis_paths list
        if preds_available:
            # print("pred avaliable")
            graph_found = False
            for pred in main_graph.predecessors(current_node):
                current_path.append(pred)
            path_graph = main_graph.subgraph(current_path).copy()
            gfound_nodes = set(path_graph.nodes) 
            # print(f'set found nodes: {set(path_graph.nodes)}')
            for i, g in enumerate(synthesis_paths):
                g_nodes = set(g.nodes) 
                # len_gfound = len(gfound_nodes)
                # len_g = len(g_nodes)
                # print(f'set g nodes: {set(g.nodes)}')
                # print(len(f'found: {len_gfound}'))
                # print(len(f'g: {len_g}'))
                # ture if g is within path_graph
                if gfound_nodes < g_nodes and gfound_nodes != g_nodes:
                    # print( f'found < g: len found:{len_gfound}, len g: {len_g}' )
                    graph_found = True
                    break
                elif g_nodes < gfound_nodes and g_nodes != gfound_nodes:
                    # print(f'found > g: len found:{len_gfound}, len g: {len_g}' )
                    # print(proof_path(G, path_graph, amino_acid))
                    synthesis_paths[i] = path_graph
                    graph_found = True
                    break
                elif gfound_nodes == g_nodes:
                    # print(f'found == g: len found:{len_gfound}, len g: {len_g}' )
                    graph_found = True
                    break
            # If the path_graph is not already in the synthesis_paths list, add it
            if not graph_found:
                # print("graph not found")
                # print(proof_path(G, path_graph, amino_acid))
                synthesis_paths.append(path_graph)
        else: # If not all preds are available, continue the search by enqueueing the preds of the current node
            # print(list(G.predecessors(current_node)))
            for pred in main_graph.predecessors(current_node):
                new_path = current_path.copy()
                new_path.add_node(pred)
                new_path.add_edge(pred, current_node)
                if main_graph.nodes[pred]['available'] or pred in current_path.nodes(): # or pred in current_path 
                    continue
                queue.put([pred, new_path, main_graph.copy()])
    # Return the list of all paths that synthesize the amino acid
    return [synthesis_paths]


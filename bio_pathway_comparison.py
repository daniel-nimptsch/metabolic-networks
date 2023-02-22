import networkx as nx
import parse_smiles
import traversal
import constants
# import util

for organism_name in constants.path_dict:
    G = parse_smiles(organism_name, add_leucine=True)
    print(f'{organism_name} full:', len(list(G.nodes)))
    H = nx.DiGraph()
    P = traversal.get_subgraph_glucose_S(G)
    for amino_acid in constants.amino_acid_list:
        H = nx.compose(H, traversal.breadth_first_search_reverse_amino_acid(G, amino_acid))
    R = P.copy()
    R.remove_nodes_from(n for n in P if n not in H)
    R.remove_edges_from(e for e in P.edges if e not in H.edges)
    print(f'{organism_name} AS-bio-pw:', len(list(R.nodes)))
    list_of_all_products = traversal.breadth_first_search(R, "D-glucose")
    for amino_acid in constants.amino_acid_list:
        if amino_acid in list_of_all_products:
            print(amino_acid + " is there")
        else:
            print(amino_acid + " is NOT there")


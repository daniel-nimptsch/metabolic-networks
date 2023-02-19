import networkx as nx
import traversal
import constants
import util

def main():
    G = parse_smiles("ecoli_cimIV")

    # # --- WP1 -> 3. breadth first search (all amino acids need to be found) -------------
    # # --- to find L-Leucine the data needs to be updated (see moodle)
    # list_of_all_products = traversal.breadth_first_search(G, "D-glucose")
    # print(list_of_all_products)
    # for amino_acid in constants.amino_acid_list:
    #     if amino_acid in list_of_all_products:
    #         print(amino_acid + " is there")
    #     else:
    #         print(amino_acid + " is NOT there")

    # # --- WP1 -> 4. subgraph based on Glucose --------------------------------------------
    # print(traversal.get_subgraph_glucose_S(G))

    # # --- WP1 -> 5. reverse traversal ---------------------------------------------------- 
    # Subgraphs include nearly every node - seems wrong
    # for amino_acid in constants.amino_acid_list:
    #     H = traversal.breadth_first_search_reverse_amino_acid(G, amino_acid)
    #     print(amino_acid + ":")
    #     print(H)
    #     print("---")

    # # --- WP1 -> 6. Create the full Amino Acid Biosynthesis Pathway ----------------------
    H = nx.DiGraph()
    P = traversal.get_subgraph_glucose_S(G)
    for amino_acid in constants.amino_acid_list:
        H = nx.compose(H, traversal.breadth_first_search_reverse_amino_acid(G, amino_acid))
    H = nx.intersection(H, P)

    # # --- WP1 -> different reconstructions based on cultivation media? (adam/cimIV)
    # impossible paths!
    for path in nx.all_shortest_paths(G, "D-glucose", "L-arginine"):
        util.path_to_string(G, path)
        input()


def parse_smiles(organism: str):
    #TODO parse every file
    entry = []
    with open(str(constants.path_dict[organism]), "r") as file:
        entry = file.read().split("\n\n")

    G = nx.DiGraph()

    counter = 0
    for i in entry:
        try:
            # Line of reaction
            react = i.split('\n')[2]
            # First line
            first_line = i.split('\n')[0]

            # get big id from first line
            try:r_attr_BiggID = first_line.split(" ")[2]
            except IndexError:
                pass
            # get metanetxid from first line
            try:r_attr_MetaNetXId = first_line.split(" ")[4]
            except IndexError:
                pass
            # get reversible from first line
            try:r_attr_Reversible = first_line.split(" ")[6]
            except IndexError:
                pass
            # get smiles line
            r_attr_Smiles = i.split('\n')[3]

            react = i.split('\n')[2]
            left, right =react.split("=")
            counter += 1

            reaction = "r" + str(counter)
        
            # create the ceaction node
            G.add_node(reaction)
            G.nodes[reaction]['BiggId'] = str(r_attr_BiggID)
            G.nodes[reaction]['MetaNetXId'] = str(r_attr_MetaNetXId)
            G.nodes[reaction]['Reversible'] = str(r_attr_Reversible)
            G.nodes[reaction]['Smiles'] = str(r_attr_Smiles)
            G.nodes[reaction]['available'] = False
            G.nodes[reaction]['reaction'] = True

            for educt in left.split(' + '):
                G.add_edge(educt.strip(), reaction)
                G.add_node(educt.strip(), available=False, reaction=False)
            for product in right.split(' + '):
                G.add_edge(reaction, product.strip())
                G.add_node(product.strip(), available=False, reaction=False)

            # add all reactions that are reversible
            if G.nodes[reaction]['Reversible'] == 'True':
                reaction = "r" + str(counter) + "reverse"
                G.add_node(reaction)
                G.nodes[reaction]['BiggId'] = str(r_attr_BiggID)
                G.nodes[reaction]['MetaNetXId'] = str(r_attr_MetaNetXId)
                G.nodes[reaction]['Reversible'] = str(r_attr_Reversible)
                G.nodes[reaction]['Smiles'] = str(r_attr_Smiles)
                G.nodes[reaction]['available'] = False
                G.nodes[reaction]['reaction'] = True

                for educt in left.split(' + '):
                    G.add_edge(reaction, educt.strip())
                for product in right.split(' + '):
                    G.add_edge(product.strip(), reaction)

        except ValueError:
            pass #TODO ignore None
        
        #add availiable edcuts
        for educt in constants.cofactors:
            G.add_node(educt, available=True, reaction=False)

    return G


if __name__ == "__main__":
    main()
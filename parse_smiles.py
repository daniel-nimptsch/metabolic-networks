import networkx as nx
import traversal
import constants
import util

def main():
    G = parse_smiles("ecoli_cimIV", add_leucine=True)

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
    # H = nx.intersection(H, P)
    R = P.copy()
    R.remove_nodes_from(n for n in P if n not in H)
    R.remove_edges_from(e for e in P.edges if e not in H.edges)

    # # --- WP1 -> different reconstructions based on cultivation media? (adam/cimIV)
    # impossible paths!
    # for path in nx.all_shortest_paths(G, "D-glucose", "L-arginine"):
    #     util.path_to_string(G, path)
    #     input()
    return R

def parse_smiles(organism: str, add_leucine: bool=False):
    #TODO parse every file
    entry = []
    with open(str(constants.path_dict[organism]), "r") as file:
        entry = file.read().split("\n\n")

    if add_leucine:
        with open(str(constants.path_dict['leucine']), "r") as file:
            entry.append(file.read())

    Biggs = []
    G = nx.DiGraph()
    for i in entry:
        i = i.lstrip() # some lines have \n at the begining
        i = i.rstrip()
        try:
            first_line = i.split('\n')[0] # First line
            try:r_attr_BiggID = str(first_line.split(" ")[2]) # get big id from first line
            except IndexError:
                print("Error no BiggId")
            Biggs.append(r_attr_BiggID)  
            try:r_attr_MetaNetXId = str(first_line.split(" ")[4]) # get metanetxid from first line
            except IndexError:
                print("Error no MetaNetXId")
                pass
            try:r_attr_Reversible = str(first_line.split(" ")[6]) # get reversible from first line
            except IndexError:
                print("Error no Reversible")
            try:r_attr_Smiles = str(i.split('\n')[3]) # get smiles line
            except IndexError:
                print("Error no Smiles")
                pass

            react = i.split('\n')[2] # reaction line
            educts, products = react.split(" = ") # split in educts and products

            educts = educts.strip().split(" + ") 
            educts_counts = {}  # determine the counts of educts if they occure twice
            for educt in educts:
                if not educts_counts:
                    educts_counts[educt] = 1
                else:
                    for key in list(educts_counts.keys()):
                        if educt == key:
                            educts_counts[educt] += 1
                        else:
                            educts_counts[educt] = 1
            
            products = products.strip().split(" + ")
            products_counts = {}
            for product in products:
                if not products_counts:
                    products_counts[product] = 1
                else:
                    for key in list(products_counts.keys()):
                        if product == key:
                            products_counts[product] += 1
                        else:
                            products_counts[product] = 1

            reaction = r_attr_BiggID
        
            # create the ceaction node
            G.add_node(reaction)
            G.nodes[reaction]['BiggId'] = r_attr_BiggID
            G.nodes[reaction]['MetaNetXId'] = r_attr_MetaNetXId
            G.nodes[reaction]['Reversible'] = r_attr_Reversible
            G.nodes[reaction]['Smiles'] = r_attr_Smiles
            G.nodes[reaction]['available'] = False
            G.nodes[reaction]['reaction'] = True
            G.nodes[reaction]['educts'] = educts_counts # add dict of educt counts to reaction node. important for flux 
            G.nodes[reaction]['products'] = products_counts # same for products

            for educt in educts_counts:
                G.add_edge(educt, reaction)
                G.add_node(educt, available=False, reaction=False)
            for product in products_counts:
                G.add_edge(reaction, product)
                G.add_node(product, available=False, reaction=False)

            # add all reactions that are reversible
            if G.nodes[reaction]['Reversible'] == 'True':
                reaction = r_attr_BiggID + "_reverse"
                G.add_node(reaction)
                G.nodes[reaction]['BiggId'] = r_attr_BiggID
                G.nodes[reaction]['MetaNetXId'] = r_attr_MetaNetXId
                G.nodes[reaction]['Reversible'] = r_attr_Reversible
                G.nodes[reaction]['Smiles'] = r_attr_Smiles
                G.nodes[reaction]['available'] = False
                G.nodes[reaction]['reaction'] = True
                G.nodes[reaction]['educts'] = products_counts
                G.nodes[reaction]['products'] = educts_counts

                for educt in educts_counts:
                    G.add_edge(reaction, educt)
                for product in products_counts:
                    G.add_edge(product, reaction)

        except ValueError:
            pass #TODO ignore None
        
        #add availiable edcuts
        for educt in constants.cofactors:
            G.add_node(educt, available=True, reaction=False)

    if (len(Biggs) != len(set(Biggs))):
        print("Warning: duplicated BiggIds!")

    return G

if __name__ == "__main__":
    main()

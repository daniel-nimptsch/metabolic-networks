import networkx as nx
import constants
from queue import SimpleQueue

queue = SimpleQueue()

def breadth_first_search(G: nx.DiGraph, metabolite: str, show_reaction_nodes: bool=True):
    """This function works as described in the course description in moodle"""
    queue.__init__()
    G.add_node(metabolite, available=True) # add the given educt
    available_nodes = [x for x,y in G.nodes(data=True) if y['available']==True] # find all available nodes
    for educt in available_nodes:
        for reaction_to_check in G.neighbors(educt): # find all reactions that are possible
            all_available = True
            for req_educt in G.predecessors(reaction_to_check): # check if every educt is available
                if not G.nodes[req_educt]['available'] == True:
                    all_available = False
            if all_available:
                queue.put(reaction_to_check) # add reaction to queue

    while not queue.empty():
        current_reaction = queue.get()
        if(G.nodes[current_reaction]['available'] == True): # skip already added reactions
            continue
        new_products = []
        for new_educt in G.neighbors(current_reaction): # add educts of successfull reactions
            G.add_node(new_educt, available=True)
            G.add_node(current_reaction, available=True) # reactions that are "available" don't need to be visited again
            new_products.append(new_educt)

        for educt in new_products: 
            for new_reaction in G.neighbors(educt): # add possible reactions
                all_available = True
                for req_educt in G.predecessors(new_reaction): # check if all educts are available
                    if G.nodes[req_educt]['available'] == False:
                        all_available = False
                if all_available:
                    if G.nodes[new_reaction]['available'] == False: #educts of reaction were already added
                        queue.put(new_reaction)

    ret_list  = []
    for ret in G.nodes:
        if show_reaction_nodes:
            if G.nodes[ret]['available'] == True:
                ret_list.append(ret)
        else:
            if G.nodes[ret]['available'] == True and G.nodes[ret]['reaction'] == False:
                ret_list.append(ret)
    return ret_list

def get_subgraph_glucose_S(G: nx.DiGraph):
    H = G.subgraph(breadth_first_search(G, "D-glucose"))
    return H

def breadth_first_search_reverse_amino_acid(G: nx.DiGraph, amino_acid: str="all"): #TODO work in progress
    queue.__init__()
    # Mark all nodes as not available
    for node in G.nodes:
        G.nodes[node]['available'] = False

    if(amino_acid == "all"):
        for amino_acid in constants.amino_acid_list:
            G.nodes[amino_acid]['available'] = True
    else:
        G.nodes[amino_acid]['available'] = True
    available_nodes = [x for x,y in G.nodes(data=True) if y['available']==True] # find all available nodes
    # Add cofactors
    for cofactor in constants.cofactors:
        G.nodes[cofactor]['available'] = False

    for product in available_nodes:
        for reaction_to_check in G.predecessors(product): # find all reactions that are needed
            queue.put(reaction_to_check) # add reaction to queue

    while not queue.empty():
        current_reaction = queue.get()
        if(G.nodes[current_reaction]['available'] == True): # skip already added reactions
            continue
        for educt_needed in G.predecessors(current_reaction): # add educts of successfull reactions
            if educt_needed not in constants.cofactors:
                G.nodes[educt_needed]['available'] = True
                G.nodes[current_reaction]['available'] = True
                for reaction_needed in G.predecessors(educt_needed):
                    queue.put(reaction_needed)

    ret_list  = []
    for ret in G.nodes:
        if G.nodes[ret]['available'] == True:
            ret_list.append(ret)
    return G.subgraph(ret_list).copy()

def amino_acid_reverse_alternative(G: nx.DiGraph, amino_acid: str="all"):
    for cofactor in constants.cofactors:
        G.add_node(cofactor, available=False)
    if(amino_acid == "all"):
        for amino_acid in constants.amino_acid_list:
            G.add_node(amino_acid, available=True)
    H = nx.DiGraph()
    try:
        for path in nx.all_shortest_paths(G, "D-glucose", amino_acid):
            print(path)
            nx.add_path(H, path)
    except nx.NetworkXNoPath:
        pass #amino_acid cannot be reached
    return H

import networkx as nx

def path_to_string(G: nx.DiGraph, path: list):
    str_print = ""
    for i in path:
        if G.nodes[i]['reaction'] == True:
            for educt in G.predecessors(i):
                str_print += educt + " + "
            str_print = str_print[:-3] + " --> "
            for product in G.neighbors(i):
                str_print += product + " + "
            str_print = str_print[:-3] +  "\n"
    print(str_print)
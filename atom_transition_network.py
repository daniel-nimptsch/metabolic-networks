import enum
from queue import SimpleQueue
import networkx as nx
import pyvis
import constants
import random
from pyvis.network import Network

@enum.unique
class TransitionType(enum.IntEnum):
    """Possible SMILES token types"""
    NO_TRANSITION = 0
    SYMMETRY = 1
    REACTION = 2
    HYDROGEN_GROUP = 3
    HYDROGEN_REACTION = 4    
    HYDROGEN_FREE = 5

def main():
    ATN_blongum = nx.read_gml("../cleaned_networks/blongum_adam_cleaned.gml")
    ATN_btheta = nx.read_gml("../cleaned_networks/btheta_adam_cleaned.gml")
    ATN_ecoli = nx.read_gml("../cleaned_networks/ecoli_adam_cleaned.gml")
    ATN_eramosum = nx.read_gml("../cleaned_networks/eramosum_adam_cleaned.gml")
    
    ATN = delete_NO_TRANSITION_edges(ATN_blongum)
    print("blongum_adam:")
    print_statistics(ATN)
    
    # ATN = delete_NO_TRANSITION_edges(ATN_btheta)
    # print("btheta_adam:")
    # print_statistics(ATN)
    
    # ATN = delete_NO_TRANSITION_edges(ATN_ecoli)
    # print("ecoli_adam:")
    # print_statistics(ATN)
    
    # ATN = delete_NO_TRANSITION_edges(ATN_eramosum)
    # print("eramosum_adam:")
    # print_statistics(ATN)
    
    source = "D-glucose"
    target = ['L-arginine']

    # count_paths(delete_NO_TRANSITION_edges(ATN_blongum),source,"C",target[0],"C") # not working (too many paths)

    atn_name = ATN_eramosum
    print("ATN: betheta")
    print("all")
    all = BFS_glucose_C(atn_name, "all")
    print("noCO2")
    noCO2 = BFS_glucose_C(remove_molecule(atn_name, "CO2"), "noCO2")
    print("noAMP")
    noAMP = BFS_glucose_C(remove_molecule(atn_name, "AMP"), "noAMP")


    compound_names, path= find_path_compound_names(ATN_ecoli, source, target[0])
    
    ATN = delete_others(ATN_ecoli, set(compound_names))
    
    # draw_Graph(ATN, source, target, compound_names, "a", path=path)

    
    # ATN = nx.read_gml("../cleaned_networks/ecoli_adam_cleaned.gml")
    # ATN = remove_CO2(ATN)
    
    # compound_names = []
    # for protein in constants.amino_acid_list:
    #     compound_names += find_path_compound_names(ATN, source, protein)
    
    
    # ATN = delete_others(ATN, compound_names)

    # draw_Graph(ATN, source, target, compound_names, "b", True)

def draw_Graph(ATN: nx.Graph, source: str=None, target: list=[], compound_names: list=[], name: str="nx", path=[]):
    draw = nx.Graph()
    draw.add_nodes_from(ATN.nodes())
    draw.add_edges_from(ATN.edges())

    for n in ATN.nodes():
        draw.nodes[n]['title'] = ATN.nodes[n]['compound_name']
        draw.nodes[n]['label'] = ATN.nodes[n]['element']
        if 'hcount' in ATN.nodes[n]:
            if ATN.nodes[n]['hcount'] != 0:
                draw.nodes[n]['label'] += str(ATN.nodes[n]['hcount'])+'H'
    if len(compound_names) != 0:
        color = 0
        color_target = 0
        for compound in compound_names:
            if compound in target:
                for n in ATN.nodes():
                    if ATN.nodes[n]['compound_name'] == compound:
                        draw.nodes[n]['color'] = constants.COLOR_LIST_TARGET[color_target]
                color_target = (color_target + 1) % len(constants.COLOR_LIST_TARGET)
            elif compound == source:
                for n in ATN.nodes():
                    if ATN.nodes[n]['compound_name'] == compound:
                        draw.nodes[n]['color'] = "lime"
            else:
                for n in ATN.nodes():
                    if ATN.nodes[n]['compound_name'] == compound:
                        draw.nodes[n]['color'] = constants.COLOR_LIST_COMPOUNDS[color]
                color = (color + 1) % len(constants.COLOR_LIST_COMPOUNDS)

    for e in ATN.edges():
        if ATN.edges[e]['transition'] == 'TransitionType.SYMMETRY':
            draw.edges[e]['color'] = "green"
        elif ATN.edges[e]['transition'] == 'TransitionType.REACTION':
            draw.edges[e]['color'] = "red"
        elif ATN.edges[e]['transition'] == 'TransitionType.HYDROGEN_GROUP' or ATN.edges[e]['transition'] == 'TransitionType.HYDROGEN_REACTION':
            draw.edges[e]['color'] = "LightSkyBlue"
        elif ATN.edges[e]['transition'] == 'TransitionType.HYDROGEN_FREE':
            draw.edges[e]['color'] = "DarkBlue"
        else:
            # draw.edges[e]['label'] = str(ATN.edges[e]['order'])
            draw.edges[e]['color'] = "black"
            draw.edges[e]['value'] = 70


    # for n in ATN.nodes():
    #     draw.nodes[n]['size'] = nx.degree(ATN,n)
    if len(path) > 0:
        del_list = []
        for e in draw.edges():
            if ATN.edges[e]["transition"] != "TransitionType.NO_TRANSITION":
                del_list.append(e)
        for e in del_list:
            a,b = e
            draw.remove_edge(a,b)
        lol = 0
        a = None
        b = None
        for e in path:
            if lol == 0:
                a = e
            else:
                b = e
                draw.add_edge(a,b,color="magenta")
                a = b
            lol += 1

    nt = Network('750px', '750px')
    nt.show_buttons()
    nt.from_nx(draw)
    nt.show(name + '.html')

def remove_molecule(ATN: nx.Graph, molecule: str):
    S = ATN.copy()
    del_list = []
    for n in S.nodes():
        if S.nodes[n]['compound_name'] == molecule:
            del_list.append(n)
            
    for del_node in del_list:
        S.remove_node(del_node)
    return S

def print_statistics(ATN: nx.Graph):
    print(ATN)
    print("Number of connected components: " + str(nx.number_connected_components(ATN)))
    print("Density: " + str(nx.density(ATN)))
    S = [ATN.subgraph(c).copy() for c in nx.connected_components(ATN)]
    # component_counter = 1
    # for subgraph in S:
    #     print("  Component " + str(component_counter))
    #     print("    " + str(subgraph))
    #     print("    Density: " + str(nx.density(subgraph)))
    #     component_counter += 1

def find_path_compound_names(ATN: nx.Graph, source: str, target: str):
    glucose = None
    protein = None
    for n in ATN.nodes:
        if ATN.nodes[n]['compound_name'] == source and ATN.nodes[n]['element'] == "C":
            glucose = n
        if ATN.nodes[n]['compound_name'] == target and ATN.nodes[n]['element'] == "C":
            protein = n

    compound_names = []
    S = delete_NO_TRANSITION_edges(ATN)
    paths = nx.all_simple_paths(S, glucose, protein, cutoff=10)
    path = None
    counter = 0
    rand = random.randint(1,100)
    for p in paths:
        path = p
        if(counter == rand):
            break
        counter += 1
    for p in path:
        compound_names.append(ATN.nodes[p]['compound_name'])
    compounds = []
    for n in ATN.nodes():
        if ATN.nodes[n]['compound_name'] in compound_names:
            compounds.append(n)
    return [compound_names, path]

def delete_others(ATN: nx.Graph, compound_names: list):
    del_list = []
    for n in ATN.nodes():
        if not ATN.nodes[n]['compound_name'] in compound_names:
            del_list.append(n)
    
    for del_node in del_list:
        ATN.remove_node(del_node)
    return ATN

def delete_NO_TRANSITION_edges(ATN: nx.Graph):
    S = ATN.copy()
    del_list = []
    for e in S.edges():
        if S.edges[e]['transition'] == 'TransitionType.NO_TRANSITION':
            del_list.append(e)
    for del_item in del_list:
        a, b = del_item
        try:
            S.remove_edge(a,b)
        except nx.NetworkXError:
            pass # edge not in Graph
    return S

def count_paths(ATN: nx.Graph, source: str, source_atom: str, target: str, target_atom: str):
    start = None
    end = None
    for n in ATN.nodes:
        if ATN.nodes[n]['compound_name'] == source and ATN.nodes[n]['element'] == source_atom:
            start = n
            break
    for n in ATN.nodes:
        if ATN.nodes[n]['compound_name'] == target and ATN.nodes[n]['element'] == target_atom:
            end = n
            break
    paths = nx.all_simple_paths(ATN, start, end, cutoff=50)
    print("Paths between " + source + "(" + source_atom + ") and " + target + "(" + target_atom + ")")
    counter = 0 
    for i in paths:
        counter += 1
    print("Amount: " + str(counter))

def BFS_glucose_C(ATN: nx.Graph, name:str):
    S = nx.Graph()
    G = delete_NO_TRANSITION_edges(ATN)
    queue = SimpleQueue()
    queue.__init__()
    visited = []
    for n in ATN.nodes:
        if ATN.nodes[n]['compound_name'] == "D-glucose" and ATN.nodes[n]['element'] == "C":
            queue.put(n)
    while not queue.empty():
        next_node = queue.get()
        for connected_node in G.neighbors(next_node):
            if not connected_node in visited:
                S.add_node(next_node, compound_name=ATN.nodes[next_node]['compound_name'])
                S.add_node(next_node, element=ATN.nodes[next_node]['element'])
                S.add_node(connected_node, compound_name=ATN.nodes[connected_node]['compound_name'])
                S.add_node(connected_node, element=ATN.nodes[connected_node]['element'])
                S.add_edge(next_node,connected_node,transition="TransitionType.REACTION")
                visited.append(connected_node)
                queue.put(connected_node)
    compound_names = []
    for i in visited:
        compound_names.append(ATN.nodes[i]['compound_name'])
    compound_names = set(compound_names)
    print("Molecule reached: " + str(len(compound_names)-1))
    counter = 0
    for i in constants.amino_acid_list:
        if i in compound_names:
            counter += 1
    print("Amino Acids reached: " + str(counter) + "/20")
    show_whole_molecule = ["D-glucose", "CO2", "AMP", "pyruvate", "L-glutamate", "2-oxoglutarate", "(6S)-5,6,7,8-tetrahydrofolate", "D-glucose 6-phosphate", "prephenate"]
    high_degree = []
    for compound in compound_names:
        degree_counter = 0
        for n in S.nodes():
            if S.nodes[n]['compound_name'] == compound:
                degree_counter += nx.degree(S, n)
        if degree_counter > 35 or compound in show_whole_molecule:
            high_degree.append((compound, degree_counter))
    print("Molecule with highest degree:")
    print(high_degree)
    print("Density: " + str(nx.density(S)))
    print(nx.degree_histogram(S))
    for e in ATN.edges():
        a, b = e
        for compound in show_whole_molecule:
            if ATN.nodes[a]['compound_name'] == compound and ATN.nodes[b]['compound_name'] == compound:
                S.add_node(a, compound_name=ATN.nodes[a]['compound_name'])
                S.add_node(a, element=ATN.nodes[a]['element'])
                S.add_node(b, compound_name=ATN.nodes[b]['compound_name'])
                S.add_node(b, element=ATN.nodes[b]['element'])
                S.add_edge(a,b,transition="TransitionType.NO_TRANSITION")
    draw_Graph(S,"D-glucose",compound_names=compound_names, name=name)
    return compound_names


main()
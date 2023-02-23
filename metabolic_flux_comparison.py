import parse_smiles
import flux_analysis
import constants
import fastaparser as fp
import os
import pandas as pd
import networkx as nx
from pyvis.network import Network
import time

def analyze():
    for organism_name in constants.path_dict:
        try:
            analyze_organism(organism_name)
        except:
            pass


def visualize_organisms():
    for organism_name in constants.path_dict:
        with open(constants.proteom_dict[organism_name]) as fasta_file:
            parser = fp.Reader(fasta_file)
            counter = 0
            for seq in parser:
                # seq is a FastaSequence object
                protein = seq.sequence_as_string()
                protein_name = seq.id
                counter += 1
                if counter > 0:
                    break
            try:
                G = parse_smiles.AS_bio_pathway_for_organism(organism_name)
                G = flux_analysis.flux_analysis(G.copy(), protein).copy()
                H = filter_low_flux(G.copy()).copy()
                H = add_size_attr_graph(H.copy()).copy()
                # visualization
                outdir = f'output/{organism_name}/'
                outpath = f"visualization_{organism_name}_{protein_name}.html"
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                outpath = os.path.join(outdir, outpath)
                nt = Network('1000px', '1000px', select_menu=True, filter_menu=True)
                nt.from_nx(H)
                nt.save_graph(f"visualization_{organism_name}_{protein_name}.html")
            except:
                pass

def analyze_organism(organism_name):
    # organism_name = "ecoli_cimIV"
    G = parse_smiles.AS_bio_pathway_for_organism(organism_name)
    outdir = f'output/{organism_name}/'
    # Create a dictionary of the information
    graph_info = {
        'Number_of_nodes': 'Number_of_nodes',
        'Number_of_edges': 'Number_of_edges',
        'Density': 'Density',
        'Protein': 'Protein',
        'Protein_description': 'Protein_description',
    }
    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(graph_info, index=[0])
    outpath = f'flux_graph_info_{organism_name}.csv'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = os.path.join(outdir, outpath)
    df.to_csv(outpath, index=False, header=False)
    # Create a dictionary of the information
    node_info = {
        'Name':'Name',
        'Flux':'Flux',
        'Reaction':'Reaction',
        'Protein':'Protein',
    }
    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(node_info, index=[0])
    outpath = f'flux_node_info_{organism_name}.csv'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = os.path.join(outdir, outpath)
    df.to_csv(outpath, index=False, header=False)
    with open(constants.proteom_dict[organism_name]) as fasta_file:
        parser = fp.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            protein = seq.sequence_as_string()
            protein_name = seq.id
            protein_description = seq.description
            analyze_metabolix_flux_prot(G, outdir, protein, protein_name, organism_name, protein_description)

def add_size_attr_graph(H: nx.DiGraph):
    flux_values = nx.get_node_attributes(H, 'flux').values()
    min_flux = min(flux_values)
    max_flux = max(flux_values)
    nodes = [x for x in H.nodes()]
    for node in nodes:
        flux = H.nodes[node]['flux']
        size = (((flux - min_flux)/(max_flux - min_flux)) * 100) + 10   # adjust size based on flux
        H.nodes[node]['size'] = size
    return(H)

def filter_low_flux(G: nx.DiGraph):
    # Collect all metabolites and reaction node names
    metabolites = [x for x,y in G.nodes(data=True) if y['reaction']==False]
    reactions = [x for x,y in G.nodes(data=True) if y['reaction']==True]
    # Exclusion_range = range(-0.1, 0.1)
    def check_inclusion_range(value):
        if -0.01 <= value <= 0.01:
            return False
        else:
            return True
    nodes_to_keep= []
    for metabolite in metabolites:
        flux = G.nodes[metabolite]['flux']
        if 'input' not in metabolite:
            if check_inclusion_range(flux):
                nodes_to_keep.append(metabolite)
    for reaction in reactions:
        flux = G.nodes[reaction]['flux']
        # if 'input' not in reaction:
        if check_inclusion_range(flux):
            nodes_to_keep.append(reaction)
    nodes_to_keep = list(set(nodes_to_keep))
    H = nx.subgraph(G, nodes_to_keep)
    H = nx.Graph(H) # unfreeze
    isolates = list(nx.isolates(H))
    H.remove_nodes_from(isolates)
    return(H)

def analyze_metabolix_flux_prot(G: nx.DiGraph, outdir:str, protein: str, protein_name: str, organism_name: str, protein_description: str="no_description"):
    G = flux_analysis.flux_analysis(G, protein)
    H = filter_low_flux(G.copy()).copy()
    # Get graph information
    num_nodes = len(H.nodes)
    num_edges = len(H.edges)
    density = nx.density(H)
    # Create a dictionary of the information
    graph_info = {
        'Number_of_nodes': num_nodes,
        'Number_of_edges': num_edges,
        'Density': density,
        'Protein': protein_name,
        'Protein_description': protein_description,
    }
    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(graph_info, index=[0])
    # Save the DataFrame to a CSV file
    outpath = f'flux_graph_info_{organism_name}.csv'
    outpath = os.path.join(outdir, outpath)
    df.to_csv(outpath, mode='a', index=False, header=False)
    # get node attributes as a dictionary
    protein_list = [protein_name] * len(list(H.nodes))
    node_info = {
        'Name': list(H.nodes),
        'Flux': list(nx.get_node_attributes(H, 'flux').values()),
        'Reaction': list(nx.get_node_attributes(H, 'reaction').values()),
        'Protein': list(protein_list),
    }
    # convert to pandas DataFrame
    df = pd.DataFrame(node_info)
    # save to CSV file
    outpath = f'flux_node_info_{organism_name}.csv'
    outpath = os.path.join(outdir, outpath)
    df.to_csv(outpath, mode='a', index=False, header=False)


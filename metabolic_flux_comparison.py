import parse_smiles
import flux_analysis
import constants
import fastaparser as fp
import os
import pandas as pd
import networkx as nx
# from pyvis.network import Network

def analyze_organism(organism_name):
    G = parse_smiles.AS_bio_pathway_for_organism(organism_name)
   
    with open(constants.proteom_dict[organism_name]) as fasta_file:
        parser = fp.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            # print('ID:', seq.id)
            # print('Description:', seq.description)
            # print('Sequence:', seq.sequence_as_string())
            # print()
            protein = seq.sequence_as_string()
            protein_name = seq.id
            protein_description = seq.description
            # analyze_metabolix_flux_prot(G, protein, protein_name, organism_name)

def analyze_metabolix_flux_prot(G, protein, protein_name, organism_name):

    G, S = flux_analysis.flux_analysis(G, protein)

    # Collect all metabolites and reaction node names
    metabolites = [x for x,y in G.nodes(data=True) if y['reaction']==False]
    reactions = [x for x,y in G.nodes(data=True) if y['reaction']==True]

    # Exclusion_range = range(-0.1, 0.1)
    def check_inclusion_range(value):
        if -0.1 <= value <= 0.1:
            return False
        else:
            return True

    nodes_to_keep= []
    for metabolite in metabolites:
        nodes_to_keep.append(metabolite)

    for reaction in reactions:
        flux = G.nodes[reaction]['flux']
        if 'input' not in reaction:
            if check_inclusion_range(flux):
                nodes_to_keep.append(reaction)

    nodes_to_keep = list(set(nodes_to_keep))
    H = nx.subgraph(G, nodes_to_keep)
    H = nx.Graph(H) # unfreeze

    isolates = list(nx.isolates(H))
    H.remove_nodes_from(isolates)

    nodes = [x for x in H.nodes()]
    for node in nodes:
        flux = H.nodes[node]['flux']
        size = int(abs(flux * 0.05) + 10)  # adjust size based on flux
        H.nodes[node]['size'] = size

    # # visualization
    # nt = Network('1000px', '1000px', select_menu=True, filter_menu=True)
    # nt.from_nx(H)
    # nt.show('nx.html')

    # Get graph information
    num_nodes = len(H.nodes)
    num_edges = len(H.edges)
    density = nx.density(H)

    # Create a dictionary of the information
    graph_info = {
        'Number of nodes': num_nodes,
        'Number of edges': num_edges,
        'Density': density,
        'Protein': protein_name
    }

    # Convert the dictionary to a pandas DataFrame
    df = pd.DataFrame(graph_info, index=[0])

    # Save the DataFrame to a CSV file
    outdir = f'output/{organism_name}/'
    outpath = f'flux_graph_info_{organism_name}.csv'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = os.path.join(outdir, outpath)
    df.to_csv(outpath, mode='a', index=False, header=False)

    # # get node attributes as a dictionary
    # node_attrs = nx.get_node_attributes(G, 'flux')
    # node_attrs.update(nx.get_node_attributes(G, 'reaction'))
    #
    # # convert to pandas DataFrame
    # df = pd.DataFrame.from_dict(node_attrs, orient='index', columns=['flux', 'reaction'])
    #
    # # save to CSV file
    # df.to_csv(f'output/{organism_name}_{protein_name}.csv')


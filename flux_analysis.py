from pulp import LpMaximize, LpProblem, LpStatus, lpSum, LpVariable
import pandas as pd
import networkx as nx
import numpy as np
import constants
import parse_smiles

def main():
    G = parse_smiles.main()
    protein = "MIVQKELVAIYDYEVPVPEDPFSFRLEIHKCSELFTGSVYRLERFRLRPTFHQRDREDADPLINDALIYIRDECIDERKLRGESPETVIAIFNRELQNIFNQEIE"
    flux_analysis(G, protein)

def flux_analysis(G, protein): 

    # Amino acid dict with counts for protein
    as_counts = {} # create empty dict for amino acid (AS/as) names and counts
    for amino_acid in protein:
        as_current = constants.amino_acid_dict.get(amino_acid)  # as_current will be the AS name corresponding to letter
        as_current = as_current.get("name")
        if as_current in as_counts:
            as_counts[as_current] += 1   # add counts to as_counts
        else:
            as_counts[as_current] = 1

    # Amino acid dict with fractions for protein
    total_counts = sum(as_counts.values())   # determine total counts
    fractions = [count / total_counts for count in as_counts.values()]   # determine fractions from counts and the total count
    as_fractions = as_counts.copy()
    for amino_acid in as_fractions: # replace the counts in the dict with fractions
        index = list(as_fractions.keys()).index(amino_acid)
        as_fractions[amino_acid] = fractions[index]

    # Add biomass node biomass reaction node and edges to the graph 
    as_set_biomass = as_counts.keys()   # AS names required for the protein
    R_biomass = "R_biomass"  # reaction name for node
    biomass = "biomass" # biomass name for node
    G.add_node(R_biomass, reaction=True, group =1) # add AS edges to r_biomass node and from it an edge to biomass node
    for amino_acid in as_set_biomass:
        G.add_edge(amino_acid, R_biomass)
    G.nodes[R_biomass]['educts'] = as_counts
    G.nodes[R_biomass]['products'] = {biomass: 1} 
    G.add_node(biomass, reaction=False, group =0)
    G.add_edge(R_biomass, biomass)

    # Add glucose input reaction and energy input reactions
    G.add_node("input_node", reaction=False, group =0)
    input_reaction = "R_input_D-glucose"
    G.add_node(input_reaction, reaction=True, group =1)
    G.add_edge("input_node", input_reaction)
    G.add_edge(input_reaction, "D-glucose")
    G.nodes[input_reaction]['educts'] = {"input_node": 1}
    G.nodes[input_reaction]['products'] = {"D-glucose": 1}

    for cofactor in constants.cofactors:
        input_reaction = "R_input_" + cofactor
        G.add_node(input_reaction, reaction=True, group =1)
        G.add_edge("input_node", input_reaction)
        G.add_edge(input_reaction, cofactor)
        G.nodes[input_reaction]['educts'] = {"input_node": 1}
        G.nodes[input_reaction]['products'] = {str(cofactor): 1}

    # Collect all metabolites and reaction node names
    metabolites = [x for x,y in G.nodes(data=True) if y['reaction']==False]
    reactions = [x for x,y in G.nodes(data=True) if y['reaction']==True]

    # Construct the stochiometric matrix
    S = np.zeros((len(metabolites), len(reactions)))
    for j, reaction in enumerate(reactions):    # j: reaction index
        educts = G.nodes[reaction]['educts']
        products = G.nodes[reaction]['products']
        for i, metabolite in enumerate(metabolites):    # i: metaboliet index
            if metabolite in educts:
                S[i, j] = -educts[metabolite]
            elif metabolite in products:
                S[i, j] = products[metabolite]

    # Upper and lower bounds
    lb = {}
    ub = {}
    for r in reactions:
        if "input_D-glucose" in r:
            lb.update({r: -10})
            ub.update({r: 1000})
        elif "input_" in r:
            lb.update({r: -1000})
            ub.update({r: 1000})
        else:
            lb.update({r: 0})
            ub.update({r: 1000})

    # Define the LP problem
    model = LpProblem('maximize_protein_production', LpMaximize)

    # Define the decision variables: reactions with lower and upper bounds
    reaction_variables = {}
    for r in reactions:
        reaction_variables[r] = LpVariable(r, lowBound=lb[r], upBound=ub[r])

    # target metabolite: protein
    target_product = "R_biomass"
    # Define the objective function
    objective_variable = reaction_variables[target_product]
    model += objective_variable

    # Define the constraints based on the stoichiometric matrix
    # Idea: for each metabolite (rows in stochiometric matrix) define a constraint 
    # that has the rule that the sum of each reaction flux multiplied by the stochiometry should be greater than 0
    # eg -3 H2O*flux + 2 H2O*flux + 4 H2O*flux >= 0
    for i in range(S.shape[0]):
        flux = lpSum(S[i, j] * reaction_variables[reactions[j]] for j in range(S.shape[1]))
        model += flux >= 0

    # Print model
    # print(model)

    # Solve the LP problem
    model.solve()
    # status = model.solve()
    # print(status)

    # Add flux attribute to reactions nodes
    for j, reaction in enumerate(reactions):
        flux = reaction_variables[reaction].varValue 
        G.nodes[reaction]['flux'] = flux   
    
    # # Convert the stoichiometric matrix to pd df 
    # S_df = pd.DataFrame(S, index=metabolites, columns=reactions)
    #
    # # Create a dictionary to store the flux values for reaction nodes
    # flux_dict = nx.get_node_attributes(G, 'flux')
    #
    # # Convert the flux dictionary to a pandas series
    # flux_series = pd.Series(flux_dict)
    #
    # # Select the flux values for reaction nodes
    # reaction_flux = flux_series[flux_series.index.isin(reactions)]
    #
    # # Calculate the sum of the absolute values of the flux for each metabolite
    # abs_flux_sum = S_df.abs().dot(reaction_flux.abs())
    #
    # # Calculate the relative flux for each metabolite
    # relative_flux = (S_df * reaction_flux).sum(axis=1) / abs_flux_sum
    # relative_flux = relative_flux.fillna(0)

    for metabolite in metabolites:
        G.nodes[metabolite]['flux'] = 0
        # G.nodes[metabolite]['flux'] = relative_flux[metabolite]

    # return the graph and the stochiometric matrix
    return[G, S]
    

from pulp import LpMaximize, LpProblem, LpStatus, lpSum, LpVariable
import networkx
import numpy as np
import constants
import parse_smiles
from pyvis.network import Network

G = parse_smiles.main()
# TODO make a definition/function for this script that takes a graph and a protein sequence as inputs
# and that returns the results as outputs
# Test protein
# >ABD18687.1 ecoli.faa
protein = "MIVQKELVAIYDYEVPVPEDPFSFRLEIHKCSELFTGSVYRLERFRLRPTFHQRDREDADPLINDALIYIRDECIDERKLRGESPETVIAIFNRELQNIFNQEIE"

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
G.add_node(R_biomass, reaction=True) # add AS edges to r_biomass node and from it an edge to biomass node
for amino_acid in as_set_biomass:
    G.add_edge(amino_acid, R_biomass)
G.nodes[R_biomass]['educts'] = as_counts
G.nodes[R_biomass]['products'] = {biomass: 1} 
G.add_node(biomass, reaction=False)
G.add_edge(R_biomass, biomass)

# Add glucose input reaction and energy input reactions
G.add_node("input_node", reaction=False)
input_reaction = "R_input_D-glucose"
G.add_node(input_reaction, reaction=True)
G.add_edge("input_node", input_reaction)
G.add_edge(input_reaction, "D-glucose")
G.nodes[input_reaction]['educts'] = {"input_node": 1}
G.nodes[input_reaction]['products'] = {"D-glucose": 1}

for cofactor in constants.cofactors:
    input_reaction = "R_input_" + cofactor
    # print(input_reaction)
    G.add_node(input_reaction, reaction=True)
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
    # print("---")
    for i, metabolite in enumerate(metabolites):    # i: metaboliet index
        if metabolite in educts:
            S[i, j] = -educts[metabolite]
            # print(-educts[metabolite])
        elif metabolite in products:
            S[i, j] = products[metabolite]
            # print(products[metabolite])

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
for i,j in enumerate(reactions):
    if j == target_product:
        target_product_index = i

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
print(model)

# Solve the LP problem
status = model.solve()
print(status)

# Print the optimal reaction_variables
for reaction in reactions:
    if reaction_variables[reaction].varValue > 0:
        print(f"{reaction}: {reaction_variables[reaction].varValue}")
   
H = G.copy()
stoichiometric_flux = np.copy(S)
for i, metabolite in enumerate(metabolites):
     for j, reaction in enumerate(reactions):
        flux = reaction_variables[reaction].varValue 
        x = stoichiometric_flux[i,j] * flux
        stoichiometric_flux[i, j] = x
        if x <= 0:
            H.remove_nodes_from(metabolite)
        else:
            H.nodes[metabolite]['stoichiometric_flux'] = x

for j, reaction in enumerate(reactions):
    flux = reaction_variables[reaction].varValue 
    if flux <= 0:
        H.remove_nodes_from(reaction)
    else:
        H.nodes[reaction]['flux'] = flux   

# print(S.shape[0])
# print(S.shape[1])
# print(len(metabolites))
# print(len(reactions))

len(list(G.nodes))
len(list(H.nodes))
H.nodes.data("stoichiometric_flux")
H.nodes.data("flux")

# visualization
nt = Network('500px', '500px')
nt.from_nx(H)
nt.show('nx.html')


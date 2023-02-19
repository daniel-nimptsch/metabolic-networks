import pulp
import networkx
import numpy as np
import constants
import parse_smiles

G = parse_smiles.main()

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
    print(input_reaction)
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
    print("---")
    for i, metabolite in enumerate(metabolites):    # i: metaboliet index
        if metabolite in educts:
            S[i, j] = -educts[metabolite]
            print(-educts[metabolite])
        elif metabolite in products:
            S[i, j] = products[metabolite]
            print(products[metabolite])

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
prob = LpProblem('maximize_protein_production', LpMaximize)

# Define the decision variables
fluxes = {}
for r in reactions:
    fluxes[r] = pulp.LpVariable(r, lowBound=lb[r], upBound=ub[r])

# Define the objective function
obj_func = pulp.lpSum(fluxes[reaction] for reaction in reactions)
prob += obj_func

# Define the constraints based on the stoichiometric matrix
for i in range(S.shape[0]):
    flux = pulp.lpSum(S[i, j] * fluxes[reactions[j]] for j in range(S.shape[1]))
    prob += flux == 0

# Solve the LP problem
prob.solve()

# Print the optimal fluxes
for reaction in reactions:
    print(f"{reaction}: {fluxes[reaction].varValue}")

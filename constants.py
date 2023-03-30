# LIST OF ALL PROTEINOGENIC AMINO ACIDS

amino_acid_list = ["L-arginine", "L-valine", "L-methionine", "L-glutamate", "L-glutamine", "L-tyrosine", "L-tryptophan",
                   "L-proline", "L-cysteine", "L-histidine", "L-asparagine", "L-aspartate", "L-phenylalanine",
                   "L-threonine", "L-lysine", "L-serine", "L-isoleucine", "glycine", "L-alanine", "L-leucine"]

# DICTIONARY MAPPING AMINO ACID ABBREVIATION

amino_acid_dict = {
    "R": {"name": "L-arginine"},
    "V": {"name": "L-valine"},
    "M": {"name": "L-methionine"},
    "E": {"name": "L-glutamate"},
    "Q": {"name": "L-glutamine"},
    "Y": {"name": "L-tyrosine"},
    "W": {"name": "L-tryptophan"},
    "P": {"name": "L-proline"},
    "C": {"name": "L-cysteine"},
    "H": {"name": "L-histidine"},
    "N": {"name": "L-asparagine"},
    "D": {"name": "L-aspartate"},
    "F": {"name": "L-phenylalanine"},
    "T": {"name": "L-threonine"},
    "K": {"name": "L-lysine"},
    "S": {"name": "L-serine"},
    "I": {"name": "L-isoleucine"},
    "G": {"name": "glycine"},
    "A": {"name": "L-alanine"},
    "L": {"name": "L-leucine"}
}

# LIST OF COFACTORS

# cofactors = ["AMP", "ADP", "ATP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "CTP", "CoA", "H2O", "NH4(+)", "hydrogen sulfide"]
# cofactors = ["AMP", "ATP", "ADP", "GDP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "FAD", "FADH2", "UTP", "CTP", "heme b", "CoA", "FMN", "H2O", "NH4(+)", "phosphate", "CO2", "hydrogen sulfide"]
cofactors = ["AMP", "ATP", "ADP", "GDP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "FAD", "FADH2", "UTP", "CTP", "heme b", "CoA", "FMN", "H2O", "NH4(+)", "phosphate", "hydrogen sulfide"]
all_cofactors = ["AMP", "ATP", "ADP", "GDP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "FAD", "FADH2", "UTP", "CTP", "heme b", "CoA", "FMN", "H2O", "NH4(+)", "phosphate", "CO2", "hydrogen sulfide"]

# LIST TO SMILES_LIST

path_dict = {
    "acacae_adam": "../sihumix/acacae_adam/acacae_adam.smiles_list",
    "acacae_cimIV": "../sihumix/acacae_cimIV/acacae_cimIV.smiles_list",
    "blongum_adam": "../sihumix/blongum_adam/blongum_adam.smiles_list",
    "blongum_cimIV": "../sihumix/blongum_cimIV/blongum_cimIV.smiles_list",
    "bproducta_adam": "../sihumix/bproducta_adam/bproducta_adam.smiles_list",
    "bproducta_cimIV": "../sihumix/bproducta_cimIV/bproducta_cimIV.smiles_list",
    "btheta_adam": "../sihumix/btheta_adam/btheta_adam.smiles_list",
    "btheta_cimIV": "../sihumix/btheta_cimIV/btheta_cimIV.smiles_list",
    "cbuty_adam": "../sihumix/cbuty_adam/cbuty_adam.smiles_list",
    "cbuty_cimIV": "../sihumix/cbuty_cimIV/cbuty_cimIV.smiles_list",
    "ecoli_adam": "../sihumix/ecoli_adam/ecoli_adam.smiles_list",
    "ecoli_cimIV": "../sihumix/ecoli_cimIV/ecoli_cimIV.smiles_list",
    "eramosum_adam": "../sihumix/eramosum_adam/eramosum_adam.smiles_list",
    "eramosum_cimIV": "../sihumix/eramosum_cimIV/eramosum_cimIV.smiles_list",
    "lplantarum_adam": "../sihumix/lplantarum_adam/lplantarum_adam.smiles_list",
    "lplantarum_cimIV": "../sihumix/lplantarum_cimIV/lplantarum_cimIV.smiles_list",
}

leucine_path = { 
    "leucine": "../L_leucine.smiles_list"
}

proteom_dict = {
    "acacae_adam": "./data/Proteom/Anaerostipes_caccae.faa",
    "acacae_cimIV": "./data/Proteom/Anaerostipes_caccae.faa",
    "blongum_adam": "./data/Proteom/Bifidobacterium_longum.faa",
    "blongum_cimIV": "./data/Proteom/Bifidobacterium_longum.faa",
    "bproducta_adam": "./data/Proteom/Blautia_producta.faa",
    "bproducta_cimIV": "./data/Proteom/Blautia_producta.faa",
    "btheta_adam": "./data/Proteom/Bacteroides_thetaiotaomicron.faa",
    "btheta_cimIV": "./data/Proteom/Bacteroides_thetaiotaomicron.faa",
    "cbuty_adam": "./data/Proteom/Clostridium_butyricum.faa",
    "cbuty_cimIV": "./data/Proteom/Clostridium_butyricum.faa",
    "ecoli_adam": "./data/Proteom/Escherichia_coli.faa",
    "ecoli_cimIV": "./data/Proteom/Escherichia_coli.faa",
    "eramosum_adam": "./data/Proteom/Erysipelatoclostridium_ramosum.faa",
    "eramosum_cimIV": "./data/Proteom/Erysipelatoclostridium_ramosum.faa",
    "lplantarum_adam": "./data/Proteom/Lactobacillus_plantarum.faa",
    "lplantarum_cimIV": "./data/Proteom/Lactobacillus_plantarum.faa",
}

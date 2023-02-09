from typing import Literal, List, TypedDict


AminoAcid = Literal[
    "L-arginine", 
    "L-valine", 
    "L-methionine", 
    "L-glutamate", 
    "L-glutamine", 
    "L-tyrosine", 
    "L-tryptophan",
    "L-proline", 
    "L-cysteine", 
    "L-histidine", 
    "L-asparagine", 
    "L-aspartate", 
    "L-phenylalanine",
    "L-threonine", 
    "L-lysine", 
    "L-serine", 
    "L-isoleucine", 
    "glycine", 
    "L-alanine", 
    "L-Leucine"
]

SpeciesMediumCombination = Literal[
    "acacae_adam",
    "acacae_cimIV",
    "blongum_adam",
    "blongum_cimIV",
    "bproducta_adam",
    "bproducta_cimIV",
    "btheta_adam",
    "btheta_cimIV",
    "cbuty_adam",
    "cbuty_cimIV",
    "ecoli_adam",
    "ecoli_cimIV",
    "eramosum_adam",
    "eramosum_cimIV",
    "lplantarum_adam",
    "lplantarum_cimIV"
]

Protein = Literal[None]

class Reaction(TypedDict):
    bigg_id: str
    metanetx_id: str
    reversible: bool
    educts: List[str]
    products: List[str]
    smiles_educts: List[str]
    smiles_products: List[str]
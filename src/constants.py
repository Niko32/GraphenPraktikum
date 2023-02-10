from typing import get_args, List

from custom_types import AminoAcid, Protein, SpeciesMediumCombination


AMINO_ACIDS: List[AminoAcid] = get_args(AminoAcid)
PROTEINS: List[Protein] = get_args(Protein)
SEPCIES_MEDIUM_COMBINATIONS: List[SpeciesMediumCombination] = get_args(SpeciesMediumCombination)
COFACTORS = ["AMP", "ADP", "ATP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "CTP", "CoA", "H2O", "NH4(+)", "hydrogen sulfide"] \
    + ["GDP", "FAD", "FADH2", "UTP", "heme b", "FMN", "phosphate", "CO2"]
AMINO_ACID_DICT = {
    "R": "L-arginine",
    "V": "L-valine",
    "M": "L-methionine",
    "E": "L-glutamate",
    "Q": "L-glutamine",
    "Y": "L-tyrosine",
    "W": "L-tryptophan",
    "P": "L-proline",
    "C": "L-cysteine",
    "H": "L-histidine",
    "N": "L-asparagine",
    "D": "L-aspartate",
    "F": "L-phenylalanine",
    "T": "L-threonine",
    "K": "L-lysine",
    "S": "L-serine",
    "I": "L-isoleucine",
    "G": "glycine",
    "A": "L-alanine",
    "L": "L-leucine"
}
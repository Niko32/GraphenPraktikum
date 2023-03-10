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

SPECIES_DICT = {
    "acacae": "Anaerostipes_caccae",
    "blongum": "Bifidobacterium_longum",
    "bproducta": "Blautia_producta",
    "btheta": "Bacteroides_thetaiotaomicron",
    "cbuty": "Clostridium_butyricum",
    "ecoli": "Escherichia_coli",
    "eramosum": "Erysipelatoclostridium_ramosum",
    "lplantarum": "Lactobacillus_plantarum"
}

LILA = {0: "#F0AFDF", 1: "#D937B0", 2: "#AC25C7", 3: "#7B37D9"}
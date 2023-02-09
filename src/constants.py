from typing import get_args

from custom_types import AminoAcid, Protein, SpeciesMediumCombination


AMINO_ACIDS = get_args(AminoAcid)
PROTEINS = get_args(Protein)
SEPCIES_MEDIUM_COMBINATIONS = get_args(SpeciesMediumCombination)
COFACTORS = ["AMP", "ADP", "ATP", "NAD(+)", "NADH", "NADP(+)", "NADPH", "CTP", "CoA", "H2O", "NH4(+)", "hydrogen sulfide"] \
    + ["GDP", "FAD", "FADH2", "UTP", "heme b", "FMN", "phosphate", "CO2"]
# IsotopeMappingToolKit

The IsotopeMappingToolKit enables the conversion of metabolic reactions summarized in SMBL format to an atomic transition networkbased on MetaNetX.

# Prerequisites

Download the MetaNetX database into the according subfolder.

The following python libraries need to be available:
* BeautifulSoup
* RXNMapper
* rdkit
* networkx
* pyvis (optional)

# How to Use

Run as

```bash

./01_bigg_to_smiles_reactions.py [SMBL Xml] [SMILES]
./02_atommap_smiles_reactions.py [SMILES] [Mapped SMILES]
./03_generate_ATN.py [Mapped SMILES] [ATN in GML format]

```


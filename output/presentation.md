# Presentation

## WP1/BFS/Network
- Which organism can produce which AAs?
- How large are 'subgraphs' (i.e. Glucose/AAs)?

## WP2/FBA
- Implementation specifies (e.g. Input/Export, Constraints, Protein Assembly, Exploration of 'settings')
- Flux distribution (few reactions with strong flow? many non-zero flows?)
- Did you play around with removing network part?
- Did the flow change over multiple runs? (any patterns between different runs?)

## WP3/ATNs
- How large, how many (weakly) connected components? (excluding 'NO-TRANSITION' edges)
- How many paths between 2 molecules? (e.g. Glucose C_1 -> pyruvate C_x) (What is 'between' 2 molecules, i.e. now do those paths relate to eachother)
- density/connectivity measures
    - for example Page rank etc. 
    - you can explore networkx algorithms for that and think about / evaluate usefulness
- Distribution of Glucose C-atoms with and without CO2/AMP present in the ATN
    - other source atoms than Glucose also interesting
- different endpoints for BFS starting at Glucose Carbon atoms?
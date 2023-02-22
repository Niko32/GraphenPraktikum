# Presentation

## WP1/BFS/Network
Q: Which organism can produce which AAs?
A: [('acacae_adam', 20), ('acacae_cimIV', 20), ('blongum_adam', 20), ('blongum_cimIV', 20), ('bproducta_adam', 19), ('bproducta_cimIV', 19), ('btheta_adam', 20), ('btheta_cimIV', 20), ('cbuty_adam', 20), ('cbuty_cimIV', 20), ('ecoli_adam', 20), ('ecoli_cimIV', 20), ('eramosum_adam', 20), ('eramosum_cimIV', 20), ('lplantarum_adam', 0), ('lplantarum_cimIV', 0)]
A: ![Bar Plot](output/plots/species_medium_aa_bar_plot.png)
A: ![Heatmap](output/plots/species_medium_aa_heatmap.png)
Q: How large are 'subgraphs' (i.e. Glucose/AAs)?
A: TODO: wp1_analyzer.plot_sizes (Nicola)

## WP2/FBA
Q: Implementation specifies (e.g. Input/Export, Constraints, Protein Assembly, Exploration of 'settings')
A: there are about 15% non zero flux values
Q: Flux distribution (few reactions with strong flow? many non-zero flows?)
A: no changes for removing ATP
Q: Did you play around with removing network part?
A: no changes, at least objective value does not change
Q: Did the flow change over multiple runs? (any patterns between different runs?)

## WP3/ATNs
Q: How large, how many (weakly) connected components? (excluding 'NO-TRANSITION' edges)
A: 1. barplot: x: species_medium_combination, y: no. of connected components
A: 2. barplot: x: species_medium_combination --> four bars each (min, max, average, mean), y: size
Q: How many paths between 2 molecules? (e.g. Glucose C_1 -> pyruvate C_x) (What is 'between' 2 molecules, i.e. now do those paths relate to eachother)
Q: density/connectivity measures
    - for example Page rank etc. 
    - you can explore networkx algorithms for that and think about / evaluate usefulness
Q: Distribution of Glucose C-atoms with and without CO2/AMP present in the ATN
    - other source atoms than Glucose also interesting
Q: different endpoints for BFS starting at Glucose Carbon atoms?
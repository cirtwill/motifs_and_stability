# motif_role_stability
Simulation project testing for links between species roles and their extinction risk


# Order of operations:

- Steps 1 and 4 take quite a long time and should only be run once.
- Step 7 is also slow.

1) generate initial networks with code/test_parameters.jl
2) check numbers of components and, if only one, collect path lengths and edge lists with code/check_components.R
3) calculate initial roles (and motif participation) with code/calculate_motifs_roles_triads/calculate_motifs_and_roles.py and code/calculate_motifs_roles_triads/calculate_motif_participation.py
4) read in initial webs and simulate removals of each species in turn with code/simulate_removals.jl
5) extract extinction orders from removals using code/extinction_orders_and_biomasses.R
 *This produces the extinction order and biomass files that are on Git*
6) compile extinction orders and initial roles into one file per network (all removals) with code/roles_to_extorder.py: gives the roles/matched_to_extorder/ files. Also compiles motif participation (counts and normalised), degrees, and TL's (shortest TL). 
7) test roles! Test mean extinction order vs. role with code/stat_analyses/roles_vs_extorder_removal.R *slow*
8) When all permutations are finished, collect results using code/stat_analyses/roles_vs_extorder_collector.py. *This produces the permanova files that are on Git*
9) Display permanova results using code/figure_creation/roles_vs_extorder_permanova_summary.py and make table with code/figure_creation/create_tables/automatically.py

[[Done to here, now brain is required]]
10) Test whether PCA positions are consistent across webs with code/stat_analyses/testing_role_consistency.R 
11) Get PCA axes of roles using code/stat_analyses/role_PCAs_forplots.R Plot using code/figure_creation/mean_loadings_PCAaxes.tsv 
12) Check correlations of extinction orders using code/stat_analyses/mean_correlation_extorder_tests.R Plot using code/figure_creation/extinction_order_correlations.py
13) Fit lmers for mean time to extinction vs. PCA positions and participation in key (one-way, except loop) motifs using code/stat_analyses/extorder_vs_PCApositions.R and code/stat_analyses/motif_participation_vs_extorder.R

 

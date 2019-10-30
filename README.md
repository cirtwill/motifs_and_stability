# motif_role_stability
Simulation project testing for links between species roles and their extinction risk


# Order of operations:

- Steps 1 and 4 take quite a long time and should only be run once.
- Step 7 is also slow.

1) generate initial networks with code/test_parameters.jl
1A) if getting initial networks from Git, run initial_network_extractor.jl to generate .json files with the initial networks.
2) check numbers of components and, if only one, collect path lengths and edge lists with code/check_components.R
3) calculate initial roles (and motif participation) with code/calculate_motifs_roles_triads/calculate_motifs_and_roles.py 
4) read in initial webs and simulate removals of each species in turn with code/simulate_removals.jl
5) extract extinction orders from removals using code/extinction_orders_and_biomasses.py
 *This produces the extinction order and biomass files that are on Git*
6) compile extinction orders and initial roles into one file per network (all removals) with code/roles_to_extorder.py: gives the roles/matched_to_extorder/ files. Also compiles motif participation (counts and normalised), degrees, and TL's (shortest TL).
7) test roles! Test mean extinction order vs. role with code/stat_analyses/roles_vs_extorder_removal.R
8) When all permutations are finished, collect results using code/stat_analyses/roles_vs_extorder_collector.py. *This produces the permanova files that are on Git*
9) Display permanova results using code/figure_creation/roles_vs_extorder_permanova_summary.py
10) Test whether extinction order is related to dissimilarity between focal and removed species. Doing this with a mantel and an lm because I'm unsure. code/stat_analyses/rolediff_vs_orders_Mantel.R and code/stat_analyses/rolediff_vs_orders_lmer.R
11) Fit lmers for mean time to extinction vs. PCA positions and participation in key (one-way, except loop) motifs using code/stat_analyses/extorder_vs_PCApositions.R and code/stat_analyses/motif_participation_vs_extorder.R

 

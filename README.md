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


- Fit partial least squares regressions of degree-normalized and network-normalized motifs (plus S, C, S:C, degree, STL) against persistence using code/stat_analysis/partial_least_squares_regression.R
- plot coefficients from above across all optimum axes using code/figure_generation/PLS_coefficients.py (produces manuscript/figures/PLS/total_coefficients.eps)


7) test roles! Test mean extinction order vs. role with code/stat_analyses/roles_vs_extorder_removal.R *slow* [[done]]
8) When all permutations are finished, collect results using code/stat_analyses/roles_vs_extorder_collector.py. *This produces the permanova files that are on Git* [[done]]
9) Display permanova results using code/figure_creation/roles_vs_extorder_permanova_summary.py and make table with code/figure_creation/create_tables_automatically.py [[done]]


- Test motifs (counts and proportions) vs. degree and STL using code/stat_analysis/motifs_vs_degree.R
- plot these using code/figure_generation/motif_correlations.py


12) Check correlations of extinction orders using code/stat_analyses/mean_correlation_extorder_tests.R Plot using code/figure_creation/extinction_order_correlations.py [[done both]]

13) Fit lmers for mean time to extinction vs. participation in motifs using [[likely unnecessary code/stat_analyses/extorder_vs_PCApositions.R]] and code/stat_analyses/motif_participation_vs_extorder.R [[done both]]


# Haven't necessarily updated all figures yet


# Figures:
create_tables_automatically.py : makes large table of PERMANOVA results. 
display_motifs.py : displays the motifs and descriptions. Fig. 1
extinction_order_correlations.py : shows similarity of extinction order across groups of nets. Fig. 2
full_motif_lm.py : shows coefficients and ranges of effects for lm of persistence ~ motifs. Fig. 5
plot_motif_correlations.py : shows correlations between all motifs and degree. Fig. 6
PLS_coefficients.py : shows aggregate coefficients of each motif. Fig. 3. Needs uncertainties.
roles_vs_extorder_permanova_summary : shows all of the PERMANOVA results, F and p-values. Fig. 4



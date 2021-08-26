## matchMS
#   calc_pepmasses2 (ms1_precursors.R)
#   bin_ms1masses (bin_masses.R)
#   load_mgfs (mgf.R)
#   ms2match (msmsmatches2.R)
#   calc_pepscores (scores.R)
#   calc_peploc
#   add_prot_acc (quant.R)
#   calc_protfdr (scores.R)
#   add_rptrs
#   try_psmC2Q
#     psmC2Q
#       grp_prots (quant2.R)
# 
# ======================================
# bin_masses.R
#   - bin_ms1masses
#    - binTheoPeps
#     - bin_theopeps
#       - find_ms1_cutpoints
#    - cbind_theopepes
# 
# mgf.R
#   - readMGF
#    - read_mgf_chunks
#     - proc_mgf_chunks_i
#      - proc_mgf_chunks
#       - proc_mgfs
#     - proc_mgfs
#    - find_ms1_interval
# 
# calc_tmtint (quant.R)
#   find_reporter_ints
# 
# 

# ms1_precursors.R: 
#   - calc_pepmasses2
#     - calc_aamasses
#       - add_fixvar_masses
#       - parse_aamasses
#     - split_fastaseqs
#     - distri_fpeps
#       - make_fastapeps0
#     - ms1masses_bare
#       - ms1masses_noterm
#         - calcms1mass_noterm
#           - calcms1mass_noterm_byprot
#             - calcms1mass_noterm_bypep
#     - distri_peps
#       - subpeps_by_vmods
#     - add_term_mass
#     - helpers below
# 
# helpers at sets of realized modifications: 
#   (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
#   - ms1_a0_fnl1_byprot
#       ms1_a0_fnl1_bypep
#         delta_ms1_a0_fnl1
#   
#  (7-8) "amods+ tmod- vnl- fnl-", "amods+ tmod+ vnl- fnl-"
#    (9-10) "amods+ tmod- vnl+ fnl-", "amods+ tmod+ vnl+ fnl-"
#    (11-12) "amods+ tmod- vnl- fnl+", "amods+ tmod+ vnl- fnl+"
#    (13-14) "amods+ tmod- vnl+ fnl+", "amods+ tmod+ vnl+ fnl+"
#   - ms1_a1_vnl0_fnl0_byprot
#       ms1_a1_vnl0_fnl0_bypep
#         [by nested combinatorial conditions; no explicit functions]
# 
# 
## ms2match (msmsmatches2.R)
# 
# ms2base.R: (1, 2) "amods- tmod+ vnl- fnl-", "amods- tmod- vnl- fnl-"
#   ms2match_base 
#     purge_search_space (utils_engine.R)
#     hms2_base (helper)
#       frames_adv_base (frame-advancing)
#         gen_ms2ions_base (for specific pep_seq)
#           ms2ions_by_type (ion_ladder.R)
#             byions, czions, axions
#         search_mgf2
#           find_ms2_bypep
#     post_ms2match (utils_engine.R)
# 
# ms2_a0_vnl0_fnl1.R: (5, 6) "amods- tmod+ vnl- fnl+", "amods- tmod- vnl- fnl+"
#   ms2match_a0_vnl0_fnl1 
#     purge_search_space
#     hms2_a0_vnl0_fnl1
#       frames_adv_a0_vnl0_fnl1
#         gen_ms2ions_a0_fnl1
#           ms2ions_by_type (ion_ladder.R)
#             byions, czions, axions
# 
# ms2_a1_vnl0_fnl0.R: (7, 8) "amods+ tmod+ vnl- fnl-", "amods+ tmod- vnl- fnl-"
#   ms2match_a1_vnl0_fnl0 
#     purge_search_space
#     hms2_a1_vnl0_fnl0
#       frames_adv_a1_vnl0_fnl0
#         gen_ms2ions_a1_vnl0_fnl0
#           combi_mvmods2
#             combi_vmods2
#           find_intercombi_p2
#           check_ms1_mass_vmods2
#           calc_ms2ions_a1_vnl0_fnl0
#             ms2ions_by_type
#               byions, czions, axions
#           add_hexcodes (sitecombi.R)
# 
# ms2_a1_vnl1_fnl0.R: (9, 10) "amods+ tmod+ vnl+ fnl-", "amods+ tmod- vnl+ fnl-"
#   ms2match_a1_vnl1_fnl0 
#     purge_search_space
#     hms2_a1_vnl1_fnl0
#       frames_adv_a1_vnl1_fnl0
#         gen_ms2ions_a1_vnl1_fnl0
#           combi_mvmods2
#             combi_vmods2
#           find_intercombi_p2
#           check_ms1_mass_vmods2 (ms2_a1_vnl0_fnl0.R)
#           calc_ms2ions_a1_vnl1_fnl0
#             ms2ions_by_type
#               byions, czions, axions
#           add_hexcodes_vnl2
# 
# ms2_a1_vnl0_fnl1.R: (11, 12) "amods+ tmod+ vnl- fnl+", "amods+ tmod- vnl- fnl+"
#   ms2match_a1_vnl0_fnl1 
#     purge_search_space
#     hms2_a1_vnl0_fnl1
#       frames_adv_a1_vnl0_fnl1
#         gen_ms2ions_a1_vnl0_fnl1
#           combi_mvmods2
#             combi_vmods2
#           find_intercombi_p2
#           check_ms1_mass_vmods2
#           calc_ms2ions_a1_vnl0_fnl1
#             ms2ions_by_type
#               byions, czions, axions
#           add_hexcodes_fnl2

## calc_pepscores (scores.R)
#   scalc_pepprobs
#     calc_probi
#       calc_probi_bypep
#         calc_probi_byvmods
#           add_seions
#           find_ppm_outer_bycombi (msmsmatches.R)
# 
## calc_protfdr (scores.R)
#   calc_protfdr_i
#   fit_protfdr

## grp_prots (quant2.R)
#   groupProts2
#     map_pepprot2
#     cut_protgrps2
#       as_lgldist
#     greedysetcover3
#     




## (Tentative) same-site rules: no additive varmods
# 
# No additive terminal mods (fixed/fixed; var/var; fixed/var)!!!
# 
# (1) No more than one fixedmod on the same residue or N/C terminal.
#   [N] fixed "Oxidation (M)" + fixed "Carbamidomethyl (M)"
#   [N] fixed "TMT6plex (N-term)" + fixed "Acetyl (Protein N-term)"
# 
#   # otherwise additive
#   [Y] fixed "Oxidation (M)" + fixed "TMT6plex (N-term)" with M on the N-term
# (2) OK different variable mods to the same residue at different sites
#   [Y] variable "Oxidation (M) @3" + variable "Carbamidomethyl (M) @4"
#   [Y] variable "Oxidation (M)" + variable "TMT6plex (N-term)" with M on the N-term
# (3) no conflict between fixedmods and variable mods (additive)
#   [Y] fixed "Oxidation (M)" + variable "TMT6plex (N-term)" with M on the N-term
#   [N] fixed "Oxidation (M) @3" + variable "Carbamidomethyl (M) @3"
#   [N] fixed "TMT6plex (N-term)" + variable "Acetyl (Protein N-term)"




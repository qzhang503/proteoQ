# $cluego_api.R
# [1] "cluego"
# 
# $corrplots.R
#  [1] "plotCorr"           "my_custom_cor"      "  percent_of_range" "plot_corr_sub"      "  my_fn"            "  lm_with_cor"     
#  [7] "  panel_cor"        "  panel_hist"       "  panel_lm"         "  my_lower"         "  my_lower_no_sm"   "  my_diag"         
# [13] "pepCorr_logFC"      "pepCorr_logInt"     "prnCorr_logFC"      "prnCorr_logInt"    
# 
# $datasets.R
# character(0)
# 
# $dbs.R
# [1] "read_fasta"  "write_fasta" "load_fasta"  "calc_avgpep"
# 
# $funs.R
# [1] "#   [2] \"find_mascot_psmraws <-function(filelist = NULL, dat_dir = NULL)\""
# 
# $global.R
# [1] "set_dat_dir"    "get_gl_dat_dir"
# 
# $go.R
#  [1] "create_db_path"        "proc_obo"              "proc_gaf"              "annot_from_to"         "get_full_entrez"      
#  [6] "find_human_orthologs"  "sp_lookup_go"          "find_abbr_species"     "set_db_outname"        "dl_msig"              
# [11] "proc_gmt"              "prepGO"                "prepMSig"              "map_to_entrez"         "Uni2Entrez"           
# [16] "Ref2Entrez"            "map_to_entrez_os_name" "create_os_lookup"      "  convert_default_os" 
# 
# $gsea.R
# [1] "make_gct" "make_cls" "prnGSEA"  "fml_gsea"
# 
# $gspa.R
#  [1] "pepGSPA"            "prnGSPA"            "gspaTest"           "fml_gspa"           "ok_min_size"        "gspa_summary_mean" 
#  [7] "gspa_summary_limma" "lm_gspa"            "prep_gspa"          "map_essential"      "prnGSPAHM"          "gspaHM"            
# [13] "byfml_gspahm"       "byfile_gspahm"      "gspa_colAnnot"      "greedysetcover"    
# 
# $gsva.R
# [1] "prnGSVA"  "gsvaTest"
# 
# $histo.R
# [1] "plotHisto" "pepHist"   "prnHist"  
# 
# $hm.R
# [1] "my_pheatmap" "plotHM"      "pepHM"       "prnHM"      
# 
# $informatics.R
# [1] "info_anal"     "find_pri_df"   "find_sec_df"   "vararg_secmsg"
# 
# $kin.R
# [1] "annot_KinSub" "KinSubTest"   "anal_KinSub" 
# 
# $lda.R
# [1] "plotLDA"  "scoreLDA" "pepLDA"   "prnLDA"  
# 
# $mds.R
#  [1] "geom_lower_text"   "plotMDS"           "plotEucDist"       "scoreMDS"          "prep_folded_tdata" "scoreEucDist"     
#  [7] "pepMDS"            "prnMDS"            "pepEucDist"        "prnEucDist"       
# 
# $metadata.R
#  [1] "prep_label_scheme"          "prep_fraction_scheme"       "update_frac_smry"           "make_frac_smry"            
#  [5] "check_frac_multipsms"       "check_exptfrac_raws"        "read_metadata"              "write_metadata"            
#  [9] "check_unnamed_cols"         "check_required_cols"        "check_optional_cols"        "check_empty_rows"          
# [13] "check_tmt126_row"           "check_tmt_nc"               "check_tmt_plex"             "check_channel_prefix"      
# [17] "check_tmt_ref"              "rm_fully_empty_tmt_sets"    "replace_empty_with_na"      "replace_na_with_empty"     
# [21] "check_metadata_integers"    "syn_empty_channels"         "triv_optcols_at_empty"      "recheck_tmt_channels"      
# [25] "check_tmt_rawfiles"         "check_tmt_sampleids"        "check_tmt_emptyids"         "check_tmt_chans_vs_levs"   
# [29] "check_tmt_mixplexes"        "check_dups_at_lcms_and_sid" "check_lfq_exptraws"         "load_dbs"                  
# [33] "load_expts"                 "reload_expts"               "channelInfo"                "n_TMT_sets"                
# [37] "n_LCMS"                     "TMT_plex"                   "TMT_plex2"                  "TMT_levels"                
# [41] "simple_label_scheme"        "find_mascot_tmtplex"       
# 
# $nmf.R
#  [1] "analNMF"         "plotNMFCon"      "plotNMFCoef"     "plotNMFmeta"     "anal_pepNMF"     "anal_prnNMF"     "plot_pepNMFCon" 
#  [8] "plot_prnNMFCon"  "plot_pepNMFCoef" "plot_prnNMFCoef" "plot_metaNMF"   
# 
# $normalization.R
#  [1] "find_fit_nms"     "hfind_fit_nms"    "calc_sd_fcts"     "center_df"        "add_mean_dev"     "find_n_comp"     
#  [7] "my_which_max"     "ok_file_ncomp"    "spline_coefs"     "normMulGau"       "dblTrim"          "sumdnorm"        
# [13] "  wt_dnorm"       "normSD"           "fitKernelDensity" "nmix_params"     
# 
# $outliers.R
#  [1] "locate_outliers" "dixon_test"      "q_dixon"         "p_dixon"         "q_table"         "Dixon_outliers"  "Grubbs_outliers"
#  [8] "grubbs_test"     "q_grubbs"        "          f"     "p_grubbs"        "Rosner_outliers"
# 
# $pca.R
# [1] "plotPCA"  "scorePCA" "pepPCA"   "prnPCA"  
# 
# $peptable.R
#  [1] "newColnames"           "hsubColnames"          "use_mq_peptable"       "use_mf_peptable"       "use_mf_prottable"     
#  [6] "use_mq_prottable"      "fillMFprnnums"         "fillMQprnnums"         "check_mq_df"           "calc_mq_log2r"        
# [11] "extract_mq_ints"       "pep_mq_lfq"            "pep_mf_lfq"            "calc_lfq_log2r"        "na_single_lfq"        
# [16] "calcLFQPepNums"        "impRefNA"              "find_species_centers"  "find_species_centers2" "point_in_polygon"     
# [21] "calcLFQPepInts"        "spreadPepNums"         "countSpecs"            "add_spec_counts"       "aggrNumLCMS"          
# [26] "pad_grp_samples"       "rep_ls_groups"         "normPep"               "find_ms1filepath"      "sep_mbrfiles"         
# [31] "extract_mbry"          "groupMZPepZ"           "groupMZPepZ2"          "med_summarise_keys"    "load_prior"           
# [36] "fmt_num_cols"          "mergePep"              "standPep"              "Pep2Prn"               "map_peps_prots"       
# [41] "find_prot_family_rows" "my_sum_n"              "calcLFQprnnums"        "calc_tmt_prnnums"      "pep_to_prn"           
# [46] "assign_duppeps"        "replace_by_rowone"     "impPepNA"              "find_mbr_ms1files"    
# 
# $prntable.R
# [1] "standPrn"
# 
# $proteoQ.R
# character(0)
# 
# $psmtable.R
#   [1] "extract_raws"                                                   
#   [2] "find_mascot_psmraws <-function(filelist = NULL, dat_dir = NULL)"
#   [3] "find_mq_psmraws"                                                
#   [4] "find_sm_psmraws"                                                
#   [5] "find_mf_psmraws"                                                
#   [6] "find_pq_psmraws"                                                
#   [7] "extract_psm_raws"                                               
#   [8] "set_mascot_colnms"                                              
#   [9] "batchPSMheader"                                                 
#  [10] "rmPSMHeaders"                                                   
#  [11] "add_mod_conf"                                                   
#  [12] "find_mascot_vmods"                                              
#  [13] "add_mascot_pepseqmod"                                           
#  [14] "add_empai"                                                      
#  [15] "add_quality_cols"                                               
#  [16] "psm_msplit"                                                     
#  [17] "find_shared_prots"                                              
#  [18] "hfind_shared_prots"                                             
#  [19] "add_shared_genes"                                               
#  [20] "add_shared_prot_accs_mf"                                        
#  [21] "add_shared_sm_genes"                                            
#  [22] "procPSMs"                                                       
#  [23] "check_raws"                                                     
#  [24] "splitPSM_ma"                                                    
#  [25] "compile_hdr_mods"                                               
#  [26] "find_mascot_plex"                                               
#  [27] "make_mascot_ratios"                                             
#  [28] "make_mascot_intensities"                                        
#  [29] "add_mascot_raw"                                                 
#  [30] "pad_mascot_channels"                                            
#  [31] "pad_mascot_fields"                                              
#  [32] "psm_mcleanup"                                                   
#  [33] "cleanupPSM"                                                     
#  [34] "mcPSM"                                                          
#  [35] "annotPSM"                                                       
#  [36] "normPSM"                                                        
#  [37] "rm_cols_mqpsm"                                                  
#  [38] "groupMZPSM1"                                                    
#  [39] "groupMZPSM0"                                                    
#  [40] "hrm_lowIntPSMs"                                                 
#  [41] "rm_lowIntPSMs"                                                  
#  [42] "rm_psm_spikes"                                                  
#  [43] "rm_apexrt_outliers"                                             
#  [44] "keep_psm_bestint"                                               
#  [45] "preMZpepLFQ"                                                    
#  [46] "postMZpepLFQ"                                                   
#  [47] "checkmzLFQRT"                                                   
#  [48] "calcLFQPeptide"                                                 
#  [49] "calcTMTLFQPeptide"                                              
#  [50] "psm_to_pep"                                                     
#  [51] "PSM2Pep"                                                        
#  [52] "add_pep_n_apexes"                                               
#  [53] "cleanRTChim"                                                    
#  [54] "alignPSMApexs"                                                  
#  [55] "rm_RToutliers"                                                  
#  [56] "makeLFQXYTS"                                                    
#  [57] "cleanPSMYTS"                                                    
#  [58] "rm_2RToutliers"                                                 
#  [59] "rm_nRToutliers"                                                 
#  [60] "updatePSMsd"                                                    
#  [61] "fitMS1RT"                                                       
#  [62] "hsub_by_bestz"                                                  
#  [63] "sub_by_bestz"                                                   
#  [64] "calibPSMApexs2"                                                 
#  [65] "my_tolower"                                                     
#  [66] "my_upper"                                                       
#  [67] "add_maxquant_pepseqmod"                                         
#  [68] "add_msfragger_pepseqmod"                                        
#  [69] "add_cols_at"                                                    
#  [70] "replace_cols_at"                                                
#  [71] "reloc_col"                                                      
#  [72] "reloc_col_after"                                                
#  [73] "reloc_col_after_last"                                           
#  [74] "reloc_col_after_first"                                          
#  [75] "reloc_col_before"                                               
#  [76] "reloc_col_before_last"                                          
#  [77] "reloc_col_before_first"                                         
#  [78] "find_preceding_colnm"                                           
#  [79] "order_psm_cols"                                                 
#  [80] "order_mascot_psm_cols"                                          
#  [81] "pad_mq_channels"                                                
#  [82] "splitPSM_mq"                                                    
#  [83] "pad_sm_channels"                                                
#  [84] "add_sm_pepseqmod"                                               
#  [85] "pad_mf_channels"                                                
#  [86] "splitPSM_mf"                                                    
#  [87] "round_msfmod_decimals"                                          
#  [88] "procMSGFmods"                                                   
#  [89] "fml2mass"                                                       
#  [90] "join_mgfs"                                                      
#  [91] "find_padding_pos"                                               
#  [92] "pad_tmt_channels"                                               
#  [93] "pad_psm_fields"                                                 
#  [94] "add_pepseqmod"                                                  
#  [95] "add_msgf_pepseqmod"                                             
#  [96] "check_dup_unimods"                                              
#  [97] "splitPSM_mz"                                                    
#  [98] "add_pep_retsd"                                                  
#  [99] "calc_pep_retsd"                                                 
# [100] "add_n_pepexpz"                                                  
# 
# $purge.R
# [1] "lgl_cleanup" "purge_by_cv" "purge_by_qt" "purge_by_n"  "psm_mpurge"  "purgePSM"    "purgePep"   
# 
# $qlfq.R
#  [1] "haddApexRTs_allsets"    "haddApexRTs_oneset"     "addApexRTs"             "addApexes"              "hpepLFQ"               
#  [6] "addMBRpeps"             "makeLFQSpecies"         "fillMBRnaRTs"           "pepLFQ"                 "mpepLFQ"               
# [11] "combine_lfq123"         "find_best_lfqvec"       "collapseSTY"            "mergeSTY"               "reduce_lfq_rtstep"     
# [16] "calc_mat_log2s_to_refs" "calc_vec_log2s_to_refs" "calc_mat_sds"           "saddMBR"               
# 
# $roadmap.R
# character(0)
# 
# $sigtests.R
# [1] "rowVars"          "filterData"       "prepFml"          "my_padj"          "lm_summary"       "model_onechannel"
# [7] "sigTest"          "pepSig"           "prnSig"          
# 
# $simulData.R
# [1] "simulUniprotPSM"         "mascot_refseq2uniprot"   "maxquant_refseq2uniprot"
# 
# $specialty.R
# [1] "labEffPSM"   "procTMT0"    "calcTMTLabs" "proteo_hm"  
# 
# $string.R
# [1] "annot_stringdb" "stringTest"     "anal_prnString" "prepString"     "load_stringdbs"
# 
# $trends.R
# [1] "checkdots_analTrend" "analTrend"           "plotTrend"           "anal_prnTrend"       "plot_prnTrend"      
# 
# $utils.R
#   [1] "prepDM"                      "find_NorZ"                   "replace_NorZ_names"          "`names_pos<-`"              
#   [5] "reorder_files2"              "reorderCols"                 "reorderCols2"                "ins_cols_after"             
#   [9] "sort_tmt_lcms"               "na_zeroIntensity"            "not_allzero_rows"            "aggrNums"                   
#  [13] "aggrTopn"                    "aggrLFQs"                    "tmt_wtmean"                  "  my_trim"                  
#  [17] "  my_wtmean"                 "TMT_top_n"                   "not_all_zero"                "not_all_NA"                 
#  [21] "not_all_nan"                 "is_all_nan"                  "is_trivial_col"              "is_trivial_dbl"             
#  [25] "replace_trivial_with_na"     "colAnnot"                    "setHMColor"                  "setHMlims"                  
#  [29] "ratio_toCtrl"                "imputeNA"                    "  my_mice"                   "  handleNA"                 
#  [33] "pepImp"                      "prnImp"                      "sp_lookup"                   "taxid_lookup"               
#  [37] "taxid_lookup_rev"            "sp_lookup_Ul"                "add_refseq_gene"             "add_entrez"                 
#  [41] "add_custom_entrez"           "parse_acc"                   "parse_fasta"                 "hparse_fasta"               
#  [45] "na_genes_by_acc"             "na_species_by_org"           "add_fasta_attrs"             "add_prot_desc"              
#  [49] "annotPrn"                    "  na_fasta_name_by_prot_acc" "  na_acc_type_to_other"      "annotKin"                   
#  [53] "save_call"                   "match_fmls"                  "match_call_arg"              "match_gset_nms"             
#  [57] "#' foo"                      "match_logi_gv"               "match_prnSig_scale_log2r"    "match_pepSig_scale_log2r"   
#  [61] "replace_na_genes"            "na_genes_by_acc"             "find_pep_pos"                "annotPeppos"                
#  [65] "subset_fasta"                "add_prot_icover"             "calc_cover"                  "to_linfc"                   
#  [69] "rm_sglval_cols"              "cmbn_meta"                   "gg_imgname"                  "rm_pval_whitespace"         
#  [73] "filters_in_call"             "arrangers_in_call"           "calc_sd_fcts_psm"            "calcSD_Splex"               
#  [77] "sd_violin"                   "rptr_violin"                 "my_geomean"                  "count_phosphopeps"          
#  [81] "count_pepmiss"               "contain_str"                 "contain_chars_in"            "not_contain_str"            
#  [85] "not_contain_chars_in"        "start_with_str"              "end_with_str"                "start_with_chars_in"        
#  [89] "ends_with_chars_in"          "rows_are_all"                "rows_are_not_all"            "concat_fml_dots"            
#  [93] "gn_rollup"                   "identical_dots"              "my_complete_cases"           "my_union"                   
#  [97] "find_fn_bases"               "find_fn_exts"                "check_dots"                  "check_depreciated_args"     
# [101] "to_complete_cases"           "check_gset_nms"              "ok_existing_params"          "env_where"                  
# [105] "set_cutpoints"               "set_cutpoints2"              "load_craps"                  "find_int_cols"              
# [109] "find_ratio_cols"             "make_mq_meta"                "make_mq_meta2"               "update_nodes"               
# [113] "pre_mq_meta"                 "add_entry_ids"               "check_formalArgs"            "check_formalArgs2"          
# [117] "set_ggsave_dots"             "find_search_engine"          "check_ggplot_aes"            "prep_range"                 
# [121] "find_delim"                  "recur_flatten"               "subset_keepattr"             "list_to_dataframe"          
# [125] "find_psmQ_files"             "get_col_types"               "write_excel_wb"              "load_ls_group"              
# [129] "parse_col_select"            "parse_filename"              "subsetBrukerMGF"             "msubsetBrukerMGF"           
# [133] "procBrukerMGF_v1"            "procBrukerMGF"               "mprocBrukerMGF"              "fixBrukerMGF"               
# [137] "mfixBrukerMGF"               "check_aes_length"            "add_col_rawfile"             "load_psmC"                  
# [141] "add_scan_id"                 "my_split1"                   "my_split2"                  
# 
# $utils_mzion.R
# [1] "find_mod_indexesQ"   "deduce_mod_indexes"  "hdeduce_mod_indexes" "find_path_ms1"      
# 
# $volcanos.R
#  [1] "plotVolcano"         "pval_complete_cases" "byfml_volcano"       "byfile_plotVolcano"  "subset_volcano_df"  
#  [6] "fullVolcano"         "summ_venn"           "plot_venn"           "gsVolcano"           "pepVol"             
# [11] "prnVol"              "gspaMap"             "geom_table"          "to_csv_"            
# 
# $wrappers.R
# [1] "my_dist" "cos_sim"
# 
# $zzz.R
# [1] ".onAttach"
# 

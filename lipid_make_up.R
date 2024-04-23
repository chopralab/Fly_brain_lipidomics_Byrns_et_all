
rm(list=ls())

library(plyr)
library(tidyverse)

library(readxl)

library(limma)

getwd()


setwd("C:/Users/Connor Beveridge/Box/connor_beveridge/New/age_flies/")

read_lipids <- function(prefix, marker, exp_names, types, add_Blank = NULL) {
  lapply(names(types), function(type){
    dir_name = paste0(prefix, "/")
    # dir_name = paste0(prefix, "/", type, "/")
    # suffix = paste0("_", types[[type]], "_results.xlsx")
    suffix =  "_12082022_results.xlsx"
    all_names <- lapply(exp_names, function(n) {
      file_name = paste0(dir_name,types[[type]],"_",n,  suffix)
      
      read_xlsx(file_name) %>%
        rename(lipid = `Lipid name`) %>%
        rename(!!n := `Total Intensity`) %>%
        mutate(type = type) %>%
        na.omit
    })
    
    if (!is.null(add_Blank)) {
      Blank = read_xlsx(paste0(dir_name, add_Blank, suffix)) %>%
        rename(lipid = `Lipid name`) %>%
        rename(!!add_Blank := `Total Intensity`) %>%
        mutate(type = type) %>%
        na.omit
      
      all_names = append(all_names, list(Blank))
    }
    
    Reduce(function(x, y) merge(x, y), all_names)
  }) %>% bind_rows() %>% as_tibble() %>% filter(!grepl("STD_", lipid))
}




lipid_types <- c( "FFA" = "FFA","AC"="AC",
                  "CE"="CE","CER"="CER","PCandSM"="PCandSM","PG"="PG","PE"="PE","PI"="PI", "TAG1" = "TAG1", "TAG2" = "TAG2") #Need to add PS#,"PS"="PS")






###5 vs 40 N
experiment_1_names <- c('40DN1','40DN2','40DN3','40DN4',"40DN5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/40DN.csv")



experiment_1_names <- c('5DN1','5DN2','5DN3','5DN4',"5DN5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/5DN.csv")


experiment_1_names <- c('40GFP1','40GFP2','40GFP3','40GFP4',"40GFP5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/40GFP.csv")

experiment_1_names <- c('5GFP1','5GFP2','5GFP3','5GFP4',"5GFP5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/5GFP.csv")


experiment_1_names <- c('42C1','42C2','42C3','42C4',"42C5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/42C.csv")

experiment_1_names <- c('7C1','7C2','7C3','7C4',"7C5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/7C.csv")

experiment_1_names <- c('40DP1','40DP2','40DP3','40DP4',"40DP5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/40DP.csv")



experiment_1_names <- c('42E1','42E2','42E3','42E4',"42E5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/42E.csv")

experiment_1_names <- c('7E1','7E2','7E3','7E4',"7E5","SolventBlank1")

cells_lipid_expr_1 = read_lipids(
  "fly_brain_12_8_22/results", "MEDIA",
  experiment_1_names, lipid_types
)
write_csv(cells_lipid_expr_1,"lipid_make_up/7E.csv")






# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # ##40 DP vs DN
# # experiment_4_names <- c('40DP1','40DP2','40DP3','40DP4',"40DP5",'40DN1','40DN2','40DN3','40DN4',"40DN5",  "SolventBlank1")
# # ##40GFP vs DN
# experiment_4_names <- c('40GFP1','40GFP2','40GFP3','40GFP4',"40GFP5",'40DN1','40DN2','40DN3','40DN4',"40DN5","SolventBlank1")
# 
# 
# ##40 GFP vs DN
# experiment_5_names <- c('40DN1','40DN2','40DN3','40DN4',"40DN5",'40GFP1','40GFP2','40GFP3','40GFP4',"40GFP5",  "SolventBlank1")
# 
# 
# #40 DP vs 5DP
# # experiment_6_names <- c('40DP1','40DP2','40DP3','40DP4',"40DP5",'5DP1','5DP2','5DP3','5DP4',"5DP5",  "SolventBlank1")
# 
# ##40 e vs 40 C
# experiment_7_names <- c('42E1','42E2','42E3','42E4',"42E5", '42C1','42C2','42C3','42C4',"42C5",  "SolventBlank1")
# 
# 
# 
# ##42E vs 7 E
# experiment_8_names <- c('42E1','42E2','42E3','42E4',"42E5",'7E1','7E2','7E3','7E4',"7E5",  "SolventBlank1")
# 
# ###7E vs 7C
# experiment_9_names <- c('7E1','7E2','7E3','7E4',"7E5",'7C1','7C2','7C3','7C4',"7C5",  "SolventBlank1")
# 
# ###5GFP vs 5DN
# experiment_10_names <- c('5GFP1','5GFP2','5GFP3','5GFP4',"5GFP5",'5DN1','5DN2','5DN3','5DN4',"5DN5",  "SolventBlank1")
# 
# ###7E vs 40C
# experiment_11_names <- c('7E1','7E2','7E3','7E4',"7E5",'42C1','42C2','42C3','42C4',"42C5",  "SolventBlank1")
# 
# 
# ###7E vs 40C
# experiment_11_names <- c('7E1','7E2','7E3','7E4',"7E5",'42C1','42C2','42C3','42C4',"42C5",  "SolventBlank1")
# 
# 
# ##40 DP vs GFP
# experiment_12_names <- c('40DP1','40DP2','40DP3','40DP4',"40DP5",'40GFP1','40GFP2','40GFP3','40GFP4',"40GFP5",  "SolventBlank1")
# 
# ##40 DP vs DN
# experiment_13_names <- c('40DP1','40DP2','40DP3','40DP4',"40DP5",'40DN1','40DN2','40DN3','40DN4',"40DN5", "SolventBlank1")
# 
# 
# 
# experiment_14_names <- c('40DP1','40DP2','40DP3','40DP4',"40DP5",'5GFP1','5GFP2','5GFP3','5GFP4',"5GFP5", "SolventBlank1")
# 
# 
# 
# title_for_plot_1 <- "40DN_vs_5DN_other_lipids"
# title_for_plot_2 <- "40GFP_vs_5GFP_other_lipids"
# title_for_plot_3 <- "42C_vs_7C_other_lipids"
# 
# title_for_plot_4 <- "40GFP_vs_40DN_other_lipids"
# title_for_plot_5 <- "40DN_vs_40GFP_other_lipids"
# title_for_plot_6 <- "40DP_vs_5DP_other_lipids"
# title_for_plot_7 <- "42E_vs_42C_other_lipids"
# title_for_plot_8 <- "42E_vs_7E_other_lipids"
# title_for_plot_9 <- "7E_vs_7C_other_lipids"
# title_for_plot_10 <- "5GFP_vs_5DN_other_lipids"
# title_for_plot_11 <- "7E_vs_42C_other_lipids"
# title_for_plot_12 <- "40DP_vs_40GFP_other_lipids"
# title_for_plot_13 <- "40DP_vs_40DN_other_lipids"
# 
# title_for_plot_14 <- "40DP_vs_5GFP_other_lipids"
# 
# 
# 
# 
# 
# 
# 
# cells_lipid_expr_1 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_1_names, lipid_types
# )
# 
# 
# cells_lipid_expr_2 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_2_names, lipid_types
# )
# 
# 
# 
# cells_lipid_expr_3 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_3_names, lipid_types
# )
# 
# 
# cells_lipid_expr_4 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_4_names, lipid_types
# )
# 
# cells_lipid_expr_5 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_5_names, lipid_types
# )
# 
# 
# # cells_lipid_expr_6 = read_lipids(
# #   "fly_brain_12_8_22/results", "MEDIA",
# #   experiment_6_names, lipid_types
# # )
# 
# cells_lipid_expr_7 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_7_names, lipid_types
# )
# # 
# cells_lipid_expr_8 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_8_names, lipid_types
# )
# 
# cells_lipid_expr_9 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_9_names, lipid_types
# )
# 
# cells_lipid_expr_10 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_10_names, lipid_types
# )
# 
# cells_lipid_expr_11 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_11_names, lipid_types
# )
# 
# cells_lipid_expr_12 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_12_names, lipid_types
# )
# 
# cells_lipid_expr_13 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_13_names, lipid_types
# )
# 
# 
# cells_lipid_expr_14 = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   experiment_14_names, lipid_types
# )
# 
# 
# 
# PCA_names = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, lipid_types
# )
# PCA_names_to_parse_No_Neurons
# 
# PCA_names_No_Neurons = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse_No_Neurons, lipid_types
# )
# 
# 
# TAG = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("TAG1" = "TAG1", "TAG2" = "TAG2")
# )
# write_csv(TAG,"heatmap_data_files/TAG_part4.csv")
# 
# FFA = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("FFA"="FFA")
# )
# write_csv(FFA,"heatmap_data_files/FFA_part3.csv")
# 
# 
# 
# CE = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("CE"="CE")
# )
# write_csv(CE,"heatmap_data_files/CE_part3.csv")
# 
# 
# CER = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("CER"="CER")
# )
# write_csv(CER,"heatmap_data_files/CER_part3.csv")
# 
# 
# 
# 
# PCandSM = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("PCandSM"="PCandSM")
# )
# write_csv(PCandSM,"heatmap_data_files/PCandSM_part3.csv")
# 
# 
# PG = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("PG"="PG")
# )
# write_csv(PG,"heatmap_data_files/PG_part3.csv")
# 
# PE = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("PE"="PE")
# )
# write_csv(PE,"heatmap_data_files/PE_part3.csv")
# 
# PI = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("PI"="PI")
# )
# write_csv(PI,"heatmap_data_files/PI_part3.csv")
# 
# 
# AC= read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("AC"="AC")
# )
# write_csv(AC,"heatmap_data_files/AC_part3.csv")


# PS = read_lipids(
#   "fly_brain_12_8_22/results", "MEDIA",
#   PCA_names_to_parse, c("PS"="PS")
# )
# write_csv(PS,"heatmap_data_files/PS_part3.csv")





# library(writexl)
# # 
# # write_xlsx(x = blank_names_data, path = "Blanks_8_2_22.xlsx", col_names = TRUE)
# 
# 
# write_xlsx(x = PCA_names, path = "Fly_Brain_Intensities_other_lipids.xlsx", col_names = TRUE)
# 
# write_xlsx(x = PCA_names_No_Neurons, path = "Fly_Brain_Intensities_other_lipids_No_Neuron.xlsx", col_names = TRUE)
# 
# library(edgeR)
# 
# # 
# 
# 
# 
# # Groups for experiment 1
# gr_expr_5 = c("GR1","GR1","GR1", "GR1","GR1","GR2","GR2","GR2", "GR2","GR2",
#               "SolventBlank1") %>%
#   factor(levels = c("SolventBlank1", "GR1", "GR2"))
# 
# 
# 
# design_expr_5 = model.matrix(~gr_expr_5)
# 
# contrasts_expr_5 = makeContrasts(
#   H = gr_expr_5GR1 - gr_expr_5GR2,
#   levels = design_expr_5
# )
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # gr_expr_1 = c("5DN","5DN","5DN","5DN","5DN", "40DN","40DN","40DN","40DN","40DN",
# #               "SolventBlank1") %>%
# #   factor(levels = c("SolventBlank1", "5DN", "40DN"))
# # 
# # 
# # 
# # design_expr_1 = model.matrix(~gr_expr_1)
# # 
# # contrasts_expr_1 = makeContrasts(
# #   H = gr_expr_15DN - gr_expr_140DN,
# #   levels = design_expr_1
# # )
# # # contrasts_expr_1 = makeContrasts(
# # #   H = gr_expr_1Control3x - gr_expr_1APIblock,
# # #   levels = design_expr_1
# # # )
# # 
# # 
# # # Groups for experiment 1
# # gr_expr_2 = c("5GFP","5GFP","5GFP","5GFP","5GFP","40GFP","40GFP","40GFP","40GFP","40GFP",
# #               "SolventBlank1") %>%
# #   factor(levels = c("SolventBlank1", "5GFP", "40GFP"))
# # 
# # 
# # 
# # design_expr_2 = model.matrix(~gr_expr_2)
# # 
# # contrasts_expr_2 = makeContrasts(
# #   H = gr_expr_25GFP - gr_expr_240GFP,
# #   levels = design_expr_2
# # )
# # 
# # 
# # 
# # # Groups for experiment 1
# # gr_expr_3 = c("DP","DP","DP", "DP","DP","GFP","GFP","GFP", "GFP","GFP",
# #               "SolventBlank1") %>%
# #   factor(levels = c("SolventBlank1", "DP", "GFP"))
# # 
# # 
# # 
# # design_expr_3 = model.matrix(~gr_expr_3)
# # 
# # contrasts_expr_3 = makeContrasts(
# #   H = gr_expr_340DP - gr_expr_340GFP,
# #   levels = design_expr_3
# # )
# # 
# # 
# # 
# # # Groups for experiment 1
# # gr_expr_4 = c("DP","DP","DP", "DP","DP","DN","DN","DN", "DN","DN",
# #               "SolventBlank1") %>%
# #   factor(levels = c("SolventBlank1", "DP", "DN"))
# # 
# # 
# # 
# # design_expr_4 = model.matrix(~gr_expr_4)
# # 
# # contrasts_expr_4 = makeContrasts(
# #   H = gr_expr_4DP - gr_expr_4DN,
# #   levels = design_expr_4
# # )
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # Groups for experiment 1
# # gr_expr_5 = c("GFP","GFP","GFP", "GFP","GFP","DN","DN","DN", "DN","DN",
# #               "SolventBlank1") %>%
# #   factor(levels = c("SolventBlank1", "GFP", "DN"))
# # 
# # 
# # 
# # design_expr_5 = model.matrix(~gr_expr_5)
# # 
# # contrasts_expr_5 = makeContrasts(
# #   H = gr_expr_5GFP - gr_expr_5DN,
# #   levels = design_expr_5
# # )
# # 
# 
# 
# 
# 
# perform_analysis_raw <- function(counts, design_mat, gr) {
#   
#   data.edgeR <- DGEList(counts = counts %>%
#                           na.omit %>%
#                           mutate(lipid = make.unique(lipid)) %>%
#                           select(-Transition, -type) %>%
#                           column_to_rownames("lipid"),
#                         group = gr
#   )
#   
#   data.edgeR <- calcNormFactors(data.edgeR, method="TMM")
#   data.edgeR <- estimateCommonDisp(data.edgeR, design=design_mat)
#   data.edgeR
# }
# 
# calculate_significance <- function(dge, contrast) {
#   dge %>%
#     glmFit() %>%
#     glmLRT(contrast = contrast)
# }
# 
# 
# 
# 
# 
# 
# 
# # 
# # experiment_1_helper <- function(df) {
# #   df %>%
# #     perform_analysis_raw(design_expr_1, gr_expr_1) %>%
# #     calculate_significance(contrasts_expr_1) %>%
# #     topTags(1500000) %>%
# #     as.data.frame() %>%
# #     rownames_to_column("lipid") %>%
# #     as_tibble()
# #   
# # }
# 
# 
# 
# # 
# # experiment_2_helper <- function(df) {
# #   df %>%
# #     perform_analysis_raw(design_expr_2, gr_expr_2) %>%
# #     calculate_significance(contrasts_expr_2) %>%
# #     topTags(1500000) %>%
# #     as.data.frame() %>%
# #     rownames_to_column("lipid") %>%
# #     as_tibble()
# #   
# # }
# # 
# # 
# # 
# # 
# # experiment_3_helper <- function(df) {
# #   df %>%
# #     perform_analysis_raw(design_expr_3, gr_expr_3) %>%
# #     calculate_significance(contrasts_expr_3) %>%
# #     topTags(1500000) %>%
# #     as.data.frame() %>%
# #     rownames_to_column("lipid") %>%
# #     as_tibble()
# #   
# # }
# # 
# # 
# # 
# # 
# # experiment_4_helper <- function(df) {
# #   df %>%
# #     perform_analysis_raw(design_expr_4, gr_expr_4) %>%
# #     calculate_significance(contrasts_expr_4) %>%
# #     topTags(1500000) %>%
# #     as.data.frame() %>%
# #     rownames_to_column("lipid") %>%
# #     as_tibble()
# #   
# # }
# # 
# 
# 
# 
# experiment_5_helper <- function(df) {
#   df %>%
#     perform_analysis_raw(design_expr_5, gr_expr_5) %>%
#     calculate_significance(contrasts_expr_5) %>%
#     topTags(1500000) %>%
#     as.data.frame() %>%
#     rownames_to_column("lipid") %>%
#     as_tibble()
#   
# }
# 
# 
# # perform_analysis_raw(design_expr_3, gr_expr_3)
# # 
# # 
# # calculate_significance(contrasts_expr_2)
# # 
# 
# 
# cl_e1_tbl <-
#   cells_lipid_expr_1 %>%
#   experiment_5_helper
# cl_e1_tbl
# 
# 
# 
# cl_e2_tbl <-
#   cells_lipid_expr_2 %>%
#   experiment_5_helper
# cl_e2_tbl
# # 
# cl_e3_tbl <-
#   cells_lipid_expr_3 %>%
#   experiment_5_helper
# cl_e3_tbl
# # 
# # 
# # 
# cl_e4_tbl <-
#   cells_lipid_expr_4 %>%
#   experiment_5_helper
# cl_e4_tbl
# 
# 
# 
# cl_e5_tbl <-
#   cells_lipid_expr_5 %>%
#   experiment_5_helper
# cl_e5_tbl
# 
# 
# # cl_e6_tbl <-
# #   cells_lipid_expr_6 %>%
# #   experiment_5_helper
# 
# 
# 
# cl_e7_tbl <-
#   cells_lipid_expr_7 %>%
#   experiment_5_helper
# 
# cl_e8_tbl <-
#   cells_lipid_expr_8 %>%
#   experiment_5_helper
# 
# cl_e9_tbl <-
#   cells_lipid_expr_9 %>%
#   experiment_5_helper
# 
# cl_e10_tbl <-
#   cells_lipid_expr_10 %>%
#   experiment_5_helper
# 
# cl_e11_tbl <-
#   cells_lipid_expr_11 %>%
#   experiment_5_helper
# 
# 
# cl_e12_tbl <-
#   cells_lipid_expr_12 %>%
#   experiment_5_helper
# 
# 
# 
# cl_e13_tbl <-
#   cells_lipid_expr_13 %>%
#   experiment_5_helper
# 
# 
# cl_e14_tbl <-
#   cells_lipid_expr_14 %>%
#   experiment_5_helper
# 
# 
# 
# 
# make_volcano_plot <- function(df, title) {
#   df %>%
#     mutate(sig = factor(FDR < 0.10)) %>%
#     ggplot(aes(logFC, -log10(FDR), color = sig)) +
#     geom_point() +
#     scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#     guides(color = F) +
#     ggtitle(title)
# }
# 
# cl_e1_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_1,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_1,".png",sep=''))
# 
# 
# cl_e2_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_2,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_2,".png",sep=''))
# 
# # 
# cl_e3_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_3,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_3,".png",sep=''))
# 
# 
# cl_e4_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_4,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_4,".png",sep=''))
# 
# 
# 
# cl_e5_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_5,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_5,".png",sep=''))
# # 
# # cl_e6_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_6,sep=''))
# # ggsave(paste("plots/Volcano_Plot_",title_for_plot_6,".png",sep=''))
# # 
# 
# cl_e7_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_7,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_7,".png",sep=''))
# 
# cl_e8_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_8,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_8,".png",sep=''))
# 
# cl_e9_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_9,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_9,".png",sep=''))
# cl_e10_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_10,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_10,".png",sep=''))
# cl_e11_tbl %>% make_volcano_plot(paste("Volcano_Plot_",title_for_plot_11,sep=''))
# ggsave(paste("plots/Volcano_Plot_",title_for_plot_11,".png",sep=''))
# 
# write_summary_and_results <- function(tbl, df, name) {
#   
#   tbl %>%
#     merge(df) %>%
#     as_tibble() %>%
#     arrange(FDR) -> results
#   
#   write_csv(results, paste0("results/", name, "_full.csv"))
#   
#   results %>%
#     group_by(type) %>%
#     summarise(ab = sum(logFC < 0 & FDR < 0.1),
#               ve = sum(logFC > 0 & FDR < 0.1)) %>%
#     write_csv(paste0("results/", name, "_summary.csv"))
# }
# 
# dir.create("results", F)
# 
# cl_e1_tbl %>% write_summary_and_results(cells_lipid_expr_1, title_for_plot_1)
# 
# 
# cl_e2_tbl %>% write_summary_and_results(cells_lipid_expr_2, title_for_plot_2)
# # 
# cl_e3_tbl %>% write_summary_and_results(cells_lipid_expr_3, title_for_plot_3)
# # 
# cl_e4_tbl %>% write_summary_and_results(cells_lipid_expr_4, title_for_plot_4)
# cl_e5_tbl %>% write_summary_and_results(cells_lipid_expr_5, title_for_plot_5)
# 
# 
# 
# # 
# # 
# # cl_e6_tbl %>% write_summary_and_results(cells_lipid_expr_6, title_for_plot_6)
# 
# cl_e7_tbl %>% write_summary_and_results(cells_lipid_expr_7, title_for_plot_7)
# cl_e8_tbl %>% write_summary_and_results(cells_lipid_expr_8, title_for_plot_8)
# cl_e9_tbl %>% write_summary_and_results(cells_lipid_expr_9, title_for_plot_9)
# cl_e10_tbl %>% write_summary_and_results(cells_lipid_expr_10, title_for_plot_10)
# cl_e11_tbl %>% write_summary_and_results(cells_lipid_expr_11, title_for_plot_11)
# 
# cl_e12_tbl %>% write_summary_and_results(cells_lipid_expr_12, title_for_plot_12)
# 
# cl_e13_tbl %>% write_summary_and_results(cells_lipid_expr_13, title_for_plot_13)
# 
# 
# 
# 
# 
# 
# cl_e14_tbl %>% write_summary_and_results(cells_lipid_expr_14, title_for_plot_14)
# 
# 
# 
# 
# 
# 
# 
# source("ggbiplot.R")
# 
# 
# get_DE_lipids <- function(counts, design_mat, gr, contrasts, p.value = 0.1) {
#   dls <-
#     counts %>%
#     perform_analysis_raw(design_mat, gr) %>%
#     calculate_significance(contrasts) %>%
#     decideTestsDGE(p.value = p.value)
#   
#   rownames(dls)[dls %>% as.logical()]
# }
# 
# make_pca_plot <- function(tp, design_mat, gr, contrasts,
#                           title = "PCA plot",
#                           ellipse = T, var.axes = F,
#                           labels = T) {
#   
#   print(labels)
#   print(labels)
#   tp_edger <-
#     tp %>%
#     get_DE_lipids(design_mat, gr, contrasts)
#   
#   if(length(tp_edger) == 0) {
#     cat("No significant lipids for ", title)
#     return()
#   }
#   
#   if(length(tp_edger) == 1) {
#     cat("Single significant lipids for ", title, " is ", tp_edger[1])
#     return()
#   }
#   
#   tp %>%
#     na.omit %>%
#     mutate(lipid = make.unique(lipid)) %>%
#     filter(lipid %in% tp_edger) %>%
#     select(-Transition, -type, -SolventBlank1) %>%
#     column_to_rownames("lipid") %>%
#     as.matrix() %>%
#     t %>%
#     prcomp(center = T, scale = T) ->
#     prcomp_data
#   
#   groups = NULL
#   
#   tp %>%
#     select(-Transition, -type, -SolventBlank1, -lipid) %>%
#     colnames() ->
#     labels.tmp
#   
#   groups = substr(labels.tmp, 1, 2)
#   print(groups)
#   print(substr(labels.tmp, 1, 4))
#   if (!is.null(labels)) {
#     labels = labels.tmp
#   }
#   
#   prcomp_data %>%
#     ggbiplot(ellipse = ellipse,
#              labels = labels,
#              groups = groups,
#              var.axes = var.axes
#     ) +
#     ggtitle(title) +
#     cowplot::theme_cowplot()
# }
# 
# 
# 
# 
# 
# make_pca_plot2 <- function(tp, design_mat, gr, contrasts,
#                           title = "PCA plot",
#                           ellipse = T, var.axes = F,
#                           labels = T) {
#   
#   print(labels)
#   print(labels)
#   tp_edger <-
#     tp %>%
#     get_DE_lipids(design_mat, gr, contrasts)
#   
#   if(length(tp_edger) == 0) {
#     cat("No significant lipids for ", title)
#     return()
#   }
#   
#   if(length(tp_edger) == 1) {
#     cat("Single significant lipids for ", title, " is ", tp_edger[1])
#     return()
#   }
#   
#   tp %>%
#     na.omit %>%
#     mutate(lipid = make.unique(lipid)) %>%
#     filter(lipid %in% tp_edger) %>%
#     select(-Transition, -type, -SolventBlank1) %>%
#     column_to_rownames("lipid") %>%
#     as.matrix() %>%
#     t %>%
#     prcomp(center = T, scale = T) ->
#     prcomp_data
#   
#   groups = NULL
#   
#   tp %>%
#     select(-Transition, -type, -SolventBlank1, -lipid) %>%
#     colnames() ->
#     labels.tmp
#   
#   groups = substr(labels.tmp, 1, 4)
#   print(groups)
#   print(substr(labels.tmp, 1, 5))
#   if (!is.null(labels)) {
#     labels = labels.tmp
#   }
#   
#   prcomp_data %>%
#     ggbiplot(ellipse = ellipse,
#              labels = labels,
#              groups = groups,
#              var.axes = var.axes
#     ) +
#     ggtitle(title) +
#     cowplot::theme_cowplot()
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# cells_lipid_expr_1 %>%
#   make_pca_plot(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_1,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_1,".svg",sep=''))
# 
# cells_lipid_expr_2 %>%
#   make_pca_plot(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_2,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_2,".svg",sep=''))
# # 
# 
# cells_lipid_expr_3 %>%
#   make_pca_plot(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_3,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_3,".svg",sep=''))
# # 
# # 
# # 
# cells_lipid_expr_4 %>%
#   make_pca_plot2(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_4,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_4,".svg",sep=''))
# 
# 
# cells_lipid_expr_5 %>%
#   make_pca_plot2(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_5,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_5,".svg",sep=''))
# 
# 
# # cells_lipid_expr_6 %>%
# #   make_pca_plot(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_6,sep=''))
# # ggsave(paste("plots/PCA_",title_for_plot_6,".svg",sep=''))
# # 
# 
# 
# 
# cells_lipid_expr_7 %>%
#   make_pca_plot2(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_7,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_7,".svg",sep=''))
# 
# 
# 
# cells_lipid_expr_8 %>%
#   make_pca_plot(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_8,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_8,".svg",sep=''))
# 
# 
# 
# cells_lipid_expr_9 %>%
#   make_pca_plot(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_9,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_9,".svg",sep=''))
# 
# 
# cells_lipid_expr_10 %>%
#   make_pca_plot2(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_10,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_10,".svg",sep=''))
# 
# 
# 
# cells_lipid_expr_11 %>%
#   make_pca_plot(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_11,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_11,".svg",sep=''))
# 
# 
# 
# cells_lipid_expr_12 %>%
#   make_pca_plot2(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_12,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_12,".svg",sep=''))
# 
# cells_lipid_expr_13 %>%
#   make_pca_plot2(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_13,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_13,".svg",sep=''))
# 
# 
# cells_lipid_expr_14 %>%
#   make_pca_plot2(design_expr_5, gr_expr_5, contrasts_expr_5, paste("PCA_Plot_",title_for_plot_14,sep=''))
# ggsave(paste("plots/PCA_",title_for_plot_14,".svg",sep=''))
# 
# 
# 
# 
# get_DE_lipids2 <- function(counts, design_mat, gr, contrasts, p.value = .01) {
#   dls <-
#     counts %>%
#     perform_analysis_raw(design_mat, gr) %>%
#     calculate_significance(contrasts) %>%
#     decideTestsDGE(p.value = p.value)
#   
#   rownames(dls)[dls %>% as.logical()]
# }
# 
# 
# make_heatmap <- function(tp, design_mat, gr, contrasts, title = "Heat-map") {
#   
#   DElist <-
#     tp %>%
#     get_DE_lipids2(design_mat, gr, contrasts)
#   
#   if(length(DElist) == 0) {
#     cat("No significant lipids for ", title)
#     return()
#   }
#   
#   if(length(DElist) == 1) {
#     cat("Single significant lipids for ", title, " is ", DElist[1])
#     return()
#   }
#   
#   tp %>%
#     mutate(lipid = make.unique(lipid)) %>%
#     filter(lipid %in% DElist) %>%
#     select(-Transition, -type) %>%
#     # mutate_if(is.numeric, log2) %>%
#     # mutate_if(is.numeric, list(~ . - SolventBlank1)) %>%
#     mutate_if(is.numeric, log2) %>%
#     select(-SolventBlank1) %>%
#     column_to_rownames("lipid") %>%
#     as.matrix() %>%
#     pheatmap::pheatmap(main = title)
# }
# 
# cells_lipid_expr_1 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_1)
# ggsave(paste("plots/Heat_map_",title_for_plot_1,".svg",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# 
# 
# 
# cells_lipid_expr_2 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_2)
# ggsave(paste("plots/Heat_map_",title_for_plot_2,".svg",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# # 
# # 
# cells_lipid_expr_3 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_3)
# ggsave(paste("plots/Heat_map_",title_for_plot_3,".svg",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# cells_lipid_expr_4 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_4)
# ggsave(paste("plots/Heat_map_",title_for_plot_4,".svg",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# 
# 
# cells_lipid_expr_5 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_5)
# ggsave(paste("plots/Heat_map_",title_for_plot_5,".pdf",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# 
# 
# 
# # cells_lipid_expr_6 %>%
# #   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_6)
# # 
# # 
# # ggsave(paste("plots/Heat_map_",title_for_plot_5,".pdf",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# # 
# # 
# 
# 
# cells_lipid_expr_7 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_7)
# ggsave(paste("plots/Heat_map_",title_for_plot_5,".pdf",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# 
# cells_lipid_expr_8 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_8)
# ggsave(paste("plots/Heat_map_",title_for_plot_5,".pdf",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# cells_lipid_expr_9 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_9)
# ggsave(paste("plots/Heat_map_",title_for_plot_5,".pdf",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# cells_lipid_expr_10 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_10)
# ggsave(paste("plots/Heat_map_",title_for_plot_5,".pdf",sep=''),width = 50, height = 25, dpi = 750,limitsize = FALSE)
# 
# cells_lipid_expr_11 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_11)
# ggsave(paste("plots/Heat_map_",title_for_plot_5,".pdf",sep=''),width = 50, height = 200, dpi = 750,limitsize = FALSE)
# 
# 
# 
# 
# cells_lipid_expr_12 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_12)
# ggsave(paste("plots/Heat_map_",title_for_plot_12,".pdf",sep=''),width = 50, height = 200, dpi = 750,limitsize = FALSE)
# 
# cells_lipid_expr_13 %>%
#   make_heatmap(design_expr_5, gr_expr_5, contrasts_expr_5, title_for_plot_13)
# ggsave(paste("plots/Heat_map_",title_for_plot_13,".pdf",sep=''),width = 50, height = 200, dpi = 750,limitsize = FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###Probably need to change
# # 
# # 
# # make_heatmap_jf_e1 <- function(tp, design_mat, gr, contrasts, title = "Heat-map", FDR = 0.10) {
# #   
# #   DElist <-
# #     tp %>%
# #     get_DE_lipids(design_mat, gr, contrasts, FDR)
# #   
# #   
# #   DElist
# #   
# #   if(length(DElist) == 0) {
# #     cat("No significant lipids for ", title)
# #     return()
# #   }
# #   
# #   if(length(DElist) == 1) {
# #     cat("Single significant lipids for ", title, " is ", DElist[1])
# #     return()
# #   }
# #   
# #   
# #   tp %>%
# #     mutate(lipid = make.unique(lipid)) %>%
# #     filter(lipid %in% DElist) %>%
# #     select(-Transition, -type) %>%
# #     rowwise() %>%
# #     mutate(mean_C = mean(c(KP365_C1_500nM, KP365_C2_500nM))) %>%
# #     
# #     mutate(mean_T = mean(c(KP365_T1_500nM, KP365_T2_500nM))) %>%
# #     
# #     mutate(mean = mean(c(mean_T, mean_C))) %>%
# #     mutate_if(is.numeric, log10) %>%
# #     mutate_if(is.numeric, list(~ . - mean)) %>%
# #     mutate(sd = sd(c(mean_T, mean_C))) %>%
# #     mutate_if(is.numeric, list(~ . / sd))
# #     select(-Blank, -mean) ->
# #       vals
# #     
# #     
# #     # mutate(mean = mean(c(KP365_C1_500nM, KP365_C2_500nM, KP365_T1_500nM, KP365_T2_500nM))) %>%
# #     # mutate_if(is.numeric, log2) %>%
# #     # mutate_if(is.numeric, list(~ . - mean)) %>%
# #     # select(-Blank, -mean) ->
# #     # vals
# #     # 
# #   vals %>%
# #     column_to_rownames("lipid") %>%
# #     as.matrix() %>%
# #     pheatmap::pheatmap(main = title)
# # }
# # 
# # cells_lipid_expr_1 %>%
# #   make_heatmap_jf_e1(design_expr_1, gr_expr_1, contrasts_expr_1, title_for_plot, 0.1)
# # 
# 
# 


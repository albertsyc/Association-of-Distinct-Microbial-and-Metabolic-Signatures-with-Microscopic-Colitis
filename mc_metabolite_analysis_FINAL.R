rm(list=ls())

library(readxl)
library(readr)
library(openxlsx)
library(Maaslin2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(fossil)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(stringr)
library(RColorBrewer)
library(forcats)
library(writexl)
library(tidyverse)

main.dir<- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R"
setwd(main.dir)

####---------------  Figure PCoA_cross-sectional  ----------------------#####
load("mc_metabolite_input_2024-04-25.RData") 
load("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/match_2024-05-13.RData")
load("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mc_input_2024-07-08.RData")

# Only keep rows with both microbiome and metabolomics data
mc_hc$mc_all <- as.factor(mc_hc$mc_all)
mc_cd_hc <- full_join(mc_cd, mc_hc)


## Change compund names
index <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/index.xlsx")
chem_to_comp <- index[, c(1, 3)]
# Rename columns
metabolite <- metabolite %>%
  rename_at(vars(as.character(chem_to_comp$CHEM_ID)),
            ~ as.character(chem_to_comp$COMP_ID)) 


metabolite$id <- gsub("_1$", "A", metabolite$id)
metabolite$id <- gsub("_2$", "B", metabolite$id)
rownames(metabolite) <- metabolite$id




metabolite <- metabolite %>% inner_join(select(mc_cd_hc, id), by = c("ID" = "id"))
rownames(metabolite) <- metabolite$id

metabolite_1 <- metabolite[grepl("^G", metabolite$id), ]
metabolite_2 <- metabolite[grepl("^MC", metabolite$id), ]
metabolite_2 <- metabolite_2 %>% filter(symptoms=="active_diarrhea")

to_remove <- metabolite_2 %>%
  arrange(ID) %>%
  filter(duplicated(ID) | duplicated(ID, fromLast = TRUE)) # 11 MC samples are repeated

metabolite_2 <- metabolite_2 %>%
  arrange(id) %>%
  distinct(ID, .keep_all = TRUE) # number of samples = 230 (11 repeated MC removed)

metabolite_2 %>% filter(group=="MC") %>% nrow() # active MC = 92
metabolite_2 %>% filter(group=="Diarrhea") %>% nrow() # active MC = 138
metabolite_active <- rbind(metabolite_1, metabolite_2)

subset.metabolite <- metabolite_active[, !colnames(metabolite_active) %in% c("id","ID","group","disease","symptoms")]
rownames(subset.metabolite) <- metabolite_active$id

bray <- vegdist(subset.metabolite, "bray")  # pairwise dissmilarity: 0=same; 1=maximally dissimilar
pcoa <- cmdscale(bray, eig = T)
pcoap <- data.frame(pcoa$points)
pcoap$ID <- rownames(pcoap)
pcoap <- pcoap[order(pcoap$ID),]

# percentages of variation explained by PCO1 & 2
eigs <- pcoa$eig
pc1.pct <- eigs[1]/sum(eigs) 
pc1.pct # pc1.pct = 15.3%
pc2.pct <- eigs[2]/sum(eigs) 
pc2.pct # pc2.pct = 4.4%

metadata$id <- gsub("_1$", "A", metadata$id)
metadata$id <- gsub("_2$", "B", metadata$id)
rownames(metadata) <- metadata$id

metadata.pcoa <- left_join(pcoap, metadata, by = c("ID" = "id"))
metadata.pcoa <- metadata.pcoa %>% rename("id" = "ID.y")

metadata.pcoa$group <- factor(metadata.pcoa$group, levels = c("Diarrhea", "MC", "Control"),
                              labels = c("Chronic diarrhea", "MC", "Control without diarrhea"))

fig_PCoA_metabolite <- ggscatter(metadata.pcoa, y = "X2", x = "X1",
                                 ylab = paste0('PCo2 (', round(pc2.pct * 100, 1), '%)'),
                                 xlab = paste0('PCo1 (', round(pc1.pct * 100, 1), '%)'),
                                 font.x = 12, font.y = 12, font.tickslab = 11,
                                 color = "group", 
                                 size = 2,
                                 legend = "bottom", legend.title = "Legend Title") +
  stat_ellipse(aes(color = group), linetype = 2) +
  guides(color = guide_legend(title = "Disease type")) +
  scale_color_manual(values = c("Chronic diarrhea" = "#F9D937", 
                                "MC" = "#900C3F",
                                "Control without diarrhea" = "#189AF9"),
                     labels = c("Chronic diarrhea" = "Chronic diarrhea", 
                                "MC" = "MC", 
                                "Control without diarrhea" = "Control without diarrhea"))

fig_PCoA_metabolite

# ggsave(filename = "fig_PCoA_metabolite.pdf",
#        plot=fig_PCoA_metabolite, units = "in", height = 5, width = 7, dpi=1200)


####---------------  Spearman correlation between the first axes of ordination for metabolite and microbiome (Cross-sectional)  ----------------------#####
## Cross-sectional
mgx.pcoa.cross <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/alpha_pcoa_cross.xlsx")
corr_cross_1 <- metadata.pcoa %>% select(ID, X1)
corr_cross_1 <- corr_cross_1 %>% rename("mbx_X1" = "X1")

corr_cross_2 <- mgx.pcoa.cross %>% select(ID, X1)
corr_cross_2 <- corr_cross_2 %>% rename("mgx_X1" = "X1")
corr_cross_filled <- left_join(corr_cross_1, corr_cross_2) %>% filter(!is.na(mbx_X1) & !is.na(mgx_X1))
corr_cross_unfilled <- left_join(corr_cross_1, corr_cross_2) %>% filter(is.na(mbx_X1) | is.na(mgx_X1)) 
# manually add mgx_X1 data
corr_cross_unfilled$mgx_X1 <- c("-0.1754828752", "0.4465455316", "0.0153313721", "0.2582425685", "0.1458070235", "-0.2089990123", "0.1845611016", "0.1956138555")

corr_cross <- rbind(corr_cross_filled, corr_cross_unfilled)
corr_cross$mgx_X1 <- as.numeric(corr_cross$mgx_X1)

# Calculate Spearman correlation 
spearman_corr_cross <- cor.test(corr_cross$mbx_X1, corr_cross$mgx_X1, method = "spearman")
spearman_corr_cross # rho = 0.579, p-value < 2.2e-16





####---------------  MaAsLin: Cross-sectional (all)  ----------------------#####
index <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/index.xlsx")
chem_to_comp <- index[, c(1, 3)]

comp_metabolite <- log_metabolite
# Rename columns
comp_metabolite <- comp_metabolite %>%
  rename_at(vars(as.character(chem_to_comp$CHEM_ID)),
            ~ as.character(chem_to_comp$COMP_ID)) 


comp_metabolite$id <- gsub("_1$", "A", comp_metabolite$id)
comp_metabolite$id <- gsub("_2$", "B", comp_metabolite$id)
rownames(comp_metabolite) <- comp_metabolite$id

comp_metabolite_maaslin <- comp_metabolite %>% select(6:1780)


## Active MC at baseline (with both mgx and mbx)
activeMC_baseline <- metadata %>% filter(group=="MC" & symptoms=="active_diarrhea") # Active MC at baseline with both mgx and mbx: n=92
activeMC_baseline <- activeMC_baseline %>% inner_join(select(mc_cd_hc, id), by = c("ID" = "id"))

activeMC_baseline <- activeMC_baseline %>%
  arrange(id) %>%
  distinct(ID, .keep_all = TRUE) # 11 MC samples are repeated --> remove

## HC at baseline (with both mgx and mbx)
HC_baseline <- metadata %>% filter(disease=="Control" & symptoms=="no_diarrhea") 
HC_baseline <- HC_baseline %>% inner_join(select(mc_cd_hc, id)) # HC at baseline with both mgx and mbx: n=105


## Active CD at baseline (with both mgx and mbx)
CD_baseline <- metadata %>% filter(disease=="Diarrhea" & symptoms=="active_diarrhea") 
CD_baseline <- CD_baseline %>% inner_join(select(mc_cd_hc, id)) # CD at baseline with both mgx and mbx: n=138

cross_baseline <- bind_rows(activeMC_baseline, HC_baseline, CD_baseline)
cross_baseline <- cross_baseline %>% left_join(select(mc_cd_hc, id, bmi), by = c("ID" = "id"))
cross_baseline <- as.data.frame(cross_baseline) 
rownames(cross_baseline) <- cross_baseline$id


# met_cross_all <- Maaslin2(input_data = comp_metabolite_maaslin,
#                           input_metadata = cross_baseline,
#                           output = "met_cross_all",
#                           min_abundance = 0.001,
#                           min_prevalence = 0.1,
#                           normalization = "NONE",
#                           transform = "None",
#                           analysis_method = "LM",
#                           max_significance = 0.25,
#                           fixed_effects = c("group","age","sex","bmi"),
#                           reference = c("group,MC","sex,Female"))




####---------------  Pathway Enrichment analysis: Cross-sectional_aMC_HC  ----------------------#####
sig <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/significant_results.tsv")
sig <- sig %>% filter(value=="Control") %>% select(1,4,5,6,8,9) 
sig$coef <- -(sig$coef)
sig <- sig %>% mutate(t = coef/stderr)
sig$feature<- substr(sig$feature, 2,6)

pathway <- index[, c(3,6)]
pathway$COMP_ID <- as.character(pathway$COMP_ID)
chemical <- index[, c(3,11)]
chemical$COMP_ID <- as.character(chemical$COMP_ID)

sig <- sig %>% left_join(pathway, by = c("feature" = "COMP_ID"))
sig <- sig %>% left_join(chemical, by = c("feature" = "COMP_ID"))
sig <- sig %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

# Only keep pathways with 3 or more metabolites
sig <- sig %>% group_by(pathway) %>% filter(n()>2) 

path <- sig[, c(1,7,8)] %>% filter(!is.na(pathway))

# Calculate the median t-statistics for each pathway
path_medians <- path %>%
  group_by(pathway) %>%
  summarise(median_t = median(t))

path_medians$positive <- ifelse(path_medians$median_t > 0, 1, 0)
path_medians$positive <- as.factor(path_medians$positive)

# Reorder the levels of the "pathway" factor based on median t-statistics
path$order <- factor(path$pathway, levels = path_medians$pathway[order(path_medians$median_t, decreasing = FALSE)])
path <- left_join(path, path_medians, by = "pathway")

path_counts <- path %>% 
  group_by(pathway) %>% 
  summarize(n = n())
path <- left_join(path, path_counts, by = "pathway")

path <- path %>%
  mutate(pathway_n = paste0(pathway, " (n = ", n, ")"))

path <- path %>% arrange(median_t)

path %>% filter(positive==0) %>% distinct(pathway)

# boxplot
bx_path_aMC_HC <- ggboxplot(path, y = "t", x = "pathway_n",
                            xlab = "Metabolite class",
                            ylab = "t-statistics", 
                            fill = "positive",
                            font.x = 12, font.y = 12, font.tickslab = 11, 
                            legend.position = "bottom") +  
  theme(axis.text.x = element_text(hjust = 0.5), axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 1.75, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -1.75, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("0" = "#189AF9", "1" = "#900C3F"),
                    labels = c("Enriched in Control without diarrhea", "Enriched in MC"))+
  coord_flip()

bx_path_aMC_HC


# ggsave(filename = "fig_bx_path_aMC_HC.pdf",
#        plot=bx_path_aMC_HC, units = "in", height = 12, width = 10, dpi=1200)



####--------------- Pathway Enrichment analysis: Cross-sectional_aMC_HC (All)  ----------------------#####
all <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/all_results.tsv")
all <- all %>% filter(value=="Control") %>% select(1,4,5,6,8,9) 
all$coef <- -(all$coef)
all <- all %>% mutate(t = coef/stderr)
all$feature<- substr(all$feature, 2,6)

all <- all %>% left_join(pathway, by = c("feature" = "COMP_ID"))
all <- all %>% left_join(chemical, by = c("feature" = "COMP_ID"))
all <- all %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

# Only keep pathways with 3 or more metabolites
all <- all %>% group_by(pathway) %>% filter(n()>2) 

path_all <- all[, c(1,7,8)] %>% filter(!is.na(pathway))

# Calculate the median t-statistics for each pathway
path_medians_all <- path_all %>%
  group_by(pathway) %>%
  summarise(median_t = median(t))

path_medians_all$positive <- ifelse(path_medians_all$median_t > 0, 1, 0)
path_medians_all$positive <- as.factor(path_medians_all$positive)

# Reorder the levels of the "pathway" factor based on median t-statistics
path_all$order <- factor(path_all$pathway, levels = path_medians_all$pathway[order(path_medians_all$median_t, decreasing = FALSE)])
path_all <- left_join(path_all, path_medians_all, by = "pathway")

path_counts_all <- path_all %>% 
  group_by(pathway) %>% 
  summarize(n = n())
path_all <- left_join(path_all, path_counts_all, by = "pathway")

path_all <- path_all %>%
  mutate(pathway_n = paste0(pathway, " (n = ", n, ")"))

path_all <- path_all %>% arrange(median_t)

# boxplot
bx_path_all_aMC_HC <- ggboxplot(path_all, y = "t", x = "pathway_n",
                                xlab = "Metabolite class",
                                ylab = "t-statistics", 
                                fill = "positive",
                                font.x = 12, font.y = 12, font.tickslab = 11, 
                                legend.position = "bottom") +  # Position legend at the bottom
  theme(axis.text.x = element_text(angle = 80, hjust = 1), axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 1.75, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -1.75, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("0" = "#189AF9", "1" = "#900C3F"),
                    labels = c("Enriched in Control without diarrhea", "Enriched in MC"))

bx_path_all_aMC_HC


# ggsave(filename = "fig_bx_path_all_aMC_HC.png",
#        plot=bx_path_all_aMC_HC, units = "in", height = 10, width = 20, dpi=300)




## For multi-omics association
met_cross_residuals <- readRDS("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/fits/residuals.rds") %>% 
  as.data.frame() %>% t()
# Remove "X" from column names
colnames(met_cross_residuals) <- gsub("^X", "", colnames(met_cross_residuals))

# Only select metabolites with at least abundance >= 0.001, prevalence >= 0.1
all_1 <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/all_results.tsv")
all_1 <- all_1 %>% filter(value=="Control") %>% select(1,4,5,6,8,9) # 1377 metabolites have at least abundance >= 0.001, prevalence >= 0.1
all_1$feature<- substr(all_1$feature, 2,6)

all_1 <- all_1 %>% left_join(pathway, by = c("feature" = "COMP_ID"))
all_1 <- all_1 %>% left_join(chemical, by = c("feature" = "COMP_ID"))
all_1 <- all_1 %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

all_1 <- all_1 %>% filter(!is.na(pathway)) # 1042 known metabolites have at least abundance >= 0.001, prevalence >= 0.1
features <- all_1$feature
met_cross_residuals <- met_cross_residuals[, colnames(met_cross_residuals) %in% features]

# Create a named vector for mapping
name_map <- setNames(all_1$metabolite, all_1$feature)

# Match and change the column names
colnames(met_cross_residuals) <- name_map[colnames(met_cross_residuals)]
met_cross_residuals <- met_cross_residuals %>% as.data.frame()

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_residuals.xlsx"
# write.xlsx(met_cross_residuals, file_path, rowNames = TRUE)


# met_cross_residuals$id <- rownames(met_cross_residuals)
# 
# MC <- cross_baseline %>% filter(group=="MC")
# CD <- cross_baseline %>% filter(group=="Diarrhea")
# HC <- cross_baseline %>% filter(group=="Control")
# 
# met_MC_halla <- met_cross_residuals %>% inner_join(select(MC, id))
# rownames(met_MC_halla) <- met_MC_halla$id
# met_MC_halla$id <- NULL
# # file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_MC_halla.xlsx"
# # write.xlsx(met_MC_halla, file_path, rowNames = TRUE)
# 
# met_CD_halla <- met_cross_residuals %>% inner_join(select(CD, id))
# rownames(met_CD_halla) <- met_CD_halla$id
# met_CD_halla$id <- NULL
# # file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_CD_halla.xlsx"
# # write.xlsx(met_CD_halla, file_path, rowNames = TRUE)
# 
# met_HC_halla <- met_cross_residuals %>% inner_join(select(HC, id))
# rownames(met_HC_halla) <- met_HC_halla$id
# met_HC_halla$id <- NULL
# # file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_HC_halla.xlsx"
# # write.xlsx(met_HC_halla, file_path, rowNames = TRUE)

# met_cross_residuals$id <- NULL


####---------------  HAllA for differentially abundant metabolite and species in MC compared to both CD and HC  ----------------------#####
target <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/halla_target.xlsx")
target <- target$metabolite

met_halla_residual <- met_cross_residuals %>% select(all_of(target))
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_halla_residual.xlsx"
# write.xlsx(met_halla_residual, file_path, rowNames = TRUE)

# met_MC_halla_target <- met_MC_halla %>% select(all_of(target))
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_MC_halla_target.xlsx"
# write.xlsx(met_MC_halla_target, file_path, rowNames = TRUE)

# met_CD_halla_target <- met_CD_halla %>% select(all_of(target))
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_CD_halla_target.xlsx"
# write.xlsx(met_CD_halla_target, file_path, rowNames = TRUE)

# met_HC_halla_target <- met_HC_halla %>% select(all_of(target))
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_HC_halla_target.xlsx"
# write.xlsx(met_HC_halla_target, file_path, rowNames = TRUE)




####---------------  Pathway Enrichment analysis: Cross-sectional_aMC_CD  ----------------------#####
sig_aMC_CD <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/significant_results.tsv")
sig_aMC_CD <- sig_aMC_CD %>% filter(value=="Diarrhea") %>% select(1,4,5,6,8,9) 
sig_aMC_CD$coef <- -(sig_aMC_CD$coef)
sig_aMC_CD <- sig_aMC_CD %>% mutate(t = coef/stderr)
sig_aMC_CD$feature<- substr(sig_aMC_CD$feature, 2,6)

sig_aMC_CD <- sig_aMC_CD %>% left_join(pathway, by = c("feature" = "COMP_ID"))
sig_aMC_CD <- sig_aMC_CD %>% left_join(chemical, by = c("feature" = "COMP_ID"))
sig_aMC_CD <- sig_aMC_CD %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

# Only keep pathways with 3 or more metabolites
sig_aMC_CD <- sig_aMC_CD %>% group_by(pathway) %>% filter(n()>2) 

path_aMC_CD <- sig_aMC_CD[, c(1,7,8)] %>% filter(!is.na(pathway))

# Calculate the median t-statistics for each pathway
path_medians_aMC_CD <- path_aMC_CD %>%
  group_by(pathway) %>%
  summarize(median_t = median(t))

path_medians_aMC_CD$positive <- ifelse(path_medians_aMC_CD$median_t > 0, 1, 0)
path_medians_aMC_CD$positive <- as.factor(path_medians_aMC_CD$positive)

# Reorder the levels of the "pathway" factor based on median t-statistics
path_aMC_CD$order <- factor(path_aMC_CD$pathway, levels = path_medians_aMC_CD$pathway[order(path_medians_aMC_CD$median_t, decreasing = FALSE)])
path_aMC_CD <- left_join(path_aMC_CD, path_medians_aMC_CD, by = "pathway")

path_counts_aMC_CD <- path_aMC_CD %>% 
  group_by(pathway) %>% 
  summarize(n = n())
path_aMC_CD <- left_join(path_aMC_CD, path_counts_aMC_CD, by = "pathway")

path_aMC_CD <- path_aMC_CD %>%
  mutate(pathway_n = paste0(pathway, " (n = ", n, ")"))

path_aMC_CD <- path_aMC_CD %>% arrange(median_t)

path_aMC_CD %>% filter(positive==0) %>% distinct(pathway)

# boxplot
bx_path_aMC_CD <- ggboxplot(path_aMC_CD, y = "t", x = "pathway_n",
                            xlab = "Metabolite class",
                            ylab = "t-statistics", 
                            fill = "positive",
                            font.x = 12, font.y = 12, font.tickslab = 11, 
                            legend.position = "bottom") +  # Position legend at the bottom
  theme(axis.text.x = element_text(hjust = 0.5), axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 1.75, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -1.75, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("0" = "#F9D937", "1" = "#900C3F"),
                    labels = c("Enriched in Chronic diarrhea", "Enriched in MC"))+
  coord_flip()

bx_path_aMC_CD


# ggsave(filename = "fig_bx_path_aMC_CD.pdf",
#        plot=bx_path_aMC_CD, units = "in", height = 8, width = 8, dpi=1200)




####--------------- Pathway Enrichment analysis: Cross-sectional_aMC_CD (All)  ----------------------#####
all_aMC_CD <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/all_results.tsv")
all_aMC_CD <- all_aMC_CD %>% filter(value=="Diarrhea") %>% select(1,4,5,6,8,9) 
all_aMC_CD$coef <- -(all_aMC_CD$coef)
all_aMC_CD <- all_aMC_CD %>% mutate(t = coef/stderr)
all_aMC_CD$feature<- substr(all_aMC_CD$feature, 2,6)

all_aMC_CD <- all_aMC_CD %>% left_join(pathway, by = c("feature" = "COMP_ID"))
all_aMC_CD <- all_aMC_CD %>% left_join(chemical, by = c("feature" = "COMP_ID"))
all_aMC_CD <- all_aMC_CD %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

# Only keep pathways with 3 or more metabolites
all_aMC_CD <- all_aMC_CD %>% group_by(pathway) %>% filter(n()>2) 

path_all_aMC_CD <- all_aMC_CD[, c(1,7,8)] %>% filter(!is.na(pathway))

# Calculate the median t-statistics for each pathway
path_medians_all_aMC_CD <- path_all_aMC_CD %>%
  group_by(pathway) %>%
  summarise(median_t = median(t))

path_medians_all_aMC_CD$positive <- ifelse(path_medians_all_aMC_CD$median_t > 0, 1, 0)
path_medians_all_aMC_CD$positive <- as.factor(path_medians_all_aMC_CD$positive)

# Reorder the levels of the "pathway" factor based on median t-statistics
path_all_aMC_CD$order <- factor(path_all_aMC_CD$pathway, levels = path_medians_all_aMC_CD$pathway[order(path_medians_all_aMC_CD$median_t, decreasing = FALSE)])
path_all_aMC_CD <- left_join(path_all_aMC_CD, path_medians_all_aMC_CD, by = "pathway")

path_counts_all_aMC_CD <- path_all_aMC_CD %>% 
  group_by(pathway) %>% 
  summarize(n = n())
path_all_aMC_CD <- left_join(path_all_aMC_CD, path_counts_all_aMC_CD, by = "pathway")

path_all_aMC_CD <- path_all_aMC_CD %>%
  mutate(pathway_n = paste0(pathway, " (n = ", n, ")"))

path_all_aMC_CD <- path_all_aMC_CD %>% arrange(median_t)

# boxplot
bx_path_all_aMC_CD <- ggboxplot(path_all_aMC_CD, y = "t", x = "pathway_n",
                                xlab = "Metabolite class",
                                ylab = "t-statistics", 
                                fill = "positive",
                                font.x = 12, font.y = 12, font.tickslab = 11, 
                                legend.position = "bottom") +  # Position legend at the bottom
  theme(axis.text.x = element_text(angle = 80, hjust = 1), axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 1.75, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -1.75, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("0" = "#F9D937", "1" = "#900C3F"),
                    labels = c("Enriched in Chronic diarrhea", "Enriched in MC"))

bx_path_all_aMC_CD


# ggsave(filename = "fig_bx_path_all_aMC_CD.png",
#        plot=bx_path_all_aMC_CD, units = "in", height = 10, width = 20, dpi=300)





####---------------  Find metabolite pathways different in both MC vs HC & MC vs CD ----------------------#####
more_mc.vs.hc <- path_medians %>% 
  filter(positive==1) %>% 
  rename("median_t_mc.vs.hc" = "median_t")

less_mc.vs.hc <- path_medians %>% 
  filter(positive==0) %>% 
  rename("median_t_mc.vs.hc" = "median_t")

more_mc.vs.cd <- path_medians_aMC_CD %>% 
  filter(positive==1) %>% 
  rename("median_t_mc.vs.cd" = "median_t")

less_mc.vs.cd <- path_medians_aMC_CD %>% 
  filter(positive==0) %>% 
  rename("median_t_mc.vs.cd" = "median_t")

## Increased species in both MC vs HC & MC vs CD
more_both <- inner_join(more_mc.vs.hc, more_mc.vs.cd)

## Decreased species in both MC vs HC & MC vs CD
less_both <- inner_join(less_mc.vs.hc, less_mc.vs.cd)

list_class <- c(more_both$pathway, less_both$pathway)

## Box plot for MC vs HC
path_both <- path %>% filter(pathway %in% list_class)

bx_aMC_HC <- ggboxplot(path_both, y = "t", x = "pathway_n",
                       xlab = "Metabolite class",
                       ylab = "t-statistics", 
                       fill = "positive",
                       font.x = 12, font.y = 12, font.tickslab = 11, 
                       legend.position = "bottom") +  
  theme(axis.text.x = element_text(hjust = 0.5), axis.text.y = element_text(size = 10)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 1.75, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -1.75, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("0" = "#189AF9", "1" = "#900C3F"),
                    labels = c("Enriched in Control without diarrhea", "Enriched in MC"))+
  coord_flip()

bx_aMC_HC


# ggsave(filename = "bx_aMC_HC.pdf",
#        plot=bx_aMC_HC, units = "in", height = 6, width = 6.5, dpi=1200)


## Box plot for MC vs CD
path_both_CD <- path_aMC_CD %>% filter(pathway %in% list_class)
path_both_CD$pathway_n <- factor(path_both_CD$pathway_n, levels = c(
  "Fatty Acid, Dicarboxylate (n = 5)", 
  "Sterol (n = 3)", 
  "Androgenic Steroids (n = 15)",
  "Tryptophan Metabolism (n = 5)",
  "Food Component/Plant (n = 12)",
  "Methionine, Cysteine, SAM and Taurine Metabolism (n = 3)",
  "Tyrosine Metabolism (n = 4)",
  "Pyrimidine Metabolism, Uracil containing (n = 3)",
  "Phospholipid Metabolism (n = 3)",
  "Monoacylglycerol (n = 3)",
  "Dipeptide (n = 5)",
  "Lysophospholipid (n = 15)",
  "Ceramides (n = 5)",
  "Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated) (n = 3)",
  "Lysoplasmalogen (n = 8)",
  "Long Chain Polyunsaturated Fatty Acid (n3 and n6) (n = 7)",
  "Lactosylceramides (LCER) (n = 3)"
))

bx_aMC_CD <- ggboxplot(path_both_CD, y = "t", x = "pathway_n",
                       xlab = "Metabolite class",
                       ylab = "t-statistics", 
                       fill = "positive",
                       font.x = 12, font.y = 12, font.tickslab = 11, 
                       legend.position = "bottom") + 
  theme(axis.text.x = element_text(hjust = 0.5), axis.text.y = element_text(size = 10)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 1.75, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -1.75, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("0" = "#F9D937", "1" = "#900C3F"),
                    labels = c("Enriched in Chronic diarrhea", "Enriched in MC"))+
  coord_flip()

bx_aMC_CD


# ggsave(filename = "bx_aMC_CD.pdf",
#        plot=bx_aMC_CD, units = "in", height = 6, width = 6.5, dpi=1200)





####---------------  Figure PCoA_longitudinal  ----------------------#####
# Select MC with both active and remission data
metabolite_long <- metabolite %>% filter(group=="MC") %>% group_by(ID) %>% filter(n()>1)
metabolite_long <- metabolite_long %>%
  group_by(ID) %>%
  filter(
    any(symptoms == "active_diarrhea") & any(symptoms == "no_diarrhea")
  )
metabolite_long <- as.data.frame(metabolite_long)
rownames(metabolite_long) <- metabolite_long$id
sum(metabolite_long$symptoms == "active_diarrhea") 
sum(metabolite_long$symptoms == "no_diarrhea") # MC with both active and remission and both mgx and mbx = 40 patients, 80 samples


### Temporarily remove MC151, MC167, MC232, MC244, MC385, MC401, MC471 for longitudinal study because they either don't have both data or have active/active
metabolite_long <- metabolite_long %>% filter(!ID=="MC151" & !ID=="MC167" & !ID=="MC232" & !ID=="MC244" & !ID=="MC385" & !ID=="MC401" & !ID=="MC471")



# Create metadata for longitudinal study
metadata_long <- metadata %>% inner_join(select(metabolite_long, id)) 
metadata_long <- as.data.frame(metadata_long)
rownames(metadata_long) <- metadata_long$id

metadata_long <- metadata_long %>% left_join(select(mc_cd_hc, id, bmi), by = c("ID" = "id"))

subset.metabolite_long <- metabolite_long[, !colnames(metabolite_long) %in% c("id", "ID", "group","disease","symptoms")]

bray_long <- vegdist(subset.metabolite_long, "bray")  # pairwise dissmilarity: 0=same; 1=maximally dissimilar
pcoa_long <- cmdscale(bray_long, eig = T)
pcoap_long <- data.frame(pcoa_long$points)
pcoap_long$ID <- rownames(pcoap_long)
pcoap_long <- pcoap_long[order(pcoap_long$ID),]

# percentages of variation explained by PCO1 & 2
eigs_long <- pcoa_long$eig
pc1.pct_long <- eigs_long[1]/sum(eigs_long) 
pc1.pct_long # pc1.pct_long = 17.4%
pc2.pct_long <- eigs_long[2]/sum(eigs_long) 
pc2.pct_long # pc2.pct_long = 8.8%

metadata_long.pcoa <- left_join(metadata_long, pcoap_long, by = c("id" = "ID"))

fig_PCoA_metabolite_long <- ggscatter(metadata_long.pcoa, y = "X2", x = "X1",
                                 ylab = paste0('PCo2 (', round(pc2.pct_long * 100, 1), '%)'),
                                 xlab = paste0('PCo1 (', round(pc1.pct_long * 100, 1), '%)'),
                                 font.x = 12, font.y = 12, font.tickslab = 11,
                                 color = "symptoms", 
                                 size = 2,
                                 legend = "bottom", legend.title = "Legend Title") +
  stat_ellipse(aes(color = symptoms), linetype = 2) +
  guides(color = guide_legend(title = "MC activity")) +
  scale_color_manual(values = c(active_diarrhea = "#900C3F", 
                                no_diarrhea = "antiquewhite4"),
                     labels = c(active_diarrhea = "active", 
                                no_diarrhea = "remission"))

fig_PCoA_metabolite_long

# ggsave(filename = "fig_PCoA_metabolite_long.pdf",
#        plot=fig_PCoA_metabolite_long, units = "in", height = 5, width = 7, dpi=1200)


####---------------  Spearman correlation between the first axes of ordination for metabolite and microbiome (Longitudinal)  ----------------------#####
## Longitudinal
mgx.pcoa.long <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/alpha.pcoa_long.xlsx")
corr_long_1 <- metadata_long.pcoa %>% select(id, X1, symptoms)
corr_long_1 <- corr_long_1 %>% rename("mbx_X1" = "X1")

corr_long_2 <- mgx.pcoa.long %>% select(ID, X1)
corr_long_2 <- corr_long_2 %>% 
  rename("id" = "ID",
    "mgx_X1" = "X1")

corr_long <- left_join(corr_long_1, corr_long_2) %>% filter(!is.na(mbx_X1) & !is.na(mgx_X1)) # 7 pairs temporarily removed because they are both "active" (according to microbiome data), wait for Jenny's data

ggscatter(corr_long, 
          x = "mgx_X1", 
          y = "mbx_X1", 
          color = "symptoms",
          palette = c("active_diarrhea" = "#900C3F", "no_diarrhea" = "antiquewhite4")) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "The first axis of ordination of the microbiome", 
       y = "The first axis of ordination of the metabolome",
       color = "MC activity") +
  scale_color_manual(labels = c("active_diarrhea" = "Active", "no_diarrhea" = "Remission"),
                     values = c("active_diarrhea" = "#900C3F", "no_diarrhea" = "antiquewhite4")) +
  theme_minimal()+
  geom_smooth(method = "lm", se = FALSE, color = "black")

# Calculate Spearman correlation 
spearman_corr_long <- cor.test(corr_long$mbx_X1, corr_long$mgx_X1, method = "spearman")
spearman_corr_long # rho = -0.550, p-value = 2.547e-06






####---------------  MaAsLin: Longitudinal  ----------------------#####
comp_metabolite_maaslin_long <- comp_metabolite %>% select(ID,6:1780)
unique_ids <- unique(metadata_long$ID)
comp_metabolite_maaslin_long <- subset(comp_metabolite_maaslin_long, ID %in% unique_ids)
all(comp_metabolite_maaslin_long$ID %in% metadata_long$ID) #TRUE
all(metadata_long$ID %in% comp_metabolite_maaslin_long$ID) #TRUE
rownames(metadata_long) <- metadata_long$id

comp_metabolite_maaslin_long <- comp_metabolite_maaslin_long %>% select(-ID)



# met_long <- Maaslin2(input_data = comp_metabolite_maaslin_long, 
#                      input_metadata = metadata_long,
#                      output = "met_long",
#                      min_abundance = 0.001,
#                      min_prevalence = 0.1,
#                      normalization = "NONE", 
#                      transform = "None",
#                      analysis_method = "LM",
#                      max_significance = 0.25,
#                      fixed_effects = c("symptoms","age","sex","bmi"),
#                      reference = c("symptoms,no_diarrhea","sex,Female"))





####--------------- Pathway Enrichment analysis: Longitudinal  ----------------------#####
sig_long <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_long/significant_results.tsv")
sig_long <- sig_long %>% filter(value=="active_diarrhea") %>% select(1,4,5,6,8,9) 
sig_long <- sig_long %>% mutate(t = coef/stderr)
sig_long$feature<- substr(sig_long$feature, 2,6)


sig_long <- sig_long %>% left_join(pathway, by = c("feature" = "COMP_ID"))
sig_long <- sig_long %>% left_join(chemical, by = c("feature" = "COMP_ID"))
sig_long <- sig_long %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

# Only keep pathways with 3 or more metabolites
sig_long <- sig_long %>% group_by(pathway) %>% filter(n()>2) 

path_long <- sig_long[, c(1,7,8)] %>% filter(!is.na(pathway))

# Calculate the median t-statistics for each pathway
path_medians_long <- path_long %>%
  group_by(pathway) %>%
  summarise(median_t = median(t))

path_medians_long$positive <- ifelse(path_medians_long$median_t > 0, 1, 0)
path_medians_long$positive <- as.factor(path_medians_long$positive)

# Reorder the levels of the "pathway" factor based on median t-statistics
path_long$order <- factor(path_long$pathway, levels = path_medians_long$pathway[order(path_medians_long$median_t, decreasing = FALSE)])
path_long <- left_join(path_long, path_medians_long, by = "pathway")

path_counts_long <- path_long %>% 
  group_by(pathway) %>% 
  summarize(n = n())
path_long <- left_join(path_long, path_counts_long, by = "pathway")

path_long <- path_long %>%
  mutate(pathway_n = paste0(pathway, " (n = ", n, ")"))

path_long <- path_long %>% arrange(median_t)

# boxplot
bx_path_long <- ggboxplot(path_long, y = "t", x = "pathway_n",
                     xlab = "Pathway",
                     ylab = "t-statistics", 
                     fill = "positive",
                     font.x = 12, font.y = 12, font.tickslab = 11, 
                     legend.position = "bottom") + 
  theme(axis.text.x = element_text(hjust = 0.5), axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 2.64, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -2.64, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("1" = "#900C3F"),
                    labels = c("Increased in active"))+
  coord_flip()

bx_path_long


# ggsave(filename = "fig_bx_path_long.png",
#        plot=bx_path_long, units = "in", height = 4, width = 9, dpi=600)

####--------------- Pathway Enrichment analysis: Longitudinal (All)  ----------------------#####
all_long <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_long/all_results.tsv")
all_long <- all_long %>% filter(value=="active_diarrhea") %>% select(1,4,5,6,8,9) 
all_long <- all_long %>% mutate(t = coef/stderr)
all_long$feature<- substr(all_long$feature, 2,6)


all_long <- all_long %>% left_join(pathway, by = c("feature" = "COMP_ID"))
all_long <- all_long %>% left_join(chemical, by = c("feature" = "COMP_ID"))
all_long <- all_long %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

# Only keep pathways with 3 or more metabolites
all_long <- all_long %>% group_by(pathway) %>% filter(n()>2) 

path_long_all <- all_long[, c(1,7,8)] %>% filter(!is.na(pathway))

# Calculate the median t-statistics for each pathway
path_medians_all_long <- path_long_all %>%
  group_by(pathway) %>%
  summarise(median_t = median(t))

path_medians_all_long$positive <- ifelse(path_medians_all_long$median_t > 0, 1, 0)
path_medians_all_long$positive <- as.factor(path_medians_all_long$positive)

# Reorder the levels of the "pathway" factor based on median t-statistics
path_long_all$order <- factor(path_long_all$pathway, levels = path_medians_all_long$pathway[order(path_medians_all_long$median_t, decreasing = FALSE)])
path_long_all <- left_join(path_long_all, path_medians_all_long, by = "pathway")

path_counts_all_long <- path_long_all %>% 
  group_by(pathway) %>% 
  summarize(n = n())
path_long_all <- left_join(path_long_all, path_counts_all_long, by = "pathway")

path_long_all <- path_long_all %>%
  mutate(pathway_n = paste0(pathway, " (n = ", n, ")"))

path_long_all <- path_long_all %>% arrange(median_t, decreasing = FALSE)

# boxplot
bx_path_long_all <- ggboxplot(path_long_all, y = "t", x = "pathway_n",
                          xlab = "Pathway",
                          ylab = "t-statistics", 
                          fill = "positive",
                          font.x = 12, font.y = 12, font.tickslab = 11, 
                          legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 80, hjust = 1), axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 2.64, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -2.64, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c("0" = "antiquewhite4","1" = "#900C3F"),
                    labels = c("Increased in remission","Increased in active"))

bx_path_long_all


# ggsave(filename = "fig_bx_path_long_all.png",
#        plot=bx_path_long_all, units = "in", height = 10, width = 20, dpi=600)





## For multi-omics association
comp_metabolite_maaslin_long$id <- rownames(comp_metabolite_maaslin_long)
met_halla_long <- comp_metabolite_maaslin_long %>% inner_join(select(metadata_long, id), by = "id")
rownames(met_halla_long) <- met_halla_long$id
# Only select metabolites with at least abundance >= 0.001, prevalence >= 0.1
all_long_1 <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_long/all_results.tsv")
all_long_1 <- all_long %>% filter(!is.na(pathway))
features_long <- all_long_1$feature
met_halla_long <- met_halla_long[, colnames(met_halla_long) %in% features_long]

# Create a named vector for mapping
name_map_long <- setNames(all_long_1$metabolite, all_long_1$feature)

# Match and change the column names in `met_halla_long`
colnames(met_halla_long) <- name_map_long[colnames(met_halla_long)]

# Remove the 'id' columns
comp_metabolite_maaslin_long$id <- NULL
met_halla_long$id <- NULL

# Add row names as a new column
met_halla_long$id <- rownames(met_halla_long)

# Convert to data frame 
met_halla_long <- as.data.frame(met_halla_long)

# Clean column names to ensure they are valid
colnames(met_halla_long) <- make.names(colnames(met_halla_long), unique = TRUE)

# Select the 'id' column and the first 985 features
met_halla_long <- met_halla_long %>% select(id, 1:985)

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_halla_long.xlsx"
# write_xlsx(met_halla_long, file_path)



####--------------- Lactosylceramides analysis: cross-sectional  ----------------------#####
altered_114 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/altered_met_114.xlsx") %>% select(4,5) %>% rename("pathway" = "pathway_HC")
altered_114 %>% distinct(pathway)

lactosylceramides <- subset(altered_114, grepl("Lactosylceramides", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
lactosylceramides_features <- lactosylceramides$feature
lactosylceramides <- metabolite %>% select(all_of(lactosylceramides_features))
lactosylceramides <- lactosylceramides %>% mutate(lactosylceramides = rowSums(across(everything())))
lactosylceramides$id <- rownames(lactosylceramides)

lactosylceramides_aMC <- lactosylceramides %>% inner_join(select(activeMC_baseline, id, disease)) %>% mutate(mc_all_label = "MC")
lactosylceramides_HC <- lactosylceramides %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control") %>% mutate(disease = "Control without diarrhea")
lactosylceramides_CD <- lactosylceramides %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea") %>% mutate(disease = "Chronic diarrhea")
lactosylceramides <- rbind(lactosylceramides_aMC,lactosylceramides_HC,lactosylceramides_CD)
lactosylceramides <- lactosylceramides %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 

lactosylceramides_filled <- lactosylceramides %>% filter(!is.na(bmi))
lactosylceramides_unfilled <- lactosylceramides %>% filter(is.na(bmi))
lactosylceramides_unfilled$ID <- substr(lactosylceramides_unfilled$id, 1, 5)
lactosylceramides_unfilled <- lactosylceramides_unfilled %>% select(-age, -sex, -bmi)
lactosylceramides_unfilled <- lactosylceramides_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
lactosylceramides <- rbind(lactosylceramides_filled, lactosylceramides_unfilled)


lactosylceramides_log <- lactosylceramides
lactosylceramides_log$lactosylceramides <- log(lactosylceramides_log$lactosylceramides)

## Box plot
lactosylceramides_log$mc_all_label <- factor(lactosylceramides_log$mc_all_label, 
                                             levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                             labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_lactosylceramides <- ggplot(lactosylceramides_log, aes(x = mc_all_label, y = lactosylceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Lactosylceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_lactosylceramides

# ggsave(filename = "bx_cross_lactosylceramides.pdf",
#        plot = bx_cross_lactosylceramides, units = "in", width=10, height=4, dpi = 1200)


## Individual metabolite
lactosylceramides_individual <- lactosylceramides %>% select(1,5:10) 
lactosylceramides_individual <- lactosylceramides_individual %>%
  mutate(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` = log(`53010`)) %>%
  select(-`53010`)

lactosylceramides_individual$mc_all_label <- factor(lactosylceramides_individual$mc_all_label, 
                                             levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                             labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_lactosylceramides_individual <- ggplot(lactosylceramides_individual, aes(x = mc_all_label, y = `lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(lactosyl-N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_lactosylceramides_individual

# ggsave(filename = "bx_cross_lactosylceramides_individual.pdf",
#        plot = bx_cross_lactosylceramides_individual, units = "in", width=10, height=4, dpi = 1200)

lactosylceramides_log$mc_all_label <- factor(lactosylceramides_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_log$sex <- factor(lactosylceramides_log$sex)
lm_lactosylceramides_log <- lm(lactosylceramides ~ mc_all_label + age + sex + bmi, data = lactosylceramides_log)
summary(lm_lactosylceramides_log)

lactosylceramides_cross <- lactosylceramides


lactosylceramides_individual$mc_all_label <- factor(lactosylceramides_individual$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_individual$sex <- factor(lactosylceramides_individual$sex)
lm_lactosylceramides_individual <- lm(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` ~ mc_all_label + age + sex + bmi, data = lactosylceramides_individual)
summary(lm_lactosylceramides_individual)


### LC and CC subgroup 
lactosylceramides_log$disease <- factor(lactosylceramides_log$disease, 
                                        levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_log_1 <- lactosylceramides_log %>% filter(!is.na(disease))

bx_cross_lactosylceramides_sub <- ggplot(lactosylceramides_log_1, aes(x = disease, y = lactosylceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Lactosylceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_lactosylceramides_sub

# ggsave(filename = "bx_cross_lactosylceramides_sub.pdf",
#        plot = bx_cross_lactosylceramides_sub, units = "in", width=10, height=5, dpi = 1200)


## Individual metabolite
lactosylceramides_individual$disease <- factor(lactosylceramides_log$disease, 
                                        levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_individual_1 <- lactosylceramides_individual %>% filter(!is.na(disease))

bx_cross_lactosylceramides_sub_individual <- ggplot(lactosylceramides_individual_1, aes(x = disease, y = `lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(lactosyl-N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_lactosylceramides_sub_individual

# ggsave(filename = "bx_cross_lactosylceramides_sub_individual.pdf",
#        plot = bx_cross_lactosylceramides_sub_individual, units = "in", width=10, height=5, dpi = 1200)


lactosylceramides_log_1$disease <- factor(lactosylceramides_log_1$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_log_1$sex <- factor(lactosylceramides_log_1$sex)
lm_lactosylceramides_log_1 <- lm(lactosylceramides ~ disease + age + sex + bmi, data = lactosylceramides_log_1)
summary(lm_lactosylceramides_log_1)

lactosylceramides_log_1$disease <- factor(lactosylceramides_log_1$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_log_1$sex <- factor(lactosylceramides_log_1$sex)
lm_lactosylceramides_log_1 <- lm(lactosylceramides ~ disease + age + sex + bmi, data = lactosylceramides_log_1)
summary(lm_lactosylceramides_log_1)


lactosylceramides_individual_1$disease <- factor(lactosylceramides_individual_1$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_individual_1$sex <- factor(lactosylceramides_individual_1$sex)
lm_lactosylceramides_individual_1 <- lm(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` ~ disease + age + sex + bmi, data = lactosylceramides_individual_1)
summary(lm_lactosylceramides_individual_1)

lactosylceramides_individual_1$disease <- factor(lactosylceramides_individual_1$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
lactosylceramides_individual_1$sex <- factor(lactosylceramides_individual_1$sex)
lm_lactosylceramides_individual_1 <- lm(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` ~ disease + age + sex + bmi, data = lactosylceramides_individual_1)
summary(lm_lactosylceramides_individual_1)




####--------------- Lactosylceramides analysis: longitudinal  ----------------------#####
lactosylceramides <- subset(altered_114, grepl("Lactosylceramides", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
lactosylceramides_features <- lactosylceramides$feature
lactosylceramides <- metabolite %>% select(all_of(lactosylceramides_features))
lactosylceramides <- lactosylceramides %>% mutate(lactosylceramides = rowSums(across(everything())))
lactosylceramides$id <- rownames(lactosylceramides)

lactosylceramides_long <- lactosylceramides %>% inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 


lactosylceramides_long_log <- lactosylceramides_long
lactosylceramides_long_log$lactosylceramides <- log(lactosylceramides_long_log$lactosylceramides)

## Box plot
lactosylceramides_long_log$symptoms <- factor(lactosylceramides_long_log$symptoms, 
                                              levels = c("active_diarrhea", "no_diarrhea"),
                                              labels = c("active", "remission"))

bx_lactosylceramides_long <- ggplot(lactosylceramides_long_log, aes(x = symptoms, y = lactosylceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(Lactosylceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_lactosylceramides_long

# ggsave(filename = "bx_lactosylceramides_long.pdf",
#        plot = bx_lactosylceramides_long, units = "in", width=10, height=4, dpi = 1200)


## Individual metabolite
lactosylceramides_individual_long <- lactosylceramides_long %>% select(1,5:9) 
lactosylceramides_individual_long <- lactosylceramides_individual_long %>%
  mutate(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` = log(`53010`)) %>%
  select(-`53010`)

lactosylceramides_individual_long$symptoms <- factor(lactosylceramides_individual_long$symptoms, 
                                              levels = c("active_diarrhea", "no_diarrhea"),
                                              labels = c("active", "remission"))

bx_lactosylceramides_individual_long <- ggplot(lactosylceramides_individual_long, aes(x = symptoms, y = `lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(lactosyl-N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_lactosylceramides_individual_long

# ggsave(filename = "bx_lactosylceramides_individual_long.pdf",
#        plot = bx_lactosylceramides_individual_long, units = "in", width=10, height=4, dpi = 1200)


lactosylceramides_long_log$symptoms <- factor(lactosylceramides_long_log$symptoms, levels = c("remission", "active"))
lm_lactosylceramides_long_log <- lm(lactosylceramides ~ symptoms + age + sex + bmi, data = lactosylceramides_long_log)
summary(lm_lactosylceramides_long_log)

lactosylceramides_individual_long$symptoms <- factor(lactosylceramides_individual_long$symptoms, levels = c("remission", "active"))
lm_lactosylceramides_individual_long_log <- lm(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` ~ symptoms + age + sex + bmi, data = lactosylceramides_individual_long)
summary(lm_lactosylceramides_individual_long_log)


### LC and CC subgroup 
lactosylceramides_long_log_1 <- lactosylceramides_long_log %>% left_join(select(lactosylceramides_log, id, disease))
lactosylceramides_long_log_1$ID <- substr(lactosylceramides_long_log_1$id, 1, 5)
lactosylceramides_long_log_1$disease <- as.character(lactosylceramides_long_log_1$disease)

lactosylceramides_long_log_1 <- lactosylceramides_long_log_1 %>%
  group_by(ID) %>%
  mutate(disease = ifelse(is.na(disease), first(na.omit(disease)), disease)) %>%
  ungroup()



lactosylceramides_long_log_1$symptoms <- factor(lactosylceramides_long_log$symptoms,
                                                levels = c("active", "remission"))
lactosylceramides_long_log_1_LC <- lactosylceramides_long_log_1 %>% filter(disease == "LC")

bx_lactosylceramides_long_LC <- ggplot(lactosylceramides_long_log_1_LC, aes(x = symptoms, y = lactosylceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(Lactosylceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_lactosylceramides_long_LC
# ggsave(filename = "bx_lactosylceramides_long_LC.pdf",
#        plot = bx_lactosylceramides_long_LC, units = "in", width=10, height=4, dpi = 1200)

lactosylceramides_long_log_1_LC$symptoms <- factor(lactosylceramides_long_log_1_LC$symptoms, levels = c("remission", "active"))
lm_lactosylceramides_long_log_1_LC <- lm(lactosylceramides ~ symptoms + age + sex + bmi, data = lactosylceramides_long_log_1_LC)
summary(lm_lactosylceramides_long_log_1_LC)



lactosylceramides_long_log_1_CC <- lactosylceramides_long_log_1 %>% filter(disease == "CC")

bx_lactosylceramides_long_CC <- ggplot(lactosylceramides_long_log_1_CC, aes(x = symptoms, y = lactosylceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(Lactosylceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_lactosylceramides_long_CC
# ggsave(filename = "bx_lactosylceramides_long_CC.pdf",
#        plot = bx_lactosylceramides_long_CC, units = "in", width=10, height=4, dpi = 1200)


## Individual metabolite
lactosylceramides_long_log_2 <- lactosylceramides_long_log_1 %>% select(1,5:11) 
lactosylceramides_long_log_2 <- lactosylceramides_long_log_2 %>%
  mutate(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` = log(`53010`)) %>%
  select(-`53010`)

lactosylceramides_long_log_2_LC <- lactosylceramides_long_log_2 %>% filter(disease == "LC")

bx_lactosylceramides_individual_long_LC <- ggplot(lactosylceramides_long_log_2_LC, aes(x = symptoms, y = `lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(lactosyl-N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_lactosylceramides_individual_long_LC
# ggsave(filename = "bx_lactosylceramides_individual_long_LC.pdf",
#        plot = bx_lactosylceramides_individual_long_LC, units = "in", width=10, height=4, dpi = 1200)

lactosylceramides_long_log_2_LC$symptoms <- factor(lactosylceramides_long_log_2_LC$symptoms, levels = c("remission", "active"))
lm_lactosylceramides_long_log_2_LC <- lm(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` ~ symptoms + age + sex + bmi, data = lactosylceramides_long_log_2_LC)
summary(lm_lactosylceramides_long_log_2_LC)



lactosylceramides_long_log_2_CC <- lactosylceramides_long_log_2 %>% filter(disease == "CC")

bx_lactosylceramides_individual_long_CC <- ggplot(lactosylceramides_long_log_2_CC, aes(x = symptoms, y = `lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(lactosyl-N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_lactosylceramides_long_CC
# ggsave(filename = "bx_lactosylceramides_individual_long_CC.pdf",
#        plot = bx_lactosylceramides_individual_long_CC, units = "in", width=10, height=4, dpi = 1200)


lactosylceramides_long_log_2_CC$symptoms <- factor(lactosylceramides_long_log_2_CC$symptoms, levels = c("remission", "active"))
lm_lactosylceramides_long_log_2_CC <- lm(`lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)` ~ symptoms + age + bmi, data = lactosylceramides_long_log_2_CC)
summary(lm_lactosylceramides_long_log_2_CC)





####--------------- Ceramides analysis: cross-sectional  ----------------------#####
ceramides <- subset(altered_114, grepl("Ceramides", pathway, ignore.case = TRUE)) 
ceramides <- subset(ceramides, !grepl("LCER|HCER|Dihydroceramides", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))

ceramides_features <- ceramides$feature
ceramides <- metabolite %>% select(all_of(ceramides_features))
ceramides <- ceramides %>% mutate(ceramides = rowSums(across(everything())))
ceramides$id <- rownames(ceramides)

ceramides_aMC <- ceramides %>% inner_join(select(activeMC_baseline, id, disease)) %>% mutate(mc_all_label = "MC")
ceramides_HC <- ceramides %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control") %>% mutate(disease = "Control without diarrhea")
ceramides_CD <- ceramides %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea") %>% mutate(disease = "Chronic diarrhea")
ceramides <- rbind(ceramides_aMC,ceramides_HC,ceramides_CD)
ceramides <- ceramides %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 

ceramides_filled <- ceramides %>% filter(!is.na(bmi))
ceramides_unfilled <- ceramides %>% filter(is.na(bmi))
ceramides_unfilled$ID <- substr(ceramides_unfilled$id, 1, 5)
ceramides_unfilled <- ceramides_unfilled %>% select(-age, -sex, -bmi)
ceramides_unfilled <- ceramides_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
ceramides <- rbind(ceramides_filled, ceramides_unfilled)

ceramides_log <- ceramides
ceramides_log$ceramides <- log(ceramides_log$ceramides)

## Box plot
ceramides_log$mc_all_label <- factor(ceramides_log$mc_all_label, 
                                     levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                     labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_ceramides <- ggplot(ceramides_log, aes(x = mc_all_label, y = ceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Ceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_ceramides

# ggsave(filename = "bx_cross_ceramides.pdf",
#        plot = bx_cross_ceramides, units = "in", width=10, height=4, dpi = 1200)



## Individual metabolite
ceramides_individual <- ceramides %>% select(1,7:12) 
ceramides_individual <- ceramides_individual %>%
  mutate(`N-palmitoyl-sphingosine (d18:1/16:0)` = log(`44877`)) %>%
  select(-`44877`)

## Box plot
ceramides_individual$mc_all_label <- factor(ceramides_individual$mc_all_label, 
                                     levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                     labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_ceramides_individual <- ggplot(ceramides_individual, aes(x = mc_all_label, y = `N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_ceramides_individual

# ggsave(filename = "bx_cross_ceramides_individual.pdf",
#        plot = bx_cross_ceramides_individual, units = "in", width=10, height=4, dpi = 1200)


ceramides_log$mc_all_label <- factor(ceramides_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_ceramides_log <- lm(ceramides ~ mc_all_label + age + sex + bmi, data = ceramides_log)
summary(lm_ceramides_log)

ceramides_cross <- ceramides

ceramides_individual$mc_all_label <- factor(ceramides_individual$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_ceramides_individual_log <- lm(`N-palmitoyl-sphingosine (d18:1/16:0)` ~ mc_all_label + age + sex + bmi, data = ceramides_individual)
summary(lm_ceramides_individual_log)




### LC and CC subgroup 
ceramides_log$disease <- factor(ceramides_log$disease, 
                                levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
ceramides_log_1 <- ceramides_log %>% filter(!is.na(disease))

bx_cross_ceramides_sub <- ggplot(ceramides_log_1, aes(x = disease, y = ceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Ceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_ceramides_sub

# ggsave(filename = "bx_cross_ceramides_sub.pdf",
#        plot = bx_cross_ceramides_sub, units = "in", width=10, height=5, dpi = 1200)


## Individual metabolite
ceramides_individual$disease <- factor(ceramides_individual$disease, 
                                levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
ceramides_individual_1 <- ceramides_individual %>% filter(!is.na(disease))

bx_cross_ceramides_individual_sub <- ggplot(ceramides_individual_1, aes(x = disease, y = `N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_ceramides_individual_sub

# ggsave(filename = "bx_cross_ceramides_individual_sub.pdf",
#        plot = bx_cross_ceramides_individual_sub, units = "in", width=10, height=5, dpi = 1200)


ceramides_log_1$disease <- factor(ceramides_log_1$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
ceramides_log_1$sex <- factor(ceramides_log_1$sex)
lm_ceramides_log_1 <- lm(ceramides ~ disease + age + sex + bmi, data = ceramides_log_1)
summary(lm_ceramides_log_1)

ceramides_log_1$disease <- factor(ceramides_log_1$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
ceramides_log_1$sex <- factor(ceramides_log_1$sex)
lm_ceramides_log_1 <- lm(ceramides ~ disease + age + sex + bmi, data = ceramides_log_1)
summary(lm_ceramides_log_1)


ceramides_individual_1$disease <- factor(ceramides_individual_1$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
ceramides_individual_1$sex <- factor(ceramides_individual_1$sex)
lm_ceramides_individual_1 <- lm(`N-palmitoyl-sphingosine (d18:1/16:0)` ~ disease + age + sex + bmi, data = ceramides_individual_1)
summary(lm_ceramides_individual_1)

ceramides_individual_1$disease <- factor(ceramides_individual_1$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
ceramides_individual_1$sex <- factor(ceramides_individual_1$sex)
lm_ceramides_individual_1 <- lm(`N-palmitoyl-sphingosine (d18:1/16:0)` ~ disease + age + sex + bmi, data = ceramides_individual_1)
summary(lm_ceramides_individual_1)





####--------------- Ceramides analysis: longitudinal  ----------------------#####
ceramides <- subset(altered_114, grepl("Ceramides", pathway, ignore.case = TRUE)) 
ceramides <- subset(ceramides, !grepl("LCER|HCER|Dihydroceramides", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))

ceramides_features <- ceramides$feature
ceramides <- metabolite %>% select(all_of(ceramides_features))
ceramides <- ceramides %>% mutate(ceramides = rowSums(across(everything())))
ceramides$id <- rownames(ceramides)

ceramides_long <- ceramides %>% inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 


ceramides_long_log <- ceramides_long
ceramides_long_log$ceramides <- log(ceramides_long_log$ceramides)

## Box plot
ceramides_long_log$symptoms <- factor(ceramides_long_log$symptoms, 
                                      levels = c("active_diarrhea", "no_diarrhea"),
                                      labels = c("active", "remission"))

bx_ceramides_long <- ggplot(ceramides_long_log, aes(x = symptoms, y = ceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(Ceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_ceramides_long

# ggsave(filename = "bx_ceramides_long.pdf",
#        plot = bx_ceramides_long, units = "in", width=10, height=4, dpi = 1200)

ceramides_long_log$symptoms <- factor(ceramides_long_log$symptoms, levels = c("remission", "active"))
lm_ceramides_long_log <- lm(ceramides ~ symptoms + age + sex + bmi, data = ceramides_long_log)
summary(lm_ceramides_long_log)


## Individual metabolite
ceramides_individual_long <- ceramides_long %>% select(1,7:11) 
ceramides_individual_long <- ceramides_individual_long %>%
  mutate(`N-palmitoyl-sphingosine (d18:1/16:0)` = log(`44877`)) %>%
  select(-`44877`)

## Box plot
ceramides_individual_long$symptoms <- factor(ceramides_individual_long$symptoms, 
                                      levels = c("active_diarrhea", "no_diarrhea"),
                                      labels = c("active", "remission"))

bx_ceramides_individual_long <- ggplot(ceramides_individual_long, aes(x = symptoms, y = `N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_ceramides_individual_long

# ggsave(filename = "bx_ceramides_individual_long.pdf",
#        plot = bx_ceramides_individual_long, units = "in", width=10, height=4, dpi = 1200)

ceramides_individual_long$symptoms <- factor(ceramides_individual_long$symptoms, levels = c("remission", "active"))
lm_ceramides_individual_long_log <- lm(`N-palmitoyl-sphingosine (d18:1/16:0)` ~ symptoms + age + sex + bmi, data = ceramides_individual_long)
summary(lm_ceramides_individual_long_log)



### LC and CC subgroup 
ceramides_long_log_1 <- ceramides_long_log %>% left_join(select(ceramides_log, id, disease))
ceramides_long_log_1$ID <- substr(ceramides_long_log_1$id, 1, 5)
ceramides_long_log_1$disease <- as.character(ceramides_long_log_1$disease)

ceramides_long_log_1 <- ceramides_long_log_1 %>%
  group_by(ID) %>%
  mutate(disease = ifelse(is.na(disease), first(na.omit(disease)), disease)) %>%
  ungroup()



ceramides_long_log_1$symptoms <- factor(ceramides_long_log$symptoms,
                                                levels = c("active", "remission"))
ceramides_long_log_1_LC <- ceramides_long_log_1 %>% filter(disease == "LC")

bx_ceramides_long_LC <- ggplot(ceramides_long_log_1_LC, aes(x = symptoms, y = ceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(Ceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_ceramides_long_LC
# ggsave(filename = "bx_ceramides_long_LC.pdf",
#        plot = bx_ceramides_long_LC, units = "in", width=10, height=4, dpi = 1200)

ceramides_long_log_1_LC$symptoms <- factor(ceramides_long_log_1_LC$symptoms, levels = c("remission", "active"))
lm_ceramides_long_log_1_LC <- lm(ceramides ~ symptoms + age + sex + bmi, data = ceramides_long_log_1_LC)
summary(lm_ceramides_long_log_1_LC)



ceramides_long_log_1_CC <- ceramides_long_log_1 %>% filter(disease == "CC")

bx_ceramides_long_CC <- ggplot(ceramides_long_log_1_CC, aes(x = symptoms, y = ceramides)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(Ceramides)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_ceramides_long_CC
# ggsave(filename = "bx_ceramides_long_CC.pdf",
#        plot = bx_ceramides_long_CC, units = "in", width=10, height=4, dpi = 1200)

ceramides_long_log_1_CC$symptoms <- factor(ceramides_long_log_1_CC$symptoms, levels = c("remission", "active"))
lm_ceramides_long_log_1_CC <- lm(ceramides ~ symptoms + age + bmi, data = ceramides_long_log_1_CC)
summary(lm_ceramides_long_log_1_CC)


## Individual metabolite
ceramides_long_log_2 <- ceramides_long_log_1 %>% select(1,7:13) 
ceramides_long_log_2 <- ceramides_long_log_2 %>%
  mutate(`N-palmitoyl-sphingosine (d18:1/16:0)` = log(`44877`)) %>%
  select(-`44877`)

ceramides_long_log_2$symptoms <- factor(ceramides_long_log_2$symptoms,
                                        levels = c("active", "remission"))
ceramides_long_log_2_LC <- ceramides_long_log_2 %>% filter(disease == "LC")

bx_ceramides_individual_long_LC <- ggplot(ceramides_long_log_2_LC, aes(x = symptoms, y = `N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_ceramides_individual_long_LC
# ggsave(filename = "bx_ceramides_individual_long_LC.pdf",
#        plot = bx_ceramides_individual_long_LC, units = "in", width=10, height=4, dpi = 1200)

ceramides_long_log_2_LC$symptoms <- factor(ceramides_long_log_2_LC$symptoms, levels = c("remission", "active"))
lm_ceramides_long_log_2_LC <- lm(`N-palmitoyl-sphingosine (d18:1/16:0)` ~ symptoms + age + sex + bmi, data = ceramides_long_log_2_LC)
summary(lm_ceramides_long_log_2_LC)



ceramides_long_log_2_CC <- ceramides_long_log_2 %>% filter(disease == "CC")

bx_ceramides_individual_long_CC <- ggplot(ceramides_long_log_2_CC, aes(x = symptoms, y = `N-palmitoyl-sphingosine (d18:1/16:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(N-palmitoyl-sphingosine (d18:1/16:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_ceramides_individual_long_CC
# ggsave(filename = "bx_ceramides_individual_long_CC.pdf",
#        plot = bx_ceramides_individual_long_CC, units = "in", width=10, height=4, dpi = 1200)

ceramides_long_log_2_CC$symptoms <- factor(ceramides_long_log_2_CC$symptoms, levels = c("remission", "active"))
lm_ceramides_long_log_2_CC <- lm(`N-palmitoyl-sphingosine (d18:1/16:0)` ~ symptoms + age + bmi, data = ceramides_long_log_2_CC)
summary(lm_ceramides_long_log_2_CC)





####--------------  Sphingosine and palmitate analysis: cross-sectional ----------------####
## Sphingosine
sphingosine <- metabolite %>% 
  select("17747")%>%
  rename("sphingosine" = "17747")%>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(lactosylceramides_log, id, mc_all_label, age, sex, bmi)) 

sphingosine_log <- sphingosine
sphingosine_log$sphingosine <- log(sphingosine_log$sphingosine)

## Box plot
sphingosine_log$mc_all_label <- factor(sphingosine_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_sphingosines <- ggplot(sphingosine_log, aes(x = mc_all_label, y = sphingosine)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Sphingosine)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_sphingosines

# ggsave(filename = "bx_cross_sphingosines.pdf",
#        plot = bx_cross_sphingosines, units = "in", width=10, height=4, dpi = 1200)


sphingosine_log$mc_all_label <- factor(sphingosine_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_sphingosine_log <- lm(sphingosine ~ mc_all_label + age + sex + bmi, data = sphingosine_log)
summary(lm_sphingosine_log)



## Palmitate
palmitate <- metabolite %>% 
  select("1336")%>%
  rename("palmitate" = "1336")%>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(lactosylceramides_log, id, mc_all_label, age, sex, bmi)) 

palmitate_log <- palmitate
palmitate_log$palmitate <- log(palmitate_log$palmitate)

## Box plot
palmitate_log$mc_all_label <- factor(palmitate_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_palmitate <- ggplot(palmitate_log, aes(x = mc_all_label, y = palmitate)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Palmitate)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_palmitate

# ggsave(filename = "bx_cross_palmitate.pdf",
#        plot = bx_cross_palmitate, units = "in", width=10, height=4, dpi = 1200)


palmitate_log$mc_all_label <- factor(palmitate_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_palmitate_log <- lm(palmitate ~ mc_all_label + age + sex + bmi, data = palmitate_log)
summary(lm_palmitate_log)




####--------------  Sphingosine and palmitate analysis: longitudinal ----------------####
## Sphingosine
sphingosine_long <- metabolite %>% 
  select("17747")%>%
  rename("sphingosine" = "17747")%>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 


sphingosine_long_log <- sphingosine_long
sphingosine_long_log$sphingosine <- log(sphingosine_long_log$sphingosine)

## Box plot
sphingosine_long_log$symptoms <- factor(sphingosine_long_log$symptoms, 
                                      levels = c("active_diarrhea", "no_diarrhea"),
                                      labels = c("active", "remission"))

bx_sphingosine_long <- ggplot(sphingosine_long_log, aes(x = symptoms, y = sphingosine)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(Sphingosine)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_sphingosine_long

# ggsave(filename = "bx_sphingosine_long.pdf",
#        plot = bx_sphingosine_long, units = "in", width=10, height=4, dpi = 1200)

sphingosine_long_log$symptoms <- factor(sphingosine_long_log$symptoms, levels = c("remission", "active"))
lm_sphingosine_long_log <- lm(sphingosine ~ symptoms + age + sex + bmi, data = sphingosine_long_log)
summary(lm_sphingosine_long_log)




## Palmitate
palmitate_long <- metabolite %>% 
  select("1336")%>%
  rename("palmitate" = "1336")%>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 

palmitate_long_log <- palmitate_long
palmitate_long_log$palmitate <- log(palmitate_long_log$palmitate)

## Box plot
palmitate_long_log$symptoms <- factor(palmitate_long_log$symptoms, 
                                      levels = c("active_diarrhea", "no_diarrhea"),
                                      labels = c("active", "remission"))

bx_palmitate_long <- ggplot(palmitate_long_log, aes(x = symptoms, y = palmitate)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(Palmitate)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_palmitate_long

# ggsave(filename = "bx_palmitate_long.pdf",
#        plot = bx_palmitate_long, units = "in", width=10, height=4, dpi = 1200)

palmitate_long_log$symptoms <- factor(palmitate_long_log$symptoms, levels = c("remission", "active"))
lm_palmitate_long_log <- lm(palmitate ~ symptoms + age + sex + bmi, data = palmitate_long_log)
summary(lm_palmitate_long_log)





####--------------- Lysophospholipids analysis: cross-sectional  ----------------------#####
Lysophospholipid <- subset(altered_114, grepl("Lysophospholipid", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
Lysophospholipid_features <- Lysophospholipid$feature
Lysophospholipid <- metabolite %>% select(all_of(Lysophospholipid_features))
Lysophospholipid <- Lysophospholipid %>% mutate(Lysophospholipid = rowSums(across(everything())))
Lysophospholipid$id <- rownames(Lysophospholipid)


Lysophospholipid_aMC <- Lysophospholipid %>% inner_join(select(activeMC_baseline, id, disease)) %>% mutate(mc_all_label = "MC")
Lysophospholipid_HC <- Lysophospholipid %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control") %>% mutate(disease = "Control without diarrhea")
Lysophospholipid_CD <- Lysophospholipid %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea") %>% mutate(disease = "Chronic diarrhea")
Lysophospholipid <- rbind(Lysophospholipid_aMC,Lysophospholipid_HC,Lysophospholipid_CD)
Lysophospholipid <- Lysophospholipid %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 

Lysophospholipid_filled <- Lysophospholipid %>% filter(!is.na(bmi))
Lysophospholipid_unfilled <- Lysophospholipid %>% filter(is.na(bmi))
Lysophospholipid_unfilled$ID <- substr(Lysophospholipid_unfilled$id, 1, 5)
Lysophospholipid_unfilled <- Lysophospholipid_unfilled %>% select(-age, -sex, -bmi)
Lysophospholipid_unfilled <- Lysophospholipid_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
Lysophospholipid <- rbind(Lysophospholipid_filled, Lysophospholipid_unfilled)

Lysophospholipid_log <- Lysophospholipid
Lysophospholipid_log$Lysophospholipid <- log(Lysophospholipid_log$Lysophospholipid)

## Box plot
Lysophospholipid_log$mc_all_label <- factor(Lysophospholipid_log$mc_all_label, 
                                            levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                            labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_Lysophospholipid <- ggplot(Lysophospholipid_log, aes(x = mc_all_label, y = Lysophospholipid)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Lysophospholipids)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysophospholipid

# ggsave(filename = "bx_cross_Lysophospholipid.pdf",
#        plot = bx_cross_Lysophospholipid, units = "in", width=10, height=4, dpi = 1200)


## Individual metabolite
Lysophospholipid_individual <- Lysophospholipid %>% select(4,16:21) 
Lysophospholipid_individual <- Lysophospholipid_individual %>%
  mutate(`1-stearoyl-GPC (18:0)` = log(`33961`)) %>%
  select(-`33961`)

## Box plot
Lysophospholipid_individual$mc_all_label <- factor(Lysophospholipid_individual$mc_all_label, 
                                            levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                            labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_Lysophospholipid_individual <- ggplot(Lysophospholipid_individual, aes(x = mc_all_label, y = `1-stearoyl-GPC (18:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(1-stearoyl-GPC (18:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysophospholipid_individual

# ggsave(filename = "bx_cross_Lysophospholipid_individual.pdf",
#        plot = bx_cross_Lysophospholipid_individual, units = "in", width=10, height=4, dpi = 1200)



Lysophospholipid_log$mc_all_label <- factor(Lysophospholipid_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_Lysophospholipid_log <- lm(Lysophospholipid ~ mc_all_label + age + sex + bmi, data = Lysophospholipid_log)
summary(lm_Lysophospholipid_log)

Lysophospholipid_cross <- Lysophospholipid

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/Lysophospholipid_halla.xlsx"
# write.xlsx(Lysophospholipid_cross, file_path, rowNames = TRUE)

Lysophospholipid_individual$mc_all_label <- factor(Lysophospholipid_individual$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_Lysophospholipid_individual_log <- lm(`1-stearoyl-GPC (18:0)` ~ mc_all_label + age + sex + bmi, data = Lysophospholipid_individual)
summary(lm_Lysophospholipid_individual_log)




### LC and CC subgroup 
Lysophospholipid_log$disease <- factor(Lysophospholipid_log$disease, 
                                       levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
Lysophospholipid_log_1 <- Lysophospholipid_log %>% filter(!is.na(disease))
Lysophospholipid_log_1$disease <- factor(Lysophospholipid_log_1$disease, 
                                         levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_Lysophospholipid_sub <- ggplot(Lysophospholipid_log_1, aes(x = disease, y = Lysophospholipid)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Lysophospholipids)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysophospholipid_sub

# ggsave(filename = "bx_cross_Lysophospholipid_sub.pdf",
#        plot = bx_cross_Lysophospholipid_sub, units = "in", width=10, height=5, dpi = 1200)


Lysophospholipid_log_1$disease <- factor(Lysophospholipid_log_1$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
Lysophospholipid_log_1$sex <- factor(Lysophospholipid_log_1$sex)
lm_Lysophospholipid_log_1 <- lm(Lysophospholipid ~ disease + age + sex + bmi, data = Lysophospholipid_log_1)
summary(lm_Lysophospholipid_log_1)

Lysophospholipid_log_1$disease <- factor(Lysophospholipid_log_1$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
Lysophospholipid_log_1$sex <- factor(Lysophospholipid_log_1$sex)
lm_Lysophospholipid_log_1 <- lm(Lysophospholipid ~ disease + age + sex + bmi, data = Lysophospholipid_log_1)
summary(lm_Lysophospholipid_log_1)


## Individual metabolite
Lysophospholipid_log_2 <- Lysophospholipid_log_1 %>% select(4,16:21) 
Lysophospholipid_log_2 <- Lysophospholipid_log_2 %>%
  mutate(`1-stearoyl-GPC (18:0)` = log(`33961`)) %>%
  select(-`33961`)

Lysophospholipid_log_2$disease <- factor(Lysophospholipid_log_2$disease, 
                                       levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_Lysophospholipid_individual_sub <- ggplot(Lysophospholipid_log_2, aes(x = disease, y = `1-stearoyl-GPC (18:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(1-stearoyl-GPC (18:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysophospholipid_individual_sub

# ggsave(filename = "bx_cross_Lysophospholipid_individual_sub.pdf",
#        plot = bx_cross_Lysophospholipid_individual_sub, units = "in", width=10, height=5, dpi = 1200)

Lysophospholipid_log_2$disease <- factor(Lysophospholipid_log_2$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
Lysophospholipid_log_2$sex <- factor(Lysophospholipid_log_2$sex)
lm_Lysophospholipid_log_2 <- lm(`1-stearoyl-GPC (18:0)` ~ disease + age + sex + bmi, data = Lysophospholipid_log_2)
summary(lm_Lysophospholipid_log_2)

Lysophospholipid_log_2$disease <- factor(Lysophospholipid_log_2$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
Lysophospholipid_log_2$sex <- factor(Lysophospholipid_log_2$sex)
lm_Lysophospholipid_log_2 <- lm(`1-stearoyl-GPC (18:0)` ~ disease + age + sex + bmi, data = Lysophospholipid_log_2)
summary(lm_Lysophospholipid_log_2)



####--------------- Lysophospholipids analysis: longitudinal  ----------------------#####
Lysophospholipid <- subset(altered_114, grepl("Lysophospholipid", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
Lysophospholipid_features <- Lysophospholipid$feature
Lysophospholipid <- metabolite %>% select(all_of(Lysophospholipid_features))
Lysophospholipid <- Lysophospholipid %>% mutate(Lysophospholipid = rowSums(across(everything())))
Lysophospholipid$id <- rownames(Lysophospholipid)

Lysophospholipid_long <- Lysophospholipid %>% inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 


Lysophospholipid_long_log <- Lysophospholipid_long
Lysophospholipid_long_log$Lysophospholipid <- log(Lysophospholipid_long_log$Lysophospholipid)

## Box plot
Lysophospholipid_long_log$symptoms <- factor(Lysophospholipid_long_log$symptoms, 
                                             levels = c("active_diarrhea", "no_diarrhea"),
                                             labels = c("active", "remission"))

bx_Lysophospholipid_long <- ggplot(Lysophospholipid_long_log, aes(x = symptoms, y = Lysophospholipid)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(Lysophospholipids)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysophospholipid_long

# ggsave(filename = "bx_Lysophospholipid_long.pdf",
#        plot = bx_Lysophospholipid_long, units = "in", width=10, height=4, dpi = 1200)

Lysophospholipid_long_log$symptoms <- factor(Lysophospholipid_long_log$symptoms, levels = c("remission", "active"))
lm_Lysophospholipid_long_log <- lm(Lysophospholipid ~ symptoms + age + sex + bmi, data = Lysophospholipid_long_log)
summary(lm_Lysophospholipid_long_log)


## Individual metabolite
Lysophospholipid_individual_long <- Lysophospholipid_long %>% select(4,16:20) 
Lysophospholipid_individual_long <- Lysophospholipid_individual_long %>%
  mutate(`1-stearoyl-GPC (18:0)` = log(`33961`)) %>%
  select(-`33961`)

## Box plot
Lysophospholipid_individual_long$symptoms <- factor(Lysophospholipid_individual_long$symptoms, 
                                             levels = c("active_diarrhea", "no_diarrhea"),
                                             labels = c("active", "remission"))

bx_Lysophospholipid_individual_long <- ggplot(Lysophospholipid_individual_long, aes(x = symptoms, y = `1-stearoyl-GPC (18:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(1-stearoyl-GPC (18:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysophospholipid_individual_long

# ggsave(filename = "bx_Lysophospholipid_individual_long.pdf",
#        plot = bx_Lysophospholipid_individual_long, units = "in", width=10, height=4, dpi = 1200)

Lysophospholipid_individual_long$symptoms <- factor(Lysophospholipid_individual_long$symptoms, levels = c("remission", "active"))
lm_Lysophospholipid_individual_long_log <- lm(`1-stearoyl-GPC (18:0)` ~ symptoms + age + sex + bmi, data = Lysophospholipid_individual_long)
summary(lm_Lysophospholipid_individual_long_log)




### LC and CC subgroup 
Lysophospholipid_long_log_1 <- Lysophospholipid_long_log %>% left_join(select(Lysophospholipid_log, id, disease))
Lysophospholipid_long_log_1$ID <- substr(Lysophospholipid_long_log_1$id, 1, 5)
Lysophospholipid_long_log_1$disease <- as.character(Lysophospholipid_long_log_1$disease)

Lysophospholipid_long_log_1 <- Lysophospholipid_long_log_1 %>%
  group_by(ID) %>%
  mutate(disease = ifelse(is.na(disease), first(na.omit(disease)), disease)) %>%
  ungroup()



Lysophospholipid_long_log_1$symptoms <- factor(Lysophospholipid_long_log$symptoms,
                                        levels = c("active", "remission"))
Lysophospholipid_long_log_1_LC <- Lysophospholipid_long_log_1 %>% filter(disease == "LC")

bx_Lysophospholipid_long_LC <- ggplot(Lysophospholipid_long_log_1_LC, aes(x = symptoms, y = Lysophospholipid)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(Lysophospholipids)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysophospholipid_long_LC
# ggsave(filename = "bx_Lysophospholipid_long_LC.pdf",
#        plot = bx_Lysophospholipid_long_LC, units = "in", width=10, height=4, dpi = 1200)

Lysophospholipid_long_log_1_LC$symptoms <- factor(Lysophospholipid_long_log_1_LC$symptoms, levels = c("remission", "active"))
lm_Lysophospholipid_long_log_1_LC <- lm(Lysophospholipid ~ symptoms + age + sex + bmi, data = Lysophospholipid_long_log_1_LC)
summary(lm_Lysophospholipid_long_log_1_LC)



Lysophospholipid_long_log_1_CC <- Lysophospholipid_long_log_1 %>% filter(disease == "CC")

bx_Lysophospholipid_long_CC <- ggplot(Lysophospholipid_long_log_1_CC, aes(x = symptoms, y = Lysophospholipid)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(Lysophospholipids)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysophospholipid_long_CC
# ggsave(filename = "bx_Lysophospholipid_long_CC.pdf",
#        plot = bx_Lysophospholipid_long_CC, units = "in", width=10, height=4, dpi = 1200)

Lysophospholipid_long_log_1_CC$symptoms <- factor(Lysophospholipid_long_log_1_CC$symptoms, levels = c("remission", "active"))
lm_Lysophospholipid_long_log_1_CC <- lm(Lysophospholipid ~ symptoms + age + bmi, data = Lysophospholipid_long_log_1_CC)
summary(lm_Lysophospholipid_long_log_1_CC)


## Individual metabolite
Lysophospholipid_long_log_2 <- Lysophospholipid_long_log_1 %>% select(4, 16:22)
Lysophospholipid_long_log_2 <- Lysophospholipid_long_log_2 %>%
  mutate(`1-stearoyl-GPC (18:0)` = log(`33961`)) %>%
  select(-`33961`)

Lysophospholipid_long_log_2_LC <- Lysophospholipid_long_log_2 %>% filter(disease == "LC")

bx_Lysophospholipid_individual_long_LC <- ggplot(Lysophospholipid_long_log_2_LC, aes(x = symptoms, y = `1-stearoyl-GPC (18:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(1-stearoyl-GPC (18:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysophospholipid_individual_long_LC
# ggsave(filename = "bx_Lysophospholipid_individual_long_LC.pdf",
#        plot = bx_Lysophospholipid_individual_long_LC, units = "in", width=10, height=4, dpi = 1200)

Lysophospholipid_long_log_2_LC$symptoms <- factor(Lysophospholipid_long_log_2_LC$symptoms, levels = c("remission", "active"))
lm_Lysophospholipid_long_log_2_LC <- lm(`1-stearoyl-GPC (18:0)` ~ symptoms + age + sex + bmi, data = Lysophospholipid_long_log_2_LC)
summary(lm_Lysophospholipid_long_log_2_LC)



Lysophospholipid_long_log_2_CC <- Lysophospholipid_long_log_2 %>% filter(disease == "CC")

bx_Lysophospholipid_individual_long_CC <- ggplot(Lysophospholipid_long_log_2_CC, aes(x = symptoms, y = `1-stearoyl-GPC (18:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(1-stearoyl-GPC (18:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysophospholipid_individual_long_CC
# ggsave(filename = "bx_Lysophospholipid_individual_long_CC.pdf",
#        plot = bx_Lysophospholipid_individual_long_CC, units = "in", width=10, height=4, dpi = 1200)

Lysophospholipid_long_log_2_CC$symptoms <- factor(Lysophospholipid_long_log_2_CC$symptoms, levels = c("remission", "active"))
lm_Lysophospholipid_long_log_2_CC <- lm(`1-stearoyl-GPC (18:0)` ~ symptoms + age + bmi, data = Lysophospholipid_long_log_2_CC)
summary(lm_Lysophospholipid_long_log_2_CC)





####--------------- Lysoplasmalogen analysis: cross-sectional  ----------------------#####
Lysoplasmalogen <- subset(altered_114, grepl("Lysoplasmalogen", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
Lysoplasmalogen_features <- Lysoplasmalogen$feature
Lysoplasmalogen <- metabolite %>% select(all_of(Lysoplasmalogen_features))
Lysoplasmalogen <- Lysoplasmalogen %>% mutate(Lysoplasmalogen = rowSums(across(everything())))
Lysoplasmalogen$id <- rownames(Lysoplasmalogen)

Lysoplasmalogen_aMC <- Lysoplasmalogen %>% inner_join(select(activeMC_baseline, id, disease)) %>% mutate(mc_all_label = "MC")
Lysoplasmalogen_HC <- Lysoplasmalogen %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control") %>% mutate(disease = "Control without diarrhea")
Lysoplasmalogen_CD <- Lysoplasmalogen %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea") %>% mutate(disease = "Chronic diarrhea")
Lysoplasmalogen <- rbind(Lysoplasmalogen_aMC,Lysoplasmalogen_HC,Lysoplasmalogen_CD)
Lysoplasmalogen <- Lysoplasmalogen %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 

Lysoplasmalogen_filled <- Lysoplasmalogen %>% filter(!is.na(bmi))
Lysoplasmalogen_unfilled <- Lysoplasmalogen %>% filter(is.na(bmi))
Lysoplasmalogen_unfilled$ID <- substr(Lysoplasmalogen_unfilled$id, 1, 5)
Lysoplasmalogen_unfilled <- Lysoplasmalogen_unfilled %>% select(-age, -sex, -bmi)
Lysoplasmalogen_unfilled <- Lysoplasmalogen_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
Lysoplasmalogen <- rbind(Lysoplasmalogen_filled, Lysoplasmalogen_unfilled)

Lysoplasmalogen_log <- Lysoplasmalogen
Lysoplasmalogen_log$Lysoplasmalogen <- log(Lysoplasmalogen_log$Lysoplasmalogen)

## Box plot
Lysoplasmalogen_log$mc_all_label <- factor(Lysoplasmalogen_log$mc_all_label, 
                                           levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                           labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_Lysoplasmalogen <- ggplot(Lysoplasmalogen_log, aes(x = mc_all_label, y = Lysoplasmalogen)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Lysoplasmalogens)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysoplasmalogen

# ggsave(filename = "bx_cross_Lysoplasmalogen.pdf",
#        plot = bx_cross_Lysoplasmalogen, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_log$mc_all_label <- factor(Lysoplasmalogen_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_Lysoplasmalogen_log <- lm(Lysoplasmalogen ~ mc_all_label + age + sex + bmi, data = Lysoplasmalogen_log)
summary(lm_Lysoplasmalogen_log)

Lysoplasmalogen_cross <- Lysoplasmalogen

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/Lysoplasmalogen_halla.xlsx"
# write.xlsx(Lysoplasmalogen_cross, file_path, rowNames = TRUE)


## Individual metabolite
Lysoplasmalogen_individual <- Lysoplasmalogen %>% select(2,10:15) 
Lysoplasmalogen_individual <- Lysoplasmalogen_individual %>%
  mutate(`1-stearyl-GPE (O-18:0)*` = log(`61501`)) %>%
  select(-`61501`)

## Box plot
Lysoplasmalogen_individual$mc_all_label <- factor(Lysoplasmalogen_individual$mc_all_label, 
                                           levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                           labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_Lysoplasmalogen_individual <- ggplot(Lysoplasmalogen_individual, aes(x = mc_all_label, y = `1-stearyl-GPE (O-18:0)*`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(1-stearyl-GPE (O-18:0)*)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysoplasmalogen_individual

# ggsave(filename = "bx_cross_Lysoplasmalogen_individual.pdf",
#        plot = bx_cross_Lysoplasmalogen_individual, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_individual$mc_all_label <- factor(Lysoplasmalogen_individual$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_Lysoplasmalogen_individual <- lm(`1-stearyl-GPE (O-18:0)*` ~ mc_all_label + age + sex + bmi, data = Lysoplasmalogen_individual)
summary(lm_Lysoplasmalogen_individual)




### LC and CC subgroup 
Lysoplasmalogen_log$disease <- factor(Lysoplasmalogen_log$disease, 
                                       levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
Lysoplasmalogen_log_1 <- Lysoplasmalogen_log %>% filter(!is.na(disease))

bx_cross_Lysoplasmalogen_sub <- ggplot(Lysoplasmalogen_log_1, aes(x = disease, y = Lysoplasmalogen)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(Lysoplasmalogens)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysoplasmalogen_sub

# ggsave(filename = "bx_cross_Lysoplasmalogen_sub.pdf",
#        plot = bx_cross_Lysoplasmalogen_sub, units = "in", width=10, height=5, dpi = 1200)


Lysoplasmalogen_log_1$disease <- factor(Lysoplasmalogen_log_1$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
Lysoplasmalogen_log_1$sex <- factor(Lysoplasmalogen_log_1$sex)
lm_Lysoplasmalogen_log_1 <- lm(Lysoplasmalogen ~ disease + age + sex + bmi, data = Lysoplasmalogen_log_1)
summary(lm_Lysoplasmalogen_log_1)

Lysoplasmalogen_log_1$disease <- factor(Lysoplasmalogen_log_1$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
Lysoplasmalogen_log_1$sex <- factor(Lysoplasmalogen_log_1$sex)
lm_Lysoplasmalogen_log_1 <- lm(Lysoplasmalogen ~ disease + age + sex + bmi, data = Lysoplasmalogen_log_1)
summary(lm_Lysoplasmalogen_log_1)


## Individual metabolite
Lysoplasmalogen_log_2 <- Lysoplasmalogen_log_1 %>% select(2, 10:15)
Lysoplasmalogen_log_2 <- Lysoplasmalogen_log_2 %>%
  mutate(`1-stearyl-GPE (O-18:0)*` = log(`61501`)) %>%
  select(-`61501`)

bx_cross_Lysoplasmalogen_individual_sub <- ggplot(Lysoplasmalogen_log_2, aes(x = disease, y = `1-stearyl-GPE (O-18:0)*`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = disease), width = 0.2, size = 2) +  
  geom_point(aes(color = disease), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(1-stearyl-GPE (O-18:0)*)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_Lysoplasmalogen_individual_sub

# ggsave(filename = "bx_cross_Lysoplasmalogen_individual_sub.pdf",
#        plot = bx_cross_Lysoplasmalogen_individual_sub, units = "in", width=10, height=5, dpi = 1200)


Lysoplasmalogen_log_2$disease <- factor(Lysoplasmalogen_log_2$disease, levels = c("CC", "LC", "Chronic diarrhea", "Control without diarrhea"))
Lysoplasmalogen_log_2$sex <- factor(Lysoplasmalogen_log_2$sex)
lm_Lysoplasmalogen_log_2 <- lm(`1-stearyl-GPE (O-18:0)*` ~ disease + age + sex + bmi, data = Lysoplasmalogen_log_2)
summary(lm_Lysoplasmalogen_log_2)

Lysoplasmalogen_log_2$disease <- factor(Lysoplasmalogen_log_2$disease, levels = c("LC", "CC", "Chronic diarrhea", "Control without diarrhea"))
Lysoplasmalogen_log_2$sex <- factor(Lysoplasmalogen_log_2$sex)
lm_Lysoplasmalogen_log_2 <- lm(`1-stearyl-GPE (O-18:0)*` ~ disease + age + sex + bmi, data = Lysoplasmalogen_log_2)
summary(lm_Lysoplasmalogen_log_2)




####--------------- Lysoplasmalogen analysis: longitudinal  ----------------------#####
Lysoplasmalogen <- subset(altered_114, grepl("Lysoplasmalogen", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
Lysoplasmalogen_features <- Lysoplasmalogen$feature
Lysoplasmalogen <- metabolite %>% select(all_of(Lysoplasmalogen_features))
Lysoplasmalogen <- Lysoplasmalogen %>% mutate(Lysoplasmalogen = rowSums(across(everything())))
Lysoplasmalogen$id <- rownames(Lysoplasmalogen)

Lysoplasmalogen_long <- Lysoplasmalogen %>% inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 


Lysoplasmalogen_long_log <- Lysoplasmalogen_long
Lysoplasmalogen_long_log$Lysoplasmalogen <- log(Lysoplasmalogen_long_log$Lysoplasmalogen)

## Box plot
Lysoplasmalogen_long_log$symptoms <- factor(Lysoplasmalogen_long_log$symptoms, 
                                            levels = c("active_diarrhea", "no_diarrhea"),
                                            labels = c("active", "remission"))

bx_Lysoplasmalogen_long <- ggplot(Lysoplasmalogen_long_log, aes(x = symptoms, y = Lysoplasmalogen)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(Lysoplasmalogens)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysoplasmalogen_long

# ggsave(filename = "bx_Lysoplasmalogen_long.pdf",
#        plot = bx_Lysoplasmalogen_long, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_long_log$symptoms <- factor(Lysoplasmalogen_long_log$symptoms, levels = c("remission", "active"))
lm_Lysoplasmalogen_long_log <- lm(Lysoplasmalogen ~ symptoms + age + sex + bmi, data = Lysoplasmalogen_long_log)
summary(lm_Lysoplasmalogen_long_log)


## Individual metabolite
Lysoplasmalogen_individual_long <- Lysoplasmalogen_long %>% select(2,10:14) 
Lysoplasmalogen_individual_long <- Lysoplasmalogen_individual_long %>%
  mutate(`1-stearyl-GPE (O-18:0)*` = log(`61501`)) %>%
  select(-`61501`)

## Box plot
Lysoplasmalogen_individual_long$symptoms <- factor(Lysoplasmalogen_individual_long$symptoms, 
                                            levels = c("active_diarrhea", "no_diarrhea"),
                                            labels = c("active", "remission"))

bx_Lysoplasmalogen_individual_long <- ggplot(Lysoplasmalogen_individual_long, aes(x = symptoms, y = `1-stearyl-GPE (O-18:0)*`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(1-stearyl-GPE (O-18:0)*)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysoplasmalogen_individual_long

# ggsave(filename = "bx_Lysoplasmalogen_individual_long.pdf",
#        plot = bx_Lysoplasmalogen_individual_long, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_individual_long$symptoms <- factor(Lysoplasmalogen_individual_long$symptoms, levels = c("remission", "active"))
lm_Lysoplasmalogen_individual_long_log <- lm(`1-stearyl-GPE (O-18:0)*` ~ symptoms + age + sex + bmi, data = Lysoplasmalogen_individual_long)
summary(lm_Lysoplasmalogen_individual_long_log)




### LC and CC subgroup 
Lysoplasmalogen_long_log_1 <- Lysoplasmalogen_long_log %>% left_join(select(Lysoplasmalogen_log, id, disease))
Lysoplasmalogen_long_log_1$ID <- substr(Lysoplasmalogen_long_log_1$id, 1, 5)
Lysoplasmalogen_long_log_1$disease <- as.character(Lysoplasmalogen_long_log_1$disease)

Lysoplasmalogen_long_log_1 <- Lysoplasmalogen_long_log_1 %>%
  group_by(ID) %>%
  mutate(disease = ifelse(is.na(disease), first(na.omit(disease)), disease)) %>%
  ungroup()



Lysoplasmalogen_long_log_1$symptoms <- factor(Lysoplasmalogen_long_log$symptoms,
                                               levels = c("active", "remission"))
Lysoplasmalogen_long_log_1_LC <- Lysoplasmalogen_long_log_1 %>% filter(disease == "LC")

bx_Lysoplasmalogen_long_LC <- ggplot(Lysoplasmalogen_long_log_1_LC, aes(x = symptoms, y = Lysoplasmalogen)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(Lysoplasmalogens)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysoplasmalogen_long_LC
# ggsave(filename = "bx_Lysoplasmalogen_long_LC.pdf",
#        plot = bx_Lysoplasmalogen_long_LC, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_long_log_1_LC$symptoms <- factor(Lysoplasmalogen_long_log_1_LC$symptoms, levels = c("remission", "active"))
lm_Lysoplasmalogen_long_log_1_LC <- lm(Lysoplasmalogen ~ symptoms + age + sex + bmi, data = Lysoplasmalogen_long_log_1_LC)
summary(lm_Lysoplasmalogen_long_log_1_LC)



Lysoplasmalogen_long_log_1_CC <- Lysoplasmalogen_long_log_1 %>% filter(disease == "CC")

bx_Lysoplasmalogen_long_CC <- ggplot(Lysoplasmalogen_long_log_1_CC, aes(x = symptoms, y = Lysoplasmalogen)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(Lysoplasmalogens)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysoplasmalogen_long_CC
# ggsave(filename = "bx_Lysoplasmalogen_long_CC.pdf",
#        plot = bx_Lysoplasmalogen_long_CC, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_long_log_1_CC$symptoms <- factor(Lysoplasmalogen_long_log_1_CC$symptoms, levels = c("remission", "active"))
lm_Lysoplasmalogen_long_log_1_CC <- lm(Lysoplasmalogen ~ symptoms + age + bmi, data = Lysoplasmalogen_long_log_1_CC)
summary(lm_Lysoplasmalogen_long_log_1_CC)



## Individual metabolite
Lysoplasmalogen_long_log_2 <- Lysoplasmalogen_long_log_1 %>% select(2,10:16) 
Lysoplasmalogen_long_log_2 <- Lysoplasmalogen_long_log_2 %>%
  mutate(`1-stearyl-GPE (O-18:0)*` = log(`61501`)) %>%
  select(-`61501`)


Lysoplasmalogen_long_log_2$symptoms <- factor(Lysoplasmalogen_long_log$symptoms,
                                              levels = c("active", "remission"))
Lysoplasmalogen_long_log_2_LC <- Lysoplasmalogen_long_log_2 %>% filter(disease == "LC")

bx_Lysoplasmalogen_individual_long_LC <- ggplot(Lysoplasmalogen_long_log_2_LC, aes(x = symptoms, y = `1-stearyl-GPE (O-18:0)*`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "LC activity", y = "Log(1-stearyl-GPE (O-18:0)*)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysoplasmalogen_individual_long_LC
# ggsave(filename = "bx_Lysoplasmalogen_individual_long_LC.pdf",
#        plot = bx_Lysoplasmalogen_individual_long_LC, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_long_log_2_LC$symptoms <- factor(Lysoplasmalogen_long_log_2_LC$symptoms, levels = c("remission", "active"))
lm_Lysoplasmalogen_long_log_2_LC <- lm(`1-stearyl-GPE (O-18:0)*` ~ symptoms + age + sex + bmi, data = Lysoplasmalogen_long_log_2_LC)
summary(lm_Lysoplasmalogen_long_log_2_LC)



Lysoplasmalogen_long_log_2_CC <- Lysoplasmalogen_long_log_2 %>% filter(disease == "CC")

bx_Lysoplasmalogen_individual_long_CC <- ggplot(Lysoplasmalogen_long_log_2_CC, aes(x = symptoms, y = `1-stearyl-GPE (O-18:0)*`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "CC activity", y = "Log(1-stearyl-GPE (O-18:0)*)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_Lysoplasmalogen_individual_long_CC
# ggsave(filename = "bx_Lysoplasmalogen_individual_long_CC.pdf",
#        plot = bx_Lysoplasmalogen_individual_long_CC, units = "in", width=10, height=4, dpi = 1200)

Lysoplasmalogen_long_log_2_CC$symptoms <- factor(Lysoplasmalogen_long_log_2_CC$symptoms, levels = c("remission", "active"))
lm_Lysoplasmalogen_long_log_2_CC <- lm(`1-stearyl-GPE (O-18:0)*` ~ symptoms + age + bmi, data = Lysoplasmalogen_long_log_2_CC)
summary(lm_Lysoplasmalogen_long_log_2_CC)




####--------------- Scatter plot: microbe and metabolite read in ----------------------#####
ceramides_mbx <- ceramides_cross %>% 
  select("57432", "57434", "65035", "57443", "44877", "id", "mc_all_label") %>%
  rename("ceramide (d18:1/14:0, d16:1/16:0)*" = "57432",
         "ceramide (d18:1/17:0, d17:1/18:0)*" = "57434",
         "ceramide (d18:2/16:0, d18:1/16:1, d16:1/18:1)*" = "65035",
         "ceramide (d18:2/24:1, d18:1/24:2)*" = "57443",
         "N-palmitoyl-sphingosine (d18:1/16:0)" = "44877")%>%
  mutate(id = case_when(
    id == "MC058" ~ "MC058A",
    id == "MC062" ~ "MC062A",
    id == "MC063" ~ "MC063A",
    id == "MC075" ~ "MC075A",
    id == "MC139A" ~ "MC139",
    id == "MC167A" ~ "MC167",
    id == "MC244A" ~ "MC244",
    id == "MC408A" ~ "MC408B",
    TRUE ~ id  # Keep the original value if no match
  ))




dihydroceramides_mbx <- metabolite %>% 
  select("52604")%>%
  rename("N-palmitoyl-sphinganine (d18:0/16:0)" = "52604") %>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(ceramides_cross, id, mc_all_label)) %>%
  mutate(id = case_when(
    id == "MC058" ~ "MC058A",
    id == "MC062" ~ "MC062A",
    id == "MC063" ~ "MC063A",
    id == "MC075" ~ "MC075A",
    id == "MC139A" ~ "MC139",
    id == "MC167A" ~ "MC167",
    id == "MC244A" ~ "MC244",
    id == "MC408A" ~ "MC408B",
    TRUE ~ id  # Keep the original value if no match
  ))




dihydrosphingomyelins_mbx <- metabolite %>% 
  select("52434")%>%
  rename("palmitoyl dihydrosphingomyelin (d18:0/16:0)*" = "52434") %>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(ceramides_cross, id, mc_all_label)) %>%
  mutate(id = case_when(
    id == "MC058" ~ "MC058A",
    id == "MC062" ~ "MC062A",
    id == "MC063" ~ "MC063A",
    id == "MC075" ~ "MC075A",
    id == "MC139A" ~ "MC139",
    id == "MC167A" ~ "MC167",
    id == "MC244A" ~ "MC244",
    id == "MC408A" ~ "MC408B",
    TRUE ~ id  # Keep the original value if no match
  ))



glucosylceramides_mbx <- metabolite %>% 
  select("53013")%>%
  rename("glycosyl-N-palmitoyl-sphingosine (d18:1/16:0)" = "53013") %>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(ceramides_cross, id, mc_all_label)) %>%
  mutate(id = case_when(
    id == "MC058" ~ "MC058A",
    id == "MC062" ~ "MC062A",
    id == "MC063" ~ "MC063A",
    id == "MC075" ~ "MC075A",
    id == "MC139A" ~ "MC139",
    id == "MC167A" ~ "MC167",
    id == "MC244A" ~ "MC244",
    id == "MC408A" ~ "MC408B",
    TRUE ~ id  # Keep the original value if no match
  )) 



lactosylceramides_mbx <- lactosylceramides_cross %>% 
  select("53010", "57370", "57422", "id", "mc_all_label") %>%
  rename("lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)" = "53010",
         "lactosyl-N-nervonoyl-sphingosine (d18:1/24:1)*" = "57370",
         "lactosyl-N-behenoyl-sphingosine (d18:1/22:0)*" = "57422")%>%
  mutate(id = case_when(
    id == "MC058" ~ "MC058A",
    id == "MC062" ~ "MC062A",
    id == "MC063" ~ "MC063A",
    id == "MC075" ~ "MC075A",
    id == "MC139A" ~ "MC139",
    id == "MC167A" ~ "MC167",
    id == "MC244A" ~ "MC244",
    id == "MC408A" ~ "MC408B",
    TRUE ~ id  # Keep the original value if no match
  )) 



sphingomyelins_mbx <- metabolite %>% 
  select("37506")%>%
  rename("palmitoyl sphingomyelin (d18:1/16:0)" = "37506") %>% 
  mutate(id = rownames(metabolite)) %>% 
  inner_join(select(ceramides_cross, id, mc_all_label)) %>%
  mutate(id = case_when(
    id == "MC058" ~ "MC058A",
    id == "MC062" ~ "MC062A",
    id == "MC063" ~ "MC063A",
    id == "MC075" ~ "MC075A",
    id == "MC139A" ~ "MC139",
    id == "MC167A" ~ "MC167",
    id == "MC244A" ~ "MC244",
    id == "MC408A" ~ "MC408B",
    TRUE ~ id  # Keep the original value if no match
  )) 



veillonella_parvula_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/veillonella_parvula_mgx.xlsx") %>%
  rename("id" = "...1")

veillonella_rogosae_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/veillonella_rogosae_mgx.xlsx") %>%
  rename("id" = "...1")

intestinibacter_bartlettii_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/intestinibacter_bartlettii_mgx.xlsx") %>%
  rename("id" = "...1")

methylobacterium_SGB15164_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/methylobacterium_SGB15164_mgx.xlsx") %>%
  rename("id" = "...1")

collinsella_SGB4121_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/collinsella_SGB4121_mgx.xlsx") %>%
  rename("id" = "...1")

mediterraneibacter_butyricigenes_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/mediterraneibacter_butyricigenes_mgx.xlsx") %>%
  rename("id" = "...1")

haemophilus_parainfluenzae_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/haemophilus_parainfluenzae_mgx.xlsx") %>%
  rename("id" = "...1")







####--------------- Scatter plot: microbe and metabolite analysis ----------------------#####
## veillonella_parvula & ceramides scatter plot
veillonella_parvula_ceramide <- veillonella_parvula_mgx %>% left_join(ceramides_mbx)

veillonella_parvula_ceramide$mc_all_label <- factor(veillonella_parvula_ceramide$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
                                                    labels = c("MC", "Control without diarrhea", "Chronic diarrhea"))

veillonella_parvula_ceramide_1 <- ggplot(veillonella_parvula_ceramide, aes(x = `ceramide (d18:2/16:0, d18:1/16:1, d16:1/18:1)*`, y = Veillonella_parvula, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of ceramide (d18:2/16:0, d18:1/16:1, d16:1/18:1)*", y = "Relative abundance of Veillonella parvula", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  theme_minimal()+
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 2)) + 
  annotate("text", x = 10, y = 1.5, label = "Spearman rho = 0.177, FDR = 0.010", size = 4, color = "black")
veillonella_parvula_ceramide_1
# ggsave(filename = "veillonella_parvula_ceramide_1.pdf",
#        plot = veillonella_parvula_ceramide_1, units = "in", width=7, height=5, dpi = 1200)



## veillonella_parvula & dihydroceramides scatter plot
veillonella_parvula_dihydroceramides <- veillonella_parvula_mgx %>% left_join(dihydroceramides_mbx)

veillonella_parvula_dihydroceramides$mc_all_label <- factor(veillonella_parvula_dihydroceramides$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
                                                            labels = c("MC", "HC", "CD"))

veillonella_parvula_dihydroceramides_1 <- ggplot(veillonella_parvula_dihydroceramides, aes(x = `N-palmitoyl-sphinganine (d18:0/16:0)`, y = Veillonella_parvula, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of N-palmitoyl-sphinganine (d18:0/16:0)", y = "Relative abundance of Veillonella parvula", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  theme_minimal()+
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 2)) + 
  annotate("text", x = 10, y = 1.5, label = "Spearman rho = 0.122, FDR = 0.072", size = 4, color = "black")
veillonella_parvula_dihydroceramides_1
# ggsave(filename = "veillonella_parvula_dihydroceramides_1.pdf",
#        plot = veillonella_parvula_dihydroceramides_1, units = "in", width=7, height=5, dpi = 1200)



## intestinibacter_bartlettii & dihydrosphingomyelins scatter plot
intestinibacter_bartlettii_dihydrosphingomyelins <- intestinibacter_bartlettii_mgx %>% left_join(dihydrosphingomyelins_mbx)

intestinibacter_bartlettii_dihydrosphingomyelins$mc_all_label <- factor(intestinibacter_bartlettii_dihydrosphingomyelins$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
                                                                        labels = c("MC", "HC", "CD"))

intestinibacter_bartlettii_dihydrosphingomyelins_1 <- ggplot(intestinibacter_bartlettii_dihydrosphingomyelins, aes(x = `palmitoyl dihydrosphingomyelin (d18:0/16:0)*`, y = Intestinibacter_bartlettii, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of palmitoyl dihydrosphingomyelin (d18:0/16:0)*", y = "Relative abundance of Intestinibacter bartlettii", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.075)) + 
  annotate("text", x = 10, y = 0.065, label = "Spearman rho = 0.117, FDR = 0.087", size = 4, color = "black")
# ggsave(filename = "intestinibacter_bartlettii_dihydrosphingomyelins_1.png",
#        plot = intestinibacter_bartlettii_dihydrosphingomyelins_1, units = "in", width=7, height=5, dpi = 300)



## Methylobacterium_SGB15164 & glucosylceramides scatter plot
methylobacterium_SGB15164_glucosylceramides <- methylobacterium_SGB15164_mgx %>% left_join(glucosylceramides_mbx)

methylobacterium_SGB15164_glucosylceramides$mc_all_label <- factor(methylobacterium_SGB15164_glucosylceramides$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
                                                                   labels = c("MC", "Control without diarrhea", "Chronic diarrhea"))

methylobacterium_SGB15164_glucosylceramides_1 <- ggplot(methylobacterium_SGB15164_glucosylceramides, aes(x = `glycosyl-N-palmitoyl-sphingosine (d18:1/16:0)`, y = Methylobacterium_SGB15164, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of glycosyl-N-palmitoyl-sphingosine (d18:1/16:0)", y = "Relative abundance of Methylobacterium_SGB15164", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  theme_minimal()+
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.20)) + 
  annotate("text", x = 10, y = 0.18, label = "Spearman rho = -0.153, FDR = 0.033", size = 4, color = "black")
methylobacterium_SGB15164_glucosylceramides_1
# ggsave(filename = "methylobacterium_SGB15164_glucosylceramides_1.pdf",
#        plot = methylobacterium_SGB15164_glucosylceramides_1, units = "in", width=7, height=5, dpi = 1200)



## Mediterraneibacter_butyricigenes & lactosylceramides scatter plot
mediterraneibacter_butyricigenes_lactosylceramides <- mediterraneibacter_butyricigenes_mgx %>% left_join(lactosylceramides_mbx)

mediterraneibacter_butyricigenes_lactosylceramides$mc_all_label <- factor(mediterraneibacter_butyricigenes_lactosylceramides$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
                                                                          labels = c("MC", "HC", "CD"))

mediterraneibacter_butyricigenes_lactosylceramides_1 <- ggplot(mediterraneibacter_butyricigenes_lactosylceramides, aes(x = `lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)`, y = Mediterraneibacter_butyricigenes, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)", y = "Relative abundance of Mediterraneibacter butyricigenes", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  theme_minimal()+
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 0.04)) + 
  annotate("text", x = 20, y = 0.035, label = "Spearman rho = -0.269, FDR < 0.001", size = 4, color = "black")
# ggsave(filename = "mediterraneibacter_butyricigenes_lactosylceramides_1.png",
#        plot = mediterraneibacter_butyricigenes_lactosylceramides_1, units = "in", width=7, height=5, dpi = 300)


## Haemophilus_parainfluenzae & sphingomyelins scatter plot
haemophilus_parainfluenzae_sphingomyelins <- haemophilus_parainfluenzae_mgx %>% left_join(sphingomyelins_mbx)

haemophilus_parainfluenzae_sphingomyelins$mc_all_label <- factor(haemophilus_parainfluenzae_sphingomyelins$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
                                                                 labels = c("MC", "HC", "CD"))

haemophilus_parainfluenzae_sphingomyelins_1 <- ggplot(haemophilus_parainfluenzae_sphingomyelins, aes(x = `palmitoyl sphingomyelin (d18:1/16:0)`, y = Haemophilus_parainfluenzae, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of palmitoyl sphingomyelin (d18:1/16:0)", y = "Relative abundance of Haemophilus parainfluenzae", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  theme_minimal()+
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 0.3)) + 
  annotate("text", x = 40, y = 0.27, label = "Spearman rho = 0.252, FDR < 0.001", size = 4, color = "black")
# ggsave(filename = "haemophilus_parainfluenzae_sphingomyelins_1.png",
#        plot = haemophilus_parainfluenzae_sphingomyelins_1, units = "in", width=7, height=5, dpi = 300)



####--------------- Scatter plot: enzyme and metabolite read in ----------------------#####
enzyme_ceramide_degradation_mgx <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_enzyme/enzyme_ceramide_degradation_mgx.xlsx") %>%
  rename("id" = "...1")

enzyme_ceramides_mbx <- ceramides_cross %>% 
  select("ceramides", "id", "mc_all_label") %>%
  mutate(id = case_when(
    id == "MC058" ~ "MC058A",
    id == "MC062" ~ "MC062A",
    id == "MC063" ~ "MC063A",
    id == "MC075" ~ "MC075A",
    id == "MC139A" ~ "MC139",
    id == "MC167A" ~ "MC167",
    id == "MC244A" ~ "MC244",
    id == "MC408A" ~ "MC408B",
    TRUE ~ id  # Keep the original value if no match
  ))




####--------------- Scatter plot: enzyme and metabolite analysis ----------------------#####
## enzyme_ceramide_degradation & ceramides scatter plot
enzyme_ceramide_degradation_scatter <- enzyme_ceramide_degradation_mgx %>% left_join(enzyme_ceramides_mbx)

enzyme_ceramide_degradation_scatter$mc_all_label <- factor(enzyme_ceramide_degradation_scatter$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
                                                    labels = c("MC", "Healthy control", "Chronic diarrhea"))

enzyme_ceramide_degradation_scatter_1 <- ggplot(enzyme_ceramide_degradation_scatter, aes(x = ceramides, y = ceramide_degradation, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of ceramides", y = "Relative abundance of microbial enzyme: ceramide degradation by alpha-oxidation", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.00005))
# ggsave(filename = "veillonella_parvula_ceramide_1.png",
#        plot = veillonella_parvula_ceramide_1, units = "in", width=7, height=5, dpi = 600)



####--------------- SCFA analysis: cross-sectional  ----------------------#####
scfa <- metabolite %>% select(1,3:5, "33443", "40605")
scfa <- scfa %>% mutate(scfa = rowSums(across(5:6)))

scfa_aMC <- scfa %>% inner_join(select(activeMC_baseline, id, disease)) %>% mutate(mc_all_label = "MC")
scfa_HC <- scfa %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control") %>% mutate(disease = "Control without diarrhea")
scfa_CD <- scfa %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea") %>% mutate(disease = "Chronic diarrhea")
scfa <- rbind(scfa_aMC,scfa_HC,scfa_CD)
scfa <- scfa %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 

scfa_filled <- scfa %>% filter(!is.na(bmi))
scfa_unfilled <- scfa %>% filter(is.na(bmi))
scfa_unfilled$ID <- substr(scfa_unfilled$id, 1, 5)
scfa_unfilled <- scfa_unfilled %>% select(-age, -sex, -bmi)
scfa_unfilled <- scfa_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
scfa <- rbind(scfa_filled, scfa_unfilled)

scfa_log <- scfa
scfa_log$scfa <- log(scfa_log$scfa)

## Box plot
scfa_log$mc_all_label <- factor(scfa_log$mc_all_label, 
                                levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_scfa <- ggplot(scfa_log, aes(x = mc_all_label, y = scfa)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(scfa)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_scfa

# ggsave(filename = "bx_cross_scfa.pdf",
#        plot = bx_cross_scfa, units = "in", width=10, height=4, dpi = 1200)


## Individual metabolite: valerate (5:0)
scfa_individual <- scfa %>% select(5, 1, 8:11) 
scfa_individual <- scfa_individual %>%
  mutate(`valerate (5:0)` = log(`33443`)) %>%
  select(-`33443`)

## Box plot
scfa_individual$mc_all_label <- factor(scfa_individual$mc_all_label, 
                                       levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                       labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_scfa_individual <- ggplot(scfa_individual, aes(x = mc_all_label, y = `valerate (5:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(valerate (5:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_scfa_individual

# ggsave(filename = "bx_cross_scfa_individual.pdf",
#        plot = bx_cross_scfa_individual, units = "in", width=10, height=4, dpi = 1200)



scfa_log$mc_all_label <- factor(scfa_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_scfa_log <- lm(scfa ~ mc_all_label + age + sex + bmi, data = scfa_log)
summary(lm_scfa_log)

scfa_cross <- scfa

scfa_individual$mc_all_label <- factor(scfa_individual$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_scfa_individual_log <- lm(`valerate (5:0)` ~ mc_all_label + age + sex + bmi, data = scfa_individual)
summary(lm_scfa_individual_log)


## Individual metabolite: butyrate/isobutyrate (4:0)
scfa_individual <- scfa %>% select(6, 1, 8:11) 
scfa_individual <- scfa_individual %>%
  mutate(`butyrate/isobutyrate (4:0)` = log(`40605`)) %>%
  select(-`40605`)

## Box plot
scfa_individual$mc_all_label <- factor(scfa_individual$mc_all_label, 
                                       levels = c("MC", "Chronic diarrhea", "Healthy control"),
                                       labels = c("MC", "Chronic diarrhea", "Control without diarrhea"))

bx_cross_scfa_individual <- ggplot(scfa_individual, aes(x = mc_all_label, y = `butyrate/isobutyrate (4:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
  geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
  labs(x = "Disease type", y = "Log(butyrate/isobutyrate (4:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_cross_scfa_individual

# ggsave(filename = "bx_cross_scfa_individual.pdf",
#        plot = bx_cross_scfa_individual, units = "in", width=10, height=4, dpi = 1200)

scfa_individual$mc_all_label <- factor(scfa_individual$mc_all_label, levels = c("MC", "Chronic diarrhea", "Control without diarrhea"))
lm_scfa_individual_log <- lm(`butyrate/isobutyrate (4:0)` ~ mc_all_label + age + sex + bmi, data = scfa_individual)
summary(lm_scfa_individual_log)


####--------------- SCFA analysis: longitudinal  ----------------------#####
scfa <- metabolite %>% select(1, "33443", "40605")
scfa <- scfa %>% mutate(scfa = rowSums(across(2:3)))

scfa_long <- scfa %>% inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 


scfa_long_log <- scfa_long
scfa_long_log$scfa <- log(scfa_long_log$scfa)

## Box plot
scfa_long_log$symptoms <- factor(scfa_long_log$symptoms, 
                                             levels = c("active_diarrhea", "no_diarrhea"),
                                             labels = c("active", "remission"))

bx_scfa_long <- ggplot(scfa_long_log, aes(x = symptoms, y = scfa)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(scfa)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_scfa_long

# ggsave(filename = "bx_scfa_long.pdf",
#        plot = bx_scfa_long, units = "in", width=10, height=4, dpi = 1200)

scfa_long_log$symptoms <- factor(scfa_long_log$symptoms, levels = c("remission", "active"))
lm_scfa_long_log <- lm(scfa ~ symptoms + age + sex + bmi, data = scfa_long_log)
summary(lm_scfa_long_log)


## Individual metabolite: valerate (5:0)
scfa_individual_long <- scfa_long %>% select(1,2,5:8) 
scfa_individual_long <- scfa_individual_long %>%
  mutate(`valerate (5:0)` = log(`33443`)) %>%
  select(-`33443`)

## Box plot
scfa_individual_long$symptoms <- factor(scfa_individual_long$symptoms, 
                                        levels = c("active_diarrhea", "no_diarrhea"),
                                        labels = c("active", "remission"))

bx_scfa_individual_long <- ggplot(scfa_individual_long, aes(x = symptoms, y = `valerate (5:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(valerate (5:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_scfa_individual_long

# ggsave(filename = "bx_scfa_individual_long.pdf",
#        plot = bx_scfa_individual_long, units = "in", width=10, height=4, dpi = 1200)

scfa_individual_long$symptoms <- factor(scfa_individual_long$symptoms, levels = c("remission", "active"))
lm_scfa_individual_long_log <- lm(`valerate (5:0)` ~ symptoms + age + sex + bmi, data = scfa_individual_long)
summary(lm_scfa_individual_long_log)



## Individual metabolite: valerate (5:0)
scfa_individual_long <- scfa_long %>% select(1,3,5:8) 
scfa_individual_long <- scfa_individual_long %>%
  mutate(`butyrate/isobutyrate (4:0)` = log(`40605`)) %>%
  select(-`40605`)

## Box plot
scfa_individual_long$symptoms <- factor(scfa_individual_long$symptoms, 
                                        levels = c("active_diarrhea", "no_diarrhea"),
                                        labels = c("active", "remission"))

bx_scfa_individual_long <- ggplot(scfa_individual_long, aes(x = symptoms, y = `butyrate/isobutyrate (4:0)`)) +
  geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
  geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
  geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
  scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
  labs(x = "MC activity", y = "Log(butyrate/isobutyrate (4:0))") + 
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "none") +
  coord_flip()

bx_scfa_individual_long

# ggsave(filename = "bx_scfa_individual_long.pdf",
#        plot = bx_scfa_individual_long, units = "in", width=10, height=4, dpi = 1200)

scfa_individual_long$symptoms <- factor(scfa_individual_long$symptoms, levels = c("remission", "active"))
lm_scfa_individual_long_log <- lm(`butyrate/isobutyrate (4:0)` ~ symptoms + age + sex + bmi, data = scfa_individual_long)
summary(lm_scfa_individual_long_log)



# ####-------------- X Microbial EC: serine palmitoyltransferase (Cross-sectional) --------------####
# ## Correlation between serine_palmitoyltransferase and Dihydroceramides
# serine_palmitoyltransferase <- EC_full %>% select(matches("palmitoyltransferase"))
# serine_palmitoyltransferase$id <- rownames(serine_palmitoyltransferase)
# serine_palmitoyltransferase <- serine_palmitoyltransferase %>%
#   mutate(id = case_when(
#     id == "MC058A" ~ "MC058",
#     id == "MC062A" ~ "MC062",
#     id == "MC063A" ~ "MC063",
#     id == "MC075A" ~ "MC075",
#     id == "MC139" ~ "MC139A",
#     id == "MC167" ~ "MC167A",
#     id == "MC244" ~ "MC244A",
#     id == "MC408B" ~ "MC408A",
#     TRUE ~ id  # Keep the original value if no match
#   ))
# rownames(serine_palmitoyltransferase) <- serine_palmitoyltransferase$id
# 
# 
# Dihydroceramides <- metabolite %>% select("52604")
# Dihydroceramides$id <- rownames(Dihydroceramides)
# 
# Dihydroceramides_aMC <- Dihydroceramides %>% inner_join(select(activeMC_baseline, id)) %>% mutate(mc_all_label = "MC")
# Dihydroceramides_HC <- Dihydroceramides %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control")
# Dihydroceramides_CD <- Dihydroceramides %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea")
# Dihydroceramides <- rbind(Dihydroceramides_aMC,Dihydroceramides_HC,Dihydroceramides_CD)
# Dihydroceramides <- Dihydroceramides %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 
# 
# Dihydroceramides_filled <- Dihydroceramides %>% filter(!is.na(bmi))
# Dihydroceramides_unfilled <- Dihydroceramides %>% filter(is.na(bmi))
# Dihydroceramides_unfilled$ID <- substr(Dihydroceramides_unfilled$id, 1, 5)
# Dihydroceramides_unfilled <- Dihydroceramides_unfilled %>% select(-age, -sex, -bmi)
# Dihydroceramides_unfilled <- Dihydroceramides_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
# Dihydroceramides <- rbind(Dihydroceramides_filled, Dihydroceramides_unfilled)
# Dihydroceramides <- Dihydroceramides %>% rename("Dihydroceramides" = "52604")
# 
# 
# ## serine_palmitoyltransferase & Dihydroceramides scatter plot
# serine_palmitoyltransferase_Dihydroceramides <- inner_join(serine_palmitoyltransferase, Dihydroceramides)
# 
# serine_palmitoyltransferase_Dihydroceramides$mc_all_label <- factor(serine_palmitoyltransferase_Dihydroceramides$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"),
#                                                     labels = c("MC", "Healthy control", "Chronic diarrhea"))
# 
# serine_palmitoyltransferase_Dihydroceramides_1 <- ggplot(serine_palmitoyltransferase_Dihydroceramides, aes(x = Dihydroceramides, y = `2.3.1.50: Serine C-palmitoyltransferase`, color = mc_all_label)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   labs(x = "Relative abundance of N-palmitoyl-sphinganine (d18:0/16:0)", y = "EC 2.3.1.50: Serine C-palmitoyltransferase", color = "Disease type") +
#   theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
#   scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
#   coord_cartesian(xlim = c(0, 35), ylim = c(0, 0.0008)) 
# # annotate("text", x = 10, y = 4.5, label = "Spearman rho = 0.192, FDR = 0.003", size = 4, color = "black")
# # ggsave(filename = "veillonella_parvula_ceramide_1.png",
# #        plot = veillonella_parvula_ceramide_1, units = "in", width=7, height=5, dpi = 600)
# 
# serine_palmitoyltransferase_Dihydroceramides_1
# 
# serine_palmitoyltransferase_Dihydroceramides$mc_all_label <- factor(serine_palmitoyltransferase_Dihydroceramides$mc_all_label, levels = c("MC","Healthy control", "Chronic diarrhea"))
# lm_serine_palmitoyltransferase_Dihydroceramides <- lm(`2.3.1.50: Serine C-palmitoyltransferase` ~ Dihydroceramides, data = serine_palmitoyltransferase_Dihydroceramides)
# summary(lm_serine_palmitoyltransferase_Dihydroceramides)


# ####-------------- Microbial EC: phospholipase A1 & A2 (Cross-sectional) --------------####
# phospholipase_A1 <- EC_full %>% select(matches("phospholipase A"))
# phospholipase_A1$id <- rownames(phospholipase_A1)
# phospholipase_A1 <- phospholipase_A1%>%
#   mutate(id = case_when(
#     id == "MC058A" ~ "MC058",
#     id == "MC062A" ~ "MC062",
#     id == "MC063A" ~ "MC063",
#     id == "MC075A" ~ "MC075",
#     id == "MC139" ~ "MC139A",
#     id == "MC167" ~ "MC167A",
#     id == "MC244" ~ "MC244A",
#     id == "MC408B" ~ "MC408A",
#     TRUE ~ id  # Keep the original value if no match
#   ))
# rownames(phospholipase_A1) <- phospholipase_A1$id
# 
# 
# # Impute half minimum for 0 values
# phospholipase_A1 <- phospholipase_A1 %>%
#   mutate(`3.1.1.32: Phospholipase A(1)` = ifelse(`3.1.1.32: Phospholipase A(1)` == 0, 
#                                                  min(`3.1.1.32: Phospholipase A(1)`[`3.1.1.32: Phospholipase A(1)` > 0]) / 2, 
#                                                  `3.1.1.32: Phospholipase A(1)`))
# 
# 
# phospholipase_A1_aMC <- phospholipase_A1 %>% inner_join(select(activeMC_baseline, id)) %>% mutate(mc_all_label = "MC")
# phospholipase_A1_HC <- phospholipase_A1 %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control")
# phospholipase_A1_CD <- phospholipase_A1 %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea")
# phospholipase_A1 <- rbind(phospholipase_A1_aMC,phospholipase_A1_HC,phospholipase_A1_CD)
# phospholipase_A1 <- phospholipase_A1 %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 
# 
# phospholipase_A1_filled <- phospholipase_A1 %>% filter(!is.na(bmi))
# phospholipase_A1_unfilled <- phospholipase_A1 %>% filter(is.na(bmi))
# phospholipase_A1_unfilled$ID <- substr(phospholipase_A1_unfilled$id, 1, 5)
# phospholipase_A1_unfilled <- phospholipase_A1_unfilled %>% select(-age, -sex, -bmi)
# phospholipase_A1_unfilled <- phospholipase_A1_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
# phospholipase_A1 <- rbind(phospholipase_A1_filled, phospholipase_A1_unfilled)
# 
# 
# 
# 
# phospholipase_A1_log <- phospholipase_A1
# phospholipase_A1_log$`3.1.1.32: Phospholipase A(1)` <- log(phospholipase_A1_log$`3.1.1.32: Phospholipase A(1)`)
# 
# 
# 
# ## Box plot
# phospholipase_A1_log$mc_all_label <- factor(phospholipase_A1_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Healthy control"))
# 
# bx_cross_phospholipase_A1 <- ggplot(phospholipase_A1_log, aes(x = mc_all_label, y = `3.1.1.32: Phospholipase A(1)`)) +
#   geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
#   geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
#   geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
#   scale_color_manual(values = c("#900C3F", "#F9D937", "#189AF9")) + 
#   labs(x = "Disease type", y = "Log(3.1.1.32: Phospholipase A(1))") + 
#   theme_minimal() +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.position = "none") +
#   coord_flip()
# 
# bx_cross_phospholipase_A1
# 
# # ggsave(filename = "bx_cross_phospholipase_A1.png",
# #        plot = bx_cross_phospholipase_A1, units = "in", width=10, height=4, dpi = 600)
# 
# 
# phospholipase_A1_log$mc_all_label <- factor(phospholipase_A1_log$mc_all_label, levels = c("MC","Healthy control", "Chronic diarrhea"))
# lm_phospholipase_A1_log <- lm(`3.1.1.32: Phospholipase A(1)` ~ mc_all_label + age + sex + bmi, data = phospholipase_A1_log)
# summary(lm_phospholipase_A1_log)
# 
# phospholipase_A1_cross <- phospholipase_A1
# 
# 
# 
# 
# 
# ####--------------- Microbial EC: phospholipase A1 (longitudinal)  ----------------------#####
# phospholipase_A1 <- EC_full %>% select(matches("phospholipase A"))
# phospholipase_A1$id <- rownames(phospholipase_A1)
# phospholipase_A1 <- phospholipase_A1%>%
#   mutate(id = case_when(
#     id == "MC058A" ~ "MC058",
#     id == "MC062A" ~ "MC062",
#     id == "MC063A" ~ "MC063",
#     id == "MC075A" ~ "MC075",
#     id == "MC139" ~ "MC139A",
#     id == "MC167" ~ "MC167A",
#     id == "MC244" ~ "MC244A",
#     id == "MC408B" ~ "MC408A",
#     TRUE ~ id  # Keep the original value if no match
#   ))
# rownames(phospholipase_A1) <- phospholipase_A1$id
# 
# # Impute half minimun for 0 values
# phospholipase_A1 <- phospholipase_A1 %>%
#   mutate(`3.1.1.32: Phospholipase A(1)` = ifelse(`3.1.1.32: Phospholipase A(1)` == 0, 
#                                                  min(`3.1.1.32: Phospholipase A(1)`[`3.1.1.32: Phospholipase A(1)` > 0]) / 2, 
#                                                  `3.1.1.32: Phospholipase A(1)`))
# 
# 
# phospholipase_A1_long <- phospholipase_A1 %>% inner_join(select(metadata_long, id, age, sex, bmi, symptoms))
# 
# phospholipase_A1_long_log <- phospholipase_A1_long
# phospholipase_A1_long_log$`3.1.1.32: Phospholipase A(1)` <- log(phospholipase_A1_long_log$`3.1.1.32: Phospholipase A(1)`)
# 
# 
# ## Box plot
# phospholipase_A1_long_log$symptoms <- factor(phospholipase_A1_long_log$symptoms, levels = c("active_diarrhea", "no_diarrhea"))
# 
# bx_phospholipase_A1_long <- ggplot(phospholipase_A1_long_log, aes(x = symptoms, y = `3.1.1.32: Phospholipase A(1)`)) +
#   geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
#   geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
#   geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
#   scale_color_manual(values = c("#900C3F", "antiquewhite4")) + 
#   labs(x = "Disease type", y = "log(3.1.1.32: Phospholipase A(1))") + 
#   theme_minimal() +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.position = "none") +
#   coord_flip()
# 
# bx_phospholipase_A1_long
# 
# # ggsave(filename = "bx_phospholipase_A1_long.png",
# #        plot = bx_phospholipase_A1_long, units = "in", width=10, height=4, dpi = 600)
# 
# phospholipase_A1_long_log$symptoms <- factor(phospholipase_A1_long_log$symptoms, levels = c("no_diarrhea", "active_diarrhea"))
# lm_phospholipase_A1_long_log <- lm(`3.1.1.32: Phospholipase A(1)` ~ symptoms + age + sex + bmi, data = phospholipase_A1_long_log)
# summary(lm_phospholipase_A1_long_log)
# 
# 
# 
# 
# 
# 
# ####--------------- LCFA analysis: cross-sectional  ----------------------#####
# lcfa <- subset(all, grepl("Long Chain Polyunsaturated", pathway, ignore.case = TRUE))
# 
# n3lcfa <- subset(lcfa, grepl("n3", metabolite, ignore.case = TRUE)) 
# n3lcfa <- subset(n3lcfa, !grepl("or", metabolite, ignore.case = TRUE))
# 
# n6lcfa <- subset(lcfa, grepl("n6", metabolite, ignore.case = TRUE)) 
# n6lcfa <- subset(n6lcfa, !grepl("or", metabolite, ignore.case = TRUE))
# 
# 
# n3lcfa <- metabolite %>% select("32504", "44675", "18467", "57651")
# n3lcfa <- n3lcfa %>% mutate(n3lcfa = `32504` + `44675` + `18467` + `57651`) 
# n3lcfa$id <- rownames(n3lcfa)
# 
# n6lcfa <- metabolite %>% select("1110", "32980", "17805", "37478", "57652", "1105", "32415")
# n6lcfa <- n6lcfa %>% mutate(n6lcfa = `1110` + `32980` + `17805` + `37478` + `57652` + `1105` + `32415`) 
# n6lcfa$id <- rownames(n6lcfa)
# 
# 
# n6_n3lcfa <- as.data.frame(n6lcfa[, 8] / n3lcfa[, 5])
# rownames(n6_n3lcfa) <- rownames(n3lcfa)
# n6_n3lcfa <- n6_n3lcfa %>% rename("n6_n3lcfa" = "n6lcfa[, 8]/n3lcfa[, 5]")
# n6_n3lcfa$id <- rownames(n6_n3lcfa)
# 
# n6_n3lcfa <- n6_n3lcfa %>% inner_join(select(BAratios_cross, id, age, sex, bmi, mc_all_label))
# n6_n3lcfa_log <- n6_n3lcfa
# n6_n3lcfa_log$n6_n3lcfa <- log(n6_n3lcfa_log$n6_n3lcfa)
# 
# ## Box plot
# n6_n3lcfa_log$mc_all_label <- factor(n6_n3lcfa_log$mc_all_label, levels = c("Healthy control", "MC", "Chronic diarrhea"))
# 
# bx_cross_n6_n3lcfa <- ggplot(n6_n3lcfa_log, aes(x = mc_all_label, y = n6_n3lcfa)) +
#   geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
#   geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
#   geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
#   scale_color_manual(values = c("#189AF9", "#900C3F", "#F9D937")) + 
#   labs(x = "Disease type", y = "Log(Omega-6 PUFAs/Omega-3 PUFAs)") + 
#   theme_minimal() +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.position = "none") +
#   coord_flip()
# 
# bx_cross_n6_n3lcfa
# 
# # ggsave(filename = "bx_cross_n6_n3lcfa.png",
# #        plot = bx_cross_n6_n3lcfa, units = "in", width=10, height=4, dpi = 600)
# 
# 
# n6_n3lcfa_log$mc_all_label <- factor(n6_n3lcfa_log$mc_all_label, levels = c("MC","Healthy control", "Chronic diarrhea"))
# lm_n6_n3lcfa_log <- lm(n6_n3lcfa ~ mc_all_label + age + sex + bmi, data = n6_n3lcfa_log)
# summary(lm_n6_n3lcfa_log)
# 
# 
# 
# ####--------------- LCFA analysis: longitudinal  ----------------------#####
# lcfa <- subset(all, grepl("Long Chain Polyunsaturated", pathway, ignore.case = TRUE))
# 
# n3lcfa <- subset(lcfa, grepl("n3", metabolite, ignore.case = TRUE)) 
# n3lcfa <- subset(n3lcfa, !grepl("or", metabolite, ignore.case = TRUE))
# 
# n6lcfa <- subset(lcfa, grepl("n6", metabolite, ignore.case = TRUE)) 
# n6lcfa <- subset(n6lcfa, !grepl("or", metabolite, ignore.case = TRUE))
# 
# 
# n3lcfa <- metabolite %>% select("32504", "44675", "18467", "57651")
# n3lcfa <- n3lcfa %>% mutate(n3lcfa = `32504` + `44675` + `18467` + `57651`) 
# n3lcfa$id <- rownames(n3lcfa)
# 
# n6lcfa <- metabolite %>% select("1110", "32980", "17805", "37478", "57652", "1105", "32415")
# n6lcfa <- n6lcfa %>% mutate(n6lcfa = `1110` + `32980` + `17805` + `37478` + `57652` + `1105` + `32415`) 
# n6lcfa$id <- rownames(n6lcfa)
# 
# 
# n6_n3lcfa <- as.data.frame(n6lcfa[, 8] / n3lcfa[, 5])
# rownames(n6_n3lcfa) <- rownames(n3lcfa)
# n6_n3lcfa <- n6_n3lcfa %>% rename("n6_n3lcfa" = "n6lcfa[, 8]/n3lcfa[, 5]")
# n6_n3lcfa$id <- rownames(n6_n3lcfa)
# 
# n6_n3lcfa_long <- n6_n3lcfa %>% inner_join(select(metadata_long, id, age, sex, bmi, symptoms)) 
# 
# 
# n6_n3lcfa_long_log <- n6_n3lcfa_long
# n6_n3lcfa_long_log$n6_n3lcfa <- log(n6_n3lcfa_long_log$n6_n3lcfa)
# 
# ## Box plot
# n6_n3lcfa_long_log$symptoms <- factor(n6_n3lcfa_long_log$symptoms, levels = c("no_diarrhea", "active_diarrhea"))
# 
# bx_n6_n3lcfa_long <- ggplot(n6_n3lcfa_long_log, aes(x = symptoms, y = n6_n3lcfa)) +
#   geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
#   geom_jitter(aes(color = symptoms), width = 0.2, size = 2) +  
#   geom_point(aes(color = symptoms), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
#   scale_color_manual(values = c("antiquewhite4", "#900C3F")) + 
#   labs(x = "Disease type", y = "Log(Omega-6 PUFAs/Omega-3 PUFAs)") + 
#   theme_minimal() +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.position = "none") +
#   coord_flip()
# 
# bx_n6_n3lcfa_long
# 
# # ggsave(filename = "bx_n6_n3lcfa_long.png",
# #        plot = bx_n6_n3lcfa_long, units = "in", width=10, height=4, dpi = 600)
# 
# n6_n3lcfa_long_log$symptoms <- factor(n6_n3lcfa_long_log$symptoms, levels = c("no_diarrhea", "active_diarrhea"))
# lm_n6_n3lcfa_long_log <- lm(n6_n3lcfa ~ symptoms + age + sex + bmi, data = n6_n3lcfa_long_log)
# summary(lm_n6_n3lcfa_long_log)




# ####--------------- Bile acids analysis: Cross-sectional  ----------------------#####
# primaryBA <- subset(all, grepl("Primary Bile Acid", pathway, ignore.case = TRUE))
# secondaryBA <- subset(all, grepl("Secondary Bile Acid", pathway, ignore.case = TRUE))
# BA <- metabolite
# 
# ## Primary:Secondary: cholate/deoxycholate
# cholate <- BA %>% select("22842", "52973")
# cholate <- cholate %>% mutate(cholate = `22842` + `52973`) # Primary BA
# 
# deoxycholate <- BA %>% select("1114", "63688", "52969") 
# deoxycholate <- deoxycholate %>% mutate(deoxycholate = `1114` + `63688` + `52969`) # Secondary BA
# 
# cholate_deoxycholate <- as.data.frame(cholate[, 3] / deoxycholate[, 4])
# rownames(cholate_deoxycholate) <- rownames(cholate)
# cholate_deoxycholate <- cholate_deoxycholate %>% rename("cholate_deoxycholate" = "cholate[, 3]/deoxycholate[, 4]")
# cholate_deoxycholate$id <- rownames(cholate_deoxycholate)
# 
# 
# ## Primary:Secondary: chenodeoxycholate/lithocholate
# # [chenodeoxycholic acid sulfate (1) + chenodeoxycholic acid sulfate (2)]/[lithocholate sulfate (1) + lithocholic acid sulfate (2)]
# chenodeoxycholate <- BA %>% select("63607", "63608") %>% mutate(chenodeoxycholate = `63607` + `63608`) # Primary BA
# lithocholate <- BA %>% select("1483", "62526", "62527") %>% mutate(lithocholate = `1483` + `62526` + `62527`) # Secondary BA
# 
# chenodeoxycholate_lithocholate <- as.data.frame(chenodeoxycholate[, 3] / lithocholate[, 4])
# rownames(chenodeoxycholate_lithocholate) <- rownames(cholate)
# chenodeoxycholate_lithocholate <- chenodeoxycholate_lithocholate %>% rename("chenodeoxycholate_lithocholate" = "chenodeoxycholate[, 3]/lithocholate[, 4]")
# chenodeoxycholate_lithocholate$id <- rownames(chenodeoxycholate_lithocholate)
# 
# 
# ## Primary:Secondary: chenodeoxycholate/ursodeoxycholate
# # [chenodeoxycholic acid sulfate (1) + chenodeoxycholic acid sulfate (2)]/[ursodeoxycholate + ursodeoxycholate sulfate (1)]
# chenodeoxycholate <- BA %>% select("63607", "63608") %>% mutate(chenodeoxycholate = `63607` + `63608`) # Primary BA
# ursodeoxycholate <- BA %>% select("1605", "52970") %>% mutate(ursodeoxycholate = `1605` + `52970`) # Secondary BA
# 
# chenodeoxycholate_ursodeoxycholate <- as.data.frame(chenodeoxycholate[, 3] / ursodeoxycholate[, 3])
# rownames(chenodeoxycholate_ursodeoxycholate) <- rownames(cholate)
# chenodeoxycholate_ursodeoxycholate <- chenodeoxycholate_ursodeoxycholate %>% rename("chenodeoxycholate_ursodeoxycholate" = "chenodeoxycholate[, 3]/ursodeoxycholate[, 3]")
# chenodeoxycholate_ursodeoxycholate$id <- rownames(chenodeoxycholate_ursodeoxycholate)
# 
# 
# ### Combining all Primary:Secondary
# BAratios <- left_join(cholate_deoxycholate, chenodeoxycholate_lithocholate)
# BAratios <- BAratios %>% left_join(chenodeoxycholate_ursodeoxycholate)
# rownames(BAratios) <- BAratios$id
# BAratios <- BAratios %>% select(id, 1,3,4)
# 
# ## Cross-sectional
# BAratios_cross_aMC <- BAratios %>% inner_join(select(activeMC_baseline, id)) %>% mutate(mc_all_label = "MC")
# BAratios_cross_HC <- BAratios %>% inner_join(select(HC_baseline, id)) %>% mutate(mc_all_label = "Healthy control")
# BAratios_cross_CD <- BAratios %>% inner_join(select(CD_baseline, id)) %>% mutate(mc_all_label = "Chronic diarrhea")
# BAratios_cross <- rbind(BAratios_cross_aMC,BAratios_cross_HC,BAratios_cross_CD)
# BAratios_cross <- BAratios_cross %>% left_join(select(alpha_with_disease_status, ID, age, sex, bmi), by = c("id" = "ID")) 
# BAratios_cross_filled <- BAratios_cross %>% filter(!is.na(BAratios_cross$bmi))
# BAratios_cross_unfilled <- BAratios_cross %>% filter(is.na(BAratios_cross$bmi))
# BAratios_cross_unfilled$ID <- substr(BAratios_cross_unfilled$id, 1, 5)
# BAratios_cross_unfilled <- BAratios_cross_unfilled %>% select(-age, -sex, -bmi)
# BAratios_cross_unfilled <- BAratios_cross_unfilled %>% left_join(select(cross_baseline, ID, age, sex, bmi)) %>% select(-ID)
# BAratios_cross <- rbind(BAratios_cross_filled, BAratios_cross_unfilled)
# 
# 
# BAratios_cross_log <- BAratios_cross
# BAratios_cross_log$cholate_deoxycholate <- log(BAratios_cross_log$cholate_deoxycholate)
# BAratios_cross_log$chenodeoxycholate_lithocholate <- log(BAratios_cross_log$chenodeoxycholate_lithocholate)
# BAratios_cross_log$chenodeoxycholate_ursodeoxycholate <- log(BAratios_cross_log$chenodeoxycholate_ursodeoxycholate)
# 
# ### Perform the linear regression
# BAratios_cross_log$mc_all_label <- factor(BAratios_cross_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Healthy control"))
# BAratios_cross_log$sex <- factor(BAratios_cross_log$sex, levels = c("Female", "Male"))
# 
# ## cholate_deoxycholate
# lm_cholate_deoxycholate <- lm(cholate_deoxycholate ~ mc_all_label + age + sex + bmi, data = BAratios_cross_log)
# summary(lm_cholate_deoxycholate)
# 
# ## chenodeoxycholate/lithocholate
# lm_chenodeoxycholate_lithocholate <- lm(chenodeoxycholate_lithocholate ~ mc_all_label + age + sex + bmi, data = BAratios_cross_log)
# summary(lm_chenodeoxycholate_lithocholate)
# 
# ## chenodeoxycholate/ursodeoxycholate
# lm_chenodeoxycholate_ursodeoxycholate <- lm(chenodeoxycholate_ursodeoxycholate ~ mc_all_label + age + sex + bmi, data = BAratios_cross_log)
# summary(lm_chenodeoxycholate_ursodeoxycholate)
# 
# 
# 
# ## Primary:Secondary: All
# primaryBA_features <- primaryBA$feature
# primary <- BA %>% select(all_of(primaryBA_features))
# primary <- primary %>% mutate(primary = rowSums(across(everything())))
# 
# secondaryBA_features <- secondaryBA$feature
# secondary <- BA %>% select(all_of(secondaryBA_features))
# secondary <- secondary %>% mutate(secondary = rowSums(across(everything())))
# 
# primary_secondary <- as.data.frame(primary[, 12] / secondary[, 39])
# rownames(primary_secondary) <- rownames(primary)
# primary_secondary <- primary_secondary %>% rename("primary_secondary" = "primary[, 12]/secondary[, 39]")
# primary_secondary$id <- rownames(primary_secondary)
# 
# primary_secondary <- primary_secondary %>% inner_join(select(BAratios_cross, id, age, sex, bmi, mc_all_label))
# primary_secondary_log <- primary_secondary
# primary_secondary_log$primary_secondary <- log(primary_secondary_log$primary_secondary)
# 
# primary_secondary_log$mc_all_label <- factor(primary_secondary_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Healthy control"))
# primary_secondary_log$sex <- factor(primary_secondary_log$sex, levels = c("Female", "Male"))
# 
# primary_secondary_log$mc_all_label <- factor(primary_secondary_log$mc_all_label, levels = c("MC", "Healthy control", "Chronic diarrhea"))
# lm_primary_secondary <- lm(primary_secondary ~ mc_all_label + age + sex + bmi, data = primary_secondary_log)
# summary(lm_primary_secondary)
# 
# ## Box plot
# primary_secondary_log$mc_all_label <- factor(primary_secondary_log$mc_all_label, levels = c("Healthy control", "MC", "Chronic diarrhea"))
# 
# bx_cross_primary_secondary <- ggplot(primary_secondary_log, aes(x = mc_all_label, y = primary_secondary)) +
#   geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
#   geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
#   geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
#   scale_color_manual(values = c("#189AF9", "#900C3F", "#F9D937")) + 
#   labs(x = "Disease type", y = "Log(Primary BAs / Secondary BAs)") + 
#   theme_minimal() +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.position = "none") +
#   coord_flip()
# 
# bx_cross_primary_secondary
# 
# # ggsave(filename = "bx_cross_primary_secondary.png",
# #        plot = bx_cross_primary_secondary, units = "in", width=10, height=4, dpi = 600)
# 
# 
# 
# 
# 
# ### Conjugated/unconjugated: cholate/deoxycholate
# conjugatedBA <- subset(all, grepl("Bile Acid", pathway, ignore.case = TRUE)) 
# conjugatedBA <- subset(conjugatedBA, grepl("glyco|tauro|sulf", metabolite, ignore.case = TRUE)) 
# 
# unconjugatedBA <- subset(all, grepl("Bile Acid", pathway, ignore.case = TRUE)) 
# unconjugatedBA <- subset(unconjugatedBA, !grepl("glyco|tauro|sulf", metabolite, ignore.case = TRUE))
# 
# ## Conjugated/unconjugated: All
# conjugatedBA_features <- conjugatedBA$feature
# conjugated <- BA %>% select(all_of(conjugatedBA_features))
# conjugated <- conjugated %>% mutate(conjugated = rowSums(across(everything())))
# 
# unconjugatedBA_features <- unconjugatedBA$feature
# unconjugated <- BA %>% select(all_of(unconjugatedBA_features))
# unconjugated <- unconjugated %>% mutate(unconjugated = rowSums(across(everything())))
# 
# conjugated_unconjugated <- as.data.frame(conjugated[, 31] / unconjugated[, 20])
# rownames(conjugated_unconjugated) <- rownames(conjugated)
# conjugated_unconjugated <- conjugated_unconjugated %>% rename("conjugated_unconjugated" = "conjugated[, 31]/unconjugated[, 20]")
# conjugated_unconjugated$id <- rownames(conjugated_unconjugated)
# 
# conjugated_unconjugated <- conjugated_unconjugated %>% inner_join(select(BAratios_cross, id, age, sex, bmi, mc_all_label))
# conjugated_unconjugated_log <- conjugated_unconjugated
# conjugated_unconjugated_log$conjugated_unconjugated <- log(conjugated_unconjugated_log$conjugated_unconjugated)
# 
# conjugated_unconjugated_log$mc_all_label <- factor(conjugated_unconjugated_log$mc_all_label, levels = c("MC", "Chronic diarrhea", "Healthy control"))
# conjugated_unconjugated_log$sex <- factor(conjugated_unconjugated_log$sex, levels = c("Female", "Male"))
# 
# lm_conjugated_unconjugated_log <- lm(conjugated_unconjugated ~ mc_all_label + age + sex + bmi, data = conjugated_unconjugated_log)
# summary(lm_conjugated_unconjugated_log)
# 
# ## Box plot
# conjugated_unconjugated_log$mc_all_label <- factor(conjugated_unconjugated_log$mc_all_label, levels = c("Healthy control", "MC", "Chronic diarrhea"))
# 
# bx_cross_conjugated_unconjugated <- ggplot(conjugated_unconjugated_log, aes(x = mc_all_label, y = conjugated_unconjugated)) +
#   geom_boxplot(color = "black", fill = NA, outlier.shape = NA) +  
#   geom_jitter(aes(color = mc_all_label), width = 0.2, size = 2) +  
#   geom_point(aes(color = mc_all_label), position = position_jitter(width = 0.2), size = 2, shape = 16) + 
#   scale_color_manual(values = c("#189AF9", "#900C3F", "#F9D937")) + 
#   labs(x = "Disease type", y = "Log(Conjugated BAs / Uncojugated BAs)") + 
#   theme_minimal() +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.position = "none") +
#   coord_flip()
# 
# bx_cross_conjugated_unconjugated
# 
# # ggsave(filename = "bx_cross_conjugated_unconjugated.png",
# #        plot = bx_cross_conjugated_unconjugated, units = "in", width=10, height=4, dpi = 600)


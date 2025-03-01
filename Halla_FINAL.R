rm(list=ls())

setwd("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R")

library(readxl)
library(dplyr)
library(stringr)
library(lmerTest)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lemon)
library(writexl)
library(openxlsx)
library(pheatmap)
library(grid)
library(gridExtra)
library("RColorBrewer")
library(gtable)
library(ggrepel)

#### Differentially abundant metabolite in MC vs HC ####
index <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/index.xlsx")
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
sig <- sig %>% 
  filter(!is.na(pathway)) %>% 
  select(coef, qval, pathway, metabolite) %>%
  rename("coef_HC" = "coef",
         "qval_HC" = "qval",
         "pathway_HC" = "pathway")


#### Differentially abundant metabolite in MC vs CD ####
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
sig_aMC_CD <- sig_aMC_CD %>% 
  filter(!is.na(pathway)) %>% 
  select(coef, qval, pathway, metabolite) %>%
  rename("coef_CD" = "coef",
         "qval_CD" = "qval",
         "pathway_CD" = "pathway")


#### Differentially abundant metabolite in MC compared to both CD and HC ####
r_diff_1 <- sig %>% inner_join(select(sig_aMC_CD, coef_CD, qval_CD, metabolite), by = "metabolite")

# Select rows where coef_HC and coef_CD have the same sign
r_diff_same_sign <- r_diff_1 %>% 
  filter((coef_HC > 0 & coef_CD > 0) | (coef_HC < 0 & coef_CD < 0))
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/halla_target.xlsx"
# write.xlsx(r_diff_same_sign, file_path, rowNames = TRUE)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/altered_met_114.xlsx"
# write.xlsx(r_diff_same_sign, file_path, rowNames = TRUE)


#### Volcano plot for MC vs HC
index <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/index.xlsx")
pathway <- index[, c(3,6)]
pathway$COMP_ID <- as.character(pathway$COMP_ID)
chemical <- index[, c(3,11)]
chemical$COMP_ID <- as.character(chemical$COMP_ID)

volcano_cross <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/all_results.tsv")
volcano_cross_hc <- volcano_cross %>% filter(value=="Control") %>% select(1,4,5,6,8,9) 
volcano_cross_hc$coef <- -(volcano_cross_hc$coef)
volcano_cross_hc$feature<- substr(volcano_cross_hc$feature, 2,6)

volcano_cross_hc <- volcano_cross_hc %>% left_join(pathway, by = c("feature" = "COMP_ID"))
volcano_cross_hc <- volcano_cross_hc %>% left_join(chemical, by = c("feature" = "COMP_ID"))
volcano_cross_hc <- volcano_cross_hc %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

volcano_cross_hc <- volcano_cross_hc %>% 
  mutate(diffexpressed = case_when(
    coef > 0 & qval < 0.25 ~ "up",
    coef < 0 & qval < 0.25 ~ "down",
    TRUE ~ "no"
  )) %>% 
  filter(!is.na(pathway))


volcano_cross_hc <- as.data.frame(volcano_cross_hc)

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially abundant species (NA in case they are not)
volcano_cross_hc$delabel <- ifelse(volcano_cross_hc$metabolite %in% head(volcano_cross_hc[order(volcano_cross_hc$qval), "metabolite"], 20), 
                                   volcano_cross_hc$metabolite, 
                                   NA)


fig_volcano_cross_hc <- 
  ggplot(data = volcano_cross_hc, aes(x = coef, y = -log10(qval), label = delabel)) +
  geom_vline(xintercept = c(0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.25), col = "gray", linetype = 'dashed') +
  geom_point(aes(color = diffexpressed), size = 2, alpha = ifelse(volcano_cross_hc$diffexpressed == "no", 0.5, 1)) + 
  scale_color_manual(values = c("#96DED1", "#808080", "#E35335"),
                     labels = c("Less abundant in MC", "Not significant", "More abundant in MC")) +
  labs(color = 'Change in abundance', 
       x = expression("coefficient"), y = expression("-log"[10]*"(q-value)")) + 
  scale_x_continuous(limits = c(-2.7, 2.7), breaks = seq(-1.5,1.5)) + # Adjust x-axis limits and intervals
  geom_text_repel(max.overlaps = Inf, size = 4) + # To show all labels
  theme_set(theme_classic(base_size = 12) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,5,0,0), size = rel(1.1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(5,0,0,0), size = rel(1.1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))

fig_volcano_cross_hc

# ggsave(filename = "fig_met_volcano_cross_hc.pdf",
#        plot = fig_volcano_cross_hc, units = "in", width=12, height=9, dpi = 1200)




#### Volcano plot for MC vs CD
volcano_cross_cd <- volcano_cross %>% filter(value=="Diarrhea") %>% select(1,4,5,6,8,9) 
volcano_cross_cd$coef <- -(volcano_cross_cd$coef)
volcano_cross_cd$feature<- substr(volcano_cross_cd$feature, 2,6)

volcano_cross_cd <- volcano_cross_cd %>% left_join(pathway, by = c("feature" = "COMP_ID"))
volcano_cross_cd <- volcano_cross_cd %>% left_join(chemical, by = c("feature" = "COMP_ID"))
volcano_cross_cd <- volcano_cross_cd %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

volcano_cross_cd <- volcano_cross_cd %>% 
  mutate(diffexpressed = case_when(
    coef > 0 & qval < 0.25 ~ "up",
    coef < 0 & qval < 0.25 ~ "down",
    TRUE ~ "no"
  )) %>% 
  filter(!is.na(pathway))


volcano_cross_cd <- as.data.frame(volcano_cross_cd)

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially abundant species (NA in case they are not)
volcano_cross_cd$delabel <- ifelse(volcano_cross_cd$metabolite %in% head(volcano_cross_cd[order(volcano_cross_cd$qval), "metabolite"], 20), 
                                   volcano_cross_cd$metabolite, 
                                   NA)


fig_volcano_cross_cd <- 
  ggplot(data = volcano_cross_cd, aes(x = coef, y = -log10(qval), label = delabel)) +
  geom_vline(xintercept = c(0), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.25), col = "gray", linetype = 'dashed') +
  geom_point(aes(color = diffexpressed), size = 2, alpha = ifelse(volcano_cross_cd$diffexpressed == "no", 0.5, 1)) + 
  scale_color_manual(values = c("#96DED1", "#808080", "#E35335"),
                     labels = c("Less abundant in MC", "Not significant", "More abundant in MC")) +
  labs(color = 'Change in abundance', 
       x = expression("coefficient"), y = expression("-log"[10]*"(q-value)")) + 
  scale_x_continuous(limits = c(-1.5, 1.5), breaks = seq(-0.75,0.75)) + # Adjust x-axis limits and intervals
  geom_text_repel(max.overlaps = Inf, size = 4) + # To show all labels
  theme_set(theme_classic(base_size = 12) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,5,0,0), size = rel(1.1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(5,0,0,0), size = rel(1.1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))

fig_volcano_cross_cd

# ggsave(filename = "fig_met_volcano_cross_cd.pdf",
#        plot = fig_volcano_cross_cd, units = "in", width=12, height=9, dpi = 1200)


### Heatmap for the 114 altered metabolites 
altered_heat <- r_diff_same_sign %>% 
  rename("coef_MC.vs.Control_without_diarrhea" = "coef_HC",
         "qval_MC.vs.Control_without_diarrhea" = "qval_HC",
         "coef_MC.vs.Chronic_diarrhea" = "coef_CD",
         "qval_MC.vs.Chronic_diarrhea" = "qval_CD") %>%
  as.data.frame()

# altered_heat <- altered_heat %>% select(4,3,1,2,5,6)
# altered_heat <- altered_heat %>%
#   rename("Metabolite name" = "metabolite",
#          "Metabolite class" = "pathway_HC",
#          "coef_MC.vs.Controls_without_diarrhea" = "coef_MC.vs.Control_without_diarrhea",
#          "qval_MC.vs.Controls_without_diarrhea" = "qval_MC.vs.Control_without_diarrhea") %>%
#   arrange(desc(coef_MC.vs.Controls_without_diarrhea))
# 
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/online_table_1.xlsx"
# write.xlsx(altered_heat, file_path, rowNames = TRUE)
  

rownames(altered_heat) <- altered_heat$metabolite
altered_heat$metabolite <- NULL
altered_heat$pathway_HC <- NULL
b_altered <- altered_heat %>% select(1,3)
custom_order <- rownames(b_altered[order(-b_altered$coef_MC.vs.Control_without_diarrhea), ])
b_altered <- b_altered[custom_order,] %>% as.matrix()

b_altered_q <- altered_heat %>% select(2,4) 
b_altered_q <- b_altered_q[custom_order,] %>% as.matrix()

asterisks_matrix_altered <- matrix("", nrow = nrow(b_altered_q), ncol = ncol(b_altered_q))
asterisks_matrix_altered[b_altered_q < 0.05] <- "***"
asterisks_matrix_altered[b_altered_q >= 0.05 & b_altered_q < 0.1] <- "**"
asterisks_matrix_altered[b_altered_q >= 0.1 & b_altered_q < 0.25] <- "*"


# Create breaks with 0 in the middle
breaks <- seq(-1.8, 1.8, length.out = 51)

# Ensure the breaks include 0
mid <- which.min(abs(breaks - 0))

# Create a custom color palette
my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))

# Plot the heatmap with reordered rows
heatmap_altered_met <- pheatmap(b_altered,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                display_numbers = asterisks_matrix_altered,
                                fontsize_number = 10,
                                color = my_color_palette,
                                breaks = breaks,
                                annotation_legend = TRUE,
                                angle_col = 90)

# ggsave("heatmap_altered_met.pdf", heatmap_altered_met, width = 4.8, height = 25, dpi=1200, limitsize = FALSE)





#### HAllA targeted readin ####
# read in all_associations.txt got from running HALLA
r_diff <- read.table('/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/halla_residual/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")


####### Data cleaning #######
test_reduce2 <- r_diff

##pathway class
path <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/index.xlsx")
path <- path %>% select(SUB_PATHWAY, CHEMICAL_NAME)
colnames(path) <- c('pathway', 'Y_features')

# Remove rows where pathway is missing
path <- path[!is.na(path$pathway), ]

# Assign test_reduce2 to test
test <- test_reduce2

# Merge test and path dataframes based on Y_features
test2 <- merge(test, path, by = "Y_features", all.x = TRUE)
test2$X_features_factor <- as.factor(test2$X_features)
test2$Y_features_factor <- as.factor(test2$Y_features)


####### prepare for clustered heatmap #######
mc_species <- c("Intestinibacter_bartlettii", "Veillonella_dispar", "Veillonella_parvula", 
                "Haemophilus_parainfluenzae", "Clostridium_sp_AF20_17LB", 
                "Clostridium_spiroforme", "Veillonella_rogosae", "GGB3612_SGB4882")
test2 <- test2 %>%
  mutate(species_abundance = if_else(X_features_factor %in% mc_species, "MC_enriched", "Controls_enriched"))


mc_metabolites <- r_diff_same_sign %>% filter(coef_HC>0 & coef_CD>0)
mc_metabolites <- mc_metabolites$metabolite
test2 <- test2 %>%
  mutate(metabolite_abundance = if_else(Y_features_factor %in% mc_metabolites, "MC_enriched", "Controls_enriched"))



#######clustering by feature type#######
# Subset Relevant Columns
a <- test2[, colnames(test2) %in% c("X_features", "Y_features", "association")]
a1 <- test2[, colnames(test2) %in% c("X_features", "Y_features",  "q.values")]

# Reshape data
b <- a %>% pivot_wider(names_from = "Y_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("X_features") %>% as.matrix()
q_values_matrix <- a1 %>% pivot_wider(names_from = "Y_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("X_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation <- data.frame(Disease_type = test2$species_abundance[match(rownames(b), test2$X_features)])
col_annotation <- data.frame(Disease_type = test2$metabolite_abundance[match(colnames(b), test2$Y_features)],
                             Metabolite_class = test2$pathway[match(colnames(b), test2$Y_features)])

# Set row and column names for annotations
rownames(row_annotation) <- rownames(b)
rownames(col_annotation) <- colnames(b)

# Order rows based on group annotations
ordered_rows <- rownames(b)[order(row_annotation$Disease_type)]

# Order columns based on Disease_type and Pathway
ordered_cols <- colnames(b)[order(col_annotation$Disease_type, col_annotation$Metabolite_class)]

# Reorder matrix
b <- b[ordered_rows, ordered_cols]

# Reorder annotations to match the new order
row_annotation <- row_annotation[ordered_rows, , drop = FALSE]
col_annotation <- col_annotation[ordered_cols, , drop = FALSE]


# Define colors for the pathways using a more distinct and softer palette
unique_pathways <- unique(col_annotation$Metabolite_class)
pathway_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways))
names(pathway_colors) <- unique_pathways

# Define colors for the annotations
annotation_colors <- list(
  Disease_type = c(Controls_enriched = "#189AF9", MC_enriched = "#900C3F"),
  Metabolite_class = pathway_colors
)


asterisks_matrix <- matrix("", nrow = nrow(q_values_matrix), ncol = ncol(q_values_matrix))
asterisks_matrix[q_values_matrix < 0.05] <- "***"
asterisks_matrix[q_values_matrix >= 0.05 & q_values_matrix < 0.1] <- "**"
asterisks_matrix[q_values_matrix >= 0.1 & q_values_matrix < 0.25] <- "*"

# Create symmetric breaks around 0
breaks <- seq(-0.35, 0.35, length.out = 51)

# Check if 0 is in the middle
print(breaks[which.min(abs(breaks))]) 

# Find the index where 0 is located
mid <- which.min(abs(breaks - 0))

# Create the custom palette
my_color_palette <- c(
  colorRampPalette(c("navy", "white"))(mid - 1), 
  "white", 
  colorRampPalette(c("white", "firebrick3"))(length(breaks) - mid)
)

# Create the heatmap
heatmap <- pheatmap(b,
                    annotation_row = row_annotation,
                    annotation_col = col_annotation,
                    annotation_colors = annotation_colors,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    display_numbers = asterisks_matrix,
                    fontsize_number = 10,
                    color = my_color_palette,
                    breaks = breaks,  # Explicitly add breaks
                    main = "Heatmap of Associations Grouped by Disease type",
                    annotation_legend = TRUE,
                    annotation_names_row = TRUE,
                    annotation_names_col = TRUE,
                    gaps_row = which(diff(as.numeric(factor(row_annotation$Disease_type))) != 0),
                    gaps_col = which(diff(as.numeric(factor(col_annotation$Disease_type))) != 0),
                    angle_col = 90)


# Extract the heatmap as a gtable object
heatmap_gtable <- heatmap$gtable

# Create text grobs for titles
title1 <- textGrob("Differentially abundant metabolites (n = 114)", gp = gpar(fontsize = 13), hjust = 1, vjust = -4)
title2 <- textGrob("Differentially abundant species (n = 19)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.2, vjust = 7)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot <- arrangeGrob(
  arrangeGrob(title2, heatmap_gtable, ncol = 2, widths = c(0.2, 4)),
  title1, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot)

# Save the combined plot
# ggsave("heatmap_species_metabolite_residual.pdf", combined_plot, width = 45, height = 10, dpi=1200, limitsize = FALSE)





#### Microbial pathways: HAllA readin ####
# read in all_associations.txt got from running HALLA
r_diff_enzyme <- read.table('/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/halla_path_residual/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")


####### Microbial pathways: Data cleaning #######
test_reduce2_enzyme <- r_diff_enzyme

# Assign test_reduce2_path to test_path
test_enzyme <- test_reduce2_enzyme

# Merge test and path dataframes based on Y_features
test2_enzyme <- merge(test_enzyme, path, by = "Y_features", all.x = TRUE)
test2_enzyme$X_features_factor <- as.factor(test2_enzyme$X_features)
test2_enzyme$Y_features_factor <- as.factor(test2_enzyme$Y_features)


####### Microbial pathways: Prepare for clustered heatmap #######
mc_enzymes <- c(
  "BIOTIN-BIOSYNTHESIS-PWY: biotin biosynthesis I",
  "COBALSYN-PWY: superpathway of adenosylcobalamin salvage from cobinamide I",
  "COLANSYN-PWY: colanic acid building blocks biosynthesis",
  "COMPLETE-ARO-PWY: superpathway of aromatic amino acid biosynthesis",
  "FASYN-ELONG-PWY: fatty acid elongation -- saturated",
  "FOLSYN-PWY: superpathway of tetrahydrofolate biosynthesis and salvage",
  "GALACT-GLUCUROCAT-PWY: superpathway of hexuronide and hexuronate degradation",
  "GALACTUROCAT-PWY: D-galacturonate degradation I",
  "GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation",
  "GLUCOSE1PMETAB-PWY: glucose and glucose-1-phosphate degradation",
  "GLUCUROCAT-PWY: superpathway of β-D-glucuronosides degradation",
  "GLUTORN-PWY: L-ornithine biosynthesis I",
  "METH-ACETATE-PWY: methanogenesis from acetate",
  "P41-PWY: pyruvate fermentation to acetate and (S)-lactate I",
  "POLYISOPRENSYN-PWY: polyisoprenoid biosynthesis (E. coli)",
  "PWY-2941: L-lysine biosynthesis II",
  "PWY-5100: pyruvate fermentation to acetate and lactate II",
  "PWY-5384: sucrose degradation IV (sucrose phosphorylase)",
  "PWY-5497: purine nucleobases degradation II (anaerobic)",
  "PWY-5505: L-glutamate and L-glutamine biosynthesis",
  "PWY-5659: GDP-mannose biosynthesis",
  "PWY-5667: CDP-diacylglycerol biosynthesis I",
  "PWY-5676: acetyl-CoA fermentation to butanoate II",
  "PWY-5838: superpathway of menaquinol-8 biosynthesis I",
  "PWY-5861: superpathway of demethylmenaquinol-8 biosynthesis I",
  "PWY-5897: superpathway of menaquinol-11 biosynthesis",
  "PWY-5898: superpathway of menaquinol-12 biosynthesis",
  "PWY-5899: superpathway of menaquinol-13 biosynthesis",
  "PWY-5941: glycogen degradation II",
  "PWY-5989: stearate biosynthesis II (bacteria and plants)",
  "PWY-621: sucrose degradation III (sucrose invertase)",
  "PWY-6270: isoprene biosynthesis I",
  "PWY-6282: palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)",
  "PWY-6292: superpathway of L-cysteine biosynthesis (mammalian)",
  "PWY-6317: D-galactose degradation I (Leloir pathway)",
  "PWY-6353: purine nucleotides degradation II (aerobic)",
  "PWY-6507: 4-deoxy-L-threo-hex-4-enopyranuronate degradation",
  "PWY-6519: 8-amino-7-oxononanoate biosynthesis I",
  "PWY-6527: stachyose degradation",
  "PWY-6606: guanosine nucleotides degradation II",
  "PWY-6628: superpathway of L-phenylalanine biosynthesis",
  "PWY-6629: superpathway of L-tryptophan biosynthesis",
  "PWY-6823: molybdopterin biosynthesis",
  "PWY-6859: all-trans-farnesol biosynthesis",
  "PWY-702: L-methionine biosynthesis II",
  "PWY-7237: myo-, chiro- and scyllo-inositol degradation",
  "PWY-7242: D-fructuronate degradation",
  "PWY-7323: superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis",
  "PWY-7392: taxadiene biosynthesis (engineered)",
  "PWY-7560: methylerythritol phosphate pathway II",
  "PWY-7664: oleate biosynthesis IV (anaerobic)",
  "PWY0-1296: purine ribonucleosides degradation",
  "PWY0-1319: CDP-diacylglycerol biosynthesis II",
  "PWY0-862: (5Z)-dodecenoate biosynthesis I",
  "PWY66-429: fatty acid biosynthesis initiation (mitochondria)",
  "RHAMCAT-PWY: L-rhamnose degradation I",
  "SALVADEHYPOX-PWY: adenosine nucleotides degradation II",
  "TCA: TCA cycle I (prokaryotic)",
  "TRPSYN-PWY: L-tryptophan biosynthesis"
)

test2_enzyme <- test2_enzyme %>%
  mutate(enzyme_abundance = if_else(X_features_factor %in% mc_enzymes, "MC_enriched", "Controls_enriched"))

test2_enzyme <- test2_enzyme %>%
  mutate(metabolite_abundance = if_else(Y_features_factor %in% mc_metabolites, "MC_enriched", "Controls_enriched"))



####### Microbial pathways: clustering by feature type #######
# Subset Relevant Columns
a_enzyme <- test2_enzyme[, colnames(test2_enzyme) %in% c("X_features", "Y_features", "association")]
a1_enzyme <- test2_enzyme[, colnames(test2_enzyme) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_enzyme <- a_enzyme %>% pivot_wider(names_from = "Y_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("X_features") %>% as.matrix()
q_values_matrix_enzyme <- a1_enzyme %>% pivot_wider(names_from = "Y_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("X_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation_enzyme <- data.frame(Disease_type = test2_enzyme$enzyme_abundance[match(rownames(b_enzyme), test2_enzyme$X_features)])
col_annotation_enzyme <- data.frame(Disease_type = test2_enzyme$metabolite_abundance[match(colnames(b_enzyme), test2_enzyme$Y_features)],
                             Pathway = test2_enzyme$pathway[match(colnames(b_enzyme), test2_enzyme$Y_features)])

# Set row and column names for annotations
rownames(row_annotation_enzyme) <- rownames(b_enzyme)
rownames(col_annotation_enzyme) <- colnames(b_enzyme)

# Order rows based on group annotations
ordered_rows_enzyme <- rownames(b_enzyme)[order(row_annotation_enzyme$Disease_type)]

# Order columns based on Disease_type and Pathway
ordered_cols_enzyme <- colnames(b_enzyme)[order(col_annotation_enzyme$Disease_type, col_annotation_enzyme$Pathway)]

# Reorder matrix
b_enzyme <- b_enzyme[ordered_rows_enzyme, ordered_cols_enzyme]

# Reorder annotations to match the new order
row_annotation_enzyme <- row_annotation_enzyme[ordered_rows_enzyme, , drop = FALSE]
col_annotation_enzyme <- col_annotation_enzyme[ordered_cols_enzyme, , drop = FALSE]


# Define colors for the pathways using a more distinct and softer palette
unique_pathways_enzyme <- unique(col_annotation_enzyme$Pathway)
pathway_colors_enzyme <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways_enzyme))
names(pathway_colors_enzyme) <- unique_pathways_enzyme

# Define colors for the annotations
annotation_colors_enzyme <- list(
  Disease_type = c(Controls_enriched = "#189AF9", MC_enriched = "#900C3F"),
  Pathway = pathway_colors_enzyme
)


asterisks_matrix_enzyme <- matrix("", nrow = nrow(q_values_matrix_enzyme), ncol = ncol(q_values_matrix_enzyme))
asterisks_matrix_enzyme[q_values_matrix_enzyme < 0.05] <- "***"
asterisks_matrix_enzyme[q_values_matrix_enzyme >= 0.05 & q_values_matrix_enzyme < 0.1] <- "**"
asterisks_matrix_enzyme[q_values_matrix_enzyme >= 0.1 & q_values_matrix_enzyme < 0.25] <- "*"

# Create the heatmap
heatmap_enzyme <- pheatmap(b_enzyme,
                    annotation_row = row_annotation_enzyme,
                    annotation_col = col_annotation_enzyme,
                    annotation_colors = annotation_colors_enzyme,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    display_numbers = asterisks_matrix_enzyme,
                    fontsize_number = 10,
                    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                    main = "Heatmap of associations between mircobial enzymes and metabolites grouped by disease type",
                    annotation_legend = TRUE,
                    annotation_names_row = TRUE,
                    annotation_names_col = TRUE,
                    gaps_row = which(diff(as.numeric(factor(row_annotation_enzyme$Disease_type))) != 0),
                    gaps_col = which(diff(as.numeric(factor(col_annotation_enzyme$Disease_type))) != 0),
                    angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_enzyme_gtable <- heatmap_enzyme$gtable

# Create text grobs for titles
title1_enzyme <- textGrob("Differentially abundant metabolites (n = 114)", gp = gpar(fontsize = 13), hjust = 2, vjust = 0)
title2_enzyme <- textGrob("Differentially abundant mircobial enzymes (n = 17)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.2, vjust = 7.5)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_enzyme <- arrangeGrob(
  arrangeGrob(title2_enzyme, heatmap_enzyme_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_enzyme, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_enzyme)

# Save the combined plot
# ggsave("heatmap_enzymes_metabolite.pdf", combined_plot_enzyme, width = 45, height = 10, dpi=1200, limitsize = FALSE)





#### All species vs altered metabolites: HAllA readin ####
# read in all_associations.txt got from running HALLA
r_diff_taxo_all <- read.table('/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/halla_residual_full/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")


####### All species vs altered metabolites: Data cleaning #######
test_reduce2_taxo_all <- r_diff_taxo_all
test_taxo_all  <- test_reduce2_taxo_all 

# Merge test and path dataframes based on Y_features
test2_taxo_all  <- merge(test_taxo_all, path, by = "Y_features", all.x = TRUE)
test2_taxo_all$X_features_factor <- as.factor(test2_taxo_all$X_features)
test2_taxo_all$Y_features_factor <- as.factor(test2_taxo_all$Y_features)


####### All species vs altered metabolites: clustering by feature type #######
# Subset Relevant Columns
a_taxo_all <- test2_taxo_all[, colnames(test2_taxo_all) %in% c("X_features", "Y_features", "association")]
a1_taxo_all <- test2_taxo_all[, colnames(test2_taxo_all) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_taxo_all <- a_taxo_all %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_taxo_all <- a1_taxo_all %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation_taxo_all <- data.frame(Pathway = test2_taxo_all$pathway[match(rownames(b_taxo_all), test2_taxo_all$Y_features)])

# Set column names for annotations
rownames(row_annotation_taxo_all) <- rownames(b_taxo_all)

# Order columns based on Pathway
ordered_rows_taxo_all <- rownames(b_taxo_all)[order(row_annotation_taxo_all$Pathway)]

# Reorder matrix
b_taxo_all <- b_taxo_all[ordered_rows_taxo_all, ]

# Reorder annotations to match the new order
row_annotation_taxo_all <- row_annotation_taxo_all[ordered_rows_taxo_all, , drop = FALSE]

# Define colors for the pathways
unique_pathways_taxo_all <- unique(row_annotation_taxo_all$Pathway)
pathway_colors_taxo_all <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways_taxo_all))
names(pathway_colors_taxo_all) <- unique_pathways_taxo_all

# Define colors for the annotations
annotation_colors_taxo_all <- list(Pathway = pathway_colors_taxo_all)


asterisks_matrix_taxo_all <- matrix("", nrow = nrow(q_values_matrix_taxo_all), ncol = ncol(q_values_matrix_taxo_all))
asterisks_matrix_taxo_all[q_values_matrix_taxo_all < 0.05] <- "***"
asterisks_matrix_taxo_all[q_values_matrix_taxo_all >= 0.05 & q_values_matrix_taxo_all < 0.1] <- "**"
asterisks_matrix_taxo_all[q_values_matrix_taxo_all >= 0.1 & q_values_matrix_taxo_all < 0.25] <- "*"

# Create the heatmap
heatmap_taxo_all <- pheatmap(b_taxo_all,
                             annotation_row = row_annotation_taxo_all,
                             annotation_colors = annotation_colors_taxo_all,
                             cluster_rows = FALSE,
                             cluster_cols = TRUE,
                             display_numbers = asterisks_matrix_taxo_all,
                             fontsize_number = 10,
                             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                             main = "Heatmap of associations between all species and altered metabolites",
                             annotation_legend = TRUE,
                             annotation_names_row = TRUE,
                             annotation_names_col = TRUE,
                             angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_taxo_all_gtable <- heatmap_taxo_all$gtable

# Create text grobs for titles
title1_taxo_all <- textGrob("All mircobial species (n = 467)", gp = gpar(fontsize = 13), hjust = 2, vjust = -4)
title2_taxo_all <- textGrob("Differentially abundant metabolites (n = 114)", rot = 90, gp = gpar(fontsize = 13), hjust = 0, vjust = 10)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_taxo_all <- arrangeGrob(
  arrangeGrob(title2_taxo_all, heatmap_taxo_all_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_taxo_all, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_taxo_all)

# Save the combined plot
# ggsave("heatmap_all_species_metabolite.pdf", combined_plot_taxo_all, width = 90, height = 35, dpi=1200, limitsize = FALSE)




#### All species vs Sphingolipids metabolites ####
test2_taxo_sphingolipid <- test2_taxo_all

features_to_select <- c(
  "ceramide (d18:1/14:0, d16:1/16:0)*",
  "ceramide (d18:1/17:0, d17:1/18:0)*",
  "ceramide (d18:2/16:0, d18:1/16:1, d16:1/18:1)*",
  "ceramide (d18:2/24:1, d18:1/24:2)*",
  "N-palmitoyl-sphingosine (d18:1/16:0)",
  "palmitoyl sphingomyelin (d18:1/16:0)",
  "N-palmitoyl-sphinganine (d18:0/16:0)",
  "palmitoyl dihydrosphingomyelin (d18:0/16:0)*",
  "glycosyl-N-palmitoyl-sphingosine (d18:1/16:0)",
  "lactosyl-N-behenoyl-sphingosine (d18:1/22:0)*",
  "lactosyl-N-nervonoyl-sphingosine (d18:1/24:1)*",
  "lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)"
)

test2_taxo_sphingolipid <- test2_taxo_sphingolipid %>%
  filter(Y_features %in% features_to_select)

# Subset Relevant Columns
a_taxo_sphingolipid <- test2_taxo_sphingolipid[, colnames(test2_taxo_sphingolipid) %in% c("X_features", "Y_features", "association")]
a1_taxo_sphingolipid <- test2_taxo_sphingolipid[, colnames(test2_taxo_sphingolipid) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_taxo_sphingolipid <- a_taxo_sphingolipid %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_taxo_sphingolipid <- a1_taxo_sphingolipid %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation_taxo_sphingolipid <- data.frame(Pathway = test2_taxo_sphingolipid$pathway[match(rownames(b_taxo_sphingolipid), test2_taxo_sphingolipid$Y_features)])

# Set column names for annotations
rownames(row_annotation_taxo_sphingolipid) <- rownames(b_taxo_sphingolipid)

# Order columns based on Pathway
ordered_rows_taxo_sphingolipid <- rownames(b_taxo_sphingolipid)[order(row_annotation_taxo_sphingolipid$Pathway)]

# Reorder matrix
b_taxo_sphingolipid <- b_taxo_sphingolipid[ordered_rows_taxo_sphingolipid, ]

# Reorder annotations to match the new order
row_annotation_taxo_sphingolipid <- row_annotation_taxo_sphingolipid[ordered_rows_taxo_sphingolipid, , drop = FALSE]

# Define colors for the pathways
unique_pathways_taxo_sphingolipid <- unique(row_annotation_taxo_sphingolipid$Pathway)
pathway_colors_taxo_sphingolipid <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways_taxo_sphingolipid))
names(pathway_colors_taxo_sphingolipid) <- unique_pathways_taxo_sphingolipid

# Define colors for the annotations
annotation_colors_taxo_sphingolipid <- list(Pathway = pathway_colors_taxo_sphingolipid)


asterisks_matrix_taxo_sphingolipid <- matrix("", nrow = nrow(q_values_matrix_taxo_sphingolipid), ncol = ncol(q_values_matrix_taxo_sphingolipid))
asterisks_matrix_taxo_sphingolipid[q_values_matrix_taxo_sphingolipid < 0.05] <- "***"
asterisks_matrix_taxo_sphingolipid[q_values_matrix_taxo_sphingolipid >= 0.05 & q_values_matrix_taxo_sphingolipid < 0.1] <- "**"
asterisks_matrix_taxo_sphingolipid[q_values_matrix_taxo_sphingolipid >= 0.1 & q_values_matrix_taxo_sphingolipid < 0.25] <- "*"

# Create the heatmap
heatmap_taxo_sphingolipid <- pheatmap(b_taxo_sphingolipid,
                             annotation_row = row_annotation_taxo_sphingolipid,
                             annotation_colors = annotation_colors_taxo_sphingolipid,
                             cluster_rows = FALSE,
                             cluster_cols = TRUE,
                             display_numbers = asterisks_matrix_taxo_sphingolipid,
                             fontsize_number = 10,
                             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                             main = "Heatmap of correlations between all species and sphingolipid metabolites",
                             annotation_legend = TRUE,
                             annotation_names_row = TRUE,
                             annotation_names_col = TRUE,
                             angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_taxo_sphingolipid_gtable <- heatmap_taxo_sphingolipid$gtable

# Create text grobs for titles
title1_taxo_sphingolipid <- textGrob("All mircobial species (n = 467)", gp = gpar(fontsize = 13), hjust = 2, vjust = -4)
title2_taxo_sphingolipid <- textGrob("Sphingolipid metabolites (n = 12)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 13)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_taxo_sphingolipid <- arrangeGrob(
  arrangeGrob(title2_taxo_sphingolipid, heatmap_taxo_sphingolipid_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_taxo_sphingolipid, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_taxo_sphingolipid)

# Save the combined plot
# ggsave("heatmap_all_species_sphingolipids.png", combined_plot_taxo_sphingolipid, width = 90, height = 7, dpi=300, limitsize = FALSE)




#### All species vs Sphingolipids metabolites (rho > 0.15 or < -0.15) ####
test2_taxo_sphingolipid_0.15 <- test2_taxo_sphingolipid

# Filter rows where the association is < -0.15 or > 0.15
filtered_data <- test2_taxo_sphingolipid %>%
  filter(association < -0.15 | association > 0.15)

# Get the unique X_features from the filtered data
unique_features <- filtered_data %>% select(X_features) %>% unique()
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/sphingolipids_0.15_bugs.xlsx"
# write.xlsx(unique_features, file_path, rowNames = TRUE)


# Merge the unique features with the original data to keep relevant rows
test2_taxo_sphingolipid_0.15 <- test2_taxo_sphingolipid %>%
  inner_join(unique_features, by = c("X_features"))

# Subset Relevant Columns
a_taxo_sphingolipid_0.15 <- test2_taxo_sphingolipid_0.15[, colnames(test2_taxo_sphingolipid_0.15) %in% c("X_features", "Y_features", "association")]
a1_taxo_sphingolipid_0.15 <- test2_taxo_sphingolipid_0.15[, colnames(test2_taxo_sphingolipid_0.15) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_taxo_sphingolipid_0.15 <- a_taxo_sphingolipid_0.15 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_taxo_sphingolipid_0.15 <- a1_taxo_sphingolipid_0.15 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation_taxo_sphingolipid_0.15 <- data.frame(Metabolite_class = test2_taxo_sphingolipid_0.15$pathway[match(rownames(b_taxo_sphingolipid_0.15), test2_taxo_sphingolipid_0.15$Y_features)])

# Set column names for annotations
rownames(row_annotation_taxo_sphingolipid_0.15) <- rownames(b_taxo_sphingolipid_0.15)

# Order columns based on Metabolite_class
ordered_rows_taxo_sphingolipid_0.15 <- rownames(b_taxo_sphingolipid_0.15)[order(row_annotation_taxo_sphingolipid_0.15$Metabolite_class)]

# Reorder matrix
b_taxo_sphingolipid_0.15 <- b_taxo_sphingolipid_0.15[ordered_rows_taxo_sphingolipid_0.15, ]

# Reorder annotations to match the new order
row_annotation_taxo_sphingolipid_0.15 <- row_annotation_taxo_sphingolipid_0.15[ordered_rows_taxo_sphingolipid_0.15, , drop = FALSE]

# Define colors for the pathways
unique_pathways_taxo_sphingolipid_0.15 <- unique(row_annotation_taxo_sphingolipid_0.15$Metabolite_class)
pathway_colors_taxo_sphingolipid_0.15 <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways_taxo_sphingolipid_0.15))
names(pathway_colors_taxo_sphingolipid_0.15) <- unique_pathways_taxo_sphingolipid_0.15

# Define colors for the annotations
annotation_colors_taxo_sphingolipid_0.15 <- list(Metabolite_class = pathway_colors_taxo_sphingolipid_0.15)


asterisks_matrix_taxo_sphingolipid_0.15 <- matrix("", nrow = nrow(q_values_matrix_taxo_sphingolipid_0.15), ncol = ncol(q_values_matrix_taxo_sphingolipid_0.15))
asterisks_matrix_taxo_sphingolipid_0.15[q_values_matrix_taxo_sphingolipid_0.15 < 0.05] <- "***"
asterisks_matrix_taxo_sphingolipid_0.15[q_values_matrix_taxo_sphingolipid_0.15 >= 0.05 & q_values_matrix_taxo_sphingolipid_0.15 < 0.1] <- "**"
asterisks_matrix_taxo_sphingolipid_0.15[q_values_matrix_taxo_sphingolipid_0.15 >= 0.1 & q_values_matrix_taxo_sphingolipid_0.15 < 0.25] <- "*"

# Create the heatmap
heatmap_taxo_sphingolipid_0.15 <- pheatmap(b_taxo_sphingolipid_0.15,
                                          annotation_row = row_annotation_taxo_sphingolipid_0.15,
                                          annotation_colors = annotation_colors_taxo_sphingolipid_0.15,
                                          cluster_rows = FALSE,
                                          cluster_cols = TRUE,
                                          display_numbers = asterisks_matrix_taxo_sphingolipid_0.15,
                                          fontsize_number = 10,
                                          color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                          annotation_legend = TRUE,
                                          annotation_names_row = TRUE,
                                          annotation_names_col = TRUE,
                                          angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_taxo_sphingolipid_0.15_gtable <- heatmap_taxo_sphingolipid_0.15$gtable

# Create text grobs for titles
title1_taxo_sphingolipid_0.15 <- textGrob("Selected mircobial species (n = 320)", gp = gpar(fontsize = 13), hjust = 1, vjust = -0.5)
title2_taxo_sphingolipid_0.15 <- textGrob("Sphingolipid metabolites (n = 12)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 10)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_taxo_sphingolipid_0.15 <- arrangeGrob(
  arrangeGrob(title2_taxo_sphingolipid_0.15, heatmap_taxo_sphingolipid_0.15_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_taxo_sphingolipid_0.15, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_taxo_sphingolipid_0.15)

# Save the combined plot
# ggsave("heatmap_selected_species_sphingolipids_residuals.pdf", combined_plot_taxo_sphingolipid_0.15, width = 65, height = 7, dpi=600, limitsize = FALSE)


### Simplified version
# List of species to present
species_to_present <- c(
  "Haemophilus_parainfluenzae",
  "Methylobacterium_SGB15164",
  "f_Clostridia_unclassified_GGB33586_SGB53517",
  "Clostridiales_bacterium_NSJ_32",
  "Mediterraneibacter_butyricigenes",
  "GGB9534_SGB14937",
  "Blautia_glucerasea",
  "Ruminococcus_gnavus",
  "Akkermansia_muciniphila"
)

# Filter the matrix and annotations to only include the selected species
b_taxo_sphingolipid_0.15_filtered <- b_taxo_sphingolipid_0.15[, colnames(b_taxo_sphingolipid_0.15) %in% species_to_present]
q_values_matrix_taxo_sphingolipid_0.15_filtered <- q_values_matrix_taxo_sphingolipid_0.15[, colnames(b_taxo_sphingolipid_0.15) %in% species_to_present]

# Reorder matrix
b_taxo_sphingolipid_0.15_filtered <- b_taxo_sphingolipid_0.15_filtered[ordered_rows_taxo_sphingolipid_0.15, ]

# Create the asterisks matrix for the filtered data
asterisks_matrix_taxo_sphingolipid_0.15_filtered <- matrix("", nrow = nrow(q_values_matrix_taxo_sphingolipid_0.15_filtered), ncol = ncol(q_values_matrix_taxo_sphingolipid_0.15_filtered))
asterisks_matrix_taxo_sphingolipid_0.15_filtered[q_values_matrix_taxo_sphingolipid_0.15_filtered < 0.05] <- "***"
asterisks_matrix_taxo_sphingolipid_0.15_filtered[q_values_matrix_taxo_sphingolipid_0.15_filtered >= 0.05 & q_values_matrix_taxo_sphingolipid_0.15_filtered < 0.1] <- "**"
asterisks_matrix_taxo_sphingolipid_0.15_filtered[q_values_matrix_taxo_sphingolipid_0.15_filtered >= 0.1 & q_values_matrix_taxo_sphingolipid_0.15_filtered < 0.25] <- "*"

# Create symmetric breaks around 0
breaks <- seq(-0.28, 0.28, length.out = 51)

# Check if 0 is in the middle
print(breaks[which.min(abs(breaks))]) 

# Find the index where 0 is located
mid <- which.min(abs(breaks - 0))

# Create the custom palette
my_color_palette <- c(
  colorRampPalette(c("navy", "white"))(mid - 1), 
  "white", 
  colorRampPalette(c("white", "firebrick3"))(length(breaks) - mid)
)

# Create the filtered heatmap
heatmap_taxo_sphingolipid_0.15_filtered <- pheatmap(b_taxo_sphingolipid_0.15_filtered,
                                                   annotation_row = row_annotation_taxo_sphingolipid_0.15,
                                                   annotation_colors = annotation_colors_taxo_sphingolipid_0.15,
                                                   cluster_rows = FALSE,
                                                   cluster_cols = TRUE,
                                                   display_numbers = asterisks_matrix_taxo_sphingolipid_0.15_filtered,
                                                   fontsize_number = 10,
                                                   color = my_color_palette,
                                                   breaks = breaks,
                                                   main = "Heatmap of correlations between selected species and sphingolipid metabolites",
                                                   annotation_legend = TRUE,
                                                   annotation_names_row = TRUE,
                                                   annotation_names_col = TRUE,
                                                   border_color = NA,
                                                   angle_col = 90)

# Extract the filtered heatmap as a gtable object
heatmap_taxo_sphingolipid_0.15_filtered_gtable <- heatmap_taxo_sphingolipid_0.15_filtered$gtable

# Create text grobs for titles
title1_taxo_sphingolipid_0.15_filtered <- textGrob("Selected microbial species (n = 8)", gp = gpar(fontsize = 13), hjust = 2, vjust = -0.5)
title2_taxo_sphingolipid_0.15_filtered <- textGrob("Sphingolipid metabolites (n = 12)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 8)

# Combine the filtered heatmap and titles using arrangeGrob and grid.arrange
combined_plot_taxo_sphingolipid_0.15_filtered <- arrangeGrob(
  arrangeGrob(title2_taxo_sphingolipid_0.15_filtered, heatmap_taxo_sphingolipid_0.15_filtered_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_taxo_sphingolipid_0.15_filtered, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined filtered plot
grid.newpage()
grid.draw(combined_plot_taxo_sphingolipid_0.15_filtered)

# Save the filtered combined plot
# ggsave("heatmap_selected_species_sphingolipids_filtered_residual.pdf", combined_plot_taxo_sphingolipid_0.15_filtered, width = 10, height = 7, dpi=1200, limitsize = FALSE)



#### Species-sphingolipid correlations concordant with abundance analysis ####
abundance_sphingolipids <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross/all_results.tsv")
abundance_sphingolipids <- abundance_sphingolipids %>% filter(metadata=="mc_all") %>% 
  select(1,3,4) %>%
  filter(feature %in% c("Haemophilus_parainfluenzae", "Methylobacterium_SGB15164", "GGB33586_SGB53517", 
                        "Mediterraneibacter_butyricigenes", "GGB9534_SGB14937", "Blautia_glucerasea"))

abundance_sphingolipids$coef <- -(abundance_sphingolipids$coef)

abundance_sphingolipids_avg <- abundance_sphingolipids %>%
  group_by(feature) %>%
  summarise(avg_coef = mean(coef))

function_sphingolipid <- b_taxo_sphingolipid[, colnames(b_taxo_sphingolipid) %in% c("Haemophilus_parainfluenzae", "Methylobacterium_SGB15164", "GGB33586_SGB53517", 
                                                                                    "Mediterraneibacter_butyricigenes", "GGB9534_SGB14937", "Blautia_glucerasea")]%>% t() %>% as.data.frame()

function_sphingolipid$feature <- rownames(function_sphingolipid)

concordance_sphingolipid <- abundance_sphingolipids_avg %>% left_join(function_sphingolipid)

concordance_sphingolipid <- concordance_sphingolipid %>% 
  mutate(enriched = case_when(
    avg_coef > 0 ~ "1",
    avg_coef < 0 ~ "0"
  ))
  
concordance_sphingolipid$species <- concordance_sphingolipid$feature
concordance_sphingolipid$feature <- NULL

fig_concordance_sphingolipid <- ggplot(concordance_sphingolipid, aes(x = avg_coef, y = `lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)`, color = as.factor(enriched))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Vertical line at x = 0
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point() +
  geom_text_repel(aes(label = species), size = 4) + 
  labs(x = "Species abundance in MC compared to controls (beta coeffecient)", 
       y = "Spearman rho of the species with lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)", 
       color = NULL) +
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("navy", "firebrick3"),
                     labels = c("Less abundant in MC & negatively correlated with sphingolipids", 
                                "More abundant in MC & positively correlated with sphingolipids"))

fig_concordance_sphingolipid

# ggsave("fig_concordance_sphingolipid_residual.pdf", plot = fig_concordance_sphingolipid, width = 12, height = 9, dpi = 1200)





#### All pathways vs altered metabolites: HAllA readin ####
# read in all_associations.txt got from running HALLA
r_diff_path_all <- read.table('/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/halla_path_residual/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")


####### All pathways vs altered metabolites: Data cleaning #######
test_reduce2_path_all <- r_diff_path_all
test_path_all  <- test_reduce2_path_all 

# Merge test and path dataframes based on Y_features
test2_path_all  <- merge(test_path_all, path, by = "Y_features", all.x = TRUE)
test2_path_all$X_features_factor <- as.factor(test2_path_all$X_features)
test2_path_all$Y_features_factor <- as.factor(test2_path_all$Y_features)


####### All pathways vs altered metabolites: clustering by feature type #######
# Subset Relevant Columns
a_path_all <- test2_path_all[, colnames(test2_path_all) %in% c("X_features", "Y_features", "association")]
a1_path_all <- test2_path_all[, colnames(test2_path_all) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_path_all <- a_path_all %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_path_all <- a1_path_all %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation_path_all <- data.frame(Pathway = test2_path_all$pathway[match(rownames(b_path_all), test2_path_all$Y_features)])

# Set column names for annotations
rownames(row_annotation_path_all) <- rownames(b_path_all)

# Order columns based on Pathway
ordered_rows_path_all <- rownames(b_path_all)[order(row_annotation_path_all$Pathway)]

# Reorder matrix
b_path_all <- b_path_all[ordered_rows_path_all, ]

# Reorder annotations to match the new order
row_annotation_path_all <- row_annotation_path_all[ordered_rows_path_all, , drop = FALSE]

# Define colors for the pathways
unique_pathways_path_all <- unique(row_annotation_path_all$Pathway)
pathway_colors_path_all <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways_path_all))
names(pathway_colors_path_all) <- unique_pathways_path_all

# Define colors for the annotations
annotation_colors_path_all <- list(Pathway = pathway_colors_path_all)


asterisks_matrix_path_all <- matrix("", nrow = nrow(q_values_matrix_path_all), ncol = ncol(q_values_matrix_path_all))
asterisks_matrix_path_all[q_values_matrix_path_all < 0.05] <- "***"
asterisks_matrix_path_all[q_values_matrix_path_all >= 0.05 & q_values_matrix_path_all < 0.1] <- "**"
asterisks_matrix_path_all[q_values_matrix_path_all >= 0.1 & q_values_matrix_path_all < 0.25] <- "*"

# Create the heatmap
heatmap_path_all <- pheatmap(b_path_all,
                             annotation_row = row_annotation_path_all,
                             annotation_colors = annotation_colors_path_all,
                             cluster_rows = FALSE,
                             cluster_cols = TRUE,
                             display_numbers = asterisks_matrix_path_all,
                             fontsize_number = 10,
                             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                             main = "Heatmap of correlations between all microbial enzymes and altered metabolites",
                             annotation_legend = TRUE,
                             annotation_names_row = TRUE,
                             annotation_names_col = TRUE,
                             angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_path_all_gtable <- heatmap_path_all$gtable

# Create text grobs for titles
title1_path_all <- textGrob("All mircobial pathways (n = 250)", gp = gpar(fontsize = 13), hjust = 2, vjust = -7)
title2_path_all <- textGrob("Differentially abundant metabolites (n = 114)", rot = 90, gp = gpar(fontsize = 13), hjust = 0, vjust = 10)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_path_all <- arrangeGrob(
  arrangeGrob(title2_path_all, heatmap_path_all_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_path_all, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_path_all)

# Save the combined plot
# ggsave("heatmap_all_enzymes_metabolite.png", combined_plot_path_all, width = 70, height = 35, dpi=300, limitsize = FALSE)




#### All pathways vs Sphingolipid metabolites ####
test2_sphingolipid_all <- test2_path_all

features_to_select <- c(
  "ceramide (d18:1/14:0, d16:1/16:0)*",
  "ceramide (d18:1/17:0, d17:1/18:0)*",
  "ceramide (d18:2/16:0, d18:1/16:1, d16:1/18:1)*",
  "ceramide (d18:2/24:1, d18:1/24:2)*",
  "N-palmitoyl-sphingosine (d18:1/16:0)",
  "palmitoyl sphingomyelin (d18:1/16:0)",
  "N-palmitoyl-sphinganine (d18:0/16:0)",
  "palmitoyl dihydrosphingomyelin (d18:0/16:0)*",
  "glycosyl-N-palmitoyl-sphingosine (d18:1/16:0)",
  "lactosyl-N-behenoyl-sphingosine (d18:1/22:0)*",
  "lactosyl-N-nervonoyl-sphingosine (d18:1/24:1)*",
  "lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)"
)

test2_sphingolipid_all <- test2_sphingolipid_all %>%
  filter(Y_features %in% features_to_select)


# Subset Relevant Columns
a_sphingolipid_all <- test2_sphingolipid_all[, colnames(test2_sphingolipid_all) %in% c("X_features", "Y_features", "association")]
a1_sphingolipid_all <- test2_sphingolipid_all[, colnames(test2_sphingolipid_all) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_sphingolipid_all <- a_sphingolipid_all %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_sphingolipid_all <- a1_sphingolipid_all %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation_sphingolipid_all <- data.frame(Pathway = test2_sphingolipid_all$pathway[match(rownames(b_sphingolipid_all), test2_sphingolipid_all$Y_features)])

# Set column names for annotations
rownames(row_annotation_sphingolipid_all) <- rownames(b_sphingolipid_all)

# Order columns based on Pathway
ordered_rows_sphingolipid_all <- rownames(b_sphingolipid_all)[order(row_annotation_sphingolipid_all$Pathway)]

# Reorder matrix
b_sphingolipid_all <- b_sphingolipid_all[ordered_rows_sphingolipid_all, ]

# Reorder annotations to match the new order
row_annotation_sphingolipid_all <- row_annotation_sphingolipid_all[ordered_rows_sphingolipid_all, , drop = FALSE]

# Define colors for the pathways
unique_pathways_sphingolipid_all <- unique(row_annotation_sphingolipid_all$Pathway)
pathway_colors_sphingolipid_all <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways_sphingolipid_all))
names(pathway_colors_sphingolipid_all) <- unique_pathways_sphingolipid_all

# Define colors for the annotations
annotation_colors_sphingolipid_all <- list(Pathway = pathway_colors_sphingolipid_all)


asterisks_matrix_sphingolipid_all <- matrix("", nrow = nrow(q_values_matrix_sphingolipid_all), ncol = ncol(q_values_matrix_sphingolipid_all))
asterisks_matrix_sphingolipid_all[q_values_matrix_sphingolipid_all < 0.05] <- "***"
asterisks_matrix_sphingolipid_all[q_values_matrix_sphingolipid_all >= 0.05 & q_values_matrix_sphingolipid_all < 0.1] <- "**"
asterisks_matrix_sphingolipid_all[q_values_matrix_sphingolipid_all >= 0.1 & q_values_matrix_sphingolipid_all < 0.25] <- "*"

# Create the heatmap
heatmap_sphingolipid_all <- pheatmap(b_sphingolipid_all,
                             annotation_row = row_annotation_sphingolipid_all,
                             annotation_colors = annotation_colors_sphingolipid_all,
                             cluster_rows = FALSE,
                             cluster_cols = TRUE,
                             display_numbers = asterisks_matrix_sphingolipid_all,
                             fontsize_number = 10,
                             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                             main = "Heatmap of correlations between all microbial pathways and sphingolipid metabolites ",
                             annotation_legend = TRUE,
                             annotation_names_row = TRUE,
                             annotation_names_col = TRUE,
                             angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_sphingolipid_all_gtable <- heatmap_sphingolipid_all$gtable

# Create text grobs for titles
title1_sphingolipid_all <- textGrob("All mircobial pathways (n = 250)", gp = gpar(fontsize = 13), hjust = 2, vjust = -7)
title2_sphingolipid_all <- textGrob("Sphingolipid metabolites (n = 12)", rot = 90, gp = gpar(fontsize = 13), hjust = -0.8, vjust = 10)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_sphingolipid_all <- arrangeGrob(
  arrangeGrob(title2_sphingolipid_all, heatmap_sphingolipid_all_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_sphingolipid_all, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_sphingolipid_all)

# Save the combined plot
# ggsave("heatmap_all_paths_sphingolipids_residual.pdf", combined_plot_sphingolipid_all, width = 70, height = 15, dpi=1200, limitsize = FALSE)



#### All pathways vs Sphingolipid metabolites (rho > 0.15 or < -0.15) ####
test2_sphingolipid_0.15 <- test2_sphingolipid_all

# Filter rows where the association is < -0.15 or > 0.15
filtered_data <- test2_sphingolipid_all %>%
  filter(association < -0.15 | association > 0.15)

# Get the unique X_features from the filtered data
unique_features <- filtered_data %>% select(X_features) %>% unique()
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/sphingolipids_0.15_paths.xlsx"
# write.xlsx(unique_features, file_path, rowNames = TRUE)


# Merge the unique features with the original data to keep relevant rows
test2_sphingolipid_0.15 <- test2_sphingolipid_all %>%
  inner_join(unique_features, by = c("X_features"))


# Subset Relevant Columns
a_sphingolipid_0.15 <- test2_sphingolipid_0.15[, colnames(test2_sphingolipid_0.15) %in% c("X_features", "Y_features", "association")]
a1_sphingolipid_0.15 <- test2_sphingolipid_0.15[, colnames(test2_sphingolipid_0.15) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_sphingolipid_0.15 <- a_sphingolipid_0.15 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_sphingolipid_0.15 <- a1_sphingolipid_0.15 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()

# Create annotations for the heatmap
row_annotation_sphingolipid_0.15 <- data.frame(Metabolite_class = test2_sphingolipid_0.15$pathway[match(rownames(b_sphingolipid_0.15), test2_sphingolipid_0.15$Y_features)])

# Set column names for annotations
rownames(row_annotation_sphingolipid_0.15) <- rownames(b_sphingolipid_0.15)

# Order columns based on Metabolite_class
ordered_rows_sphingolipid_0.15 <- rownames(b_sphingolipid_0.15)[order(row_annotation_sphingolipid_0.15$Metabolite_class)]

# Reorder matrix
b_sphingolipid_0.15 <- b_sphingolipid_0.15[ordered_rows_sphingolipid_0.15, ]

# Reorder annotations to match the new order
row_annotation_sphingolipid_0.15 <- row_annotation_sphingolipid_0.15[ordered_rows_sphingolipid_0.15, , drop = FALSE]

# Define colors for the Metabolite_classs
unique_pathways_sphingolipid_0.15 <- unique(row_annotation_sphingolipid_0.15$Metabolite_class)
pathway_colors_sphingolipid_0.15 <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique_pathways_sphingolipid_0.15))
names(pathway_colors_sphingolipid_0.15) <- unique_pathways_sphingolipid_0.15

# Define colors for the annotations
annotation_colors_sphingolipid_0.15 <- list(Metabolite_class = pathway_colors_sphingolipid_0.15)


asterisks_matrix_sphingolipid_0.15 <- matrix("", nrow = nrow(q_values_matrix_sphingolipid_0.15), ncol = ncol(q_values_matrix_sphingolipid_0.15))
asterisks_matrix_sphingolipid_0.15[q_values_matrix_sphingolipid_0.15 < 0.05] <- "***"
asterisks_matrix_sphingolipid_0.15[q_values_matrix_sphingolipid_0.15 >= 0.05 & q_values_matrix_sphingolipid_0.15 < 0.1] <- "**"
asterisks_matrix_sphingolipid_0.15[q_values_matrix_sphingolipid_0.15 >= 0.1 & q_values_matrix_sphingolipid_0.15 < 0.25] <- "*"

# Create the heatmap
heatmap_sphingolipid_0.15 <- pheatmap(b_sphingolipid_0.15,
                                     annotation_row = row_annotation_sphingolipid_0.15,
                                     annotation_colors = annotation_colors_sphingolipid_0.15,
                                     cluster_rows = FALSE,
                                     cluster_cols = TRUE,
                                     display_numbers = asterisks_matrix_sphingolipid_0.15,
                                     fontsize_number = 10,
                                     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                     annotation_legend = TRUE,
                                     annotation_names_row = TRUE,
                                     annotation_names_col = TRUE,
                                     angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_sphingolipid_0.15_gtable <- heatmap_sphingolipid_0.15$gtable

# Create text grobs for titles
title1_sphingolipid_0.15 <- textGrob("Selected mircobial pathways (n = 145)", gp = gpar(fontsize = 13), hjust = 1, vjust = -7)
title2_sphingolipid_0.15 <- textGrob("Sphingolipid metabolites (n = 12)", rot = 90, gp = gpar(fontsize = 13), hjust = -0.8, vjust = 7)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_sphingolipid_0.15 <- arrangeGrob(
  arrangeGrob(title2_sphingolipid_0.15, heatmap_sphingolipid_0.15_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_sphingolipid_0.15, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_sphingolipid_0.15)

# Save the combined plot
# ggsave("heatmap_sphingolipids_0.15_residual.pdf", combined_plot_sphingolipid_0.15, width = 40, height = 12, dpi=1200, limitsize = FALSE)



## Simplified version
# List of pathways to present
df <- as.data.frame(colnames(b_sphingolipid_0.15))
pathways_to_present <- c(
  "GLUCUROCAT-PWY: superpathway of &beta;-D-glucuronosides degradation",
  "PWY-7242: D-fructuronate degradation",
  "PWY-621: sucrose degradation III (sucrose invertase)",
  "PWY-5384: sucrose degradation IV (sucrose phosphorylase)",
  "PWY-7238: sucrose biosynthesis II",
  "PWY-6901: superpathway of glucose and xylose degradation",
  "GLYCOGENSYNTH-PWY: glycogen biosynthesis I (from ADP-D-Glucose)",
  "PENTOSE-P-PWY: pentose phosphate pathway",
  "PWY-841: superpathway of purine nucleotides de novo biosynthesis I",
  "PWY-7229: superpathway of adenosine nucleotides de novo biosynthesis I",
  "PWY-6126: superpathway of adenosine nucleotides de novo biosynthesis II",
  "PWY-7228: superpathway of guanosine nucleotides de novo biosynthesis I",
  "PWY-6125: superpathway of guanosine nucleotides de novo biosynthesis II",
  "PWY0-1296: purine ribonucleosides degradation",
  "PWY-7208: superpathway of pyrimidine nucleobases salvage",
  "PWY0-1586: peptidoglycan maturation (meso-diaminopimelate containing)",
  "PHOSLIPSYN-PWY: superpathway of phospholipid biosynthesis I (bacteria)"
)

# Filter the matrix and annotations to only include the selected species
b_sphingolipid_0.15_filtered <- b_sphingolipid_0.15[, colnames(b_sphingolipid_0.15) %in% pathways_to_present]
q_values_matrix_sphingolipid_0.15_filtered <- q_values_matrix_sphingolipid_0.15[, colnames(q_values_matrix_sphingolipid_0.15) %in% pathways_to_present]

# Reorder matrix
b_sphingolipid_0.15_filtered <- b_sphingolipid_0.15_filtered[ordered_rows_sphingolipid_0.15, ]

# Create the asterisks matrix for the filtered data
asterisks_matrix_sphingolipid_0.15_filtered <- matrix("", nrow = nrow(q_values_matrix_sphingolipid_0.15_filtered), ncol = ncol(q_values_matrix_sphingolipid_0.15_filtered))
asterisks_matrix_sphingolipid_0.15_filtered[q_values_matrix_sphingolipid_0.15_filtered < 0.05] <- "***"
asterisks_matrix_sphingolipid_0.15_filtered[q_values_matrix_sphingolipid_0.15_filtered >= 0.05 & q_values_matrix_sphingolipid_0.15_filtered < 0.1] <- "**"
asterisks_matrix_sphingolipid_0.15_filtered[q_values_matrix_sphingolipid_0.15_filtered >= 0.1 & q_values_matrix_sphingolipid_0.15_filtered < 0.25] <- "*"

# Create the filtered heatmap
heatmap_sphingolipid_0.15_filtered <- pheatmap(b_sphingolipid_0.15_filtered,
                                              annotation_row = row_annotation_sphingolipid_0.15,
                                              annotation_colors = annotation_colors_sphingolipid_0.15,
                                              cluster_rows = FALSE,
                                              cluster_cols = TRUE,
                                              display_numbers = asterisks_matrix_sphingolipid_0.15_filtered,
                                              fontsize_number = 10,
                                              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                              main = "Heatmap of correlations between selected pathways and sphingolipid metabolites",
                                              annotation_legend = TRUE,
                                              annotation_names_row = TRUE,
                                              annotation_names_col = TRUE,
                                              angle_col = 90)

# Extract the filtered heatmap as a gtable object
heatmap_sphingolipid_0.15_filtered_gtable <- heatmap_sphingolipid_0.15_filtered$gtable

# Create text grobs for titles
title1_sphingolipid_0.15_filtered <- textGrob("Selected microbial pathways (n = 17)", gp = gpar(fontsize = 13), hjust = 2, vjust = -0.5)
title2_sphingolipid_0.15_filtered <- textGrob("Sphingolipid metabolites (n = 12)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 0)

# Combine the filtered heatmap and titles using arrangeGrob and grid.arrange
combined_plot_sphingolipid_0.15_filtered <- arrangeGrob(
  arrangeGrob(title2_sphingolipid_0.15_filtered, heatmap_sphingolipid_0.15_filtered_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_sphingolipid_0.15_filtered, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined filtered plot
grid.newpage()
grid.draw(combined_plot_sphingolipid_0.15_filtered)

# Save the filtered combined plot
# ggsave("heatmap_selected_species_sphingolipids_filtered_residual.pdf", combined_plot_sphingolipid_0.15_filtered, width = 12, height = 10, dpi=1200, limitsize = FALSE)







#### All pathways vs Lysophospholipid metabolites ####
# Prepare input for halla
Lysophospholipid_halla <- read_excel("met_halla_residual.xlsx") %>% t() %>%
  as.data.frame() %>% 
  {setNames(.[-1, ], unlist(.[1, ]))}
Lysophospholipid_halla$id <- rownames(Lysophospholipid_halla)

# To select Lysophospholipid
altered_114 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/altered_met_114.xlsx") %>% select(4,5) %>% rename("pathway" = "pathway_HC")
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
Lysophospholipid <- subset(altered_114, grepl("Lysophospholipid", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
Lysophospholipid_metabolites <- Lysophospholipid$metabolite
Lysophospholipid_halla <- Lysophospholipid_halla %>% select(all_of(Lysophospholipid_metabolites))
Lysophospholipid_halla$id <- rownames(Lysophospholipid_halla)

# Change the ID to match microbiome ID
Lysophospholipid_halla <- Lysophospholipid_halla %>%
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
rownames(Lysophospholipid_halla) <- Lysophospholipid_halla$id
Lysophospholipid_halla$id <- NULL

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_lysophospholipid_halla_residual.xlsx"
# write.xlsx(Lysophospholipid_halla, file_path, rowNames = TRUE)

# read in all_associations.txt got from running HALLA
r_diff_lysophospholipid <- read.table('halla_lysophospholipid_path_residual/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")

# Data cleaning
test_reduce2_lysophospholipid <- r_diff_lysophospholipid
test2_lysophospholipid  <- test_reduce2_lysophospholipid

# Subset Relevant Columns
a_lysophospholipid <- test2_lysophospholipid[, colnames(test2_lysophospholipid) %in% c("X_features", "Y_features", "association")]
a1_lysophospholipid <- test2_lysophospholipid[, colnames(test2_lysophospholipid) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysophospholipid <- a_lysophospholipid %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysophospholipid <- a1_lysophospholipid %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysophospholipid <- matrix("", nrow = nrow(q_values_matrix_lysophospholipid), ncol = ncol(q_values_matrix_lysophospholipid))
asterisks_matrix_lysophospholipid[q_values_matrix_lysophospholipid < 0.05] <- "***"
asterisks_matrix_lysophospholipid[q_values_matrix_lysophospholipid >= 0.05 & q_values_matrix_lysophospholipid < 0.1] <- "**"
asterisks_matrix_lysophospholipid[q_values_matrix_lysophospholipid >= 0.1 & q_values_matrix_lysophospholipid < 0.25] <- "*"

# Create the heatmap
heatmap_lysophospholipid <- pheatmap(b_lysophospholipid,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     display_numbers = asterisks_matrix_lysophospholipid,
                                     fontsize_number = 10,
                                     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                     main = "Heatmap of correlations between all microbial pathways and lysophospholipid metabolites ",
                                     annotation_legend = TRUE,
                                     annotation_names_row = TRUE,
                                     annotation_names_col = TRUE,
                                     angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysophospholipid_gtable <- heatmap_lysophospholipid$gtable

# Create text grobs for titles
title1_lysophospholipid <- textGrob("All mircobial pathways (n = 250)", gp = gpar(fontsize = 13), hjust = 2, vjust = -7)
title2_lysophospholipid <- textGrob("Lysophospholipid metabolites (n = 14)", rot = 90, gp = gpar(fontsize = 13), hjust = -0.8, vjust = 10)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysophospholipid <- arrangeGrob(
  arrangeGrob(title2_lysophospholipid, heatmap_lysophospholipid_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysophospholipid, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysophospholipid)

# Save the combined plot
# ggsave("heatmap_all_paths_lysophospholipids.png", combined_plot_lysophospholipid, width = 70, height = 13, dpi=300, limitsize = FALSE)




#### All pathways vs Lysophospholipid metabolites (rho > 0.15 or < -0.15) ####
test2_lysophospholipid_0.15 <- test2_lysophospholipid

# Filter rows where the association is < -0.15 or > 0.15
filtered_data <- test2_lysophospholipid %>%
  filter(association < -0.15 | association > 0.15)

# Get the unique X_features from the filtered data
unique_features <- filtered_data %>% select(X_features) %>% unique()
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysophospholipids_0.15_paths.xlsx"
# write.xlsx(unique_features, file_path, rowNames = TRUE)


# Merge the unique features with the original data to keep relevant rows
test2_lysophospholipid_0.15 <- test2_lysophospholipid %>%
  inner_join(unique_features, by = c("X_features"))


# Subset Relevant Columns
a_lysophospholipid_0.15 <- test2_lysophospholipid_0.15[, colnames(test2_lysophospholipid_0.15) %in% c("X_features", "Y_features", "association")]
a1_lysophospholipid_0.15 <- test2_lysophospholipid_0.15[, colnames(test2_lysophospholipid_0.15) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysophospholipid_0.15 <- a_lysophospholipid_0.15 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysophospholipid_0.15 <- a1_lysophospholipid_0.15 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysophospholipid_0.15 <- matrix("", nrow = nrow(q_values_matrix_lysophospholipid_0.15), ncol = ncol(q_values_matrix_lysophospholipid_0.15))
asterisks_matrix_lysophospholipid_0.15[q_values_matrix_lysophospholipid_0.15 < 0.05] <- "***"
asterisks_matrix_lysophospholipid_0.15[q_values_matrix_lysophospholipid_0.15 >= 0.05 & q_values_matrix_lysophospholipid_0.15 < 0.1] <- "**"
asterisks_matrix_lysophospholipid_0.15[q_values_matrix_lysophospholipid_0.15 >= 0.1 & q_values_matrix_lysophospholipid_0.15 < 0.25] <- "*"

# Create the heatmap
heatmap_lysophospholipid_0.15 <- pheatmap(b_lysophospholipid_0.15,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         display_numbers = asterisks_matrix_lysophospholipid_0.15,
                                         fontsize_number = 10,
                                         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                         annotation_legend = TRUE,
                                         annotation_names_row = TRUE,
                                         annotation_names_col = TRUE,
                                         angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysophospholipid_0.15_gtable <- heatmap_lysophospholipid_0.15$gtable

# Create text grobs for titles
title1_lysophospholipid_0.15 <- textGrob("Selected mircobial pathways (n = 113)", gp = gpar(fontsize = 13), hjust = 0.5, vjust = -7)
title2_lysophospholipid_0.15 <- textGrob("Lysophospholipid metabolites (n = 14)", rot = 90, gp = gpar(fontsize = 13), hjust = -0.65, vjust = 3)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysophospholipid_0.15 <- arrangeGrob(
  arrangeGrob(title2_lysophospholipid_0.15, heatmap_lysophospholipid_0.15_gtable, ncol = 2, widths = c(0.15, 4)),
  title1_lysophospholipid_0.15, ncol = 1, heights = c(4, 0.15)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysophospholipid_0.15)

# Save the combined plot
# ggsave("heatmap_lysophospholipids_0.15_residual.pdf", combined_plot_lysophospholipid_0.15, width = 30, height = 13, dpi=1200, limitsize = FALSE)


# ## Simplified version
# # List of pathways to present
# df <- as.data.frame(colnames(b_lysophospholipid_0.15))
# pathways_to_present <- c(
#   "PHOSLIPSYN-PWY: superpathway of phospholipid biosynthesis I (bacteria)"
# )
# 
# # Filter the matrix and annotations to only include the selected pathways
# b_lysophospholipid_0.15_filtered <- b_lysophospholipid_0.15[, colnames(b_lysophospholipid_0.15) %in% pathways_to_present] %>% as.matrix()
# q_values_matrix_lysophospholipid_0.15_filtered <- q_values_matrix_lysophospholipid_0.15[, colnames(q_values_matrix_lysophospholipid_0.15) %in% pathways_to_present] %>% as.matrix()
# colnames(b_lysophospholipid_0.15_filtered) <- colnames(b_lysophospholipid_0.15)[colnames(b_lysophospholipid_0.15) %in% pathways_to_present]
# colnames(q_values_matrix_lysophospholipid_0.15_filtered) <- colnames(q_values_matrix_lysophospholipid_0.15)[colnames(q_values_matrix_lysophospholipid_0.15) %in% pathways_to_present]
# 
# 
# # Create the asterisks matrix for the filtered data
# asterisks_matrix_lysophospholipid_0.15_filtered <- matrix("", nrow = nrow(q_values_matrix_lysophospholipid_0.15_filtered), ncol = ncol(q_values_matrix_lysophospholipid_0.15_filtered))
# asterisks_matrix_lysophospholipid_0.15_filtered[q_values_matrix_lysophospholipid_0.15_filtered < 0.05] <- "***"
# asterisks_matrix_lysophospholipid_0.15_filtered[q_values_matrix_lysophospholipid_0.15_filtered >= 0.05 & q_values_matrix_lysophospholipid_0.15_filtered < 0.1] <- "**"
# asterisks_matrix_lysophospholipid_0.15_filtered[q_values_matrix_lysophospholipid_0.15_filtered >= 0.1 & q_values_matrix_lysophospholipid_0.15_filtered < 0.25] <- "*"
# 
# # Create the filtered heatmap
# # Create breaks with 0 in the middle
# breaks <- seq(-0.2, 0.2, length.out = 51)
# 
# # Ensure the breaks include 0
# mid <- which.min(abs(breaks + 0.03))
# 
# # Create a custom color palette
# my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))
# 
# heatmap_lysophospholipid_0.15_filtered <- pheatmap(b_lysophospholipid_0.15_filtered,
#                                               cluster_rows = TRUE,
#                                               cluster_cols = FALSE,
#                                               display_numbers = asterisks_matrix_lysophospholipid_0.15_filtered,
#                                               fontsize_number = 10,
#                                               color = my_color_palette,
#                                               breaks = breaks,
#                                               annotation_legend = TRUE,
#                                               annotation_names_row = TRUE,
#                                               annotation_names_col = TRUE,
#                                               angle_col = 90)
# 
# # Extract the filtered heatmap as a gtable object
# heatmap_lysophospholipid_0.15_filtered_gtable <- heatmap_lysophospholipid_0.15_filtered$gtable
# 
# # Create text grobs for titles
# title1_lysophospholipid_0.15_filtered <- textGrob("Selected microbial pathways (n = 1)", gp = gpar(fontsize = 13), hjust = 2, vjust = -0.5)
# title2_lysophospholipid_0.15_filtered <- textGrob("Lysophospholipid metabolites (n = 14)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 0)
# 
# # Combine the filtered heatmap and titles using arrangeGrob and grid.arrange
# combined_plot_lysophospholipid_0.15_filtered <- arrangeGrob(
#   arrangeGrob(title2_lysophospholipid_0.15_filtered, heatmap_lysophospholipid_0.15_filtered_gtable, ncol = 2, widths = c(0.2, 4)),
#   title1_lysophospholipid_0.15_filtered, ncol = 1, heights = c(4, 0.2)
# )
# 
# # Draw the combined filtered plot
# grid.newpage()
# grid.draw(combined_plot_lysophospholipid_0.15_filtered)
# 
# # Save the filtered combined plot
# # ggsave("heatmap_selected_species_lysophospholipid_filtered.pdf", combined_plot_lysophospholipid_0.15_filtered, width = 4, height = 13, dpi=1200, limitsize = FALSE)







#### All species vs Lysophospholipid metabolites ####
# read in all_associations.txt got from running HALLA
r_diff_lysophospholipid_taxo_all <- read.table('halla_lysophospholipid_species_residual/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")

# Data cleaning
test2_lysophospholipid_taxo_all <- r_diff_lysophospholipid_taxo_all

# Subset Relevant Columns
a_lysophospholipid_taxo_all <- test2_lysophospholipid_taxo_all[, colnames(test2_lysophospholipid_taxo_all) %in% c("X_features", "Y_features", "association")]
a1_lysophospholipid_taxo_all <- test2_lysophospholipid_taxo_all[, colnames(test2_lysophospholipid_taxo_all) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysophospholipid_taxo_all <- a_lysophospholipid_taxo_all %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysophospholipid_taxo_all <- a1_lysophospholipid_taxo_all %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysophospholipid_taxo_all <- matrix("", nrow = nrow(q_values_matrix_lysophospholipid_taxo_all), ncol = ncol(q_values_matrix_lysophospholipid_taxo_all))
asterisks_matrix_lysophospholipid_taxo_all[q_values_matrix_lysophospholipid_taxo_all < 0.05] <- "***"
asterisks_matrix_lysophospholipid_taxo_all[q_values_matrix_lysophospholipid_taxo_all >= 0.05 & q_values_matrix_lysophospholipid_taxo_all < 0.1] <- "**"
asterisks_matrix_lysophospholipid_taxo_all[q_values_matrix_lysophospholipid_taxo_all >= 0.1 & q_values_matrix_lysophospholipid_taxo_all < 0.25] <- "*"

# Create the heatmap
heatmap_lysophospholipid_taxo_all <- pheatmap(b_lysophospholipid_taxo_all,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     display_numbers = asterisks_matrix_lysophospholipid_taxo_all,
                                     fontsize_number = 10,
                                     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                     main = "Heatmap of correlations between all species and lysophospholipid metabolites ",
                                     annotation_legend = TRUE,
                                     annotation_names_row = TRUE,
                                     annotation_names_col = TRUE,
                                     angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysophospholipid_gtable_taxo_all <- heatmap_lysophospholipid_taxo_all$gtable

# Create text grobs for titles
title1_lysophospholipid_taxo_all <- textGrob("All species (n = 467)", gp = gpar(fontsize = 13), hjust = 1, vjust = 0)
title2_lysophospholipid_taxo_all <- textGrob("Lysophospholipid metabolites (n = 33)", rot = 90, gp = gpar(fontsize = 13), hjust = 0, vjust = 15)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysophospholipid_taxo_all <- arrangeGrob(
  arrangeGrob(title2_lysophospholipid_taxo_all, heatmap_lysophospholipid_gtable_taxo_all, ncol = 2, widths = c(0.2, 4)),
  title1_lysophospholipid_taxo_all, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysophospholipid_taxo_all)

# Save the combined plot
# ggsave("heatmap_all_species_lysophospholipids_residual.png", combined_plot_lysophospholipid_taxo_all, width = 90, height = 13, dpi=300, limitsize = FALSE)




#### All species vs Lysophospholipid metabolites (rho > 0.15 or < -0.15) ####
test2_lysophospholipid_taxo_0.15 <- test2_lysophospholipid_taxo_all

# Filter rows where the association is < -0.15 or > 0.15
filtered_data <- test2_lysophospholipid_taxo_all %>%
  filter(association < -0.15 | association > 0.15)

# Get the unique X_features from the filtered data
unique_features <- filtered_data %>% select(X_features) %>% unique()
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysophospholipids_0.15_bugs.xlsx"
# write.xlsx(unique_features, file_path, rowNames = TRUE)

# Merge the unique features with the original data to keep relevant rows
test2_lysophospholipid_taxo_0.15 <- test2_lysophospholipid_taxo_all %>%
  inner_join(unique_features, by = c("X_features"))


# Subset Relevant Columns
a_lysophospholipid_taxo_0.15 <- test2_lysophospholipid_taxo_0.15[, colnames(test2_lysophospholipid_taxo_0.15) %in% c("X_features", "Y_features", "association")]
a1_lysophospholipid_taxo_0.15 <- test2_lysophospholipid_taxo_0.15[, colnames(test2_lysophospholipid_taxo_0.15) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysophospholipid_taxo_0.15 <- a_lysophospholipid_taxo_0.15 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysophospholipid_taxo_0.15 <- a1_lysophospholipid_taxo_0.15 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysophospholipid_taxo_0.15 <- matrix("", nrow = nrow(q_values_matrix_lysophospholipid_taxo_0.15), ncol = ncol(q_values_matrix_lysophospholipid_taxo_0.15))
asterisks_matrix_lysophospholipid_taxo_0.15[q_values_matrix_lysophospholipid_taxo_0.15 < 0.05] <- "***"
asterisks_matrix_lysophospholipid_taxo_0.15[q_values_matrix_lysophospholipid_taxo_0.15 >= 0.05 & q_values_matrix_lysophospholipid_taxo_0.15 < 0.1] <- "**"
asterisks_matrix_lysophospholipid_taxo_0.15[q_values_matrix_lysophospholipid_taxo_0.15 >= 0.1 & q_values_matrix_lysophospholipid_taxo_0.15 < 0.25] <- "*"

# Create the heatmap
heatmap_lysophospholipid_taxo_0.15 <- pheatmap(b_lysophospholipid_taxo_0.15,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         display_numbers = asterisks_matrix_lysophospholipid_taxo_0.15,
                                         fontsize_number = 10,
                                         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                         annotation_legend = TRUE,
                                         annotation_names_row = TRUE,
                                         annotation_names_col = TRUE,
                                         angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysophospholipid_taxo_0.15_gtable <- heatmap_lysophospholipid_taxo_0.15$gtable

# Create text grobs for titles
title1_lysophospholipid_taxo_0.15 <- textGrob("Selected species (n = 354)", gp = gpar(fontsize = 13), hjust = 1, vjust = -3)
title2_lysophospholipid_taxo_0.15 <- textGrob("Lysophospholipid metabolites (n = 14)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.15, vjust = 9)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysophospholipid_taxo_0.15 <- arrangeGrob(
  arrangeGrob(title2_lysophospholipid_taxo_0.15, heatmap_lysophospholipid_taxo_0.15_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysophospholipid_taxo_0.15, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysophospholipid_taxo_0.15)

# Save the combined plot
# ggsave("heatmap_species_lysophospholipid_0.15_residual.pdf", combined_plot_lysophospholipid_taxo_0.15, width = 60, height = 7.5, dpi=1200, limitsize = FALSE)



### Simplified version
# List of species to present
species_to_present <- c(
  "Veillonella_dispar",
  "Veillonella_parvula",
  "Haemophilus_parainfluenzae",
  "Methylobacterium_SGB15164",
  "GGB33586_SGB53517",
  "Clostridiales_bacterium_NSJ_32",
  "Clostridiales_bacterium",
  "Collinsella_SGB4121",
  "Mediterraneibacter_butyricigenes",
  "GGB9534_SGB14937",
  "Blautia_glucerasea",
  "GGB9534_SGB14937",
  "Ruminococcus_gnavus",
  "Akkermansia_muciniphila"
)

# Filter the matrix and annotations to only include the selected species
b_lysophospholipid_taxo_0.15_filtered <- b_lysophospholipid_taxo_0.15[, colnames(b_lysophospholipid_taxo_0.15) %in% species_to_present]
q_values_matrix_lysophospholipid_taxo_0.15_filtered <- q_values_matrix_lysophospholipid_taxo_0.15[, colnames(b_lysophospholipid_taxo_0.15) %in% species_to_present]

# Create the asterisks matrix for the filtered data
asterisks_matrix_lysophospholipid_taxo_0.15_filtered <- matrix("", nrow = nrow(q_values_matrix_lysophospholipid_taxo_0.15_filtered), ncol = ncol(q_values_matrix_lysophospholipid_taxo_0.15_filtered))
asterisks_matrix_lysophospholipid_taxo_0.15_filtered[q_values_matrix_lysophospholipid_taxo_0.15_filtered < 0.05] <- "***"
asterisks_matrix_lysophospholipid_taxo_0.15_filtered[q_values_matrix_lysophospholipid_taxo_0.15_filtered >= 0.05 & q_values_matrix_lysophospholipid_taxo_0.15_filtered < 0.1] <- "**"
asterisks_matrix_lysophospholipid_taxo_0.15_filtered[q_values_matrix_lysophospholipid_taxo_0.15_filtered >= 0.1 & q_values_matrix_lysophospholipid_taxo_0.15_filtered < 0.25] <- "*"

# Create the filtered heatmap
heatmap_lysophospholipid_taxo_0.15_filtered <- pheatmap(b_lysophospholipid_taxo_0.15_filtered,
                                                   cluster_rows = TRUE,
                                                   cluster_cols = TRUE,
                                                   display_numbers = asterisks_matrix_lysophospholipid_taxo_0.15_filtered,
                                                   fontsize_number = 10,
                                                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                   main = "Heatmap of correlations between selected species and lysophospholipid metabolites",
                                                   annotation_legend = TRUE,
                                                   annotation_names_row = TRUE,
                                                   annotation_names_col = TRUE,
                                                   border_color = NA,
                                                   angle_col = 90)

# Extract the filtered heatmap as a gtable object
heatmap_lysophospholipid_taxo_0.15_filtered_gtable <- heatmap_lysophospholipid_taxo_0.15_filtered$gtable

# Create text grobs for titles
title1_lysophospholipid_taxo_0.15_filtered <- textGrob("Selected microbial species (n = 14)", gp = gpar(fontsize = 13), hjust = 2, vjust = -0.5)
title2_lysophospholipid_taxo_0.15_filtered <- textGrob("lysophospholipid metabolites (n = 14)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 0)

# Combine the filtered heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysophospholipid_taxo_0.15_filtered <- arrangeGrob(
  arrangeGrob(title2_lysophospholipid_taxo_0.15_filtered, heatmap_lysophospholipid_taxo_0.15_filtered_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysophospholipid_taxo_0.15_filtered, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined filtered plot
grid.newpage()
grid.draw(combined_plot_lysophospholipid_taxo_0.15_filtered)

# Save the filtered combined plot
# ggsave("heatmap_selected_species_lysophospholipids_filtered_residual.pdf", combined_plot_lysophospholipid_taxo_0.15_filtered, width = 8, height = 8, dpi=1200, limitsize = FALSE)



#### Species-lysophospholipid correlations concordant with abundance analysis ####
abundance_lysophospholipid <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross/all_results.tsv")
abundance_lysophospholipid <- abundance_lysophospholipid %>% filter(metadata=="mc_all") %>% 
  select(1,3,4) %>%
  filter(feature %in% c("Veillonella_dispar", "Veillonella_parvula", "Haemophilus_parainfluenzae", "Methylobacterium_SGB15164", "GGB33586_SGB53517", "Clostridiales_bacterium", 
                        "Collinsella_SGB4121", "Mediterraneibacter_butyricigenes", "GGB9534_SGB14937", "Blautia_glucerasea"))

abundance_lysophospholipid$coef <- -(abundance_lysophospholipid$coef)

abundance_lysophospholipid_avg <- abundance_lysophospholipid %>%
  group_by(feature) %>%
  summarise(avg_coef = mean(coef))

function_lysophospholipid <- b_lysophospholipid_taxo_all[, colnames(b_lysophospholipid_taxo_all) %in% c("Veillonella_dispar", "Veillonella_parvula", "Haemophilus_parainfluenzae", "Methylobacterium_SGB15164", "GGB33586_SGB53517", "Clostridiales_bacterium", 
                                                                                                        "Collinsella_SGB4121", "Mediterraneibacter_butyricigenes", "GGB9534_SGB14937", "Blautia_glucerasea")]%>% t() %>% as.data.frame()

function_lysophospholipid$feature <- rownames(function_lysophospholipid)

concordance_lysophospholipid <- abundance_lysophospholipid_avg %>% left_join(function_lysophospholipid)

concordance_lysophospholipid <- concordance_lysophospholipid %>% 
  mutate(enriched = case_when(
    avg_coef > 0 ~ "1",
    avg_coef < 0 ~ "0"
  ))

concordance_lysophospholipid$species <- concordance_lysophospholipid$feature
concordance_lysophospholipid$feature <- NULL

fig_concordance_lysophospholipid <- ggplot(concordance_lysophospholipid, aes(x = avg_coef, y = `1-linoleoyl-GPG (18:2)*`, color = as.factor(enriched))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point() +
  geom_text_repel(aes(label = species), size = 4) + 
  labs(x = "Species abundance in MC compared to controls (beta coeffecient from regression)", 
       y = "Spearman rho of the species with 1-linoleoyl-GPG (18:2)*", 
       color = NULL) +
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("navy", "firebrick3"),
                     labels = c("Less abundant in MC & negatively correlated with lysophospholipids", 
                                "More abundant in MC & positively correlated with lysophospholipids"))

fig_concordance_lysophospholipid

# ggsave("fig_concordance_lysophospholipid_residual.pdf", plot = fig_concordance_lysophospholipid, width = 12, height = 9, dpi = 1200)






#### All pathways vs Lysoplasmalogen metabolites ####
# Prepare input for halla
Lysoplasmalogen_halla <- read_excel("met_halla_residual.xlsx") %>% t() %>%
  as.data.frame() %>% 
  {setNames(.[-1, ], unlist(.[1, ]))}
Lysoplasmalogen_halla$id <- rownames(Lysoplasmalogen_halla)

# To select Lysoplasmalogen
Lysoplasmalogen <- subset(altered_114, grepl("Lysoplasmalogen", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
Lysoplasmalogen_metabolites <- Lysoplasmalogen$metabolite
Lysoplasmalogen_halla <- Lysoplasmalogen_halla %>% select(all_of(Lysoplasmalogen_metabolites))
Lysoplasmalogen_halla$id <- rownames(Lysoplasmalogen_halla)

# Change the ID to match microbiome ID
Lysoplasmalogen_halla <- Lysoplasmalogen_halla %>%
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
rownames(Lysoplasmalogen_halla) <- Lysoplasmalogen_halla$id
Lysoplasmalogen_halla$id <- NULL

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_lysoplasmalogen_halla_residual.xlsx"
# write.xlsx(Lysoplasmalogen_halla, file_path, rowNames = TRUE)

# read in all_associations.txt got from running HALLA
r_diff_lysoplasmalogen <- read.table('halla_lysoplasmalogen_path_residual/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")

# Data cleaning
test2_lysoplasmalogen <- r_diff_lysoplasmalogen
 

# Subset Relevant Columns
a_lysoplasmalogen <- test2_lysoplasmalogen[, colnames(test2_lysoplasmalogen) %in% c("X_features", "Y_features", "association")]
a1_lysoplasmalogen <- test2_lysoplasmalogen[, colnames(test2_lysoplasmalogen) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysoplasmalogen <- a_lysoplasmalogen %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysoplasmalogen <- a1_lysoplasmalogen %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysoplasmalogen <- matrix("", nrow = nrow(q_values_matrix_lysoplasmalogen), ncol = ncol(q_values_matrix_lysoplasmalogen))
asterisks_matrix_lysoplasmalogen[q_values_matrix_lysoplasmalogen < 0.05] <- "***"
asterisks_matrix_lysoplasmalogen[q_values_matrix_lysoplasmalogen >= 0.05 & q_values_matrix_lysoplasmalogen < 0.1] <- "**"
asterisks_matrix_lysoplasmalogen[q_values_matrix_lysoplasmalogen >= 0.1 & q_values_matrix_lysoplasmalogen < 0.25] <- "*"

# Create the heatmap
heatmap_lysoplasmalogen <- pheatmap(b_lysoplasmalogen,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     display_numbers = asterisks_matrix_lysoplasmalogen,
                                     fontsize_number = 10,
                                     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                     main = "Heatmap of correlations between all microbial pathways and lysoplasmalogen metabolites ",
                                     annotation_legend = TRUE,
                                     annotation_names_row = TRUE,
                                     annotation_names_col = TRUE,
                                     angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysoplasmalogen_gtable <- heatmap_lysoplasmalogen$gtable

# Create text grobs for titles
title1_lysoplasmalogen <- textGrob("All mircobial pathways (n = 250)", gp = gpar(fontsize = 13), hjust = 2, vjust = -7)
title2_lysoplasmalogen <- textGrob("Lysoplasmalogen metabolites (n = 10)", rot = 90, gp = gpar(fontsize = 13), hjust = -0.4, vjust = 10)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysoplasmalogen <- arrangeGrob(
  arrangeGrob(title2_lysoplasmalogen, heatmap_lysoplasmalogen_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysoplasmalogen, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysoplasmalogen)

# Save the combined plot
# ggsave("heatmap_all_paths_lysoplasmalogens_residual.pdf", combined_plot_lysoplasmalogen, width = 70, height = 12, dpi=1200, limitsize = FALSE)




#### All pathways vs Lysoplasmalogen metabolites (rho > 0.15 or < -0.15) ####
test2_lysoplasmalogen_0.15 <- test2_lysoplasmalogen

# Filter rows where the association is < -0.15 or > 0.15
filtered_data <- test2_lysoplasmalogen %>%
  filter(association < -0.15 | association > 0.15)

# Get the unique X_features from the filtered data
unique_features <- filtered_data %>% select(X_features) %>% unique()
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysoplasmalogens_0.15_paths.xlsx"
# write.xlsx(unique_features, file_path, rowNames = TRUE)


# Merge the unique features with the original data to keep relevant rows
test2_lysoplasmalogen_0.15 <- test2_lysoplasmalogen %>%
  inner_join(unique_features, by = c("X_features"))


# Subset Relevant Columns
a_lysoplasmalogen_0.15 <- test2_lysoplasmalogen_0.15[, colnames(test2_lysoplasmalogen_0.15) %in% c("X_features", "Y_features", "association")]
a1_lysoplasmalogen_0.15 <- test2_lysoplasmalogen_0.15[, colnames(test2_lysoplasmalogen_0.15) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysoplasmalogen_0.15 <- a_lysoplasmalogen_0.15 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysoplasmalogen_0.15 <- a1_lysoplasmalogen_0.15 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysoplasmalogen_0.15 <- matrix("", nrow = nrow(q_values_matrix_lysoplasmalogen_0.15), ncol = ncol(q_values_matrix_lysoplasmalogen_0.15))
asterisks_matrix_lysoplasmalogen_0.15[q_values_matrix_lysoplasmalogen_0.15 < 0.05] <- "***"
asterisks_matrix_lysoplasmalogen_0.15[q_values_matrix_lysoplasmalogen_0.15 >= 0.05 & q_values_matrix_lysoplasmalogen_0.15 < 0.1] <- "**"
asterisks_matrix_lysoplasmalogen_0.15[q_values_matrix_lysoplasmalogen_0.15 >= 0.1 & q_values_matrix_lysoplasmalogen_0.15 < 0.25] <- "*"

# Create the heatmap
heatmap_lysoplasmalogen_0.15 <- pheatmap(b_lysoplasmalogen_0.15,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         display_numbers = asterisks_matrix_lysoplasmalogen_0.15,
                                         fontsize_number = 10,
                                         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                         annotation_legend = TRUE,
                                         annotation_names_row = TRUE,
                                         annotation_names_col = TRUE,
                                         angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysoplasmalogen_0.15_gtable <- heatmap_lysoplasmalogen_0.15$gtable

# Create text grobs for titles
title1_lysoplasmalogen_0.15 <- textGrob("Selected mircobial pathways (n = 83)", gp = gpar(fontsize = 13), hjust = 0.8, vjust = -7)
title2_lysoplasmalogen_0.15 <- textGrob("Lysoplasmalogen metabolites (n = 8)", rot = 90, gp = gpar(fontsize = 13), hjust = -0.5, vjust = 5)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysoplasmalogen_0.15 <- arrangeGrob(
  arrangeGrob(title2_lysoplasmalogen_0.15, heatmap_lysoplasmalogen_0.15_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysoplasmalogen_0.15, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysoplasmalogen_0.15)

# Save the combined plot
# ggsave("heatmap_lysoplasmalogens_0.15_residual.pdf", combined_plot_lysoplasmalogen_0.15, width = 30, height = 10.5, dpi=1200, limitsize = FALSE)



## Simplified version
# List of pathways to present
df <- as.data.frame(colnames(b_lysoplasmalogen_0.15))
pathways_to_present <- c(
  "GLYCOGENSYNTH-PWY: glycogen biosynthesis I (from ADP-D-Glucose)"
)

# Filter the matrix and annotations to only include the selected pathways
b_lysoplasmalogen_0.15_filtered <- b_lysoplasmalogen_0.15[, colnames(b_lysoplasmalogen_0.15) %in% pathways_to_present] %>% as.matrix()
q_values_matrix_lysoplasmalogen_0.15_filtered <- q_values_matrix_lysoplasmalogen_0.15[, colnames(q_values_matrix_lysoplasmalogen_0.15) %in% pathways_to_present] %>% as.matrix()
colnames(b_lysoplasmalogen_0.15_filtered) <- colnames(b_lysoplasmalogen_0.15)[colnames(b_lysoplasmalogen_0.15) %in% pathways_to_present]
colnames(q_values_matrix_lysoplasmalogen_0.15_filtered) <- colnames(q_values_matrix_lysoplasmalogen_0.15)[colnames(q_values_matrix_lysoplasmalogen_0.15) %in% pathways_to_present]


# Create the asterisks matrix for the filtered data
asterisks_matrix_lysoplasmalogen_0.15_filtered <- matrix("", nrow = nrow(q_values_matrix_lysoplasmalogen_0.15_filtered), ncol = ncol(q_values_matrix_lysoplasmalogen_0.15_filtered))
asterisks_matrix_lysoplasmalogen_0.15_filtered[q_values_matrix_lysoplasmalogen_0.15_filtered < 0.05] <- "***"
asterisks_matrix_lysoplasmalogen_0.15_filtered[q_values_matrix_lysoplasmalogen_0.15_filtered >= 0.05 & q_values_matrix_lysoplasmalogen_0.15_filtered < 0.1] <- "**"
asterisks_matrix_lysoplasmalogen_0.15_filtered[q_values_matrix_lysoplasmalogen_0.15_filtered >= 0.1 & q_values_matrix_lysoplasmalogen_0.15_filtered < 0.25] <- "*"

# Create the filtered heatmap
# Create breaks with 0 in the middle
breaks <- seq(-0.23, 0.23, length.out = 51)

# Ensure the breaks include 0
mid <- which.min(abs(breaks-0))

# Create a custom color palette
my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))

heatmap_lysoplasmalogen_0.15_filtered <- pheatmap(b_lysoplasmalogen_0.15_filtered,
                                                   cluster_rows = TRUE,
                                                   cluster_cols = FALSE,
                                                   display_numbers = asterisks_matrix_lysoplasmalogen_0.15_filtered,
                                                   fontsize_number = 10,
                                                   color = my_color_palette,
                                                   breaks = breaks,
                                                   annotation_legend = TRUE,
                                                   annotation_names_row = TRUE,
                                                   annotation_names_col = TRUE,
                                                   angle_col = 90)

# Extract the filtered heatmap as a gtable object
heatmap_lysoplasmalogen_0.15_filtered_gtable <- heatmap_lysoplasmalogen_0.15_filtered$gtable

# Create text grobs for titles
title1_lysoplasmalogen_0.15_filtered <- textGrob("Selected microbial pathways (n = 1)", gp = gpar(fontsize = 13), hjust = 2, vjust = -0.5)
title2_lysoplasmalogen_0.15_filtered <- textGrob("lysoplasmalogen metabolites (n = 8)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 0)

# Combine the filtered heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysoplasmalogen_0.15_filtered <- arrangeGrob(
  arrangeGrob(title2_lysoplasmalogen_0.15_filtered, heatmap_lysoplasmalogen_0.15_filtered_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysoplasmalogen_0.15_filtered, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined filtered plot
grid.newpage()
grid.draw(combined_plot_lysoplasmalogen_0.15_filtered)

# Save the filtered combined plot
# ggsave("heatmap_selected_path_lysoplasmalogen_filtered_residual.pdf", combined_plot_lysoplasmalogen_0.15_filtered, width = 4.2, height = 7.5, dpi=1200, limitsize = FALSE)





#### All species vs Lysoplasmalogen metabolites ####
# read in all_associations.txt got from running HALLA
r_diff_lysoplasmalogen_taxo_all <- read.table('halla_lysoplasmalogen_species_residual/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")

# Data cleaning
test2_lysoplasmalogen_taxo_all <- r_diff_lysoplasmalogen_taxo_all

# Subset Relevant Columns
a_lysoplasmalogen_taxo_all <- test2_lysoplasmalogen_taxo_all[, colnames(test2_lysoplasmalogen_taxo_all) %in% c("X_features", "Y_features", "association")]
a1_lysoplasmalogen_taxo_all <- test2_lysoplasmalogen_taxo_all[, colnames(test2_lysoplasmalogen_taxo_all) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysoplasmalogen_taxo_all <- a_lysoplasmalogen_taxo_all %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysoplasmalogen_taxo_all <- a1_lysoplasmalogen_taxo_all %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysoplasmalogen_taxo_all <- matrix("", nrow = nrow(q_values_matrix_lysoplasmalogen_taxo_all), ncol = ncol(q_values_matrix_lysoplasmalogen_taxo_all))
asterisks_matrix_lysoplasmalogen_taxo_all[q_values_matrix_lysoplasmalogen_taxo_all < 0.05] <- "***"
asterisks_matrix_lysoplasmalogen_taxo_all[q_values_matrix_lysoplasmalogen_taxo_all >= 0.05 & q_values_matrix_lysoplasmalogen_taxo_all < 0.1] <- "**"
asterisks_matrix_lysoplasmalogen_taxo_all[q_values_matrix_lysoplasmalogen_taxo_all >= 0.1 & q_values_matrix_lysoplasmalogen_taxo_all < 0.25] <- "*"

# Create the heatmap
heatmap_lysoplasmalogen_taxo_all <- pheatmap(b_lysoplasmalogen_taxo_all,
                                              cluster_rows = TRUE,
                                              cluster_cols = TRUE,
                                              display_numbers = asterisks_matrix_lysoplasmalogen_taxo_all,
                                              fontsize_number = 10,
                                              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                              main = "Heatmap of correlations between all species and lysoplasmalogen metabolites ",
                                              annotation_legend = TRUE,
                                              annotation_names_row = TRUE,
                                              annotation_names_col = TRUE,
                                              angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysoplasmalogen_gtable_taxo_all <- heatmap_lysoplasmalogen_taxo_all$gtable

# Create text grobs for titles
title1_lysoplasmalogen_taxo_all <- textGrob("All species (n = 461)", gp = gpar(fontsize = 13), hjust = 1, vjust = 0)
title2_lysoplasmalogen_taxo_all <- textGrob("Lysoplasmalogen metabolites (n = 10)", rot = 90, gp = gpar(fontsize = 13), hjust = 0, vjust = 15)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysoplasmalogen_taxo_all <- arrangeGrob(
  arrangeGrob(title2_lysoplasmalogen_taxo_all, heatmap_lysoplasmalogen_gtable_taxo_all, ncol = 2, widths = c(0.2, 4)),
  title1_lysoplasmalogen_taxo_all, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysoplasmalogen_taxo_all)

# Save the combined plot
# ggsave("heatmap_all_species_lysoplasmalogens_residual.pdf", combined_plot_lysoplasmalogen_taxo_all, width = 90, height = 7, dpi=300, limitsize = FALSE)




#### All species vs Lysoplasmalogen metabolites (rho > 0.15 or < -0.15) ####
test2_lysoplasmalogen_taxo_0.15 <- test2_lysoplasmalogen_taxo_all

# Filter rows where the association is < -0.15 or > 0.15
filtered_data <- test2_lysoplasmalogen_taxo_all %>%
  filter(association < -0.15 | association > 0.15)

# Get the unique X_features from the filtered data
unique_features <- filtered_data %>% select(X_features) %>% unique()
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysoplasmalogens_0.15_bugs.xlsx"
# write.xlsx(unique_features, file_path, rowNames = TRUE)

# Merge the unique features with the original data to keep relevant rows
test2_lysoplasmalogen_taxo_0.15 <- test2_lysoplasmalogen_taxo_all %>%
  inner_join(unique_features, by = c("X_features"))


# Subset Relevant Columns
a_lysoplasmalogen_taxo_0.15 <- test2_lysoplasmalogen_taxo_0.15[, colnames(test2_lysoplasmalogen_taxo_0.15) %in% c("X_features", "Y_features", "association")]
a1_lysoplasmalogen_taxo_0.15 <- test2_lysoplasmalogen_taxo_0.15[, colnames(test2_lysoplasmalogen_taxo_0.15) %in% c("X_features", "Y_features", "q.values")]

# Reshape data
b_lysoplasmalogen_taxo_0.15 <- a_lysoplasmalogen_taxo_0.15 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
q_values_matrix_lysoplasmalogen_taxo_0.15 <- a1_lysoplasmalogen_taxo_0.15 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
  as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()


asterisks_matrix_lysoplasmalogen_taxo_0.15 <- matrix("", nrow = nrow(q_values_matrix_lysoplasmalogen_taxo_0.15), ncol = ncol(q_values_matrix_lysoplasmalogen_taxo_0.15))
asterisks_matrix_lysoplasmalogen_taxo_0.15[q_values_matrix_lysoplasmalogen_taxo_0.15 < 0.05] <- "***"
asterisks_matrix_lysoplasmalogen_taxo_0.15[q_values_matrix_lysoplasmalogen_taxo_0.15 >= 0.05 & q_values_matrix_lysoplasmalogen_taxo_0.15 < 0.1] <- "**"
asterisks_matrix_lysoplasmalogen_taxo_0.15[q_values_matrix_lysoplasmalogen_taxo_0.15 >= 0.1 & q_values_matrix_lysoplasmalogen_taxo_0.15 < 0.25] <- "*"

# Create the heatmap
heatmap_lysoplasmalogen_taxo_0.15 <- pheatmap(b_lysoplasmalogen_taxo_0.15,
                                              cluster_rows = TRUE,
                                              cluster_cols = TRUE,
                                              display_numbers = asterisks_matrix_lysoplasmalogen_taxo_0.15,
                                              fontsize_number = 10,
                                              color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                              annotation_legend = TRUE,
                                              annotation_names_row = TRUE,
                                              annotation_names_col = TRUE,
                                              angle_col = 90)

# Extract the heatmap as a gtable object
heatmap_lysoplasmalogen_taxo_0.15_gtable <- heatmap_lysoplasmalogen_taxo_0.15$gtable

# Create text grobs for titles
title1_lysoplasmalogen_taxo_0.15 <- textGrob("Selected species (n = 272)", gp = gpar(fontsize = 13), hjust = 0, vjust = -1)
title2_lysoplasmalogen_taxo_0.15 <- textGrob("Lysoplasmalogen metabolites (n = 8)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 7)

# Combine the heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysoplasmalogen_taxo_0.15 <- arrangeGrob(
  arrangeGrob(title2_lysoplasmalogen_taxo_0.15, heatmap_lysoplasmalogen_taxo_0.15_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysoplasmalogen_taxo_0.15, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined plot
grid.newpage()
grid.draw(combined_plot_lysoplasmalogen_taxo_0.15)

# Save the combined plot
# ggsave("heatmap_species_lysoplasmalogen_0.15_residual.pdf", combined_plot_lysoplasmalogen_taxo_0.15, width = 50, height = 7, dpi=1200, limitsize = FALSE)



### Simplified version
# List of species to present
species_to_present <- c(
  "Veillonella_parvula",
  "Haemophilus_parainfluenzae",
  "Methylobacterium_SGB15164",
  "Mediterraneibacter_butyricigenes",
  "GGB9534_SGB14937",
  "Ruminococcus_gnavus",
  "Akkermansia_muciniphila"
)

# Filter the matrix and annotations to only include the selected species
b_lysoplasmalogen_taxo_0.15_filtered <- b_lysoplasmalogen_taxo_0.15[, colnames(b_lysoplasmalogen_taxo_0.15) %in% species_to_present]
q_values_matrix_lysoplasmalogen_taxo_0.15_filtered <- q_values_matrix_lysoplasmalogen_taxo_0.15[, colnames(b_lysoplasmalogen_taxo_0.15) %in% species_to_present]

# Create the asterisks matrix for the filtered data
asterisks_matrix_lysoplasmalogen_taxo_0.15_filtered <- matrix("", nrow = nrow(q_values_matrix_lysoplasmalogen_taxo_0.15_filtered), ncol = ncol(q_values_matrix_lysoplasmalogen_taxo_0.15_filtered))
asterisks_matrix_lysoplasmalogen_taxo_0.15_filtered[q_values_matrix_lysoplasmalogen_taxo_0.15_filtered < 0.05] <- "***"
asterisks_matrix_lysoplasmalogen_taxo_0.15_filtered[q_values_matrix_lysoplasmalogen_taxo_0.15_filtered >= 0.05 & q_values_matrix_lysoplasmalogen_taxo_0.15_filtered < 0.1] <- "**"
asterisks_matrix_lysoplasmalogen_taxo_0.15_filtered[q_values_matrix_lysoplasmalogen_taxo_0.15_filtered >= 0.1 & q_values_matrix_lysoplasmalogen_taxo_0.15_filtered < 0.25] <- "*"

# Create breaks with 0 in the middle
breaks <- seq(-0.25, 0.25, length.out = 51)

# Ensure the breaks include 0
mid <- which.min(abs(breaks-0))

# Create a custom color palette
my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))

# Create the filtered heatmap
heatmap_lysoplasmalogen_taxo_0.15_filtered <- pheatmap(b_lysoplasmalogen_taxo_0.15_filtered,
                                                       cluster_rows = TRUE,
                                                       cluster_cols = TRUE,
                                                       display_numbers = asterisks_matrix_lysoplasmalogen_taxo_0.15_filtered,
                                                       fontsize_number = 10,
                                                       color = my_color_palette,
                                                       main = "Heatmap of correlations between selected species and lysoplasmalogen metabolites",
                                                       breaks = breaks,
                                                       annotation_legend = TRUE,
                                                       annotation_names_row = TRUE,
                                                       annotation_names_col = TRUE,
                                                       border_color = NA,
                                                       angle_col = 90)

# Extract the filtered heatmap as a gtable object
heatmap_lysoplasmalogen_taxo_0.15_filtered_gtable <- heatmap_lysoplasmalogen_taxo_0.15_filtered$gtable

# Create text grobs for titles
title1_lysoplasmalogen_taxo_0.15_filtered <- textGrob("Selected microbial species (n = 7)", gp = gpar(fontsize = 13), hjust = 2, vjust = -0.5)
title2_lysoplasmalogen_taxo_0.15_filtered <- textGrob("lysoplasmalogen metabolites (n = 8)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 0)

# Combine the filtered heatmap and titles using arrangeGrob and grid.arrange
combined_plot_lysoplasmalogen_taxo_0.15_filtered <- arrangeGrob(
  arrangeGrob(title2_lysoplasmalogen_taxo_0.15_filtered, heatmap_lysoplasmalogen_taxo_0.15_filtered_gtable, ncol = 2, widths = c(0.2, 4)),
  title1_lysoplasmalogen_taxo_0.15_filtered, ncol = 1, heights = c(4, 0.2)
)

# Draw the combined filtered plot
grid.newpage()
grid.draw(combined_plot_lysoplasmalogen_taxo_0.15_filtered)

# Save the filtered combined plot
# ggsave("heatmap_selected_species_lysoplasmalogens_filtered_residual.pdf", combined_plot_lysoplasmalogen_taxo_0.15_filtered, width = 8, height = 8, dpi=1200, limitsize = FALSE)



#### Species-lysoplasmalogen correlations concordant with abundance analysis ####
abundance_lysoplasmalogen <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross/all_results.tsv")
abundance_lysoplasmalogen <- abundance_lysoplasmalogen %>% filter(metadata=="mc_all") %>% 
  select(1,3,4) %>%
  filter(feature %in% c("Veillonella_parvula", "Haemophilus_parainfluenzae", "Methylobacterium_SGB15164", 
                        "Mediterraneibacter_butyricigenes", "GGB9534_SGB14937"))

abundance_lysoplasmalogen$coef <- -(abundance_lysoplasmalogen$coef)

abundance_lysoplasmalogen_avg <- abundance_lysoplasmalogen %>%
  group_by(feature) %>%
  summarise(avg_coef = mean(coef))

function_lysoplasmalogen <- b_lysoplasmalogen_taxo_all[, colnames(b_lysoplasmalogen_taxo_all) %in% c("Veillonella_parvula", "Haemophilus_parainfluenzae", "Methylobacterium_SGB15164", 
                                                                                                     "Mediterraneibacter_butyricigenes", "GGB9534_SGB14937")]%>% t() %>% as.data.frame()

function_lysoplasmalogen$feature <- rownames(function_lysoplasmalogen)

concordance_lysoplasmalogen <- abundance_lysoplasmalogen_avg %>% left_join(function_lysoplasmalogen)

concordance_lysoplasmalogen <- concordance_lysoplasmalogen %>% 
  mutate(enriched = case_when(
    avg_coef > 0 ~ "1",
    avg_coef < 0 ~ "0"
  ))

concordance_lysoplasmalogen$species <- concordance_lysoplasmalogen$feature
concordance_lysoplasmalogen$feature <- NULL

fig_concordance_lysoplasmalogen <- ggplot(concordance_lysoplasmalogen, aes(x = avg_coef, y = `1-palmityl-GPC (O-16:0)`, color = as.factor(enriched))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point() +
  geom_text_repel(aes(label = species), size = 4) + 
  labs(x = "Species abundance in MC compared to controls (beta coeffecient from regression)", 
       y = "Spearman rho of the species with 1-palmityl-GPC (O-16:0)", 
       color = NULL) +
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("navy", "firebrick3"),
                     labels = c("Less abundant in MC & negatively correlated with lysoplasmalogens", 
                                "More abundant in MC & positively correlated with lysoplasmalogens"))

fig_concordance_lysoplasmalogen

# ggsave("fig_concordance_lysoplasmalogen_residual.pdf", plot = fig_concordance_lysoplasmalogen, width = 12, height = 9, dpi = 1200)













# #### Sphingolipids 0.2: correlated pathways vs correlated species ####
# # read in all_associations.txt got from running HALLA
# r_diff_sphingolipids_bugs_0.2 <- read.table('/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/halla_sphingolipids_0.2_bugs_paths/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
# 
# # Data cleaning
# test_reduce2_sphingolipids_bugs_0.2 <- r_diff_sphingolipids_bugs_0.2
# test2_sphingolipids_bugs_0.2  <- test_reduce2_sphingolipids_bugs_0.2
# 
# # Subset Relevant Columns
# a_sphingolipids_bugs_0.2 <- test2_sphingolipids_bugs_0.2[, colnames(test2_sphingolipids_bugs_0.2) %in% c("X_features", "Y_features", "association")]
# a1_sphingolipids_bugs_0.2 <- test2_sphingolipids_bugs_0.2[, colnames(test2_sphingolipids_bugs_0.2) %in% c("X_features", "Y_features", "q.values")]
# 
# # Reshape data
# b_sphingolipids_bugs_0.2 <- a_sphingolipids_bugs_0.2 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
#   as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
# q_values_matrix_sphingolipids_bugs_0.2<- a1_sphingolipids_bugs_0.2 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
#   as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
# 
# 
# asterisks_matrix_sphingolipids_bugs_0.2 <- matrix("", nrow = nrow(q_values_matrix_sphingolipids_bugs_0.2), ncol = ncol(q_values_matrix_sphingolipids_bugs_0.2))
# asterisks_matrix_sphingolipids_bugs_0.2[q_values_matrix_sphingolipids_bugs_0.2 < 0.05] <- "***"
# asterisks_matrix_sphingolipids_bugs_0.2[q_values_matrix_sphingolipids_bugs_0.2 >= 0.05 & q_values_matrix_sphingolipids_bugs_0.2 < 0.1] <- "**"
# asterisks_matrix_sphingolipids_bugs_0.2[q_values_matrix_sphingolipids_bugs_0.2 >= 0.1 & q_values_matrix_sphingolipids_bugs_0.2 < 0.25] <- "*"
# 
# # Create the heatmap
# heatmap_sphingolipids_bugs_0.2 <- pheatmap(b_sphingolipids_bugs_0.2,
#                                            cluster_rows = TRUE,
#                                            cluster_cols = TRUE,
#                                            display_numbers = asterisks_matrix_sphingolipids_bugs_0.2,
#                                            fontsize_number = 10,
#                                            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#                                            main = "Heatmap of correlations between microbial species and pathways correlated with sphingolipids",
#                                            annotation_legend = TRUE,
#                                            annotation_names_row = TRUE,
#                                            annotation_names_col = TRUE,
#                                            angle_col = 90)
# 
# # Extract the heatmap as a gtable object
# heatmap_sphingolipids_bugs_0.2_gtable <- heatmap_sphingolipids_bugs_0.2$gtable
# 
# # Create text grobs for titles
# title1_sphingolipids_bugs_0.2 <- textGrob("Species correlated with sphingolipids (n = 160)", gp = gpar(fontsize = 13), hjust = 1, vjust = 0)
# title2_sphingolipids_bugs_0.2 <- textGrob("Mircobial pathways correlated with sphingolipids (n = 95)", rot = 90, gp = gpar(fontsize = 13), hjust = 0, vjust = 5)
# 
# # Combine the heatmap and titles using arrangeGrob and grid.arrange
# combined_plot_sphingolipids_bugs_0.2 <- arrangeGrob(
#   arrangeGrob(title2_sphingolipids_bugs_0.2, heatmap_sphingolipids_bugs_0.2_gtable, ncol = 2, widths = c(0.2, 4)),
#   title1_sphingolipids_bugs_0.2, ncol = 1, heights = c(4, 0.2)
# )
# 
# # Draw the combined plot
# grid.newpage()
# grid.draw(combined_plot_sphingolipids_bugs_0.2)
# 
# # Save the combined plot
# # ggsave("heatmap_sphingolipids_bugs_0.2.png", combined_plot_sphingolipids_bugs_0.2, width = 45, height = 25, dpi=300, limitsize = FALSE)




# #### Lysophospholipids 0.2: correlated pathways vs correlated species ####
# # read in all_associations.txt got from running HALLA
# r_diff_lysophospholipids_bugs_0.2 <- read.table('halla_lysophospholipids_0.2_bugs_paths/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
# 
# # Data cleaning
# test2_lysophospholipids_bugs_0.2 <- r_diff_lysophospholipids_bugs_0.2
# 
# 
# # Subset Relevant Columns
# a_lysophospholipids_bugs_0.2 <- test2_lysophospholipids_bugs_0.2[, colnames(test2_lysophospholipids_bugs_0.2) %in% c("X_features", "Y_features", "association")]
# a1_lysophospholipids_bugs_0.2 <- test2_lysophospholipids_bugs_0.2[, colnames(test2_lysophospholipids_bugs_0.2) %in% c("X_features", "Y_features", "q.values")]
# 
# # Reshape data
# b_lysophospholipids_bugs_0.2 <- a_lysophospholipids_bugs_0.2 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
#   as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
# q_values_matrix_lysophospholipids_bugs_0.2<- a1_lysophospholipids_bugs_0.2 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
#   as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
# 
# 
# asterisks_matrix_lysophospholipids_bugs_0.2 <- matrix("", nrow = nrow(q_values_matrix_lysophospholipids_bugs_0.2), ncol = ncol(q_values_matrix_lysophospholipids_bugs_0.2))
# asterisks_matrix_lysophospholipids_bugs_0.2[q_values_matrix_lysophospholipids_bugs_0.2 < 0.05] <- "***"
# asterisks_matrix_lysophospholipids_bugs_0.2[q_values_matrix_lysophospholipids_bugs_0.2 >= 0.05 & q_values_matrix_lysophospholipids_bugs_0.2 < 0.1] <- "**"
# asterisks_matrix_lysophospholipids_bugs_0.2[q_values_matrix_lysophospholipids_bugs_0.2 >= 0.1 & q_values_matrix_lysophospholipids_bugs_0.2 < 0.25] <- "*"
# 
# ### Create the heatmap
# # Create breaks with 0 in the middle
# breaks <- seq(-0.4, 0.4, length.out = 51)
# 
# # Ensure the breaks include 0
# mid <- which.min(abs(breaks - 0))
# 
# # Create a custom color palette
# my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))
# 
# heatmap_lysophospholipids_bugs_0.2 <- pheatmap(
#   b_lysophospholipids_bugs_0.2,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   display_numbers = asterisks_matrix_lysophospholipids_bugs_0.2,
#   fontsize_number = 10,
#   color = my_color_palette,
#   breaks = breaks,
#   main = "Heatmap of correlations between microbial species and pathways correlated with lysophospholipids",
#   annotation_legend = TRUE,
#   annotation_names_row = TRUE,
#   annotation_names_col = TRUE,
#   angle_col = 90
# )
# 
# # Extract the heatmap as a gtable object
# heatmap_lysophospholipids_bugs_0.2_gtable <- heatmap_lysophospholipids_bugs_0.2$gtable
# 
# # Create text grobs for titles
# title1_lysophospholipids_bugs_0.2 <- textGrob("Species correlated with lysophospholipids (n = 184)", gp = gpar(fontsize = 13), hjust = 1, vjust = 0)
# title2_lysophospholipids_bugs_0.2 <- textGrob("Mircobial pathways correlated with lysophospholipids (n = 98)", rot = 90, gp = gpar(fontsize = 13), hjust = 0, vjust = 5)
# 
# # Combine the heatmap and titles using arrangeGrob and grid.arrange
# combined_plot_lysophospholipids_bugs_0.2 <- arrangeGrob(
#   arrangeGrob(title2_lysophospholipids_bugs_0.2, heatmap_lysophospholipids_bugs_0.2_gtable, ncol = 2, widths = c(0.2, 4)),
#   title1_lysophospholipids_bugs_0.2, ncol = 1, heights = c(4, 0.2)
# )
# 
# # Draw the combined plot
# grid.newpage()
# grid.draw(combined_plot_lysophospholipids_bugs_0.2)
# 
# # Save the combined plot
# # ggsave("heatmap_lysophospholipids_bugs_paths_0.2.png", combined_plot_lysophospholipids_bugs_0.2, width = 45, height = 25, dpi=300, limitsize = FALSE)
# 
# 
# 
# ### Simplified version
# # List of species and pathways to present
# species_to_present <- c(
#   "Veillonella_dispar",
#   "Veillonella_parvula",
#   "Haemophilus_parainfluenzae",
#   "Methylobacterium_SGB15164",
#   "GGB33586_SGB53517",
#   "Clostridiales_bacterium_NSJ_32",
#   "Clostridiales_bacterium",
#   "Collinsella_SGB4121",
#   "Mediterraneibacter_butyricigenes",
#   "GGB9534_SGB14937",
#   "Blautia_glucerasea",
#   "GGB9534_SGB14937",
#   "Ruminococcus_gnavus",
#   "Akkermansia_muciniphila"
# )
# 
# pathways_to_present <- c(
#   "PHOSLIPSYN-PWY: superpathway of phospholipid biosynthesis I (bacteria)",
#   "PWY4FS-7: phosphatidylglycerol biosynthesis I (plastidic)",
#   "PWY4FS-8: phosphatidylglycerol biosynthesis II (non-plastidic)"
# )
# # Filter the matrix and annotations to only include the selected species
# b_lysophospholipids_bugs_0.2_filtered <- b_lysophospholipids_bugs_0.2[rownames(b_lysophospholipids_bugs_0.2)%in% pathways_to_present, colnames(b_lysophospholipids_bugs_0.2) %in% species_to_present]
# q_values_matrix_lysophospholipids_bugs_0.2_filtered <- q_values_matrix_lysophospholipids_bugs_0.2[rownames(b_lysophospholipids_bugs_0.2)%in% pathways_to_present, colnames(b_lysophospholipids_bugs_0.2) %in% species_to_present]
# 
# # Create the asterisks matrix for the filtered data
# asterisks_matrix_lysophospholipids_bugs_0.2_filtered <- matrix("", nrow = nrow(q_values_matrix_lysophospholipids_bugs_0.2_filtered), ncol = ncol(q_values_matrix_lysophospholipids_bugs_0.2_filtered))
# asterisks_matrix_lysophospholipids_bugs_0.2_filtered[q_values_matrix_lysophospholipids_bugs_0.2_filtered < 0.05] <- "***"
# asterisks_matrix_lysophospholipids_bugs_0.2_filtered[q_values_matrix_lysophospholipids_bugs_0.2_filtered >= 0.05 & q_values_matrix_lysophospholipids_bugs_0.2_filtered < 0.1] <- "**"
# asterisks_matrix_lysophospholipids_bugs_0.2_filtered[q_values_matrix_lysophospholipids_bugs_0.2_filtered >= 0.1 & q_values_matrix_lysophospholipids_bugs_0.2_filtered < 0.25] <- "*"
# 
# # Create the filtered heatmap
# heatmap_lysophospholipids_bugs_0.2_filtered <- pheatmap(b_lysophospholipids_bugs_0.2_filtered,
#                                                        cluster_rows = TRUE,
#                                                        cluster_cols = TRUE,
#                                                        display_numbers = asterisks_matrix_lysophospholipids_bugs_0.2_filtered,
#                                                        fontsize_number = 10,
#                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#                                                        main = "Heatmap of correlations between selected species and lysophospholipid metabolites",
#                                                        annotation_legend = TRUE,
#                                                        annotation_names_row = TRUE,
#                                                        annotation_names_col = TRUE,
#                                                        border_color = NA,
#                                                        angle_col = 90)
# 
# # Extract the filtered heatmap as a gtable object
# heatmap_lysophospholipids_bugs_0.2_filtered_gtable <- heatmap_lysophospholipids_bugs_0.2_filtered$gtable
# 
# # Create text grobs for titles
# title1_lysophospholipids_bugs_0.2_filtered <- textGrob("Selected microbial species (n = 6)", gp = gpar(fontsize = 13), hjust = 2, vjust = 0)
# title2_lysophospholipids_bugs_0.2_filtered <- textGrob("lysophospholipid metabolites (n = 12)", rot = 90, gp = gpar(fontsize = 13), hjust = 0.1, vjust = 0)
# 
# # Combine the filtered heatmap and titles using arrangeGrob and grid.arrange
# combined_plot_lysophospholipids_bugs_0.2_filtered <- arrangeGrob(
#   arrangeGrob(title2_lysophospholipids_bugs_0.2_filtered, heatmap_lysophospholipids_bugs_0.2_filtered_gtable, ncol = 2, widths = c(0.2, 4)),
#   title1_lysophospholipids_bugs_0.2_filtered, ncol = 1, heights = c(4, 0.2)
# )
# 
# # Draw the combined filtered plot
# grid.newpage()
# grid.draw(combined_plot_lysophospholipids_bugs_0.2_filtered)
# 
# # Save the filtered combined plot
# # ggsave("heatmap_selected_species_lysophospholipids_filtered.png", combined_plot_lysophospholipids_bugs_0.2_filtered, width = 11, height = 4.5, dpi=300, limitsize = FALSE)



# #### Lysoplasmalogens 0.2: correlated pathways vs correlated metabolites ####
# # read in all_associations.txt got from running HALLA
# r_diff_lysoplasmalogens_bugs_0.2 <- read.table('halla_lysoplasmalogens_0.2_bugs_paths/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
# 
# # Data cleaning
# test2_lysoplasmalogens_bugs_0.2 <- r_diff_lysoplasmalogens_bugs_0.2
# 
# 
# # Subset Relevant Columns
# a_lysoplasmalogens_bugs_0.2 <- test2_lysoplasmalogens_bugs_0.2[, colnames(test2_lysoplasmalogens_bugs_0.2) %in% c("X_features", "Y_features", "association")]
# a1_lysoplasmalogens_bugs_0.2 <- test2_lysoplasmalogens_bugs_0.2[, colnames(test2_lysoplasmalogens_bugs_0.2) %in% c("X_features", "Y_features", "q.values")]
# 
# # Reshape data
# b_lysoplasmalogens_bugs_0.2 <- a_lysoplasmalogens_bugs_0.2 %>% pivot_wider(names_from = "X_features", values_from = "association") %>%
#   as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
# q_values_matrix_lysoplasmalogens_bugs_0.2<- a1_lysoplasmalogens_bugs_0.2 %>% pivot_wider(names_from = "X_features", values_from = "q.values") %>%
#   as.data.frame() %>% column_to_rownames("Y_features") %>% as.matrix()
# 
# 
# asterisks_matrix_lysoplasmalogens_bugs_0.2 <- matrix("", nrow = nrow(q_values_matrix_lysoplasmalogens_bugs_0.2), ncol = ncol(q_values_matrix_lysoplasmalogens_bugs_0.2))
# asterisks_matrix_lysoplasmalogens_bugs_0.2[q_values_matrix_lysoplasmalogens_bugs_0.2 < 0.05] <- "***"
# asterisks_matrix_lysoplasmalogens_bugs_0.2[q_values_matrix_lysoplasmalogens_bugs_0.2 >= 0.05 & q_values_matrix_lysoplasmalogens_bugs_0.2 < 0.1] <- "**"
# asterisks_matrix_lysoplasmalogens_bugs_0.2[q_values_matrix_lysoplasmalogens_bugs_0.2 >= 0.1 & q_values_matrix_lysoplasmalogens_bugs_0.2 < 0.25] <- "*"
# 
# ### Create the heatmap
# # Create breaks with 0 in the middle
# breaks <- seq(-0.5, 0.5, length.out = 51)
# 
# # Ensure the breaks include 0
# mid <- which.min(abs(breaks - 0))
# 
# # Create a custom color palette
# my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))
# 
# heatmap_lysoplasmalogens_bugs_0.2 <- pheatmap(
#   b_lysoplasmalogens_bugs_0.2,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   display_numbers = asterisks_matrix_lysoplasmalogens_bugs_0.2,
#   fontsize_number = 10,
#   color = my_color_palette,
#   breaks = breaks,
#   main = "Heatmap of correlations between microbial species and pathways correlated with lysoplasmalogens",
#   annotation_legend = TRUE,
#   annotation_names_row = TRUE,
#   annotation_names_col = TRUE,
#   angle_col = 90
# )
# 
# # Extract the heatmap as a gtable object
# heatmap_lysoplasmalogens_bugs_0.2_gtable <- heatmap_lysoplasmalogens_bugs_0.2$gtable
# 
# # Create text grobs for titles
# title1_lysoplasmalogens_bugs_0.2 <- textGrob("Species correlated with lysoplasmalogens (n = 80)", gp = gpar(fontsize = 13), hjust = 1, vjust = -2)
# title2_lysoplasmalogens_bugs_0.2 <- textGrob("Mircobial pathways correlated with lysoplasmalogens (n = 46)", rot = 90, gp = gpar(fontsize = 13), hjust = 0, vjust = 3)
# 
# # Combine the heatmap and titles using arrangeGrob and grid.arrange
# combined_plot_lysoplasmalogens_bugs_0.2 <- arrangeGrob(
#   arrangeGrob(title2_lysoplasmalogens_bugs_0.2, heatmap_lysoplasmalogens_bugs_0.2_gtable, ncol = 2, widths = c(0.2, 4)),
#   title1_lysoplasmalogens_bugs_0.2, ncol = 1, heights = c(4, 0.2)
# )
# 
# # Draw the combined plot
# grid.newpage()
# grid.draw(combined_plot_lysoplasmalogens_bugs_0.2)
# 
# # Save the combined plot
# # ggsave("heatmap_lysoplasmalogens_bugs_paths_0.2.png", combined_plot_lysoplasmalogens_bugs_0.2, width = 25, height = 15, dpi=300, limitsize = FALSE)
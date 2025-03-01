rm(list=ls()) 

library(readr)
library(Maaslin2)
library(tidyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(vegan)
library(fossil)
library(ggpubr)
library(rstatix)
library(stringr)
library(readxl)
library(writexl)
library(coin)

main.dir<- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R"
setwd(main.dir)

####---------------  READIN  Input files ----------------------------####
## Metaphlan Species abundance table
all_taxo <- read.delim(file = file.path(main.dir,"all_taxo.tsv"),
                       sep = "\t", header = T, stringsAsFactors = F)
all_taxo <- all_taxo[grepl("s__", all_taxo$clade_name) & !grepl("\\|t__", all_taxo$clade_name), ]
all_taxo$clade_name <- str_extract(all_taxo$clade_name, "s__.+")

# Remove "_taxo" from each column name
current_colnames <- colnames(all_taxo)
new_colnames <- sub("_taxo", "", current_colnames)
colnames(all_taxo) <- new_colnames

# Simplify the clade_name
all_taxo$clade_name <- sub("s__", "", all_taxo$clade_name)

# Set clade_names to rownames
current_rownames <- rownames(all_taxo)
new_rownames <- all_taxo$clade_name
rownames(all_taxo) <- new_rownames
all_taxo <- all_taxo %>% select(2:1174)

all_taxo <- t(all_taxo) %>% as.data.frame()

# Count how many species in each row
number_of_species_per_sample <- data.frame(apply(all_taxo != 0, 1, sum))

# Remove patients excluded after QC
microbiome_list <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/microbiome_list.xlsx", col_names = FALSE)
microbiome_list <- microbiome_list %>%
  rename(id_microbiome = ...1)

all_taxo$id <- rownames(all_taxo)

all_taxo_hc <- all_taxo[grepl("G", all_taxo$id), ]
all_taxo_mc <- all_taxo[grepl("MC", all_taxo$id), ]

all_taxo_mc <- all_taxo_mc %>%
  filter(id %in% microbiome_list$id_microbiome)

all_taxo <- rbind(all_taxo_hc, all_taxo_mc)
all_taxo <- all_taxo %>% select(id, 1:2509)

## HuMAnn3 pathabundance 
all_path <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/all_path.xlsx")
all_path <- all_path %>% rename("path" = "...1")

# Only keep rows that contain "|"
rows_to_keep <- apply(all_path, 1, function(row) any(grepl("\\|", row)))
all_path <- all_path[rows_to_keep, ]

# Simplify column names
colnames(all_path) <- sub("_Abundance", "", colnames(all_path))
all_path <- t(all_path) %>% as.data.frame()
colnames(all_path) <- all_path[1, ]
all_path <- all_path[-1, ]

all_path$id <- rownames(all_path) 

# Remove patients excluded after QC
all_path_hc <- all_path[grepl("G", all_path$id), ]
all_path_mc <- all_path[grepl("MC", all_path$id), ]

all_path_mc <- all_path_mc %>%
  filter(id %in% microbiome_list$id_microbiome)

all_path <- rbind(all_path_hc, all_path_mc)

# Convert colnames from sample_ID to GIDER_id
df <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list_v2.xlsx")
all_path_hc <- all_path[grepl("G", all_path$id), ]
all_path_mc <- all_path[grepl("MC", all_path$id), ]

all_path_hc <- left_join(all_path_hc, df, by = c("id" = "sample_id")) # n = 732 samples
all_path_hc <- all_path_hc %>% filter(!is.na(GIDER_id)) # n = 717 samples, 15 samples that failed QC were removed
all_path_hc <- all_path_hc %>% select(GIDER_id, everything()) %>% select(-id)
all_path_hc <- all_path_hc %>% rename("id" = "GIDER_id")

# Only keep rows that are selected as HC after 1:3 matching
load("match_2024-05-13.RData") 
all_path_hc <- left_join(match3df_controls, all_path_hc, by = c("id" = "id"))
all_path_hc <- all_path_hc %>% select(-age, -sex)
all_path_hc <- as.data.frame(all_path_hc)

all_path <- rbind(all_path_hc, all_path_mc)
all_path <- all_path %>% t() %>% as.data.frame()

colnames(all_path) <- all_path[1, ]
all_path <- all_path[-1, ]
all_path$path <- rownames(all_path)
all_path <- all_path %>% select(path, everything())
rownames(all_path) <- seq_len(nrow(all_path))


####---- Find contributory species of specific pathways ----####
# GGB9534_SGB14937 <- all_path[grepl("GGB9534_SGB14937", all_path$path), ]
# GGB80011_SGB15265 <- all_path[grepl("GGB80011_SGB15265", all_path$path), ]
Bacteroides_stercoris <- all_path[grepl("Bacteroides_stercoris", all_path$path), ]
# Methylobacterium_SGB15164 <- all_path[grepl("Methylobacterium_SGB15164", all_path$path), ]
# GGB9524_SGB14923 <- all_path[grepl("GGB9524_SGB14923", all_path$path), ]
# Collinsella_SGB4121 <- all_path[grepl("Collinsella_SGB4121", all_path$path), ]
# Clostridiales_bacterium_Choco116 <-  all_path[grepl("Clostridiales_bacterium_Choco116", all_path$path), ]
Clostridiales_bacterium <- all_path[grepl("Clostridiales_bacterium", all_path$path), ]
# GGB33586_SGB53517 <- all_path[grepl("GGB33586_SGB53517", all_path$path), ]
# Blautia_glucerasea <- all_path[grepl("Blautia_glucerasea", all_path$path), ]
# Mediterraneibacter_butyricigenes <- all_path[grepl("Mediterraneibacter_butyricigenes", all_path$path), ]

# Clostridium_sp_AF20_17LB <- all_path[grepl("Clostridium_sp_AF20_17LB", all_path$path), ]
Haemophilus_parainfluenzae <- all_path[grepl("Haemophilus_parainfluenzae", all_path$path), ]
# GGB3612_SGB4882 <- all_path[grepl("GGB3612_SGB4882", all_path$path), ]
Veillonella_dispar <- all_path[grepl("Veillonella_dispar", all_path$path), ]
Clostridium_spiroforme <- all_path[grepl("Clostridium_spiroforme", all_path$path), ]
Intestinibacter_bartlettii <- all_path[grepl("Intestinibacter_bartlettii", all_path$path), ]
Veillonella_rogosae <- all_path[grepl("Veillonella_rogosae", all_path$path), ]
Veillonella_parvula <- all_path[grepl("Veillonella_parvula", all_path$path), ]


####---- Find contributory species of specific EC ----####
EC <- read.delim(file = "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/EC/ECs_relab_renamed.tsv",
                 sep = "\t", header = T, stringsAsFactors = F)
EC <- EC %>% rename("ECs" = "X..Gene.Family")

# Remove rows that contain "|"
row_to_keep <- apply(EC, 1, function(row) any(grepl("\\|", row)))
EC <- EC[row_to_keep, ]

# Simplify column names
colnames(EC) <- sub("_Abundance.RPKs", "", colnames(EC))
EC <- t(EC) %>% as.data.frame()
colnames(EC) <- EC[1, ]
EC <- EC[-1, ]

# Change the row name from "MC087" to "MC087B"
rownames(EC)[rownames(EC) == "MC087"] <- "MC087B"

EC$id <- rownames(EC)

# Remove patients excluded after QC
EC_hc <- EC[grepl("G", EC$id), ]
EC_mc <- EC[grepl("MC", EC$id), ]

EC_mc <- EC_mc %>%
  filter(id %in% microbiome_list$id_microbiome)

EC <- rbind(EC_hc, EC_mc)

EC$id <- NULL
EC <- t(EC) %>% as.data.frame()
EC$EC <- rownames(EC)

# Move 'EC' to the first column and reset row names
EC <- EC %>%
  select(EC, everything())

# Reset row names to sequential numbers
rownames(EC) <- seq_len(nrow(EC))


# GGB9534_SGB14937_EC <- EC[grepl("GGB9534_SGB14937", EC$EC), ]
# GGB80011_SGB15265_EC <- EC[grepl("GGB80011_SGB15265", EC$EC), ]
Bacteroides_stercoris_EC <- EC[grepl("Bacteroides_stercoris", EC$EC), ]
# Methylobacterium_SGB15164_EC <- EC[grepl("Methylobacterium_SGB15164", EC$EC), ]
# GGB9524_SGB14923_EC <- EC[grepl("GGB9524_SGB14923", EC$EC), ]
# Collinsella_SGB4121_EC <- EC[grepl("Collinsella_SGB4121", EC$EC), ]
# Clostridiales_bacterium_Choco116_EC <-  EC[grepl("Clostridiales_bacterium_Choco116", EC$EC), ]
Clostridiales_bacterium_EC <- EC[grepl("Clostridiales_bacterium", EC$EC), ]
# GGB33586_SGB53517_EC <- EC[grepl("GGB33586_SGB53517", EC$EC), ]
# Blautia_glucerasea_EC <- EC[grepl("Blautia_glucerasea", EC$EC), ]
# Mediterraneibacter_butyricigenes_EC <- EC[grepl("Mediterraneibacter_butyricigenes", EC$EC), ]

# Clostridium_sp_AF20_17LB_EC <- EC[grepl("Clostridium_sp_AF20_17LB", EC$EC), ]
Haemophilus_parainfluenzae_EC <- EC[grepl("Haemophilus_parainfluenzae", EC$EC), ]
# GGB3612_SGB4882_EC <- EC[grepl("GGB3612_SGB4882", EC$EC), ]
Veillonella_dispar_EC <- EC[grepl("Veillonella_dispar", EC$EC), ]
Clostridium_spiroforme_EC <- EC[grepl("Clostridium_spiroforme", EC$EC), ]
Intestinibacter_bartlettii_EC <- EC[grepl("Intestinibacter_bartlettii", EC$EC), ]
Veillonella_rogosae_EC <- EC[grepl("Veillonella_rogosae", EC$EC), ]
Veillonella_parvula_EC <- EC[grepl("Veillonella_parvula", EC$EC), ]




####---- Scatter plots for metabolites and EC abundance of specific species ----####
load("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/mc_metabolite_input_2024-04-25.RData") 
load("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/match_2024-05-13.RData")
load("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mc_input_2025-02-12.RData")

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
metabolite <- rbind(metabolite_1, metabolite_2)

altered_114 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/altered_met_114.xlsx") %>% select(4,5) %>% rename("pathway" = "pathway_HC")
all <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/met_cross_all/all_results.tsv")
all <- all %>% filter(value=="Control") %>% select(1,4,5,6,8,9) 
all$coef <- -(all$coef)
all <- all %>% mutate(t = coef/stderr)
all$feature<- substr(all$feature, 2,6)

pathway <- index[, c(3,6)]
pathway$COMP_ID <- as.character(pathway$COMP_ID)
chemical <- index[, c(3,11)]
chemical$COMP_ID <- as.character(chemical$COMP_ID)

all <- all %>% left_join(pathway, by = c("feature" = "COMP_ID"))
all <- all %>% left_join(chemical, by = c("feature" = "COMP_ID"))
all <- all %>% rename(
  'pathway' = 'SUB_PATHWAY',
  'metabolite' = 'CHEMICAL_NAME'
)

alpha_with_disease_status <- alpha_with_disease_status %>% 
  mutate(mc_all_label = case_when(
    mc_all == "1" ~ "MC",
    mc_all == "0" ~ "Chronic diarrhea",
    mc_all == "2" ~ "Controls without diarrhea")) 

alpha_with_disease_status$mc_all_label <- factor(alpha_with_disease_status$mc_all_label,
                                                 labels = c("MC", "Control without diarrhea", "Chronic diarrhea"),
                                                 levels = c("MC", "Controls without diarrhea", "Chronic diarrhea"))

## Lysophospholipids levels
Lysophospholipid <- subset(altered_114, grepl("Lysophospholipid", pathway, ignore.case = TRUE)) %>% left_join(select(all, feature, metabolite))
Lysophospholipid_features <- Lysophospholipid$feature
Lysophospholipid <- metabolite %>% select(all_of(Lysophospholipid_features))
Lysophospholipid <- Lysophospholipid %>% mutate(Lysophospholipid = rowSums(across(everything())))
Lysophospholipid$id <- rownames(Lysophospholipid)
Lysophospholipid <- Lysophospholipid%>%
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

# 1-linoleoyl-GPG (18:2)* & EC 3.1.3.27 in Veillonella_parvula
linoleoyl_GPG <- Lysophospholipid %>% select("54885", id) %>% rename("1-linoleoyl-GPG (18:2)*" = "54885")
EC_3_1_3_27 <- Veillonella_dispar_EC %>% t() %>% as.data.frame()

# Set the first row as column names
colnames(EC_3_1_3_27) <- EC_3_1_3_27[1, ]
EC_3_1_3_27 <- EC_3_1_3_27[-1, ]
EC_3_1_3_27 <- EC_3_1_3_27 %>% select("3.1.3.27: Phosphatidylglycerophosphatase|g__Veillonella.s__Veillonella_dispar")
EC_3_1_3_27$id <- rownames(EC_3_1_3_27)

linoleoyl_GPG_EC_3_1_3_27 <- linoleoyl_GPG %>%
  left_join(EC_3_1_3_27) %>%
  left_join(select(alpha_with_disease_status, ID, mc_all_label), by = c("id" = "ID")) %>%
  select(id, everything())

linoleoyl_GPG_EC_3_1_3_27$mc_all_label <- factor(linoleoyl_GPG_EC_3_1_3_27$mc_all_label, levels = c("MC", "Control without diarrhea", "Chronic diarrhea"),
                                                 labels = c("MC", "Control without diarrhea", "Chronic diarrhea"))
linoleoyl_GPG_EC_3_1_3_27[,2] <- as.numeric(linoleoyl_GPG_EC_3_1_3_27[,2])
linoleoyl_GPG_EC_3_1_3_27[,3] <- as.numeric(linoleoyl_GPG_EC_3_1_3_27[,3])

# Calculate Spearman correlation
x <- linoleoyl_GPG_EC_3_1_3_27$`1-linoleoyl-GPG (18:2)*`
y <- linoleoyl_GPG_EC_3_1_3_27$`3.1.3.27: Phosphatidylglycerophosphatase|g__Veillonella.s__Veillonella_dispar`
# Regular Spearman
spearman_result <- cor.test(
  linoleoyl_GPG_EC_3_1_3_27$`1-linoleoyl-GPG (18:2)*`, 
  linoleoyl_GPG_EC_3_1_3_27$`3.1.3.27: Phosphatidylglycerophosphatase|g__Veillonella.s__Veillonella_dispar`, 
  method = "spearman"
)
print(spearman_result)
# Spearman by permutation (consistent with HAllA)
spearman_perm <- spearman_test(x ~ y, distribution = approximate(nresample = 10000))
print(spearman_perm)


linoleoyl_GPG_EC_3_1_3_27_fig <- ggplot(linoleoyl_GPG_EC_3_1_3_27, aes(x = `1-linoleoyl-GPG (18:2)*`, y = `3.1.3.27: Phosphatidylglycerophosphatase|g__Veillonella.s__Veillonella_dispar`, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of 1-linoleoyl-GPG (18:2)*", y = "Relative abundance of EC3.1.3.27: Phosphatidylglycerophosphatase of Veillonella_parvula", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  theme_minimal()+
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 0.000025)) + 
  annotate("text", x = 30, y = 0.00002, label = "Spearman rho = 0.165, p-value = 0.002", size = 4, color = "black")

linoleoyl_GPG_EC_3_1_3_27_fig
# ggsave(filename = "linoleoyl_GPG_EC_3_1_3_27_fig.pdf",
#        plot = linoleoyl_GPG_EC_3_1_3_27_fig, units = "in", width=7, height=5, dpi = 1200)




# 1-stearoyl-GPE (18:0) & EC 4.1.1.65 in Haemophilus_parainfluenzae
stearoyl_GPE <- Lysophospholipid %>% select("42398", id) %>% rename("1-stearoyl-GPE (18:0)" = "42398")
EC_4_1_1_65 <- Haemophilus_parainfluenzae_EC %>% t() %>% as.data.frame()

# Set the first row as column names
colnames(EC_4_1_1_65) <- EC_4_1_1_65[1, ]
EC_4_1_1_65 <- EC_4_1_1_65[-1, ]
EC_4_1_1_65 <- EC_4_1_1_65 %>% select("4.1.1.65: Phosphatidylserine decarboxylase|g__Haemophilus.s__Haemophilus_parainfluenzae")
EC_4_1_1_65$id <- rownames(EC_4_1_1_65)

stearoyl_GPE_EC_4_1_1_65 <- stearoyl_GPE %>%
  left_join(EC_4_1_1_65) %>%
  left_join(select(alpha_with_disease_status, ID, mc_all_label), by = c("id" = "ID")) %>%
  select(id, everything())

stearoyl_GPE_EC_4_1_1_65$mc_all_label <- factor(stearoyl_GPE_EC_4_1_1_65$mc_all_label, levels = c("MC", "Control without diarrhea", "Chronic diarrhea"),
                                                 labels = c("MC", "Control without diarrhea", "Chronic diarrhea"))
stearoyl_GPE_EC_4_1_1_65[,2] <- as.numeric(stearoyl_GPE_EC_4_1_1_65[,2])
stearoyl_GPE_EC_4_1_1_65[,3] <- as.numeric(stearoyl_GPE_EC_4_1_1_65[,3])

# Calculate Spearman correlation
x <- stearoyl_GPE_EC_4_1_1_65$`1-stearoyl-GPE (18:0)`
y <- stearoyl_GPE_EC_4_1_1_65$`4.1.1.65: Phosphatidylserine decarboxylase|g__Haemophilus.s__Haemophilus_parainfluenzae`
# Regular Spearman
spearman_result <- cor.test(
  stearoyl_GPE_EC_4_1_1_65$`1-stearoyl-GPE (18:0)`, 
  stearoyl_GPE_EC_4_1_1_65$`4.1.1.65: Phosphatidylserine decarboxylase|g__Haemophilus.s__Haemophilus_parainfluenzae`, 
  method = "spearman"
)
print(spearman_result)
# Spearman by permutation (consistent with HAllA)
spearman_perm <- spearman_test(x ~ y, distribution = approximate(nresample = 10000))
print(spearman_perm)


stearoyl_GPE_EC_4_1_1_65_fig <- ggplot(stearoyl_GPE_EC_4_1_1_65, aes(x = `1-stearoyl-GPE (18:0)`, y = `4.1.1.65: Phosphatidylserine decarboxylase|g__Haemophilus.s__Haemophilus_parainfluenzae`, color = mc_all_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Relative abundance of 1-stearoyl-GPE (18:0)", y = "Relative abundance of EC4.1.1.65: Phosphatidylserine decarboxylase of Haemophilus_parainfluenzae", color = "Disease type") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  theme_minimal()+
  scale_color_manual(values = c("#900C3F", "#189AF9", "#F9D937")) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 0.000025)) + 
  annotate("text", x = 30, y = 0.00002, label = "Spearman rho = 0.120, p-value = 0.029", size = 4, color = "black")

stearoyl_GPE_EC_4_1_1_65_fig
# ggsave(filename = "stearoyl_GPE_EC_4_1_1_65_fig.pdf",
#        plot = stearoyl_GPE_EC_4_1_1_65_fig, units = "in", width=7, height=5, dpi = 1200)

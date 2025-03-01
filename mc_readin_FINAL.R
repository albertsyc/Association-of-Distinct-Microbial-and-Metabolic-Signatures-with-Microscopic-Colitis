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

# Remove rows that contain "|"
rows_to_remove <- apply(all_path, 1, function(row) any(grepl("\\|", row)))
all_path <- all_path[!rows_to_remove, ]

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
all_path <- all_path %>% select(id, 1:545)

# enzyme_ceramide_de_novo_synthesis <- all_path %>% select(contains("ceramide"))




# ## EC  
# EC <- read.delim(file = "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/EC/ECs_relab_renamed.tsv",
#                  sep = "\t", header = T, stringsAsFactors = F)
# EC <- EC %>% rename("ECs" = "X..Gene.Family")
# 
# # Remove rows that contain "|"
# row_to_remove <- apply(EC, 1, function(row) any(grepl("\\|", row)))
# EC <- EC[!row_to_remove, ]
# 
# # Simplify column names
# colnames(EC) <- sub("_Abundance.RPKs", "", colnames(EC))
# EC <- t(EC) %>% as.data.frame()
# colnames(EC) <- EC[1, ]
# EC <- EC[-1, ]
# 
# # Change the row name from "MC087" to "MC087B"
# rownames(EC)[rownames(EC) == "MC087"] <- "MC087B"
# 
# EC$id <- rownames(EC) 
# 
# # Remove patients excluded after QC
# EC_hc <- EC[grepl("G", EC$id), ]
# EC_mc <- EC[grepl("MC", EC$id), ]
# 
# EC_mc <- EC_mc %>%
#   filter(id %in% microbiome_list$id_microbiome)
# 
# EC <- rbind(EC_hc, EC_mc)
# EC <- EC %>% select(id, 1:2673)








####---------------  Use MaAsLin to filter the species with relative abundance thresholds of at least 0.01% and present in at least 10% of samples ----------------------------####
load("metadata_2024-05-13.RData") 
load("match_2024-05-13.RData") 
samples <- as.data.frame(rownames(all_taxo))
colnames(samples)[1] <- "ID"

# Subset rows where ID starts with "MC"
ID <- samples[grepl("^MC", samples$ID), ]
samples_mc <- as.data.frame(ID)
samples_mc$id <- substr(samples_mc$ID, 1, 5)

with_disease_status <- samples_mc %>% 
  left_join(select(baseline_microbiome, id, active_1, active_2, mc_all), by = "id")

with_disease_status$stool_time<- substr(with_disease_status$ID, 6,7)
# Treating patients with only one sample (no "A" or "B") as "A"
with_disease_status$stool_time[with_disease_status$stool_time==""] <- "A"
with_disease_status <- with_disease_status %>%
  mutate(disease_status = case_when(
    mc_all == 1 & stool_time == "A" & active_1 == "Yes" ~ "active",
    mc_all == 1 & stool_time == "B" & active_2 == "Yes" ~ "active",
    mc_all == 1 & stool_time == "A" & active_1 == "No" ~ "remission",
    mc_all == 1 & stool_time == "B" & active_2 == "No" ~ "remission",
    mc_all == 1 & id == "MC508" ~ "active",
    mc_all == 0 ~ "active",
    TRUE ~ NA_character_
  ))

with_disease_status <- with_disease_status %>%
  left_join(select(baseline_microbiome, id, mc_or_diarrhea), by = "id")

with_disease_status <- with_disease_status %>% 
  mutate(mc_all_label = case_when(
    mc_all == 1 ~ "MC",
    mc_all == 0 ~ "Chronic_diarrhea",
    TRUE ~ NA
  ))

rownames(with_disease_status) <- with_disease_status$ID

## Cleaning HC data
# Convert colnames from sample_ID to GIDER_id
df <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list_v2.xlsx")
all_taxo_hc <- all_taxo[grepl("G", all_taxo$id), ]
all_taxo_mc <- all_taxo[grepl("MC", all_taxo$id), ]

all_taxo_hc <- left_join(all_taxo_hc, df, by = c("id" = "sample_id")) # n = 732 samples
all_taxo_hc <- all_taxo_hc %>% filter(!is.na(GIDER_id)) # n = 717 samples, 15 samples that failed QC were removed
all_taxo_hc <- all_taxo_hc %>% select(GIDER_id, 1:2510) %>% select(-id)
all_taxo_hc <- all_taxo_hc %>% rename("id" = "GIDER_id")

# Only keep rows that are selected as HC after 1:3 matching
load("match_2024-05-13.RData") 
match3df_controls <- as.data.frame(match3df_controls)
rownames(match3df_controls) <- match3df_controls$id

all_taxo_hc <- left_join(match3df_controls, all_taxo_hc, by = c("id" = "id"))
all_taxo_hc <- all_taxo_hc %>% select(-age, -sex)
all_taxo_hc <- as.data.frame(all_taxo_hc)

all_taxo <- rbind(all_taxo_hc, all_taxo_mc)
rownames(all_taxo) <- all_taxo$id

# Double check if all id in with_disease_status match to id in baseline_microbiome, and vice versa
all(with_disease_status$id %in% baseline_microbiome$id) # TRUE
all(baseline_microbiome$id %in% with_disease_status$id) # TRUE

# Create a temporary metadata form maaslin filter
with_disease_status <- with_disease_status %>% left_join(select(baseline_microbiome, id, age, sex), by = c("id"="id"))
df1 <- with_disease_status %>% select(ID, age, sex)
rownames(df1) <- df1$ID
df1 <- df1 %>% rename("id" = "ID")
meta <- rbind(match3df_controls, df1)

# maaslin_filter <- Maaslin2(input_data = all_taxo,
#                            input_metadata = meta,
#                            output = "maaslin_filter",
#                            min_abundance = 0.0001,
#                            min_prevalence = 0.1,
#                            normalization = "NONE", #already relative abundance
#                            transform = "LOG",
#                            analysis_method = "LM",
#                            max_significance = 0.25,
#                            fixed_effects = c("age"))

for_filter <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_filter/all_results.tsv")

# Extract the list of species from the rows of for_filter
species_to_select <- for_filter$feature # 467 species are > 0.01% in relative abundance and present in at least 10% of samples

# Select columns from all_taxo based on the species_to_select list
all_taxo <- all_taxo[, species_to_select, drop = FALSE]
all_taxo <- as.data.frame(all_taxo)
all_taxo$id <- rownames(all_taxo)

all_taxo <- all_taxo %>% select(id,1:467)




####---------------  Use MaAsLin to filter the pathways with relative abundance thresholds of at least 0.1% and present in at least 10% of samples ----------------------------####
## Cleaning HC data
# Convert colnames from sample_ID to GIDER_id
df <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list_v2.xlsx")
all_path_hc <- all_path[grepl("G", all_path$id), ]
all_path_mc <- all_path[grepl("MC", all_path$id), ]

all_path_hc <- left_join(all_path_hc, df, by = c("id" = "sample_id")) # n = 732 samples
all_path_hc <- all_path_hc %>% filter(!is.na(GIDER_id)) # n = 717 samples, 15 samples that failed QC were removed
all_path_hc <- all_path_hc %>% select(GIDER_id, 1:546) %>% select(-id)
all_path_hc <- all_path_hc %>% rename("id" = "GIDER_id")

# Only keep rows that are selected as HC after 1:3 matching
all_path_hc <- left_join(match3df_controls, all_path_hc, by = c("id" = "id"))
all_path_hc <- all_path_hc %>% select(-age, -sex)
all_path_hc <- as.data.frame(all_path_hc)

all_path <- rbind(all_path_hc, all_path_mc)
rownames(all_path) <- all_path$id

# Convert all columns to numeric
id <- all_path$id
all_path$id <- NULL
all_path <- apply(all_path, 2, function(x) as.numeric(as.character(x)))
rownames(all_path) <- id
all_path <- as.data.frame(all_path)

all_path_full <- all_path

# Change col_names to path_1~path_545
path_col_names <- colnames(all_path)

num_cols <- ncol(all_path)
new_col_names <- paste("path", 1:num_cols, sep = "_")
colnames(all_path) <- new_col_names


# maaslin_filter_path <- Maaslin2(input_data = all_path,
#                                 input_metadata = meta,
#                                 output = "maaslin_filter_path",
#                                 min_abundance = 0.001,
#                                 min_prevalence = 0.1,
#                                 normalization = "NONE", #already relative abundance
#                                 transform = "LOG",
#                                 analysis_method = "LM",
#                                 max_significance = 0.25,
#                                 fixed_effects = c("age"))

for_filter_path <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_filter_path/all_results.tsv")

### Extract the list of pathways from the rows of for_filter_path
path_to_select <- for_filter_path$feature # 250 pathways are > 0.1% in relative abundance and present in at least 10% of samples


# Select columns from all_path based on the path_to_select list
all_path <- all_path[, path_to_select, drop = FALSE]
all_path <- as.data.frame(all_path)

# Change the colnames back to the original path_col_names
colnames(all_path) <- c(path_col_names[match(colnames(all_path), new_col_names)])
all_path$id <- rownames(all_path)

all_path <- all_path %>% select(id,1:250)






# ####---------------  Use MaAsLin to filter the ECs with relative abundance thresholds of at least 0.1% and present in at least 10% of samples ----------------------------####
# ## Cleaning HC data
# # Convert colnames from sample_ID to GIDER_id
# df <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list_v2.xlsx")
# EC_hc <- EC[grepl("G", EC$id), ]
# EC_mc <- EC[grepl("MC", EC$id), ]
# 
# EC_hc <- left_join(EC_hc, df, by = c("id" = "sample_id")) # n = 732 samples
# EC_hc <- EC_hc %>% filter(!is.na(GIDER_id)) # n = 717 samples, 15 samples that failed QC were removed
# EC_hc <- EC_hc %>% select(GIDER_id, 1:2674) %>% select(-id)
# EC_hc <- EC_hc %>% rename("id" = "GIDER_id")
# 
# # Only keep rows that are selected as HC after 1:3 matching
# EC_hc <- left_join(match3df_controls, EC_hc, by = c("id" = "id"))
# EC_hc <- EC_hc %>% select(-age, -sex)
# EC_hc <- as.data.frame(EC_hc)
# 
# EC <- rbind(EC_hc, EC_mc)
# rownames(EC) <- EC$id
# 
# # Convert all columns to numeric
# id <- EC$id
# EC$id <- NULL
# EC <- apply(EC, 2, function(x) as.numeric(as.character(x)))
# rownames(EC) <- id
# EC <- as.data.frame(EC)
# 
# EC_full <- EC
# 
# # Change col_names to EC_1~EC_2673
# path_col_names_EC <- colnames(EC)
# 
# num_cols_EC <- ncol(EC)
# new_col_names_EC <- paste("EC", 1:num_cols_EC, sep = "_")
# colnames(EC) <- new_col_names_EC
# 
# 
# # maaslin_filter_EC <- Maaslin2(input_data = EC,
# #                               input_metadata = meta,
# #                               output = "maaslin_filter_EC",
# #                               min_abundance = 0.001,
# #                               min_prevalence = 0.1,
# #                               normalization = "NONE", #already relative abundance
# #                               transform = "LOG",
# #                               analysis_method = "LM",
# #                               max_significance = 0.25,
# #                               fixed_effects = c("age"))
# 
# for_filter_EC <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_filter_EC/all_results.tsv")
# 
# ### Extract the list of ECs from the rows of for_filter_EC
# EC_to_select <- for_filter_EC$feature # 472 ECs are > 0.1% in relative abundance and present in at least 10% of samples
# 
# 
# # Select columns from EC based on the EC_to_select list
# EC <- EC[, EC_to_select, drop = FALSE]
# EC <- as.data.frame(EC)
# 
# # Change the colnames back to the original path_col_names_EC
# colnames(EC) <- c(path_col_names_EC[match(colnames(EC), new_col_names_EC)])
# EC$id <- rownames(EC)
# 
# EC <- EC %>% select(id,1:472)







####---------------  CALCULATE DIVERSITY  ----------------#####
# Use VEGAN package
# ALPHA DIVERSITY
subset.all_taxo <- all_taxo[, !colnames(all_taxo) %in% c("id")]
alpha.shannon <- vegan::diversity(subset.all_taxo, index = "shannon", base=exp(1), MARGIN=1)
hist(alpha.shannon)

# BETA DIVERSITY
bray <- vegdist(subset.all_taxo, "bray")  # pairwise dissmilarity: 0=same; 1=maximally dissimilar

# Use fossil package -- Chao1 richness measure
alpha.chao1 <- apply(subset.all_taxo, 1, function(x) chao1(x, taxa.row = T))

alpha <- data.frame(ID = names(alpha.shannon), 
                    alpha.shannon = as.numeric(alpha.shannon), 
                    alpha.chao1 = as.numeric(alpha.chao1))

# Create id from ID
alpha_hc <- alpha[grepl("G", alpha$ID), ]
alpha_mc <- alpha[grepl("MC", alpha$ID), ]

alpha_mc$id <- substr(alpha_mc$ID, 1,5)
alpha_hc$id <- alpha_hc$ID
alpha <- rbind(alpha_hc, alpha_mc)

# Create alpha_with_disease_status
load("metadata_2024-05-13.RData") 
alpha <- alpha %>%
  left_join(select(baseline_microbiome, id, active_1, active_2, mc_all, mc_or_diarrhea), by = "id")

alpha$mc_all[is.na(alpha$mc_all)] <- 2

alpha_hc <- alpha[grepl("G", alpha$ID), ]
alpha_mc <- alpha[grepl("MC", alpha$ID), ]

alpha_mc$stool_time<- substr(alpha_mc$ID, 6,7)
# Treating patients with only one sample (no "A" or "B") as "A"
alpha_mc$stool_time[alpha_mc$stool_time==""] <- "A"

alpha_hc$stool_time <- "A"
alpha <- rbind(alpha_hc, alpha_mc)


alpha <- alpha %>%
  mutate(disease_status = case_when(
    mc_all == 1 & stool_time == "A" & active_1 == "Yes" ~ "active",
    mc_all == 1 & stool_time == "B" & active_2 == "Yes" ~ "active",
    mc_all == 1 & stool_time == "A" & active_1 == "No" ~ "remission",
    mc_all == 1 & stool_time == "B" & active_2 == "No" ~ "remission",
    mc_all == 1 & id == "MC508" ~ "active",
    mc_all == 0 ~ "active",
    TRUE ~ NA_character_
  ))

alpha$disease_status[substr(alpha$ID, 1, 1) == "G"] <- "remission"
alpha$active_1[substr(alpha$ID, 1, 1) == "G"] <- "No"

alpha_with_disease_status <- alpha %>% filter(!is.na(disease_status))
alpha_with_disease_status <- alpha_with_disease_status %>% left_join(select(baseline_microbiome, id, age, sex, bmi), by = c("id" = "id"))
alpha_with_disease_status <- alpha_with_disease_status %>% left_join(select(matched_hc, id, age, sex, bmi), by = c("ID" = "id"))
alpha_with_disease_status$age.x <- ifelse(is.na(alpha_with_disease_status$age.x), alpha_with_disease_status$age.y, alpha_with_disease_status$age.x)
alpha_with_disease_status$sex.x <- ifelse(is.na(alpha_with_disease_status$sex.x), alpha_with_disease_status$sex.y, alpha_with_disease_status$sex.x)
alpha_with_disease_status$bmi.x <- ifelse(is.na(alpha_with_disease_status$bmi.x), alpha_with_disease_status$bmi.y, alpha_with_disease_status$bmi.x)
alpha_with_disease_status <- alpha_with_disease_status %>% select(-age.y, -bmi.y, -sex.y) %>% rename("age" = "age.x",
                                                                                                     "sex" = "sex.x",
                                                                                                     "bmi" = "bmi.x")
alpha_with_disease_status$mc_or_diarrhea[substr(alpha$ID, 1, 1) == "G"] <- "Healthy"


####---------------  SAVE to RDATA ---------------##########
# save(alpha_with_disease_status, all_taxo, all_path, all_path_full,
#      file=paste0("mc_input_",Sys.Date(),".RData"))



# ####--------------- (X) HC QC ----------------#####
# taxo_hc <- read.delim(file = file.path(main.dir,"hc_taxo.tsv"),
#                       sep = "\t", header = T, stringsAsFactors = F)
# taxo_hc <- taxo_hc[grepl("s__", taxo_hc$clade_name) & !grepl("\\|t__", taxo_hc$clade_name), ]
# taxo_hc$clade_name <- str_extract(taxo_hc$clade_name, "s__.+")
# 
# # Remove "_taxo" from each column name
# current_colnames <- colnames(taxo_hc)
# new_colnames <- sub("_taxo", "", current_colnames)
# colnames(taxo_hc) <- new_colnames
# 
# # Simplify the clade_name
# taxo_hc$clade_name <- sub("s__", "", taxo_hc$clade_name) 
# 
# # Set clade_names to rownames
# current_rownames <- rownames(taxo_hc)
# new_rownames <- taxo_hc$clade_name
# rownames(taxo_hc) <- new_rownames
# taxo_hc <- taxo_hc %>% select(2:733)
# 
# taxo_hc <- t(taxo_hc) %>% as.data.frame()
# 
# # Count how many species in each row
# number_of_species_per_sample <- data.frame(apply(taxo_hc != 0, 1, sum))
# # write_xlsx(number_of_species_per_sample, "number_of_species_per_sample.xlsx")
# 
# # HC: Create a QC without duplicate spreadsheet 
# qc <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/microbiome_MC/QC_hc.xlsx")
# qc <- qc %>%
#   rename("id" = "Sample Name",
#          "reads" = "# of reads post-trimming")
# 
# qc_no_duplicate <- qc %>%
#   arrange(id, desc(reads)) %>%
#   distinct(id, .keep_all = TRUE)
# # write_xlsx(qc_no_duplicate, "qc_no_duplicate.xlsx")
rm(list=ls())
library(readxl)
library(dplyr)
library(tidyr)
library(forcats)
library(tableone)
library(MatchIt)
library(openxlsx)
library(readr)


df <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list/GIDER_master_list.xlsx")
qc_717 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/HC_qced_list_717.xlsx", col_names = FALSE)
qc_717 <- qc_717 %>% rename("sample_id" = "...1")
df <- df %>% slice(460:741)
df <- left_join(qc_717, df, by = "sample_id")
df$GIDER_id <- ifelse(is.na(df$GIDER_id), df$sample_id, df$GIDER_id)
# write_xlsx(df, "GIDER_master_list_v2.xlsx")

df1 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list/GIDER_demographics_corrected.xlsx")

####---------------  Cleaning  ----------------------#####
hc <- left_join(df, df1, by = "GIDER_id")
hc <- hc %>% rename(id = GIDER_id) %>% select(id,age,sex)

# active MC microbiome samples to match
load("metadata_2024-05-13.RData") 
load("mc_input_2024-07-08.RData") 
to_match_micro_sample <- alpha_with_disease_status %>% filter(mc_all==1 & disease_status=="active") %>% select(1:4)
to_match_micro_id <- to_match_micro_sample %>% group_by(id) %>% distinct(id) # active MC patients with microbiome data to match: 131 patients; 143 samples
to_match_micro_id <- to_match_micro_id %>% left_join(select(baseline_microbiome, id, age, sex), by = c("id" = "id"))


# active MC metabolomics samples to match
load("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/mc_metabolite_input_2024-04-25.RData") 
remission_mc_metabolite_sample <- metadata %>% filter(group=="MC" & symptoms == "no_diarrhea") 
remission_mc_metabolite_id <- remission_mc_metabolite_sample %>% distinct(ID) # Remission MC with metabolomics: 82 patients; 96 samples

to_match_met_sample <- metadata %>% filter(group=="MC"& symptoms == "active_diarrhea") %>% select(1:4)
to_match_met_id <- to_match_met_sample %>% group_by(ID) %>% distinct(ID) # active MC patients with metabolite data to match: 106 patients; 118 samples
to_match_met_sample$id <- gsub("_1$", "A", to_match_met_sample$id)
to_match_met_sample$id <- gsub("_2$", "B", to_match_met_sample$id)

# active MC patients with both microbiome and metabolomics to match
both_id <- inner_join(to_match_micro_id, to_match_met_id, by = c("id" = "ID")) 
both_sample <- both_id %>% left_join(select(alpha_with_disease_status, id, ID), by = c("id" = "id")) # active MC patients with both microbiome and metabolomics: 92 patients; 143 samples

# active MC patients with either microbiome or metabolomics to match
either_id <- full_join(to_match_micro_id, to_match_met_id, by = c("id" = "ID")) # active MC patients with either microbiome or metabolomics: 145 patients

####-------- Since microbiome is the main part, we match controls to MC with microbiome data: 131 patients --------####
# Selected HC that passed QC
selected_QC <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/qc_no_duplicate.xlsx")
selected_QC <- selected_QC %>%
  rename("id" = "Sample Name",
         "reads" = "# of reads post-trimming",
         "qced_reads" = "# of reads post-trimming (excluding <7 millions)",
         "species" = "# of species...4" ,
         "qced_species" = "# of species...5") %>% select(1:5) %>% filter(!is.na(qced_reads)) %>% slice(1:717)
selected_QC <- left_join(selected_QC, df, by = c("id" = "sample_id"))
selected_QC <- selected_QC %>% filter(!is.na(GIDER_id)) %>% select(6,2:5)
selected_QC <- selected_QC %>% rename(id = GIDER_id)

hc <- selected_QC %>% left_join(hc, by = "id")

raw <- rbind(to_match_micro_id, hc)
raw$disease <- ifelse(grepl("^G", raw$id), 0, 
                      ifelse(grepl("^MC", raw$id), 1, NA))

match1 <- matchit(disease ~ age + sex, data = raw, exact = c("sex"), distance = "glm", replace = FALSE, ratio = 1)
match1
summary(match1)

match2 <- matchit(disease ~ age + sex, data = raw, exact = c("sex"), distance = "glm", replace = FALSE, ratio = 2)
match2
summary(match2)

match3 <- matchit(disease ~ age + sex, data = raw, exact = c("sex"), distance = "glm", replace = FALSE, ratio = 3)
match3
summary(match3)

match4 <- matchit(disease ~ age + sex, data = raw, exact = c("sex"), distance = "glm", replace = FALSE, ratio = 4)
match4
summary(match4) # Not all cases will get 4 matches.

# Use match3: 1 case matched to 3 controls
match3df <- match.data(match3)
match3df <- match3df %>% arrange(match3df$subclass)

match3df_controls <- match3df %>% filter(disease==0) %>% select(1:3)


## Create an excel of selected HC to calculate read depth and base pairs
QC_hc <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/microbiome_MC/QC_hc.xlsx") %>% select(1,2)
master_list <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list_v2.xlsx") 
QC_hc <- QC_hc %>% distinct(`Sample Name`, .keep_all = TRUE)
QC_hc <- left_join(QC_hc, master_list, by = c("Sample Name" = "sample_id")) %>% select(3,2) %>% rename("id" = "GIDER_id")
QC_hc <- inner_join(QC_hc, match3df_controls) %>% select(1,2)

QC_hc_bp_1 <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/microbiome_MC/total_bases_summary_batch1.tsv")
QC_hc_bp_2 <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/microbiome_MC/total_bases_summary_batch2.tsv")
QC_hc_bp_3 <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/microbiome_MC/total_bases_summary_batch3.tsv")
QC_hc_bp <- rbind(QC_hc_bp_1, QC_hc_bp_2, QC_hc_bp_3)

QC_hc_bp_clean <- QC_hc_bp %>%
  # Remove "_fastqc" from "Sample Name"
  mutate(`Sample Name` = gsub("_fastqc", "", `Sample Name`)) %>%
  # Convert "Mbp" to numeric (divide by 1000) and remove "Gbp"
  mutate(`Total bases` = ifelse(grepl("Mbp", `Total bases`),
                                as.numeric(gsub(" Mbp", "", `Total bases`)) / 1000,
                                as.numeric(gsub(" Gbp", "", `Total bases`))))

QC_hc_bp_clean <- QC_hc_bp_clean %>% distinct(`Sample Name`, .keep_all = TRUE)
QC_hc_bp_clean <- left_join(QC_hc_bp_clean, master_list, by = c("Sample Name" = "sample_id")) %>% rename("id" = "GIDER_id")
QC_hc_bp_clean <- inner_join(QC_hc_bp_clean, match3df_controls) %>% select(2,3)
QC_hc <- left_join(QC_hc, QC_hc_bp_clean)

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/microbiome_MC/QC_hc_new.xlsx"
# write.xlsx(QC_hc, file_path, rowNames = TRUE)


# Check how many selected controls have both microbiome and metabolomics data
controls_mgx_mbx <- match3df_controls %>% inner_join(select(metadata, id), by = c("id" = "id")) # Selected controls with both microbiome and metabolomics data: 105 participants
controls_met <- metadata %>% filter(group == "Control") # Controls with metabolomics data: 273 participants

# Check how many CD have both microbiome and metabolomics data
diarrhea_mgx_mbx_1 <- alpha_with_disease_status %>% filter(mc_all == 0) # Chronic diarrhea with microbiome data: 159 patients
diarrhea_mgx_mbx_2 <- metadata %>% filter(group == "Diarrhea") # Chronic diarrhea with metabolomics data: 163 patients
diarrhea_mgx_mbx <- diarrhea_mgx_mbx_1 %>% inner_join(select(diarrhea_mgx_mbx_2, id, race), by = c("id" = "id")) # Chronic diarrhea with both microbiome and metabolomics data: 138 patients

# Check how many HC have both microbiome and metabolomics data
hc_mgx_mbx_1 <- alpha_with_disease_status %>% filter(mc_all == 2) # HC with microbiome data: 393 patients
hc_mgx_mbx_2 <- metadata %>% filter(group == "Control") # HC with metabolomics data: 273 patients
hc_mgx_mbx <- hc_mgx_mbx_1 %>% inner_join(select(hc_mgx_mbx_2, id, race), by = c("id" = "id")) # HC with both microbiome and metabolomics data: 105 patients
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_selected_hc.xlsx"
# write.xlsx(hc_mgx_mbx, file_path, rowNames = TRUE)



####-------- Table 1 for patients with both microbiome and metabolomics --------####
## MC vs HC
mc_hc_1 <- both_id %>% left_join(select(baseline_microbiome, id, bmi, mc_all, active_1, active_2, hispanic, race, smoke, anti_12))
mc_hc_1 %>% filter(active_1=="Yes" | active_2=="Yes") #n=92, all MC has at least 1 active sample
mc_hc_2 <- hc_mgx_mbx %>% select(ID, age, sex, bmi, mc_all, race)
mc_hc_2 <- mc_hc_2 %>% rename("id" = "ID")
  
mc_hc <- rbind(mc_hc_1,mc_hc_2)

table_mc_hc_both <- CreateTableOne(vars = c("age", "sex", "bmi", "hispanic", "race", "smoke", "anti_12"), 
                                                  strata = "mc_all", 
                                                  data = mc_hc)
print(table_mc_hc_both, showAllLevels = TRUE, noSpaces = TRUE)



## MC vs CD
mc_cd_1 <- both_id %>% left_join(select(baseline_microbiome, id, bmi, mc_all, hispanic, race, smoke, anti_12))
mc_cd_2 <- diarrhea_mgx_mbx %>% select(ID, age, sex, bmi, mc_all, race, disease_status)
mc_cd_2 %>% filter(disease_status=="active") #n=138, all CD are active
mc_cd_2 <- mc_cd_2 %>% rename("id" = "ID")

mc_cd <- rbind(mc_cd_1, mc_cd_2)

mc_cd$mc_all <- factor(mc_cd$mc_all, levels = c("1", "0"))

table_mc_cd_both <- CreateTableOne(vars = c("age", "sex", "bmi", "hispanic", "race", "smoke", "anti_12"), 
                                   strata = "mc_all", 
                                   data = mc_cd)
print(table_mc_cd_both, showAllLevels = TRUE, noSpaces = TRUE)






####---------------  SAVE to RDATA ---------------##########
# setwd("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R")
# save(match3df, match3df_controls, mc_hc, mc_cd, 
#      file=paste0("match_",Sys.Date(),".RData"))

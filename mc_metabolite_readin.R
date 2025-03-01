rm(list=ls())
library(readxl)
library(dplyr)
library(tidyr)
library(forcats)

metabolite <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/raw_metabolite.xlsx")
log_metabolite <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/log_metabolite.xlsx")
metadata_df <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metadata_metabolite.xlsx")

# Clean metadata: raw metabolite ----------------------------------------------------------
all(metadata_df$CLIENT_IDENTIFIER == metadata_df$CLIENT_SAMPLE_ID)
metadata <- metadata_df %>% select(1,2,7,15,17,22,RACE_ETHNICITY,SUBJECT_OR_ANIMAL_ID)
full_metabolite <- left_join(metadata, metabolite, by = "PARENT_SAMPLE_NAME")
full_metabolite <- full_metabolite %>% rename(
  'id' = 'CLIENT_IDENTIFIER',
  'ID' = 'SUBJECT_OR_ANIMAL_ID',
  'age' = 'AGE',
  'sex' = 'GENDER',
  'disease' = 'GROUP_ID2',
  'symptoms' = 'MICROSCOPIC_COLITIS_ST_1',
  'race' = 'RACE_ETHNICITY'
) %>% select(2,8,3:7,9:1783)

full_metabolite <- full_metabolite %>%
  mutate(
    group = case_when(
      disease == "Control" ~ "Control",
      disease == "Diarrhea" ~ "Diarrhea",
      disease == "Active_NS" ~ "MC",
      disease == "Remission_NS" ~ "MC",
      disease == "Active_CC" ~ "MC",
      disease == "Remission_CC" ~ "MC",
      disease == "Active_LC" ~ "MC",
      disease == "Remission_LC" ~ "MC"
    )
  )


full_metabolite$disease <- fct_recode(full_metabolite$disease,
                                    "NS" = "Active_NS",
                                    "NS" = "Remission_NS",
                                    "CC" = "Active_CC",
                                    "CC" = "Remission_CC",
                                    "LC" = "Active_LC",
                                    "LC" = "Remission_LC")

# Replace NA values based on the value of disease
full_metabolite$symptoms[is.na(full_metabolite$symptoms) & full_metabolite$disease == "Control"] <- "no_diarrhea"
full_metabolite$symptoms[is.na(full_metabolite$symptoms) & full_metabolite$disease == "Diarrhea"] <- "active_diarrhea"
# Replace existing values
full_metabolite$symptoms[full_metabolite$symptoms == "Active"] <- "active_diarrhea"
full_metabolite$symptoms[full_metabolite$symptoms == "Remission"] <- "no_diarrhea"

full_metabolite <- full_metabolite %>% select(1:4,1783,5:1782)

metadata <- full_metabolite %>% select(1:8)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$id

metabolite <- full_metabolite %>% select(1,2,5:7,9:1783)
metabolite <- as.data.frame(metabolite)
rownames(metabolite) <- metabolite$id

# Clean metadata: log metabolite ----------------------------------------------------------
metadata1 <- metadata_df %>% select(1,2,7,15,17,22,RACE_ETHNICITY,SUBJECT_OR_ANIMAL_ID)
log_full_metabolite <- left_join(metadata1, log_metabolite, by = "PARENT_SAMPLE_NAME")
log_full_metabolite <- log_full_metabolite %>% rename(
  'id' = 'CLIENT_IDENTIFIER',
  'ID' = 'SUBJECT_OR_ANIMAL_ID',
  'age' = 'AGE',
  'sex' = 'GENDER',
  'disease' = 'GROUP_ID2',
  'symptoms' = 'MICROSCOPIC_COLITIS_ST_1',
  'race' = 'RACE_ETHNICITY'
) %>% select(2,8,3:7,9:1783)

log_full_metabolite <- log_full_metabolite %>%
  mutate(
    group = case_when(
      disease == "Control" ~ "Control",
      disease == "Diarrhea" ~ "Diarrhea",
      disease == "Active_NS" ~ "MC",
      disease == "Remission_NS" ~ "MC",
      disease == "Active_CC" ~ "MC",
      disease == "Remission_CC" ~ "MC",
      disease == "Active_LC" ~ "MC",
      disease == "Remission_LC" ~ "MC"
    )
  )


log_full_metabolite$disease <- fct_recode(log_full_metabolite$disease,
                                      "NS" = "Active_NS",
                                      "NS" = "Remission_NS",
                                      "CC" = "Active_CC",
                                      "CC" = "Remission_CC",
                                      "LC" = "Active_LC",
                                      "LC" = "Remission_LC")

# Replace NA values based on the value of disease
log_full_metabolite$symptoms[is.na(log_full_metabolite$symptoms) & log_full_metabolite$disease == "Control"] <- "no_diarrhea"
log_full_metabolite$symptoms[is.na(log_full_metabolite$symptoms) & log_full_metabolite$disease == "Diarrhea"] <- "active_diarrhea"
# Replace existing values
log_full_metabolite$symptoms[log_full_metabolite$symptoms == "Active"] <- "active_diarrhea"
log_full_metabolite$symptoms[log_full_metabolite$symptoms == "Remission"] <- "no_diarrhea"

log_full_metabolite <- log_full_metabolite %>% select(1:4,1783,5:1782)

log_metabolite <- log_full_metabolite %>% select(1,2,5:7,9:1783)
log_metabolite <- as.data.frame(log_metabolite)
rownames(log_metabolite) <- log_metabolite$id



####---------------  SAVE to RDATA ---------------##########
save(metadata, metabolite, log_metabolite,
     file=paste0("mc_metabolite_input_",Sys.Date(),".RData"))



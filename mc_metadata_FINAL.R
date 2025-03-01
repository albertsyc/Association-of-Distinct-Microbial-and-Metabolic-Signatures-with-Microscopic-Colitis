rm(list=ls())
library(readxl)
library(dplyr)
library(tidyr)
library(forcats)
library(tableone)
df <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mc_0229_2024.xlsx")
setwd("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R")

# Cleaning ----------------------------------------------------------------
# Rename columns
mc <- df %>% 
  rename(
    'record_id' = "Record ID:",
    'stool_time' = "Event Name",
    'id' = "Study ID", 
    'old_mc' = "Known MC?",
    'old_mc_type' = "Type MC?...5",
    'new_mc' = "New MC?",
    'new_mc_type' = "Type MC?...9",
    'chronic_diarrhea' = "Chronic Diarrhea?",
    'stool_1' = "Stool 1 collected?",
    'active_1' = "Symptomatic at time of Stool 1?",
    'stool_2' = "Stool 2 collected?", 
    'active_2' = "Symptomatic at time of Stool 2?",
    'age' = "Age",
    'sex' = "Gender",
    'inch' = "Height:",
    'lb' = "Weight:",
    'hispanic' = "Self-Reported Ethnicity: Hispanic",
    'race' = "Self-Reported Race:",
    'smoke' = "Smoking status",
    'anti_chemo_immune_12' = "In the past 12 months, have you received any antibiotics, chemotherapy treatments, or immunosuppressants (including oral corticosteroids)?",
    'anti_12' = "In the past 12 months, have you taken any oral or IV (intravenous) antibiotics?",
    'nsaid_24' = "Have you used Non-Steroidal Anti-inflammatory Drugs (NSAIDS) regularly in the past 2 years?      - includes Ibuprofen (Advil, Motrin, Nuprin) Naproxen (Aleve, Naprosyn), Sulindac (Clinoril) Indomethacin(Indocin),Nabumetone (Relafen), Piroxican (Feldene), Ketoprofen (Orudis, Oruvail)    (Do NOT include aspirin-free products such as Tylenol/acetaminophen).  (Regular use is defined as greater than twice per week.)",
    'old_mc_med' = "Prev MC meds?",
    'old_mc_med_2' = "Known MC new/current meds? Include start date",
    'new_mc_med' = "New MC new meds? Include start date.",
    'mc_med_1' = "Microscopic colitis medication name and dose:...26",
    'mc_med_2' = "Microscopic colitis medication name and dose:...27",
    'mc_med_3' = "Microscopic colitis medication name and dose:...28",
    'mc_med_4' = "Microscopic colitis medication name and dose:...29",
    'mc_med_5' = "Microscopic colitis medication name and dose:...30",
    'mc_med_6' = "Microscopic colitis medication name and dose:...31",
    'mc_med_7' = "Microscopic colitis medication name and dose:...32",
    'mc_med_8' = "Microscopic colitis medication name and dose:...33"
  ) %>% select(1:5, 8:9, 11:25, 6, 7, 10, 26:33)

mc$inch <- as.numeric(mc$inch)
mc <- mc %>% mutate(bmi = lb*0.453592/(inch*0.0254)^2) %>% select(1:14, 34, 17:33)

# create new columns for medications
mc <- mc %>%
  mutate(
    med_salicylates = case_when(
      grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(old_mc_med)) |
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(old_mc_med_2))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(new_mc_med))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_1))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_2))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_3))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_4))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_5))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_6))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_7))|
        grepl("5[- ]?aminosalicylate|mesalamine|delzicol|canasa|asacol hd|pentasa|sfrowasa|lialda|apriso|asacol|balsalazide|giazo|colazal|colazide|sulfasalazine|azulfidine|salofalk|dipentum", tolower(mc_med_8))~ 1,
      TRUE ~ 0),
    med_steroid = case_when(
      grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(old_mc_med)) |
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(old_mc_med_2))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(new_mc_med))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_1))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_2))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_3))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_4))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_5))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_6))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_7))|
        grepl("budesonide|entocort|uceris|prednisone|deltasone|prednisode|steroids|Budesonide-8/2021", tolower(mc_med_8))~ 1,
      TRUE ~ 0),
    med_immunosuppressant_biologics = case_when(
      grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(old_mc_med)) |
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(old_mc_med_2))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(new_mc_med))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_1))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_2))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_3))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_4))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_5))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_6))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_7))|
        grepl("mp|6mp|methotrexate|humira|adalimumab|azathioprine|azasan|imuran|infliximab|avsola|inflectra|ixifi|remicade|renflexis|vedolizumab|entyvio|tofacitinib|xeljanz|golimumab|simponi|certolizumab pegol|cimzia|natalizumab|tysabri|ustekinumab|stelara", tolower(mc_med_8))~ 1,
      TRUE ~ 0)
  ) %>% select(1:18, 20:21, 33:35)


mc$med_salicylates <- as.factor(mc$med_salicylates)
mc$med_steroid <- as.factor(mc$med_steroid)
mc$med_immunosuppressant_biologics <- as.factor(mc$med_immunosuppressant_biologics)

# Match sample_id based on record_id
mc <- mc %>%
  arrange(record_id) %>%
  group_by(record_id) %>%
  fill(id) %>% 
  ungroup()%>%
  select(3, 2, 4:23)

mc$stool_time <- fct_recode(mc$stool_time,
                            "1" = "Event 1",
                            "2" = "Event 2")

# Correct medication history in event 2 (patients with positive med hx in event 1 should also have positive med hx in event 2)
mc <- mc %>%
  group_by(id) %>%
  mutate(med_salicylates = if_else(stool_time == 2 & med_salicylates < lag(med_salicylates) & !is.na(lag(med_salicylates)), lag(med_salicylates), med_salicylates),
         med_steroid = if_else(stool_time == 2 & med_steroid < lag(med_steroid) & !is.na(lag(med_steroid)), lag(med_steroid), med_steroid),
         med_immunosuppressant_biologics = if_else(stool_time == 2 & med_immunosuppressant_biologics < lag(med_immunosuppressant_biologics) & !is.na(lag(med_immunosuppressant_biologics)), lag(med_immunosuppressant_biologics), med_immunosuppressant_biologics)
         ) %>%
  ungroup()

# Categorized patients into mc_lc, mc_cc, mc_ns, mc_all, diarrhea
mc <- mc %>%
  mutate(
    mc_lc = ifelse(old_mc_type == "Lymphocytic" | new_mc_type == "Lymphocytic", 1, 0),
    mc_cc = ifelse(old_mc_type == "Collagenous" | new_mc_type == "Collagenous", 1, 0),
    mc_ns = ifelse(old_mc_type == "Non-specific" | new_mc_type == "Non-specific", 1, 0),
    mc_all = ifelse(old_mc == "Yes" | !is.na(old_mc_type) | new_mc == "Yes" | !is.na(new_mc_type), 1, 0),
    diarrhea = ifelse(chronic_diarrhea == "Yes" | (old_mc == "No" & is.na(old_mc_type) & new_mc == "No" & is.na(new_mc_type)), 1, 0)
  )

# Treat NA in mc_lc, mc_cc, mc_ns, mc_all, diarrhea as 0
mc <- mc %>%
  mutate_at(vars(mc_lc, mc_cc, mc_ns, mc_all, diarrhea), ~ replace(., is.na(.), 0)) %>%
  select(1:2, 23:27, 8:22, 3:7)

# Create an mc_or_diarrhea column
mc <- mc %>%
  mutate(
    mc_or_diarrhea = case_when(
      mc_lc == 1 ~ "Lymphocytic",
      mc_cc == 1 ~ "Collagenous",
      mc_ns == 1 ~ "Non-specific",
      diarrhea == 1 ~ "diarrhea",
      TRUE ~ "0"  # Default case if none of the above conditions are met
    )
  ) %>%
  select(1:2, 28, 3:27)

# Table 1A: active MC vs active chronic diarrhea (for those with microbiome data) -----------------------------------------------------------------
baseline = mc %>% filter(stool_time==1)

microbiome_list <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/microbiome_list.xlsx", col_names = FALSE)
microbiome_list <- microbiome_list %>% 
  rename(id_microbiome = ...1)
# Remove "A" and "B" from string values
microbiome_list$id_microbiome <- gsub("[AB]", "", microbiome_list$id_microbiome)
# Remove duplicates after removing "A" and "B"
microbiome_list <- microbiome_list[!duplicated(microbiome_list$id_microbiome), ]

# Find patients with microbiome data: filter 'baseline' dataframe based on matching values in 'id' column with 'id_microbiome' column in 'microbiome_list'
baseline_microbiome <- baseline %>%
  filter(id %in% microbiome_list$id_microbiome)

# Change Unknown in smoke, hispanic and race into NA
baseline_microbiome <- baseline_microbiome %>%
  mutate_at(vars(smoke, hispanic, race), ~ ifelse(. == "Unknown", NA, .))

# Number of active MC
baseline_microbiome %>% 
  filter((mc_or_diarrhea == "Lymphocytic" | mc_or_diarrhea == "Collagenous" | mc_or_diarrhea == "Non-specific") & 
           (active_1 == "Yes"| active_2 == "Yes" | id=="MC508")) %>% # MC508 is in active phase based on Bristol = 6~7
  filter(!(id=="MC087")) # MC087 is removed because the patient only has stool 2 (Remission) # n=131



## Imputation of missing BMI with linear regression by age, sex
# Filter out rows with complete data for regression
regression <- baseline_microbiome %>% filter(!is.na(bmi) & !is.na(age) & !is.na(sex))
# Fit linear regression model
model <- lm(bmi ~ age + sex, data = regression)
# Impute missing values
baseline_microbiome$imputed_bmi <- predict(model, newdata = baseline_microbiome)
# Replace missing BMI values with imputed values
baseline_microbiome$bmi[is.na(baseline_microbiome$bmi)] <- baseline_microbiome$imputed_bmi[is.na(baseline_microbiome$bmi)]
# Remove the temporary 'imputed_bmi' column
baseline_microbiome$imputed_bmi <- NULL # TOTAL: 13 bmi imputed



# Number of active chronic diarrhea: based on Jess's info
baseline_microbiome_cross_sectional <- baseline_microbiome %>% 
  filter(((mc_or_diarrhea == "Lymphocytic" | mc_or_diarrhea == "Collagenous" | mc_or_diarrhea == "Non-specific") & 
            (active_1 == "Yes"| active_2 == "Yes"| id=="MC508"))|(mc_all == 0)& !(id=="MC277"|id=="MC282"|id=="MC319"|id=="MC358"|id=="MC409"|id=="MC473"|id=="MC479")) %>% 
  filter(!(id=="MC087")) 

baseline_microbiome_cross_sectional # n=290

table_activemc.vs.diarrhea_microbiome <- CreateTableOne(vars = c("age", "sex", "bmi", "hispanic", "race", "smoke", "active_1", "active_2", "anti_12", 
                                                           "nsaid_24", "med_salicylates", "med_steroid", "med_immunosuppressant_biologics"), 
                                                  strata = "mc_all", 
                                                  data = baseline_microbiome_cross_sectional)
print(table_activemc.vs.diarrhea_microbiome, showAllLevels = TRUE, noSpaces = TRUE)

## Find MC that were treated v.s. untreated 
MC_untreated <- baseline_microbiome_cross_sectional %>% filter(mc_all == 1) %>% filter(med_salicylates == 0 & med_steroid == 0 & med_immunosuppressant_biologics == 0) 
# MC027 is treated with citrucel
# MC088 is treated with loperamide
# MC093 is treated with fiber
# MC156 is treated with Colesevelam and loperamide
# MC188 is treated with loperamide
# MC294 is treated with loperamide
# MC299 is treated with loperamide
# MC385 is treated with loperamide
# MC437 is treated with fiber 
# MC481 is treated with Cholestyramine
# MC489 is treated with loperamide
# MC499 is treated with Cholestyramine


# Table 1C: Active vs remission MC (for those with microbiome data) -----------------------------------------------------------------
# Number of MC with both active and remission data
baseline_microbiome %>% filter(mc_or_diarrhea == "Lymphocytic" | mc_or_diarrhea == "Collagenous" | mc_or_diarrhea == "Non-specific") %>% filter(active_1 == "Yes" & active_2 == "No") # n=59
baseline_microbiome %>% filter(mc_or_diarrhea == "Lymphocytic" | mc_or_diarrhea == "Collagenous" | mc_or_diarrhea == "Non-specific") %>% filter(active_1 == "No" & active_2 == "Yes") # n=14

# Number of LC with both active and remission data
baseline_microbiome %>% filter(mc_or_diarrhea == "Lymphocytic") %>% filter(active_1 == "Yes" & active_2 == "No") # n=31
baseline_microbiome %>% filter(mc_or_diarrhea == "Lymphocytic") %>% filter(active_1 == "No" & active_2 == "Yes") # n=8

# Number of CC with both active and remission data
baseline_microbiome %>% filter(mc_or_diarrhea == "Collagenous") %>% filter(active_1 == "Yes" & active_2 == "No") # n=23
baseline_microbiome %>% filter(mc_or_diarrhea == "Collagenous") %>% filter(active_1 == "No" & active_2 == "Yes") # n=5

# Number of NS with both active and remission data
baseline_microbiome %>% filter(mc_or_diarrhea == "Non-specific") %>% filter(active_1 == "Yes" & active_2 == "No") # n=5
baseline_microbiome %>% filter(mc_or_diarrhea == "Non-specific") %>% filter(active_1 == "No" & active_2 == "Yes") # n=1

baseline_microbiome_longitudinal <- baseline_microbiome %>% 
  filter((mc_or_diarrhea == "Lymphocytic" | mc_or_diarrhea == "Collagenous" | mc_or_diarrhea == "Non-specific") & 
           ((active_1 == "Yes" & active_2 == "No") | (active_1 == "No" & active_2 == "Yes"))) %>% 
  filter(!(id=="MC087"|id=="MC167"))

table_lc.vs.cc_microbiome <- CreateTableOne(vars = c("age", "sex", "bmi", "hispanic", "race", "smoke", "anti_12", 
                                                           "nsaid_24", "med_salicylates", "med_steroid", "med_immunosuppressant_biologics"), 
                                                  strata = "mc_or_diarrhea", 
                                                  data = baseline_microbiome_longitudinal)
print(table_lc.vs.cc_microbiome, showAllLevels = TRUE, noSpaces = TRUE)



# Table 1B: Active MC vs Healthy controls (1:3 age- and sex- matched) -----------------------------------------------------------------
load("match_2024-05-13.RData") 
df1 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_master_list/GIDER_demographics_corrected.xlsx")
df1 <- df1 %>% rename(id = GIDER_id) %>% select(1,4,5)
matched_hc <- match3df_controls %>% left_join(df1, by = join_by(id)) %>% arrange(id)
matched_hc$mc_all <- 2

# Get BMI info
q1 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_metadata/GIDER_Q1_05_17_21.xlsx")
q1 <- q1 %>% select(study_id, antibiotics_lastyear, ht_in, wt_lbs, cigarettes)
q1$antibiotics_lastyear <- ifelse(is.na(q1$antibiotics_lastyear), 0, q1$antibiotics_lastyear)
q1 <- q1 %>% rename("id" = "study_id",
                    "anti_12" = "antibiotics_lastyear",
                    "ht" = "ht_in",
                    "wt" = "wt_lbs",
                    "smoke" = "cigarettes")
q1 <- q1 %>% mutate(bmi = wt*0.453592/(ht*0.0254)^2)
q1_id <- q1$id


q2 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_metadata/GIDER_Q2_5_17_21.xlsx")
q2 <- q2 %>% select(study_id, antibiotics_lastyear_cb9ef1, ht_inches, weight, cigarettes_a67d61)
q2 <- q2 %>% rename("id" = "study_id",
                    "anti_12" = "antibiotics_lastyear_cb9ef1",
                    "ht" = "ht_inches",
                    "wt" = "weight",
                    "smoke" = "cigarettes_a67d61")
q2 <- q2[!q2$id %in% q1_id, ]
df2 <- inner_join(q2,matched_hc) %>% select(-age, -sex, -hispanic, -race, -mc_all)

df2$ht <- fct_recode(df2$ht,
                     "62.4" = "5.2",
                     "69.6" = "5.8",
                     "65" = "22", # 22 inches is implausible, so I impute with median 65 inches. 
                     "64.17" = "163", # 163 is likely in cm.
                     "57.5" = "4 9 1/2",
                     "59.5" = "4 ft 11 1/2 inches",
                     "53.19" = "44325",
                     "62" = "5  2\"",
                     "70" = "5 feet 10 inches",
                     "65" = "5 foot 5 inches",
                     "60" = "5 ft",
                     "62.5" = '5\' 2.5"',
                     "62" = '5\' 2"',
                     "64" = "5' 4\"",
                     "66" = "5' 6\"",
                     "67" = "5' 7\"",
                     "62" = "5' 2\"",
                     "64" = "5' 4\"",
                     "66" = "5' 6\"",
                     "67" = "5' 7\"",
                     "64" = "5' 4\"",
                     "62" = "62 IN",
                     "64" = "64\"",
                     "66" = "66\"",
                     "67" = "67\"",
                     "69" = "69\"",
                     "74" = "74 Inches",
                     "66" = "5'6\"",
                     "67" = "5'7\"",
                     "62" = "5'2\"",
                     "63" = "5'3\"",
                     "70" = "5 10\"",
                     "70" = "5'10\"",
                     "66.5" = "5' 6.5\"",
                     "61" = "5'1\"",
                     "60.5" = "5'1/2\"",
                     "65" = "5'5",
                     "65.5" = "5'5.5\"",
                     "65" = "5'5\"",
                     "66" = "5'6",
                     "66" = "5'6''",
                     "69.75" = "5'9 3/4\"",
                     "69" = "5'9'",
                     "64" = "5'4\"",
                     "64" = "5'4'",
                     "64" = '5"4',
                     "66" = "5feet 6 inches")

df2$wt <- gsub(" lbs", "", df2$wt)
df2$wt <- gsub("lbs", "", df2$wt)

# Remove non-numeric characters from ht column
df2$ht <- gsub("[^0-9.]", "", df2$ht)

df2$ht <- as.numeric(df2$ht)
df2$wt <- as.numeric(df2$wt)
df2 <- df2 %>% mutate(bmi = wt*0.453592/(ht*0.0254)^2)
q2 <- df2
q2_id <- q2$id



q3 <- read_excel("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/GIDER_metadata/GIDER_Q3_5_17_21.xlsx")
q3 <- q3 %>% select(study_id, antibiotics_lastyr_adl, ht_pmh, weight_pmh, cigarettes_ch)
q3 <- q3 %>% rename("id" = "study_id",
                    "anti_12" = "antibiotics_lastyr_adl",
                    "ht" = "ht_pmh",
                    "wt" = "weight_pmh",
                    "smoke" = "cigarettes_ch")
q3 <- q3[!(q3$id %in% q2_id) & !(q3$id %in% q1_id), ]
df3 <- inner_join(q3,matched_hc) %>% select(-age, -sex, -hispanic, -race, -mc_all)

df3$ht <- fct_recode(df3$ht,
                     "64.5" = "5' 41/2\"",
                     "67.5" = "5'7\"",
                     "70.5" = "5ft 21/2 inches",
                     "63" = "5'3\"",
                     "63" = "5' 3\"",
                     "63.5" = "5'4 1/2",
                     "71" = "5'11\"",
                     "67" = "5' 7\"")

# Remove non-numeric characters from ht column
df3$ht <- gsub("[^0-9.]", "", df3$ht)

df3$ht <- as.numeric(df3$ht)
df3$wt <- as.numeric(df3$wt)
df3 <- df3 %>% mutate(bmi = wt*0.453592/(ht*0.0254)^2)
q3 <- df3


matched_hc <- matched_hc %>% left_join(q1, by = c("id" = "id"))
matched_hc <- matched_hc %>% left_join(q2, by = c("id" = "id"))

matched_hc$anti_12.x <- ifelse(is.na(matched_hc$anti_12.x), matched_hc$anti_12.y, matched_hc$anti_12.x)
matched_hc$bmi.x <- ifelse(is.na(matched_hc$bmi.x), matched_hc$bmi.y, matched_hc$bmi.x)
matched_hc$smoke.x <- ifelse(is.na(matched_hc$smoke.x), matched_hc$smoke.y, matched_hc$smoke.x)
matched_hc <- matched_hc %>% select(-ht.x, - wt.x, -anti_12.y, -ht.y, -wt.y, -bmi.y, -smoke.y)

matched_hc <- matched_hc %>% left_join(q3, by = c("id" = "id"))
matched_hc$anti_12 <- ifelse(is.na(matched_hc$anti_12), matched_hc$anti_12.x, matched_hc$anti_12)
matched_hc$bmi <- ifelse(is.na(matched_hc$bmi), matched_hc$bmi.x, matched_hc$bmi)
matched_hc$smoke <- ifelse(is.na(matched_hc$smoke), matched_hc$smoke.x, matched_hc$smoke)
matched_hc <- matched_hc %>% select(-anti_12.x, -ht, -wt, -bmi.x, -smoke.x)

matched_hc$anti_12 <- ifelse(matched_hc$anti_12 == 0, "No", "Yes")
matched_hc$smoke <- ifelse(matched_hc$smoke == 1, "Current",
                           ifelse(matched_hc$smoke == 2, "Past",
                                  ifelse(matched_hc$smoke == 3, "Never", matched_hc$smoke)))


## Imputation of missing BMI with linear regression by age, sex
# Filter out rows with complete data for regression
df_regression <- matched_hc %>% filter(!is.na(bmi), !is.na(age), !is.na(sex))
# Fit linear regression model
lm_model <- lm(bmi ~ age + sex, data = df_regression)
# Impute missing values
matched_hc$imputed_bmi <- predict(lm_model, newdata = matched_hc)
# Replace missing BMI values with imputed values
matched_hc$bmi[is.na(matched_hc$bmi)] <- matched_hc$imputed_bmi[is.na(matched_hc$bmi)]
# Remove the temporary 'imputed_bmi' column
matched_hc$imputed_bmi <- NULL


## Create Table 1
baseline_microbiome_cross_sectional_1 <- baseline_microbiome_cross_sectional %>% filter(mc_all==1)
hc_mc_cd <- rbind(matched_hc, baseline_microbiome_cross_sectional)

table_activemc.vs.hc_microbiome <- CreateTableOne(vars = c("age", "sex", "bmi", "hispanic", "race", "smoke", "anti_12", "med_salicylates", "med_steroid", "med_immunosuppressant_biologics"), 
                                                        strata = "mc_all", 
                                                        data = hc_mc_cd)
print(table_activemc.vs.hc_microbiome, showAllLevels = TRUE, noSpaces = TRUE)

# match3df_controls <- match3df_controls %>% arrange(id)
# write.xlsx(match3df_controls, "to_match.xlsx")

####---------------  SAVE to RDATA ---------------##########
setwd("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R")
save(baseline_microbiome, baseline_microbiome_cross_sectional, baseline_microbiome_longitudinal, matched_hc, 
     file=paste0("metadata_",Sys.Date(),".RData"))





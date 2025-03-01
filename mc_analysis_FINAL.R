#####------------------- MC Tables and Figures  ------------------######
##### What: MC vs chronic diarrhea vs healthy controls analysis for tables and figures
##### Who: Albert Chen
##### When March 2024 

rm(list=ls())


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
library(pheatmap)
library(readxl)
library(ggsignif)
library(compositions)

main.dir<- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R"
setwd(main.dir)

####---------------  Alpha diversity (longitudinal) (Chao 1) ----------------#####
load("mc_input_2024-07-08.RData") 
load("metadata_2024-05-13.RData") 
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/sample_names.xlsx"
# write.xlsx(all_taxo[1], file_path, rowNames = TRUE)

# all_taxo[, 1] <- paste0(all_taxo[, 1], ".fastq")
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/file_names.xlsx"
# write.xlsx(all_taxo[1], file_path, rowNames = TRUE)
## Chao1 index
# Longitudinal: prepare patients list (patients with both active and remission data)
wide.alpha <- alpha_with_disease_status[,c("ID","disease_status","alpha.chao1")] %>% 
  pivot_wider(names_from = disease_status, values_from = alpha.chao1) 

wide.alpha$id<- substr(wide.alpha$ID, 1,5)
wide.alpha <- wide.alpha %>%
  left_join(select(baseline_microbiome, id, mc_or_diarrhea), by = "id")

wide.alpha.longitudinal_1 <- wide.alpha %>% filter(mc_or_diarrhea == "Lymphocytic" | mc_or_diarrhea == "Collagenous" | mc_or_diarrhea == "Non-specific")
wide.alpha.longitudinal_2 <- wide.alpha.longitudinal_1 %>% group_by(id) %>% filter(n()>1)

id_with_both_active_and_remission <- wide.alpha.longitudinal_2 %>% filter((is.na(active) & !is.na(remission) & !is.na(lag(active)) & is.na(lag(remission))) |
                                       (!is.na(active) & is.na(remission) & is.na(lag(active)) & !is.na(lag(remission))))

list_id_with_both_active_and_remission <- id_with_both_active_and_remission %>% select(id)

wide.alpha.longitudinal_final <- wide.alpha.longitudinal_2 %>%
  filter(id %in% list_id_with_both_active_and_remission$id) %>% select(-ID) %>%
  group_by(id) %>%
  summarise(across(everything(), ~max(., na.rm = TRUE)))

wide.alpha.longitudinal_final

# Compare Longitudinal 
df.long <- alpha_with_disease_status[alpha_with_disease_status$id %in% wide.alpha.longitudinal_final$id, ]
stat.test_long_chao1 <- compare_means(alpha.chao1 ~ disease_status, 
                                data=df.long,
                                method = "wilcox.test", paired = T, p.adjust.method="none") 

# Reorder levels of disease_status
df.long$disease_status <- factor(df.long$disease_status, levels = c("active", "remission"))

bx.a.l <- ggboxplot(df.long, y = "alpha.chao1", x = "disease_status",
                 fill="disease_status", palette = c("#900C3F", "antiquewhite4"),
                 xlab="MC activity", ylab="Alpha diversity (Chao1 index)", 
                 ylim=c(0,400),font.x=12, font.y=12, font.tickslab=11, 
                 legend = "none")

fig_chao1_long <- bx.a.l +
  stat_pvalue_manual(stat.test_long_chao1, label = "Wilcoxon p = {p.format}*", y.position = 400) +
  annotate("text", x=1, y=15, label = "n = 66") +
  annotate("text", x=2, y=15, label = "n = 66")

fig_chao1_long

# ggsave(filename = "fig_chao1_long.pdf",
#        plot=fig_chao1_long, units = "in", height = 7, width = 8, dpi=1200)



####---------------  Alpha diversity (Cross-sectional) (Chao 1) ----------------#####
## Cross-sectional: prepare patients list (At baseline, Active MC vs Active diarrhea vs HC)
# Cross: MC vs CD
alpha.cross.sectional_active_1 <- alpha_with_disease_status %>% group_by(id) %>% 
  filter(disease_status=="active") %>% 
  distinct(id, .keep_all = TRUE)

alpha.cross.sectional_active_1 <- alpha.cross.sectional_active_1 %>% 
  mutate(mc_all_label = case_when(
    mc_all == 1 ~ "MC",
    mc_all == 0 ~ "Chronic diarrhea",
    TRUE ~ NA_character_))

alpha.cross.sectional_active_2 <- alpha_with_disease_status %>% filter(mc_all==2)
alpha.cross.sectional_active_2$mc_all_label <- "Healthy control"

alpha.cross.sectional_active <- rbind(alpha.cross.sectional_active_1, alpha.cross.sectional_active_2)

# Pivot to wide form
wide.alpha.cross.sectional_active <- alpha.cross.sectional_active[,c("ID","mc_all_label","alpha.chao1")] %>% 
  pivot_wider(names_from = mc_all_label, values_from = alpha.chao1) 

wide.alpha.cross.sectional_active_1 <- wide.alpha.cross.sectional_active[substr(wide.alpha.cross.sectional_active$ID, 1, 1) == "G", ]
wide.alpha.cross.sectional_active_2 <- wide.alpha.cross.sectional_active[substr(wide.alpha.cross.sectional_active$ID, 1, 2) == "MC", ]
wide.alpha.cross.sectional_active_1$id <- wide.alpha.cross.sectional_active_1$ID
wide.alpha.cross.sectional_active_2$id<- substr(wide.alpha.cross.sectional_active_2$ID, 1,5)

wide.alpha.cross.sectional_active <- rbind(wide.alpha.cross.sectional_active_1, wide.alpha.cross.sectional_active_2)

wide.alpha.cross.sectional_active <- wide.alpha.cross.sectional_active %>% select(5,2:4)


# Compare Cross-sectional 
df.cross <- alpha.cross.sectional_active 
df.cross$mc_all_label <- factor(df.cross$mc_all_label, 
                                levels = c("Chronic diarrhea", "MC", "Healthy control"), 
                                labels = c("Chronic diarrhea", "MC", "Control without diarrhea"))

stat.test_cross_chao1 <- compare_means(alpha.chao1 ~ mc_all_label, 
                                      data=df.cross,
                                      method = "wilcox.test", paired = F, p.adjust.method="none") 

bx.a.c <- ggboxplot(df.cross, y = "alpha.chao1", x = "mc_all_label",
                    fill = "mc_all_label", palette = c("#F9D937", "#900C3F", "#189AF9"),
                    xlab = "Disease type", ylab = "Alpha diversity (Chao1 index)", 
                    ylim = c(0, 400), font.x = 12, font.y = 12, font.tickslab = 11, 
                    legend = "none") +
  theme(axis.text.x = element_text(hjust = 0.5))

df.cross %>% filter(mc_all_label=="MC") # MC: n = 131
df.cross %>% filter(mc_all_label=="Chronic diarrhea") # Chronic diarrhea: n = 159
df.cross %>% filter(mc_all_label=="Healthy control") # Healthy control: n = 393


stat.test_cross_chao1 <- stat.test_cross_chao1 %>%
  mutate(y.position = c(345, 370, 400),
         label = c("Wilcoxon p = 0.446", "Wilcoxon p < 0.001*", "Wilcoxon p < 0.001*"))

fig_chao1_cross <- bx.a.c +
  stat_pvalue_manual(stat.test_cross_chao1)+
  annotate("text", x=1, y=320, label = "n = 159") +
  annotate("text", x=2, y=320, label = "n = 131") +
  annotate("text", x=3, y=320, label = "n = 393")

fig_chao1_cross

# ggsave(filename = "fig_chao1_cross.pdf",
#        plot=fig_chao1_cross, units = "in", height = 7, width = 8, dpi=1200)


####---------------  Alpha diversity (longitudinal) (Shannon) ----------------#####
# Longitudinal: prepare patients list (patients with both active and remission data)
wide.alpha_shannon <- alpha_with_disease_status[,c("ID","disease_status","alpha.shannon")] %>% 
  pivot_wider(names_from = disease_status, values_from = alpha.shannon) 

wide.alpha_shannon$id<- substr(wide.alpha_shannon$ID, 1,5)
wide.alpha_shannon <- wide.alpha_shannon %>%
  left_join(select(baseline_microbiome, id, mc_or_diarrhea), by = "id")

wide.alpha.longitudinal_shannon_1 <- wide.alpha_shannon %>% filter(mc_or_diarrhea == "Lymphocytic" | mc_or_diarrhea == "Collagenous" | mc_or_diarrhea == "Non-specific")
wide.alpha.longitudinal_shannon_2 <- wide.alpha.longitudinal_shannon_1 %>% group_by(id) %>% filter(n()>1)

wide.alpha.longitudinal_final_shannon <- wide.alpha.longitudinal_shannon_2 %>%
  filter(id %in% list_id_with_both_active_and_remission$id) %>% select(-ID) %>%
  group_by(id) %>%
  summarise(across(everything(), ~max(., na.rm = TRUE)))

wide.alpha.longitudinal_final_shannon

# Compare Longitudinal 
df.long_shannon <- alpha_with_disease_status[alpha_with_disease_status$id %in% wide.alpha.longitudinal_final_shannon$id, ]
stat.test_long_shannon <- compare_means(alpha.shannon ~ disease_status, 
                                      data=df.long_shannon,
                                      method = "wilcox.test", paired = T, p.adjust.method="none") 

# Reorder levels of disease_status
df.long_shannon$disease_status <- factor(df.long_shannon$disease_status, levels = c("active", "remission"))

bx.a.l_shannon <- ggboxplot(df.long_shannon, y = "alpha.shannon", x = "disease_status",
                    fill="disease_status", palette = c("#900C3F", "antiquewhite4"),
                    xlab="MC activity", ylab="Alpha diversity (Shannon index)", 
                    ylim=c(0,5),font.x=12, font.y=12, font.tickslab=11, 
                    legend = "none")

fig_shannon_long <- bx.a.l_shannon +
  stat_pvalue_manual(stat.test_long_shannon, label = "Wilcoxon p = {p.format}", y.position = 5) +
  annotate("text", x=1, y=1, label = "n = 66") +
  annotate("text", x=2, y=1, label = "n = 66")

fig_shannon_long

# ggsave(filename = "fig_shannon_long.pdf",
#        plot=fig_shannon_long, units = "in", height = 5, width = 5, dpi=1200)



####---------------  Alpha diversity (Cross-sectional) (Shannon) ----------------#####
## Cross-sectional: prepare patients list (At baseline, Active MC vs Active diarrhea vs HC)
# Cross: MC vs CD
# Pivot to wide form
wide.alpha.cross.sectional_active_shannon <- alpha.cross.sectional_active[,c("ID","mc_all_label","alpha.shannon")] %>% 
  pivot_wider(names_from = mc_all_label, values_from = alpha.shannon) 

wide.alpha.cross.sectional_active_shannon_1 <- wide.alpha.cross.sectional_active_shannon[substr(wide.alpha.cross.sectional_active_shannon$ID, 1, 1) == "G", ]
wide.alpha.cross.sectional_active_shannon_2 <- wide.alpha.cross.sectional_active_shannon[substr(wide.alpha.cross.sectional_active_shannon$ID, 1, 2) == "MC", ]
wide.alpha.cross.sectional_active_shannon_1$id <- wide.alpha.cross.sectional_active_shannon_1$ID
wide.alpha.cross.sectional_active_shannon_2$id<- substr(wide.alpha.cross.sectional_active_shannon_2$ID, 1,5)

wide.alpha.cross.sectional_active_shannon <- rbind(wide.alpha.cross.sectional_active_shannon_1, wide.alpha.cross.sectional_active_shannon_2)

wide.alpha.cross.sectional_active_shannon <- wide.alpha.cross.sectional_active_shannon %>% select(5,2:4)


# Compare Cross-sectional 
stat.test_cross_shannon <- compare_means(alpha.shannon ~ mc_all_label, 
                                       data=df.cross,
                                       method = "wilcox.test", paired = F, p.adjust.method="none") 

bx.a.s <- ggboxplot(df.cross, y = "alpha.shannon", x = "mc_all_label",
                    fill = "mc_all_label", palette = c("#F9D937", "#900C3F", "#189AF9"),
                    xlab = "Disease type", ylab = "Alpha diversity (Shannon index)", 
                    ylim = c(0, 6), font.x = 12, font.y = 12, font.tickslab = 11, 
                    legend = "none") +
  theme(axis.text.x = element_text(hjust = 0.5))


stat.test_cross_shannon <- stat.test_cross_shannon %>%
  mutate(y.position = c(5, 5.5, 6),
         label = c("Wilcoxon p = 0.130", "Wilcoxon p < 0.001*", "Wilcoxon p < 0.001*"))

fig_shannon_cross <- bx.a.s +
  stat_pvalue_manual(stat.test_cross_shannon)+
  annotate("text", x=1, y=4.6, label = "n = 159") +
  annotate("text", x=2, y=4.6, label = "n = 131") +
  annotate("text", x=3, y=4.6, label = "n = 393")

fig_shannon_cross

# ggsave(filename = "fig_shannon_cross.pdf",
#        plot=fig_shannon_cross, units = "in", height = 6, width = 6, dpi=1200)




####---------------  Alpha diversity (LC and CC) (Chao 1) ----------------#####
### Longitudinal 
list_id_with_both_active_and_remission_LC <- id_with_both_active_and_remission %>% filter(mc_or_diarrhea=="Lymphocytic") %>% select(id) # n=33
list_id_with_both_active_and_remission_CC <- id_with_both_active_and_remission %>% filter(mc_or_diarrhea=="Collagenous") %>% select(id) # n=27

wide.alpha.longitudinal_final_LC <- wide.alpha.longitudinal_2 %>%
  filter(id %in% list_id_with_both_active_and_remission_LC$id) %>% select(-ID) %>%
  group_by(id) %>%
  summarise(across(everything(), ~max(., na.rm = TRUE)))

wide.alpha.longitudinal_final_CC <- wide.alpha.longitudinal_2 %>%
  filter(id %in% list_id_with_both_active_and_remission_CC$id) %>% select(-ID) %>%
  group_by(id) %>%
  summarise(across(everything(), ~max(., na.rm = TRUE)))

# Compare Longitudinal 
df.long_LC <- alpha_with_disease_status[alpha_with_disease_status$id %in% wide.alpha.longitudinal_final_LC$id, ]
stat.test_long_chao1_LC <- compare_means(alpha.chao1 ~ disease_status, 
                                         data=df.long_LC,
                                         method = "wilcox.test", paired = T, p.adjust.method="none") 

df.long_CC <- alpha_with_disease_status[alpha_with_disease_status$id %in% wide.alpha.longitudinal_final_CC$id, ]
stat.test_long_chao1_CC <- compare_means(alpha.chao1 ~ disease_status, 
                                         data=df.long_CC,
                                         method = "wilcox.test", paired = T, p.adjust.method="none") 


# Reorder levels of disease_status
df.long_LC$disease_status <- factor(df.long_LC$disease_status, levels = c("active", "remission"))
df.long_CC$disease_status <- factor(df.long_CC$disease_status, levels = c("active", "remission"))

bx.a.l_LC <- ggboxplot(df.long_LC, y = "alpha.chao1", x = "disease_status",
                       fill="disease_status", palette = c("#900C3F", "antiquewhite4"),
                       xlab="LC activity", ylab="Alpha diversity (Chao1 index)", 
                       ylim=c(0,400),font.x=12, font.y=12, font.tickslab=11, 
                       legend = "none")

fig_chao1_long_LC <- bx.a.l_LC +
  stat_pvalue_manual(stat.test_long_chao1_LC, label = "Wilcoxon p = {p.format}", y.position = 400) +
  annotate("text", x=1, y=15, label = "n = 33") +
  annotate("text", x=2, y=15, label = "n = 33")

fig_chao1_long_LC

# ggsave(filename = "fig_chao1_long_LC.pdf",
#        plot=fig_chao1_long_LC, units = "in", height = 7, width = 8, dpi=1200)

bx.a.l_CC <- ggboxplot(df.long_CC, y = "alpha.chao1", x = "disease_status",
                       fill="disease_status", palette = c("#900C3F", "antiquewhite4"),
                       xlab="CC activity", ylab="Alpha diversity (Chao1 index)", 
                       ylim=c(0,400),font.x=12, font.y=12, font.tickslab=11, 
                       legend = "none")

fig_chao1_long_CC <- bx.a.l_CC +
  stat_pvalue_manual(stat.test_long_chao1_CC, label = "Wilcoxon p = {p.format}", y.position = 400) +
  annotate("text", x=1, y=15, label = "n = 27") +
  annotate("text", x=2, y=15, label = "n = 27")

fig_chao1_long_CC

# ggsave(filename = "fig_chao1_long_CC.pdf",
#        plot=fig_chao1_long_CC, units = "in", height = 7, width = 8, dpi=1200)


### Cross-sectional
alpha.cross.sectional_active_LC_CC <- alpha.cross.sectional_active %>%
  mutate(mc_or_diarrhea = recode(mc_or_diarrhea,
                                 "Lymphocytic" = "LC",
                                 "Collagenous" = "CC",
                                 "diarrhea" = "CD",
                                 "Non-specific" = "NS",
                                 "Healthy" = "HC")) %>%
  filter(!mc_or_diarrhea == "NS")


# Pivot to wide form
wide.alpha.cross.sectional_active_LC_CC <- alpha.cross.sectional_active_LC_CC[,c("ID","mc_or_diarrhea","alpha.chao1")] %>% 
  pivot_wider(names_from = mc_or_diarrhea, values_from = alpha.chao1) 

wide.alpha.cross.sectional_active_1_LC_CC <- wide.alpha.cross.sectional_active_LC_CC[substr(alpha.cross.sectional_active_LC_CC$ID, 1, 1) == "G", ]
wide.alpha.cross.sectional_active_2_LC_CC <- wide.alpha.cross.sectional_active_LC_CC[substr(alpha.cross.sectional_active_LC_CC$ID, 1, 2) == "MC", ]
wide.alpha.cross.sectional_active_1_LC_CC$id <- wide.alpha.cross.sectional_active_1_LC_CC$ID
wide.alpha.cross.sectional_active_2_LC_CC$id<- substr(wide.alpha.cross.sectional_active_2_LC_CC$ID, 1,5)

wide.alpha.cross.sectional_active_LC_CC <- rbind(wide.alpha.cross.sectional_active_1_LC_CC, wide.alpha.cross.sectional_active_2_LC_CC)

wide.alpha.cross.sectional_active_LC_CC <- wide.alpha.cross.sectional_active_LC_CC %>% select(6,2:5)


# Compare Cross-sectional 
df.cross_LC_CC <- alpha.cross.sectional_active_LC_CC 

df.cross_LC_CC$mc_or_diarrhea <- factor(df.cross_LC_CC$mc_or_diarrhea, levels = c("CD", "LC", "CC", "HC"), labels = c("Chronic diarrhea", "LC", "CC", "Control without diarrhea"))

stat.test_cross_chao1_LC_CC <- compare_means(alpha.chao1 ~ mc_or_diarrhea, 
                                             data=df.cross_LC_CC,
                                             method = "wilcox.test", paired = F, p.adjust.method="none") 

bx.a.c_LC_CC <- ggboxplot(df.cross_LC_CC, y = "alpha.chao1", x = "mc_or_diarrhea",
                          fill = "mc_or_diarrhea", palette = c("#F9D937", "#900C3F", "#900C3F", "#189AF9"),
                          xlab = "Disease type", ylab = "Alpha diversity (Chao1 index)", 
                          ylim = c(0, 480), font.x = 12, font.y = 12, font.tickslab = 11, 
                          legend = "none") +
  theme(axis.text.x = element_text(hjust = 0.5))

df.cross_LC_CC %>% filter(mc_or_diarrhea=="LC") # LC: n = 69
df.cross_LC_CC %>% filter(mc_or_diarrhea=="CC") # CC: n = 52
df.cross_LC_CC %>% filter(mc_or_diarrhea=="HC") # HC: n = 393
df.cross_LC_CC %>% filter(mc_or_diarrhea=="CD") # HC: n = 159


stat.test_cross_chao1_LC_CC <- stat.test_cross_chao1_LC_CC %>%
  mutate(y.position = c(325, 395, 460, 350, 430, 375),
         label = c("p = 0.103", "p = 0.751", "p < 0.001*",
                   "p = 0.136", "p < 0.001*", "p < 0.001*"))

fig_chao1_cross_LC_CC <- bx.a.c_LC_CC +
  stat_pvalue_manual(stat.test_cross_chao1_LC_CC)+
  annotate("text", x=1, y=300, label = "n = 159") +
  annotate("text", x=2, y=300, label = "n = 69") +
  annotate("text", x=3, y=300, label = "n = 52") +
  annotate("text", x=4, y=300, label = "n = 393")

fig_chao1_cross_LC_CC

# ggsave(filename = "fig_chao1_cross_LC_CC.pdf",
#        plot=fig_chao1_cross_LC_CC, units = "in", height = 7, width = 9, dpi=1200)


####---------------  Alpha diversity (LC and CC) (Shannon) ----------------#####
### Longitudinal 
wide.alpha.longitudinal_final_shannon_LC <- wide.alpha.longitudinal_shannon_2 %>%
  filter(id %in% list_id_with_both_active_and_remission_LC$id) %>% select(-ID) %>%
  group_by(id) %>%
  summarise(across(everything(), ~max(., na.rm = TRUE)))

wide.alpha.longitudinal_final_shannon_CC <- wide.alpha.longitudinal_shannon_2 %>%
  filter(id %in% list_id_with_both_active_and_remission_CC$id) %>% select(-ID) %>%
  group_by(id) %>%
  summarise(across(everything(), ~max(., na.rm = TRUE)))

# Compare Longitudinal 
stat.test_long_shannon_LC <- compare_means(alpha.shannon ~ disease_status, 
                                         data=df.long_LC,
                                         method = "wilcox.test", paired = T, p.adjust.method="none") 

stat.test_long_shannon_CC <- compare_means(alpha.shannon ~ disease_status, 
                                         data=df.long_CC,
                                         method = "wilcox.test", paired = T, p.adjust.method="none") 


bx.a.l_shannon_LC <- ggboxplot(df.long_LC, y = "alpha.shannon", x = "disease_status",
                       fill="disease_status", palette = c("#900C3F", "antiquewhite4"),
                       xlab="LC activity", ylab="Alpha diversity (Shannon index)", 
                       ylim=c(0,5),font.x=12, font.y=12, font.tickslab=11, 
                       legend = "none")

fig_shannon_long_LC <- bx.a.l_shannon_LC +
  stat_pvalue_manual(stat.test_long_shannon_LC, label = "Wilcoxon p = {p.format}", y.position = 5) +
  annotate("text", x=1, y=1, label = "n = 33") +
  annotate("text", x=2, y=1, label = "n = 33")

fig_shannon_long_LC

# ggsave(filename = "fig_shannon_long_LC.pdf",
#        plot=fig_shannon_long_LC, units = "in", height = 5, width = 5, dpi=1200)

bx.a.l_shannon_CC <- ggboxplot(df.long_CC, y = "alpha.shannon", x = "disease_status",
                       fill="disease_status", palette = c("#900C3F", "antiquewhite4"),
                       xlab="CC activity", ylab="Alpha diversity (Shannon index)", 
                       ylim=c(0,5),font.x=12, font.y=12, font.tickslab=11, 
                       legend = "none")

fig_shannon_long_CC <- bx.a.l_shannon_CC +
  stat_pvalue_manual(stat.test_long_shannon_CC, label = "Wilcoxon p = {p.format}", y.position = 5) +
  annotate("text", x=1, y=1, label = "n = 27") +
  annotate("text", x=2, y=1, label = "n = 27")

fig_shannon_long_CC

# ggsave(filename = "fig_shannon_long_CC.pdf",
#        plot=fig_shannon_long_CC, units = "in", height = 5, width = 5, dpi=1200)


### Cross-sectional
# Pivot to wide form
wide.alpha.cross.sectional_active_LC_CC_shannon <- alpha.cross.sectional_active_LC_CC[,c("ID","mc_or_diarrhea","alpha.shannon")] %>% 
  pivot_wider(names_from = mc_or_diarrhea, values_from = alpha.shannon) 

wide.alpha.cross.sectional_active_1_LC_CC_shannon <- wide.alpha.cross.sectional_active_LC_CC_shannon[substr(alpha.cross.sectional_active_LC_CC$ID, 1, 1) == "G", ]
wide.alpha.cross.sectional_active_2_LC_CC_shannon <- wide.alpha.cross.sectional_active_LC_CC_shannon[substr(alpha.cross.sectional_active_LC_CC$ID, 1, 2) == "MC", ]
wide.alpha.cross.sectional_active_1_LC_CC_shannon$id <- wide.alpha.cross.sectional_active_1_LC_CC_shannon$ID
wide.alpha.cross.sectional_active_2_LC_CC_shannon$id<- substr(wide.alpha.cross.sectional_active_2_LC_CC_shannon$ID, 1,5)

wide.alpha.cross.sectional_active_LC_CC_shannon <- rbind(wide.alpha.cross.sectional_active_1_LC_CC_shannon, wide.alpha.cross.sectional_active_2_LC_CC_shannon)

wide.alpha.cross.sectional_active_LC_CC_shannon <- wide.alpha.cross.sectional_active_LC_CC_shannon %>% select(6,2:5)


# Compare Cross-sectional 
stat.test_cross_shannon_LC_CC <- compare_means(alpha.shannon ~ mc_or_diarrhea, 
                                             data=df.cross_LC_CC,
                                             method = "wilcox.test", paired = F, p.adjust.method="none") 

bx.a.s_LC_CC <- ggboxplot(df.cross_LC_CC, y = "alpha.shannon", x = "mc_or_diarrhea",
                          fill = "mc_or_diarrhea", palette = c("#F9D937", "#900C3F", "#900C3F", "#189AF9"),
                          xlab = "Disease type", ylab = "Alpha diversity (Shannon index)", 
                          ylim = c(0, 6.8), font.x = 12, font.y = 12, font.tickslab = 11, 
                          legend = "none") +
  theme(axis.text.x = element_text(hjust = 0.5))

stat.test_cross_shannon_LC_CC <- stat.test_cross_shannon_LC_CC %>%
  mutate(y.position = c(5.0, 5.9, 6.6, 5.3, 6.3, 5.6),
         label = c("p = 0.025", "p = 0.864", "p < 0.001*",
                   "p = 0.104", "p = 0.002*", "p < 0.001*"))

fig_shannon_cross_LC_CC <- bx.a.s_LC_CC +
  stat_pvalue_manual(stat.test_cross_shannon_LC_CC)+
  annotate("text", x=1, y=4.7, label = "n = 159") +
  annotate("text", x=2, y=4.7, label = "n = 69") +
  annotate("text", x=3, y=4.7, label = "n = 52") +
  annotate("text", x=4, y=4.7, label = "n = 393")

fig_shannon_cross_LC_CC

# ggsave(filename = "fig_shannon_cross_LC_CC.pdf",
#        plot=fig_shannon_cross_LC_CC, units = "in", height = 6, width = 6, dpi=1200)





####---------------  PCoA (cross-sectional: Bray-Curtis)  ----------------------#####
pcoa_taxo_cross <- alpha.cross.sectional_active[, 1] %>% left_join(all_taxo, by = c("ID" = "id"))
pcoa_taxo_cross <- as.data.frame(pcoa_taxo_cross)
rownames(pcoa_taxo_cross) <- pcoa_taxo_cross$ID


subset.pcoa_taxo_cross <- pcoa_taxo_cross[, !colnames(pcoa_taxo_cross) %in% c("ID")]
bray_cross <- vegdist(subset.pcoa_taxo_cross, "bray")  # pairwise dissmilarity: 0=same; 1=maximally dissimilar
pcoa_cross <- cmdscale(bray_cross, eig = T)
pcoap_cross <- data.frame(pcoa_cross$points)
pcoap_cross$ID <- rownames(pcoap_cross)
pcoap_cross <- pcoap_cross[order(pcoap_cross$ID),]

# percentages of variation explained by PCO1 & 2
eigs_cross <- pcoa_cross$eig
pc1.pct_cross <- eigs_cross[1]/sum(eigs_cross) 
pc1.pct_cross # pc1.pct_cross = 10.6%
pc2.pct_cross <- eigs_cross[2]/sum(eigs_cross) 
pc2.pct_cross # pc2.pct_cross = 6.4%

alpha.pcoa_cross <- left_join(pcoap_cross , alpha_with_disease_status, by = "ID")
alpha.pcoa_cross <- alpha.pcoa_cross %>% 
  mutate(mc_all_label = case_when(
    mc_all == 2 ~ "Healthy control",
    mc_all == 1 ~ "MC",
    mc_all == 0 ~ "Chronic diarrhea",
    TRUE ~ NA_character_))

# Modify levels for disease status
alpha.pcoa_cross$disease_status <- fct_recode(alpha.pcoa_cross$disease_status,
                                              "Active" = "active",
                                              "Remission" = "remission")

alpha.pcoa_cross$mc_all_label <- factor(alpha.pcoa_cross$mc_all_label, levels = c("Chronic diarrhea", "MC", "Healthy control"), labels = c("Chronic diarrhea", "MC", "Control without diarrhea"))

fig_PCoA_cross <- ggscatter(alpha.pcoa_cross, y = "X2", x = "X1",
                            ylab = paste0('PCo2 (', round(pc2.pct_cross * 100, 1), '%)'),
                            xlab = paste0('PCo1 (', round(pc1.pct_cross * 100, 1), '%)'),
                            font.x = 12, font.y = 12, font.tickslab = 11,
                            color = "mc_all_label", 
                            palette = c("#F9D937", "#900C3F", "#189AF9"), size = 2,
                            legend = "bottom", legend.title = "Legend Title") +
  stat_ellipse(aes(color = mc_all_label), linetype = 2) +
  scale_color_manual(values = c("Chronic diarrhea" = "#F9D937", 
                                "MC" = "#900C3F",
                                "Control without diarrhea" = "#189AF9"),
                     labels = c("Chronic diarrhea" = "Chronic diarrhea", 
                                "MC" = "MC", 
                                "Control without diarrhea" = "Control without diarrhea")) +
  guides(color = guide_legend(title = "Disease type"))


fig_PCoA_cross

# ggsave(filename = "fig_PCoA_cross.pdf",
#        plot=fig_PCoA_cross, units = "in", height = 5, width = 7, dpi=1200)

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/alpha_pcoa_cross.xlsx"
# write.xlsx(alpha.pcoa_cross, file_path, rowNames = FALSE)

####---------------  PCoA (longitudinal: Bray-Curtis)  ----------------------#####
pcoa_taxo_long <- all_taxo
pcoa_taxo_long$id_number <- substr(pcoa_taxo_long$id, 1,5)
pcoa_taxo_long <- pcoa_taxo_long %>%
  inner_join(list_id_with_both_active_and_remission,  by = c("id_number" = "id"))
rownames(pcoa_taxo_long) <- pcoa_taxo_long$id

pcoa_taxo_long <- pcoa_taxo_long[, !colnames(pcoa_taxo_long) %in% c("id", "id_number")]

bray_long <- vegdist(pcoa_taxo_long, "bray")  # pairwise dissmilarity: 0=same; 1=maximally dissimilar
pcoa_long <- cmdscale(bray_long, eig = T)
pcoap_long <- data.frame(pcoa_long$points)
pcoap_long$ID <- rownames(pcoap_long)
pcoap_long <- pcoap_long[order(pcoap_long$ID),]

# percentages of variation explained by PCO1 & 2
eigs_long <- pcoa_long$eig
pc1.pct_long <- eigs_long[1]/sum(eigs_long) 
pc1.pct_long # pc1.pct_long = 11.7%
pc2.pct_long <- eigs_long[2]/sum(eigs_long) 
pc2.pct_long # pc2.pct = 7.9%

alpha.pcoa_long <- inner_join(pcoap_long, alpha_with_disease_status, by = "ID")


# Modify levels for disease status
fig_PCoA_long <- ggscatter(alpha.pcoa_long, y = "X2", x = "X1",
                           ylab = paste0('PCo2 (', round(pc2.pct_long * 100, 1), '%)'),
                           xlab = paste0('PCo1 (', round(pc1.pct_long * 100, 1), '%)'),
                           font.x = 12, font.y = 12, font.tickslab = 11,
                           color = "disease_status", 
                           palette = c("#900C3F", "antiquewhite4"), size = 2,
                           legend = "bottom", legend.title = "Legend Title") +
  stat_ellipse(aes(color = disease_status), linetype=2) +
  guides(color = guide_legend(title = "MC activity"))

fig_PCoA_long

# ggsave(filename = "fig_PCoA_long.pdf",
#        plot=fig_PCoA_long, units = "in", height = 5, width = 7, dpi=1200)

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/alpha.pcoa_long.xlsx"
# write.xlsx(alpha.pcoa_long, file_path, rowNames = FALSE)


####---------------  PCoA (Cross-sectional: Aitchison Distance)  ----------------------#####
# CLR transformation (Aitchison distance requires log-ratio transformation)
clr_data <- clr(subset.pcoa_taxo_cross)  # Centered log-ratio transformation

# Compute Aitchison distance (Euclidean distance on CLR-transformed data)
aitchison_dist <- vegdist(clr_data, method = "euclidean")

# Perform PCoA on Aitchison distance
pcoa_cross_aitchison <- cmdscale(aitchison_dist, eig = TRUE)
pcoap_cross_aitchison <- data.frame(pcoa_cross_aitchison$points)
pcoap_cross_aitchison$ID <- rownames(pcoap_cross_aitchison)
pcoap_cross_aitchison <- pcoap_cross_aitchison[order(pcoap_cross_aitchison$ID),]

# Percentages of variation explained by PCo1 & PCo2
eigs_cross_aitchison <- pcoa_cross_aitchison$eig
pc1.pct_cross_aitchison <- eigs_cross_aitchison[1] / sum(eigs_cross_aitchison)  
pc2.pct_cross_aitchison <- eigs_cross_aitchison[2] / sum(eigs_cross_aitchison)  

# Merge with disease status data
alpha.pcoa_cross_aitchison <- left_join(pcoap_cross_aitchison, alpha_with_disease_status, by = "ID") %>%
  mutate(mc_all_label = case_when(
    mc_all == 2 ~ "Healthy control",
    mc_all == 1 ~ "MC",
    mc_all == 0 ~ "Chronic diarrhea",
    TRUE ~ NA_character_
  ))

# Modify levels for disease status
alpha.pcoa_cross_aitchison$disease_status <- fct_recode(alpha.pcoa_cross_aitchison$disease_status,
                                                        "Active" = "active",
                                                        "Remission" = "remission")

alpha.pcoa_cross_aitchison$mc_all_label <- factor(alpha.pcoa_cross_aitchison$mc_all_label, 
                                                  levels = c("Chronic diarrhea", "MC", "Healthy control"), 
                                                  labels = c("Chronic diarrhea", "MC", "Control without diarrhea"))

# PCoA Plot using Aitchison Distance
fig_PCoA_cross_aitchison <- ggscatter(alpha.pcoa_cross_aitchison, y = "X2", x = "X1",
                            ylab = paste0('Aitchison Distance: PCo2 (', round(pc2.pct_cross_aitchison * 100, 1), '%)'),
                            xlab = paste0('Aitchison Distance: PCo1 (', round(pc1.pct_cross_aitchison * 100, 1), '%)'),
                            font.x = 12, font.y = 12, font.tickslab = 11,
                            color = "mc_all_label", 
                            palette = c("#F9D937", "#900C3F", "#189AF9"), size = 2,
                            legend = "bottom", legend.title = "Legend Title") +
  stat_ellipse(aes(color = mc_all_label), linetype = 2) +
  scale_color_manual(values = c("Chronic diarrhea" = "#F9D937", 
                                "MC" = "#900C3F",
                                "Control without diarrhea" = "#189AF9"),
                     labels = c("Chronic diarrhea" = "Chronic diarrhea", 
                                "MC" = "MC", 
                                "Control without diarrhea" = "Control without diarrhea")) +
  guides(color = guide_legend(title = "Disease type"))

fig_PCoA_cross_aitchison

# ggsave(filename = "fig_PCoA_cross_aitchison.pdf",
#        plot = fig_PCoA_cross_aitchison, units = "in", height = 5, width = 7, dpi = 1200)


####---------------  PCoA (longitudinal: Aitchison Distance)  ----------------------#####
# CLR transformation (Aitchison distance requires log-ratio transformation)
clr_data_long <- clr(pcoa_taxo_long)  # Centered log-ratio transformation

# Compute Aitchison distance (Euclidean distance on CLR-transformed data)
aitchison_dist_long <- vegdist(clr_data_long, method = "euclidean")

# Perform PCoA on Aitchison distance
pcoa_long_aitchison <- cmdscale(aitchison_dist_long, eig = TRUE)
pcoap_long_aitchison <- data.frame(pcoa_long_aitchison$points)
pcoap_long_aitchison$ID <- rownames(pcoap_long_aitchison)
pcoap_long_aitchison <- pcoap_long_aitchison[order(pcoap_long_aitchison$ID),]

# Percentages of variation explained by PCo1 & PCo2
eigs_long_aitchison <- pcoa_long_aitchison$eig
pc1.pct_long_aitchison <- eigs_long_aitchison[1] / sum(eigs_long_aitchison)  
pc2.pct_long_aitchison <- eigs_long_aitchison[2] / sum(eigs_long_aitchison)  

# Merge with disease status data
alpha.pcoa_long_aitchison <- inner_join(pcoap_long_aitchison, alpha_with_disease_status, by = "ID")

# PCoA Plot using Aitchison Distance
fig_PCoA_long_aitchison <- ggscatter(alpha.pcoa_long_aitchison, y = "X2", x = "X1",
                           ylab = paste0('Aitchison Distance: PCo2 (', round(pc2.pct_long_aitchison * 100, 1), '%)'),
                           xlab = paste0('Aitchison Distance: PCo1 (', round(pc1.pct_long_aitchison * 100, 1), '%)'),
                           font.x = 12, font.y = 12, font.tickslab = 11,
                           color = "disease_status", 
                           palette = c("#900C3F", "antiquewhite4"), size = 2,
                           legend = "bottom", legend.title = "Legend Title") +
  stat_ellipse(aes(color = disease_status), linetype=2) +
  guides(color = guide_legend(title = "MC activity"))

fig_PCoA_long_aitchison

# ggsave(filename = "fig_PCoA_long_aitchison.pdf",
#        plot = fig_PCoA_long_aitchison, units = "in", height = 5, width = 7, dpi = 1200)



####---------------  PERMANOVA ---------------------#####
# Species level & metadata
taxo.rel <- subset.pcoa_taxo_cross 
taxo.rel$ID <- rownames(taxo.rel)
ID <- taxo.rel$ID
taxo.rel <- taxo.rel %>% left_join(select(alpha_with_disease_status, ID, mc_all, disease_status, age, sex, bmi))
rownames(taxo.rel) <- taxo.rel$ID
taxo.rel <- taxo.rel %>% left_join(select(matched_hc, id, race, smoke), by = c("ID" = "id"))
rownames(taxo.rel) <- taxo.rel$ID
taxo.rel$id <- substr(taxo.rel$ID, 1,5)
taxo.rel <- taxo.rel %>% left_join(select(baseline_microbiome, id, race, smoke), by = c("id" = "id"))

taxo.rel$race.x <- ifelse(is.na(taxo.rel$race.x), taxo.rel$race.y, taxo.rel$race.x)
taxo.rel$smoke.x <- ifelse(is.na(taxo.rel$smoke.x), taxo.rel$smoke.y, taxo.rel$smoke.x)

taxo.rel <- taxo.rel %>% rename("race" = "race.x",
                                "smoke" = "smoke.x") %>% select(-race.y, -smoke.y, -id)
rownames(taxo.rel) <- taxo.rel$ID
taxo.rel_1 <- taxo.rel[substr(taxo.rel$ID, 1, 1) == "G", ]
taxo.rel_2 <- taxo.rel[substr(taxo.rel$ID, 1, 2) == "MC", ]

taxo.rel_1$PID <- taxo.rel_1$ID
taxo.rel_2$PID<- substr(taxo.rel_2$ID, 1,5)

taxo.rel <- rbind(taxo.rel_1, taxo.rel_2)

taxo.rel$mc_all <- as.factor(taxo.rel$mc_all)

# Count NAs in selected columns
na_counts <- colSums(is.na(taxo.rel[, c("bmi", "smoke", "race")]))
# Create a dataframe to store the NA counts
na_counts_df <- data.frame(Columns = names(na_counts), NA_Counts = na_counts)

# Remove columns with too many NA
taxo.rel <- taxo.rel %>% select(-smoke) # remove smoke due to too many NAs, which will prevent PERMANOVA calculation.
taxo.rel_clean <- na.omit(taxo.rel) # Stool sample number: 683 --> 673

#bray_curtis dissimilarity matrix
bray <- vegdist(taxo.rel_clean[,-which(colnames(taxo.rel_clean) %in% 
                                     c("ID","disease_status","mc_all","age","sex","bmi","race","PID"))], "bray")  # pairwise dissmilarity: 0=same; 1=maximally dissimilar
bray <- as.data.frame(bray)


## PERMANOVA -- microbiome ~ covariates
perm <- how(nperm = 999) 
setBlocks(perm) <- with(taxo.rel_clean, as.factor(PID))

# x_mc_all <- adonis2(formula = bray ~ mc_all, data = taxo.rel_clean, permutations = perm) 
# mc_all <- x_mc_all["mc_all", , drop = FALSE]
# 
# x_disease_status <- adonis2(formula = bray ~ disease_status, data = taxo.rel_clean, permutations = perm) 
# disease_status <- x_disease_status["disease_status", ,drop = FALSE]
# 
# x_age <- adonis2(formula = bray ~ age, data = taxo.rel_clean, permutations = perm) 
# age <- x_age["age", , drop = FALSE]
# 
# x_sex <- adonis2(formula = bray ~ sex, data = taxo.rel_clean, permutations = perm) 
# sex <- x_sex["sex", , drop = FALSE]
# 
# x_bmi <- adonis2(formula = bray ~ bmi, data = taxo.rel_clean, permutations = perm) 
# bmi <- x_bmi["bmi", , drop = FALSE]
# 
# x_race <- adonis2(formula = bray ~ race, data = taxo.rel_clean, permutations = perm) 
# race <- x_race["race", , drop = FALSE]
# 
# all_rsq <- rbind(mc_all, disease_status, age, sex, race, bmi)
# component <- c('MC/CD/HC',
#               'Active/Remission',
#               'Age at enrollment',
#               'Sex',
#               'Race',
#               'BMI')  
# all_rsq <- cbind(all_rsq, component)
# all_rsq$pct_var<- round(all_rsq$R2*100,1)
# all_rsq$stars <- cut(all_rsq$`Pr(>F)`, breaks=c(-Inf, 0.001, 0.05, Inf), label=c("**", "*", ""))
# 
# fig_permanova <- ggbarplot(all_rsq, x="component",y="pct_var",  
#                            xlab = "",ylab= "% variation explained", x.text.angle = 45, 
#                            fill="lightblue",font.x=12, font.y=12, font.tickslab=11,
#                            label=T, lab.pos ="out", 
#                            sort.val = "desc")
# 
# fig_permanova

# ggsave(filename = "fig_permanova.png",
#        plot = fig_permanova, units = "in", width=5, height=5, dpi = 300)


####---------------  Species abundance: MC vs controls (Cross-sectional)----------------#####
# Preparing dataframes
taxo.rel_maaslin_cross <- taxo.rel %>% select(1:467)
metadata_maaslin_cross <- taxo.rel %>% select(468:475)

# maaslin_cross <- Maaslin2(input_data = taxo.rel_maaslin_cross, 
#                           input_metadata = metadata_maaslin_cross,
#                           output = "maaslin_cross",
#                           min_abundance = 0.0001,
#                           min_prevalence = 0.1,
#                           normalization = "NONE", #already relative abundance
#                           transform = "LOG",
#                           analysis_method = "LM",
#                           max_significance = 0.25,
#                           fixed_effects = c("mc_all","age","bmi","sex","disease_status"),
#                           reference = c("mc_all,1","sex,Female","disease_status,remission"))


metadata_maaslin_cross_HC <- metadata_maaslin_cross %>% filter(!mc_all == 0)
metadata_maaslin_cross_CD <- metadata_maaslin_cross %>% filter(!mc_all == 2)

# maaslin_cross_separate_HC <- Maaslin2(input_data = taxo.rel_maaslin_cross,
#                                       input_metadata = metadata_maaslin_cross_HC,
#                                       output = "maaslin_cross_separate_HC",
#                                       min_abundance = 0.0001,
#                                       min_prevalence = 0.1,
#                                       normalization = "NONE", #already relative abundance
#                                       transform = "LOG",
#                                       analysis_method = "LM",
#                                       max_significance = 0.25,
#                                       fixed_effects = c("mc_all","bmi","disease_status"),
#                                       reference = c("mc_all,1","disease_status,remission"))

# maaslin_cross_separate_CD <- Maaslin2(input_data = taxo.rel_maaslin_cross,
#                                       input_metadata = metadata_maaslin_cross_CD,
#                                       output = "maaslin_cross_separate_CD",
#                                       min_abundance = 0.0001,
#                                       min_prevalence = 0.1,
#                                       normalization = "NONE", #already relative abundance
#                                       transform = "LOG",
#                                       analysis_method = "LM",
#                                       max_significance = 0.25,
#                                       fixed_effects = c("mc_all","age","bmi","sex"),
#                                       reference = c("mc_all,1","sex,Female"))

MC_HC <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross_separate_HC/all_results.tsv")
MC_HC <- MC_HC %>% filter(metadata=="mc_all") %>% select(1,3,4,5,8,9)
MC_HC$coef <- -(MC_HC$coef)

MC_CD <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross_separate_CD/all_results.tsv")
MC_CD <- MC_CD %>% filter(metadata=="mc_all") %>% select(1,3,4,5,8,9)
MC_CD$coef <- -(MC_CD$coef)

volcano_cross <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross/all_results.tsv")
volcano_cross <- volcano_cross %>% filter(metadata=="mc_all") %>% select(1,3,4,5,8,9) 

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/volcano_cross_hc.xlsx"
# write.xlsx(volcano_cross_hc, file_path, rowNames = FALSE)

### Volcano plot for MC vs HC
volcano_cross_hc <- volcano_cross %>% filter(value==2)

# Because the reference was MC and we want to see MC compare to HC, we take negative of coefficients
volcano_cross_hc$coef <- -(volcano_cross_hc$coef)

volcano_cross_hc <- volcano_cross_hc %>% 
  mutate(diffexpressed = case_when(
    coef > 0.5 & qval < 0.25 ~ "up",
    coef < (-0.5) & qval < 0.25 ~ "down",
    TRUE ~ "no"
  ))
volcano_cross_hc <- as.data.frame(volcano_cross_hc)

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially abundant species (NA in case they are not)
volcano_cross_hc$delabel <- ifelse(volcano_cross_hc$feature %in% head(volcano_cross_hc[order(volcano_cross_hc$qval), "feature"], 30), 
                                   volcano_cross_hc$feature, 
                                   NA)

fig_volcano_cross_hc <- 
  ggplot(data = volcano_cross_hc, aes(x = coef, y = -log10(qval), label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.25), col = "gray", linetype = 'dashed') +
  geom_point(aes(color = diffexpressed), size = 2) + 
  scale_color_manual(values = c("#96DED1", "#808080", "#E35335"),
                     labels = c("Less abundant in MC", "Not significant", "More abundant in MC")) +
  labs(color = 'Change in abundance', 
       x = expression("coefficient"), y = expression("-log"[10]*"(q-value)")) + 
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-2.5,0,2.5))+ # Adjust x-axis limits and intervals
  geom_text_repel(max.overlaps = Inf, size = 3)+ # To show all labels
  theme_set(theme_classic(base_size = 12) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,5,0,0), size = rel(1.1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(5,0,0,0), size = rel(1.1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))

fig_volcano_cross_hc

# ggsave(filename = "fig_volcano_cross_hc.pdf",
#        plot = fig_volcano_cross_hc, units = "in", width=12, height=8, dpi = 1200)



### Volcano plot for MC vs CD
volcano_cross_cd <- volcano_cross %>% filter(value==0)

# Because the reference was MC and we want to see MC compare to CD, we take negative of coefficients
volcano_cross_cd$coef <- -(volcano_cross_cd$coef)

volcano_cross_cd <- volcano_cross_cd %>% 
  mutate(diffexpressed = case_when(
    coef > 0.5 & qval < 0.25 ~ "up",
    coef < (-0.5) & qval < 0.25 ~ "down",
    TRUE ~ "no"
  ))
volcano_cross_cd <- as.data.frame(volcano_cross_cd)

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/volcano_cross_cd.xlsx"
# write.xlsx(volcano_cross_cd, file_path, rowNames = FALSE)


# Create a new column "delabel" to de, that will contain the name of the top 30 differentially abundant species (NA in case they are not)
volcano_cross_cd$delabel <- ifelse(volcano_cross_cd$feature %in% head(volcano_cross_cd[order(volcano_cross_cd$qval), "feature"], 30), 
                                   volcano_cross_cd$feature, 
                                   NA)

fig_volcano_cross_cd <- 
  ggplot(data = volcano_cross_cd, aes(x = coef, y = -log10(qval), label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.25), col = "gray", linetype = 'dashed') +
  geom_point(aes(color = diffexpressed), size = 2) + 
  scale_color_manual(values = c("#96DED1", "#808080", "#E35335"),
                     labels = c("Less abundant in MC", "Not significant", "More abundant in MC")) +
  labs(color = 'Change in abundance', 
       x = expression("coefficient"), y = expression("-log"[10]*"(q-value)")) + 
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5,0,2.5))+ # Adjust x-axis limits and intervals
  geom_text_repel(max.overlaps = Inf, size = 3)+ # To show all labels
  theme_set(theme_classic(base_size = 12) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,5,0,0), size = rel(1.1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(5,0,0,0), size = rel(1.1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))

fig_volcano_cross_cd

# ggsave(filename = "fig_volcano_cross_cd.pdf",
#        plot = fig_volcano_cross_cd, units = "in", width=12, height=8, dpi = 1200)




####---------------  Find species that are different in MC compared to both controls ----------------#####
more_mc.vs.hc <- volcano_cross_hc %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Control_without_diarrhea" = "coef",
         "qval_MC.vs.Control_without_diarrhea" = "qval")

less_mc.vs.hc <- volcano_cross_hc %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Control_without_diarrhea" = "coef",
         "qval_MC.vs.Control_without_diarrhea" = "qval")

more_mc.vs.cd <- volcano_cross_cd %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Chronic_diarrhea" = "coef",
         "qval_MC.vs.Chronic_diarrhea" = "qval")

less_mc.vs.cd <- volcano_cross_cd %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Chronic_diarrhea" = "coef",
         "qval_MC.vs.Chronic_diarrhea" = "qval")

## Enriched species in both MC vs HC & MC vs CD
more_both <- inner_join(more_mc.vs.hc, more_mc.vs.cd)
more_both <- more_both %>% mutate(new_species = case_when(
  feature == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
  TRUE ~ feature
))

## Depleted species in both MC vs HC & MC vs CD
less_both <- inner_join(less_mc.vs.hc, less_mc.vs.cd)
less_both <- less_both %>% mutate(new_species = case_when(
  feature == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
  feature == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
  feature == "GGB80011_SGB15265" ~ "p_Firmicutes_GGB80011_SGB15265",
  feature == "GGB9524_SGB14923" ~ "p_Firmicutes_GGB9524_SGB14923",
  TRUE ~ feature
))


### Separate models
more_MC_HC <- MC_HC %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Control_without_diarrhea" = "coef",
         "qval_MC.vs.Control_without_diarrhea" = "qval")

less_MC_HC <- MC_HC %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Control_without_diarrhea" = "coef",
         "qval_MC.vs.Control_without_diarrhea" = "qval")

more_MC_CD <- MC_CD %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Chronic_diarrhea" = "coef",
         "qval_MC.vs.Chronic_diarrhea" = "qval")

less_MC_CD <- MC_CD %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Chronic_diarrhea" = "coef",
         "qval_MC.vs.Chronic_diarrhea" = "qval")

## Enriched species in both MC vs HC & MC vs CD
more_both_spearate <- inner_join(more_MC_HC, more_MC_CD)
more_both_spearate  <- more_both_spearate %>% mutate(new_species = case_when(
  feature == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
  TRUE ~ feature
))

## Depleted species in both MC vs HC & MC vs CD
less_both_spearate <- inner_join(less_MC_HC, less_MC_CD)
less_both_spearate <- less_both_spearate %>% mutate(new_species = case_when(
  feature == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
  feature == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
  feature == "GGB80011_SGB15265" ~ "p_Firmicutes_GGB80011_SGB15265",
  feature == "GGB9524_SGB14923" ~ "p_Firmicutes_GGB9524_SGB14923",
  TRUE ~ feature
))





## Heatmap
more_heat <- more_both %>% select(-feature)
less_heat <- less_both %>% select(-feature)
altered_heat <- rbind(more_heat, less_heat)

rownames(altered_heat) <- altered_heat$new_species
altered_heat$new_species <- NULL
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
breaks <- seq(-2.6, 2.6, length.out = 51)

# Ensure the breaks include 0
mid <- which.min(abs(breaks - 0))

# Create a custom color palette
my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))

# Plot the heatmap with reordered rows
heatmap_altered <- pheatmap(b_altered,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            display_numbers = asterisks_matrix_altered,
                            fontsize_number = 10,
                            color = my_color_palette,
                            breaks = breaks,
                            annotation_legend = TRUE,
                            angle_col = 90)

# ggsave("heatmap_altered.pdf", heatmap_altered, width = 4.5, height = 11, dpi=1200, limitsize = FALSE)





### Box plots
## Enriched species in both MC vs HC & MC vs CD
box_more_both <- more_both$feature
species_box_more_both <- taxo.rel_maaslin_cross[, box_more_both, drop = FALSE]
species_box_more_both$id <- rownames(species_box_more_both)
species_box_more_both <- as.data.frame(species_box_more_both)
species_box_more_both <- species_box_more_both %>%
  left_join(select(metadata_maaslin_cross, ID, mc_all), by = c("id" = "ID")) %>%
  select(id, mc_all, 1:8)
rownames(species_box_more_both) <- species_box_more_both$id


species_box_more_both_long <- pivot_longer(species_box_more_both,
                                    cols = 3:10,
                                    names_to = "species",
                                    values_to = "relative_abundance")


species_box_more_both_long <- species_box_more_both_long %>%
  mutate(new_species = case_when(
    species == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
    TRUE ~ species
  ))

more_both_order <- c("Intestinibacter_bartlettii", "Veillonella_dispar", "Veillonella_parvula",
                     "Haemophilus_parainfluenzae", "Clostridium_sp_AF20_17LB", "Clostridium_spiroforme",
                     "Veillonella_rogosae", "p_Firmicutes_GGB3612_SGB4882")

# Arrange the data frame by the custom order
species_box_more_both_long <- species_box_more_both_long %>%
  arrange(factor(new_species, levels = more_both_order))

# Box plot
more_both_qval.to.plot <- more_both %>%
  mutate(qval.format = round(qval_MC.vs.Control_without_diarrhea, 3)) %>%
  mutate(group1 = "0") %>%
  mutate(group2 = "1") %>%
  mutate(qval_text = ifelse(qval_MC.vs.Control_without_diarrhea < 0.001, "q < 0.001", paste0("q = ", qval.format)))

more_both_qval.to.plot <- more_both_qval.to.plot %>%
  mutate(qval.format.cd = round(qval_MC.vs.Chronic_diarrhea, 3)) %>%
  mutate(group1 = "0") %>%
  mutate(group2 = "1") %>%
  mutate(qval_text_cd = ifelse(qval_MC.vs.Chronic_diarrhea < 0.001, "q < 0.001", paste0("q = ", qval.format.cd)))

more_both_qval.to.plot <- more_both_qval.to.plot %>%
  mutate(new_species = case_when(
    feature == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
    TRUE ~ feature
  ))

species_box_more_both_long <- species_box_more_both_long %>% mutate(
  mc_all_label = case_when(
    mc_all == 2 ~ "Healthy_control",
    mc_all == 1 ~ "MC",
    mc_all == 0 ~ "Chronic_diarrhea"
  )
)

species_box_more_both_long$mc_all_label <- factor(species_box_more_both_long$mc_all_label , levels = c("Healthy_control", "MC", "Chronic_diarrhea"))

# Calculate half of the minimum non-zero value
half_min_value <- min(species_box_more_both_long$relative_abundance[species_box_more_both_long$relative_abundance > 0]) / 2

# Replace zeros with the calculated half minimum value
species_box_more_both_long <- species_box_more_both_long %>%
  mutate(relative_abundance = ifelse(relative_abundance == 0, half_min_value, relative_abundance))

# Take log for relative_abundance
species_box_more_both_long$log_relative_abundance <- log(species_box_more_both_long$relative_abundance)

# # Violin plot for all enriched species
# violin.both_more <- ggplot(species_box_more_both_long, aes(factor(new_species, levels = more_both_order), y = log_relative_abundance)) +
#   geom_violin(aes(fill = mc_all_label), trim = FALSE, position = position_dodge(width = 0.8))+
#   geom_jitter(aes(fill = mc_all_label),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8), size = 1, alpha = 0.4)+
#   geom_text(data = more_both_qval.to.plot, aes(x = factor(new_species, levels = more_both_order), y = 5.4, label = paste("MC vs HC: ", qval_text)), vjust = 0, hjust = 0.5, size = 3, color = "black") +
#   geom_text(data = more_both_qval.to.plot, aes(x = factor(new_species, levels = more_both_order), y = 5, label = paste("MC vs CD: ", qval_text_cd)), vjust = 0, hjust = 0.5, size = 3, color = "black") +
#   labs(x = "Species", y = "log(relative abundance(%))", fill = "Disease type") +
#   scale_fill_manual(values = c("MC" = "#900C3F", "Chronic_diarrhea" = "#F9D937", "Healthy_control" = "#189AF9")) +
#   theme_minimal() +
#   theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1), axis.text.y = element_text(size = 8))
# 
# violin.both_more
# # ggsave(filename = "violin.both_more.png",
# #        plot = violin.both_more, units = "in", width=12, height=12, dpi = 300)

## Individual species
# V.dispar
v.dispar <- species_box_more_both_long %>% filter(new_species == "Veillonella_dispar")
v.dispar_q <- more_both_qval.to.plot %>% filter(feature=="Veillonella_dispar")

box_v.dispar <- ggplot(v.dispar, aes(factor(new_species), y = log_relative_abundance)) + 
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 1), width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = mc_all_label), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), size = 1, alpha = 0.5,show.legend = FALSE) +
  geom_signif(
    y_position = c(5, 5), xmin = c(0.65, 1.05), xmax = c(0.95, 1.35),
    annotation = c("q < 0.001", "q = 0.029"), tip_length = 0
  ) +
  labs(x = "", y = "log(relative abundance(%))", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", 
                               "Chronic_diarrhea" = "#F9D937", 
                               "Healthy_control" = "#189AF9"),
                    labels = c("MC" = "MC", 
                               "Chronic_diarrhea" = "Chronic diarrhea", 
                               "Healthy_control" = "Control without diarrhea")) +
  scale_color_manual(values = c("MC" = "#900C3F", 
                                "Chronic_diarrhea" = "#F9D937", 
                                "Healthy_control" = "#189AF9"),
                     labels = c("MC" = "MC", 
                                "Chronic_diarrhea" = "Chronic diarrhea", 
                                "Healthy_control" = "Control without diarrhea")) +
  theme_minimal() +
  theme(legend.position = "top", 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 8))

box_v.dispar


# ggsave(filename = "box_v.dispar.pdf",
#        plot = box_v.dispar, units = "in", width=6, height=6, dpi = 1200)

# V.parvula
v.parvula <- species_box_more_both_long %>% filter(new_species == "Veillonella_parvula")
v.parvula_q <- more_both_qval.to.plot %>% filter(feature=="Veillonella_parvula")

box_v.parvula <- ggplot(v.parvula, aes(factor(new_species), y = log_relative_abundance)) + 
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 1), width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = mc_all_label), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), size = 1, alpha = 0.5,show.legend = FALSE) +
  geom_signif(
    y_position = c(5, 5), xmin = c(0.65, 1.05), xmax = c(0.95, 1.35),
    annotation = c("q < 0.001", "q = 0.009"), tip_length = 0
  ) +
  labs(x = "", y = "log(relative abundance(%))", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", 
                               "Chronic_diarrhea" = "#F9D937", 
                               "Healthy_control" = "#189AF9"),
                    labels = c("MC" = "MC", 
                               "Chronic_diarrhea" = "Chronic diarrhea", 
                               "Healthy_control" = "Control without diarrhea")) +
  scale_color_manual(values = c("MC" = "#900C3F", 
                                "Chronic_diarrhea" = "#F9D937", 
                                "Healthy_control" = "#189AF9"),
                     labels = c("MC" = "MC", 
                                "Chronic_diarrhea" = "Chronic diarrhea", 
                                "Healthy_control" = "Control without diarrhea")) +
  theme_minimal() +
  theme(legend.position = "top", 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 8))

box_v.parvula

# ggsave(filename = "box_v.parvula.pdf",
#        plot = box_v.parvula, units = "in", width=6, height=6, dpi = 1200)

# H.parainfluenzae
h.parainfluenzae <- species_box_more_both_long %>% filter(new_species == "Haemophilus_parainfluenzae")
h.parainfluenzae_q <- more_both_qval.to.plot %>% filter(feature=="Haemophilus_parainfluenzae")

box_h.parainfluenzae <- ggplot(h.parainfluenzae, aes(factor(new_species), y = log_relative_abundance)) + 
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 1), width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = mc_all_label), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), size = 1, alpha = 0.5,show.legend = FALSE) +
  geom_signif(
    y_position = c(5, 5), xmin = c(0.65, 1.05), xmax = c(0.95, 1.35),
    annotation = c("q < 0.001", "q < 0.001"), tip_length = 0
  ) +
  labs(x = "", y = "log(relative abundance(%))", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", 
                               "Chronic_diarrhea" = "#F9D937", 
                               "Healthy_control" = "#189AF9"),
                    labels = c("MC" = "MC", 
                               "Chronic_diarrhea" = "Chronic diarrhea", 
                               "Healthy_control" = "Control without diarrhea")) +
  scale_color_manual(values = c("MC" = "#900C3F", 
                                "Chronic_diarrhea" = "#F9D937", 
                                "Healthy_control" = "#189AF9"),
                     labels = c("MC" = "MC", 
                                "Chronic_diarrhea" = "Chronic diarrhea", 
                                "Healthy_control" = "Control without diarrhea")) +
  theme_minimal() +
  theme(legend.position = "top", 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 8))

box_h.parainfluenzae

# ggsave(filename = "box_h.parainfluenzae.pdf",
#        plot = box_h.parainfluenzae, units = "in", width=6, height=6, dpi = 1200)


## Depleted species in both MC vs HC & MC vs CD
box_less_both <- less_both$feature
species_box_less_both <- taxo.rel_maaslin_cross[, box_less_both, drop = FALSE]
species_box_less_both$id <- rownames(species_box_less_both)
species_box_less_both <- as.data.frame(species_box_less_both)
species_box_less_both <- species_box_less_both %>%
  left_join(select(metadata_maaslin_cross, ID, mc_all), by = c("id" = "ID")) %>%
  select(id, mc_all, 1:11)
rownames(species_box_less_both) <- species_box_less_both$id



species_box_less_both_long <- pivot_longer(species_box_less_both,
                                           cols = 3:13,
                                           names_to = "species",
                                           values_to = "relative_abundance")


species_box_less_both_long <- species_box_less_both_long %>%
  mutate(new_species = case_when(
    species == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
    species == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
    species == "GGB80011_SGB15265" ~ "p_Firmicutes_GGB80011_SGB15265",
    species == "GGB9524_SGB14923" ~ "p_Firmicutes_GGB9524_SGB14923",
    TRUE ~ species
  ))

less_both_order <- c("Methylobacterium_SGB15164", "f_Clostridia_unclassified_GGB33586_SGB53517", "Mediterraneibacter_butyricigenes",
                     "f_Clostridia_unclassified_GGB9534_SGB14937", "Clostridiales_bacterium", "Collinsella_SGB4121",
                     "p_Firmicutes_GGB80011_SGB15265", "Bacteroides_stercoris", "Blautia_glucerasea", "p_Firmicutes_GGB9524_SGB14923", "Clostridiales_bacterium_Choco116")

# Arrange the data frame by the custom order
species_box_less_both_long <- species_box_less_both_long %>%
  arrange(factor(new_species, levels = less_both_order))

# Box plot
less_both_qval.to.plot <- less_both %>%
  mutate(qval.format = round(qval_MC.vs.Control_without_diarrhea, 3)) %>%
  mutate(group1 = "0") %>%
  mutate(group2 = "1") %>%
  mutate(qval_text = ifelse(qval_MC.vs.Control_without_diarrhea < 0.001, "q < 0.001", paste0("q = ", qval.format)))

less_both_qval.to.plot <- less_both_qval.to.plot %>%
  mutate(qval.format.cd = round(qval_MC.vs.Chronic_diarrhea, 3)) %>%
  mutate(group1 = "0") %>%
  mutate(group2 = "1") %>%
  mutate(qval_text_cd = ifelse(qval_MC.vs.Chronic_diarrhea < 0.001, "q < 0.001", paste0("q = ", qval.format.cd)))

less_both_qval.to.plot <- less_both_qval.to.plot %>%
  mutate(new_species = case_when(
    feature == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
    feature == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
    feature == "GGB80011_SGB15265" ~ "p_Firmicutes_GGB80011_SGB15265",
    feature == "GGB9524_SGB14923" ~ "p_Firmicutes_GGB9524_SGB14923",
    TRUE ~ feature
  ))

species_box_less_both_long <- species_box_less_both_long %>% mutate(
  mc_all_label = case_when(
    mc_all == 2 ~ "Healthy_control",
    mc_all == 1 ~ "MC",
    mc_all == 0 ~ "Chronic_diarrhea"
  )
)


species_box_less_both_long$mc_all_label <- factor(species_box_less_both_long$mc_all_label , levels = c("Healthy_control", "MC", "Chronic_diarrhea"))

# Calculate half of the minimum non-zero value
half_min_value_less <- min(species_box_less_both_long$relative_abundance[species_box_less_both_long$relative_abundance > 0]) / 2

# Replace zeros with the calculated half minimum value
species_box_less_both_long <- species_box_less_both_long %>%
  mutate(relative_abundance = ifelse(relative_abundance == 0, half_min_value_less, relative_abundance))

# Take log for relative_abundance
species_box_less_both_long$log_relative_abundance <- log(species_box_less_both_long$relative_abundance)


# Methylobacterium_SGB15164
methylobacterium <- species_box_less_both_long %>% filter(new_species == "Methylobacterium_SGB15164")
methylobacterium_q <- less_both_qval.to.plot %>% filter(feature=="Methylobacterium_SGB15164")

box_methylobacterium <- ggplot(methylobacterium, aes(factor(new_species), y = log_relative_abundance)) + 
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 1), width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(fill = mc_all_label),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), size = 1, alpha = 0.4)+
  geom_signif(
    y_position = c(5, 5), xmin = c(0.65, 1.05), xmax = c(0.95, 1.35),
    annotation = c("q < 0.001", "q = 0.21"), tip_length = 0
  ) +
  labs(x = "", y = "log(relative abundance(%))", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", 
                               "Chronic_diarrhea" = "#F9D937", 
                               "Healthy_control" = "#189AF9"),
                    labels = c("MC" = "MC", 
                               "Chronic_diarrhea" = "Chronic diarrhea", 
                               "Healthy_control" = "Control without diarrhea")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 8))
box_methylobacterium

# ggsave(filename = "box_methylobacterium.pdf",
#        plot = box_methylobacterium, units = "in", width=6, height=6, dpi = 1200)

# Mediterraneibacter_butyricigenes
m.butyricigenes <- species_box_less_both_long %>% filter(new_species == "Mediterraneibacter_butyricigenes")
m.butyricigenes_q <- less_both_qval.to.plot %>% filter(feature=="Mediterraneibacter_butyricigenes")

box_m.butyricigenes <- ggplot(m.butyricigenes, aes(factor(new_species), y = log_relative_abundance)) + 
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 1), width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(fill = mc_all_label),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), size = 1, alpha = 0.4)+
  geom_signif(
    y_position = c(5, 5), xmin = c(0.65, 1.05), xmax = c(0.95, 1.35),
    annotation = c("q < 0.001", "q = 0.176"), tip_length = 0
  ) +
  labs(x = "", y = "log(relative abundance(%))", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", 
                               "Chronic_diarrhea" = "#F9D937", 
                               "Healthy_control" = "#189AF9"),
                    labels = c("MC" = "MC", 
                               "Chronic_diarrhea" = "Chronic diarrhea", 
                               "Healthy_control" = "Control without diarrhea")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 8))
box_m.butyricigenes

# ggsave(filename = "box_m.butyricigenes.pdf",
#        plot = box_m.butyricigenes, units = "in", width=6, height=6, dpi = 1200)

# Bacteroides_stercoris
b.stercoris <- species_box_less_both_long %>% filter(new_species == "Bacteroides_stercoris")
b.stercoris_q <- less_both_qval.to.plot %>% filter(feature=="Bacteroides_stercoris")

box_b.stercoris <- ggplot(b.stercoris, aes(factor(new_species), y = log_relative_abundance)) + 
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 1), width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(fill = mc_all_label),position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1), size = 1, alpha = 0.4)+
  geom_signif(
    y_position = c(5, 5), xmin = c(0.65, 1.05), xmax = c(0.95, 1.35),
    annotation = c("q = 0.111", "q = 0.212"), tip_length = 0
  ) +
  labs(x = "", y = "log(relative abundance(%))", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", 
                               "Chronic_diarrhea" = "#F9D937", 
                               "Healthy_control" = "#189AF9"),
                    labels = c("MC" = "MC", 
                               "Chronic_diarrhea" = "Chronic diarrhea", 
                               "Healthy_control" = "Control without diarrhea")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 8))
box_b.stercoris

# ggsave(filename = "box_b.stercoris.pdf",
#        plot = box_b.stercoris, units = "in", width=6, height=6, dpi = 1200)



bx.both_more <- ggplot(species_box_more_both_long, aes(x = factor(new_species, levels = more_both_order), y = relative_abundance)) +
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 0.8), width = 0.5) +
  geom_text(data = more_both_qval.to.plot, aes(x = factor(new_species, levels = more_both_order), y = 1.05, label = paste("MC vs HC: ", qval_text)), vjust = 0, hjust = 0.5, size = 3, color = "black") +
  geom_text(data = more_both_qval.to.plot, aes(x = factor(new_species, levels = more_both_order), y = 1.03, label = paste("MC vs CD: ", qval_text_cd)), vjust = 0, hjust = 0.5, size = 3, color = "black") +
  labs(x = "Species", y = "Relative abundance (%)", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", "Chronic_diarrhea" = "#F9D937", "Healthy_control" = "#189AF9")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1), axis.text.y = element_text(size = 8)) +
  ylim(0, 1.05) 

bx.both_less <- ggplot(species_box_less_both_long, aes(x = factor(new_species, levels = less_both_order), y = relative_abundance)) +
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 0.8), width = 0.5) +
  geom_text(data = less_both_qval.to.plot, aes(x = factor(new_species, levels = less_both_order), y = 1.05, label = paste("MC vs HC: ", qval_text)), vjust = 0, hjust = 0.5, size = 3, color = "black") +
  geom_text(data = less_both_qval.to.plot, aes(x = factor(new_species, levels = less_both_order), y = 1.03, label = paste("MC vs CD: ", qval_text_cd)), vjust = 0, hjust = 0.5, size = 3, color = "black") +
  labs(x = "Species", y = "Relative abundance (%)", fill = "Disease type") +
  scale_fill_manual(values = c("MC" = "#900C3F", "Chronic_diarrhea" = "#F9D937", "Healthy_control" = "#189AF9")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1), axis.text.y = element_text(size = 8)) +
  ylim(0, 1.05) 

bx.both_less

# ggsave(filename = "bx.both_less.png",
#        plot = bx.both_less, units = "in", width=16, height=12, dpi = 600)




####---------------  HAllA for species ----------------#####
MC <- metadata_maaslin_cross %>% filter(mc_all==1)
CD <- metadata_maaslin_cross %>% filter(mc_all==0)
HC <- metadata_maaslin_cross %>% filter(mc_all==2)

more <- colnames(species_box_more_both[3:10])
less <- colnames(species_box_less_both[3:13])
more_or_less <- union(more, less)

load("match_2024-05-13.RData") 
mc_hc$mc_all <- as.factor(mc_hc$mc_all)
mc_cd_hc <- full_join(mc_cd, mc_hc)

# Residuals from maaslin_cross
res_maaslin_cross <- readRDS("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross/fits/residuals.rds") %>% 
  t() %>% as.data.frame()

res_maaslin_cross_1 <- res_maaslin_cross %>% select(all_of(more_or_less))
res_maaslin_cross_1$ID <- rownames(res_maaslin_cross_1)

res_maaslin_cross_1$id <- substr(res_maaslin_cross_1$ID, 1, 5)

## Altered species compared to both controls
mgx_halla_residual <- res_maaslin_cross_1 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(mgx_halla_residual) <- mgx_halla_residual$ID
mgx_halla_residual$ID <- NULL
mgx_halla_residual$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mgx_halla_residual.xlsx"
# write.xlsx(mgx_halla_residual, file_path, rowNames = TRUE)

# mgx_MC_halla <- taxo.rel_maaslin_cross_1 %>% inner_join(select(MC, ID))
# rownames(mgx_MC_halla) <- mgx_MC_halla$ID
# mgx_MC_halla$ID <- NULL
# # file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mgx_MC_halla.xlsx"
# # write.xlsx(mgx_MC_halla, file_path, rowNames = TRUE)
# 
# mgx_CD_halla <- taxo.rel_maaslin_cross_1 %>% inner_join(select(CD, ID))
# rownames(mgx_CD_halla) <- mgx_CD_halla$ID
# mgx_CD_halla$ID <- NULL
# # file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mgx_CD_halla.xlsx"
# # write.xlsx(mgx_CD_halla, file_path, rowNames = TRUE)
# 
# mgx_HC_halla <- taxo.rel_maaslin_cross_1 %>% inner_join(select(HC, ID))
# rownames(mgx_HC_halla) <- mgx_HC_halla$ID
# mgx_HC_halla$ID <- NULL
# # file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mgx_HC_halla.xlsx"
# # write.xlsx(mgx_HC_halla, file_path, rowNames = TRUE)



mgx_halla_residual_full <- res_maaslin_cross
matching_ids <- rownames(mgx_halla_residual)
mgx_halla_residual_full$ID <- rownames(mgx_halla_residual_full)
mgx_halla_residual_full <- mgx_halla_residual_full[mgx_halla_residual_full$ID %in% matching_ids, ]
mgx_halla_residual_full$ID <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mgx_halla_residual_full.xlsx"
# write.xlsx(mgx_halla_residual_full, file_path, rowNames = TRUE)



## Species correlated with sphingolipids (rho > 0.2 or < -0.2)
species_halla_sphingolipids_0.2 <- read_excel("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/sphingolipids_0.2_bugs.xlsx") %>%
  select(2) %>% rename("species" = "X_features")

species_halla_sphingolipids_0.2 <- species_halla_sphingolipids_0.2$species

taxo.rel_maaslin_cross_2 <- taxo.rel_maaslin_cross %>% select(all_of(species_halla_sphingolipids_0.2))
taxo.rel_maaslin_cross_2$ID <- rownames(taxo.rel_maaslin_cross_2)
taxo.rel_maaslin_cross_2$id <- substr(taxo.rel_maaslin_cross_2$ID, 1, 5)

species_halla_sphingolipids_0.2 <- taxo.rel_maaslin_cross_2 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(species_halla_sphingolipids_0.2) <- species_halla_sphingolipids_0.2$ID
species_halla_sphingolipids_0.2$ID <- NULL
species_halla_sphingolipids_0.2$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/species_halla_sphingolipids_0.2.xlsx"
# write.xlsx(species_halla_sphingolipids_0.2, file_path, rowNames = TRUE)



## Species correlated with lysophospholipids (rho > 0.2 or < -0.2)
species_halla_lysophospholipids_0.2 <- read_excel("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysophospholipids_0.2_bugs.xlsx") %>%
  select(2) %>% rename("species" = "X_features")

species_halla_lysophospholipids_0.2 <- species_halla_lysophospholipids_0.2$species

taxo.rel_maaslin_cross_3 <- taxo.rel_maaslin_cross %>% select(all_of(species_halla_lysophospholipids_0.2))
taxo.rel_maaslin_cross_3$ID <- rownames(taxo.rel_maaslin_cross_3)
taxo.rel_maaslin_cross_3$id <- substr(taxo.rel_maaslin_cross_3$ID, 1, 5)

species_halla_lysophospholipids_0.2 <- taxo.rel_maaslin_cross_3 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(species_halla_lysophospholipids_0.2) <- species_halla_lysophospholipids_0.2$ID
species_halla_lysophospholipids_0.2$ID <- NULL
species_halla_lysophospholipids_0.2$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/species_halla_lysophospholipids_0.2.xlsx"
# write.xlsx(species_halla_lysophospholipids_0.2, file_path, rowNames = TRUE)



## Species correlated with lysoplasmalogens (rho > 0.2 or < -0.2)
species_halla_lysoplasmalogens_0.2 <- read_excel("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysoplasmalogens_0.2_bugs.xlsx") %>%
  select(2) %>% rename("species" = "X_features")

species_halla_lysoplasmalogens_0.2 <- species_halla_lysoplasmalogens_0.2$species

taxo.rel_maaslin_cross_4 <- taxo.rel_maaslin_cross %>% select(all_of(species_halla_lysoplasmalogens_0.2))
taxo.rel_maaslin_cross_4$ID <- rownames(taxo.rel_maaslin_cross_4)
taxo.rel_maaslin_cross_4$id <- substr(taxo.rel_maaslin_cross_4$ID, 1, 5)

species_halla_lysoplasmalogens_0.2 <- taxo.rel_maaslin_cross_4 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(species_halla_lysoplasmalogens_0.2) <- species_halla_lysoplasmalogens_0.2$ID
species_halla_lysoplasmalogens_0.2$ID <- NULL
species_halla_lysoplasmalogens_0.2$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/species_halla_lysoplasmalogens_0.2.xlsx"
# write.xlsx(species_halla_lysoplasmalogens_0.2, file_path, rowNames = TRUE)



####---------------  Species abundance: LC vs controls (Cross-sectional)----------------#####
# Preparing dataframes
taxo.rel_LC <- taxo.rel %>% inner_join(select(alpha.cross.sectional_active_LC_CC, ID, mc_or_diarrhea)) %>% 
  filter(mc_or_diarrhea == "LC" | mc_or_diarrhea == "HC" | mc_or_diarrhea == "CD")
taxo.rel_LC <- taxo.rel_LC %>% as.data.frame()
rownames(taxo.rel_LC) <- taxo.rel_LC$ID

taxo.rel_maaslin_cross_LC <- taxo.rel_LC %>% select(1:467) 
taxo.rel_maaslin_cross_LC <- taxo.rel_maaslin_cross_LC[, more_or_less]
metadata_maaslin_cross_LC <- taxo.rel_LC %>% select(468:477)

# maaslin_cross_LC <- Maaslin2(input_data = taxo.rel_maaslin_cross_LC,
#                              input_metadata = metadata_maaslin_cross_LC,
#                              output = "maaslin_cross_LC",
#                              min_abundance = 0,
#                              min_prevalence = 0,
#                              normalization = "NONE", #already relative abundance
#                              transform = "LOG",
#                              analysis_method = "LM",
#                              max_significance = 0.25,
#                              fixed_effects = c("mc_or_diarrhea","age","bmi","sex","disease_status"),
#                              reference = c("mc_or_diarrhea,LC","sex,Female","disease_status,remission"))


volcano_cross_LC <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross_LC/all_results.tsv")
volcano_cross_LC <- volcano_cross_LC %>% filter(metadata=="mc_or_diarrhea") %>% select(1,3,4,5,8,9) 

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/volcano_cross_LC.xlsx"
# write.xlsx(volcano_cross_LC, file_path, rowNames = FALSE)

### Volcano plot for LC vs HC
volcano_cross_LC_hc <- volcano_cross_LC %>% filter(value=="HC")

# Because the reference was LC and we want to see LC compare to HC, we take negative of coefficients
volcano_cross_LC_hc$coef <- -(volcano_cross_LC_hc$coef)

### Volcano plot for LC vs CD
volcano_cross_LC_cd <- volcano_cross_LC %>% filter(value=="CD")

# Because the reference was LC and we want to see LC compare to CD, we take negative of coefficients
volcano_cross_LC_cd$coef <- -(volcano_cross_LC_cd$coef)


### Find species that are different in LC compared to both controls
more_LC.vs.hc <- volcano_cross_LC_hc %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_LC.vs.control_without_diarrhea" = "coef",
         "qval_LC.vs.control_without_diarrhea" = "qval")

less_LC.vs.hc <- volcano_cross_LC_hc %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_LC.vs.control_without_diarrhea" = "coef",
         "qval_LC.vs.control_without_diarrhea" = "qval")

more_LC.vs.cd <- volcano_cross_LC_cd %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_LC.vs.cd" = "coef",
         "qval_LC.vs.cd" = "qval")

less_LC.vs.cd <- volcano_cross_LC_cd %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_LC.vs.cd" = "coef",
         "qval_LC.vs.cd" = "qval")

## Enriched species in LC vs both controls
more_both_LC <- inner_join(more_LC.vs.hc, more_LC.vs.cd)

## Depleted species in LC vs both controls
less_both_LC <- inner_join(less_LC.vs.hc, less_LC.vs.cd)



####---------------  Species abundance: CC vs controls (Cross-sectional)----------------#####
# Preparing dataframes
taxo.rel_CC <- taxo.rel %>% inner_join(select(alpha.cross.sectional_active_LC_CC, ID, mc_or_diarrhea)) %>% 
  filter(mc_or_diarrhea == "CC" | mc_or_diarrhea == "HC" | mc_or_diarrhea == "CD")
taxo.rel_CC <- taxo.rel_CC %>% as.data.frame()
rownames(taxo.rel_CC) <- taxo.rel_CC$ID

taxo.rel_maaslin_cross_CC <- taxo.rel_CC %>% select(1:467)
taxo.rel_maaslin_cross_CC <- taxo.rel_maaslin_cross_CC[, more_or_less]
metadata_maaslin_cross_CC <- taxo.rel_CC %>% select(468:477)

# maaslin_cross_CC <- Maaslin2(input_data = taxo.rel_maaslin_cross_CC,
#                              input_metadata = metadata_maaslin_cross_CC,
#                              output = "maaslin_cross_CC",
#                              min_abundance = 0,
#                              min_prevalence = 0,
#                              normalization = "NONE", #already relative abundance
#                              transform = "LOG",
#                              analysis_method = "LM",
#                              max_significance = 0.25,
#                              fixed_effects = c("mc_or_diarrhea","age","bmi","sex","disease_status"),
#                              reference = c("mc_or_diarrhea,CC","sex,Female","disease_status,remission"))


volcano_cross_CC <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross_CC/all_results.tsv")
volcano_cross_CC <- volcano_cross_CC %>% filter(metadata=="mc_or_diarrhea") %>% select(1,3,4,5,8,9) 

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/volcano_cross_CC.xlsx"
# write.xlsx(volcano_cross_CC, file_path, rowNames = FALSE)

### Volcano plot for CC vs HC
volcano_cross_CC_hc <- volcano_cross_CC %>% filter(value=="HC")
volcano_cross_CC_hc$coef <- -(volcano_cross_CC_hc$coef)

volcano_cross_CC_cd <- volcano_cross_CC %>% filter(value=="CD")
volcano_cross_CC_cd$coef <- -(volcano_cross_CC_cd$coef)


### Find species that are different in CC compared to both controls
more_CC.vs.hc <- volcano_cross_CC_hc %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_CC.vs.control_without_diarrhea" = "coef",
         "qval_CC.vs.control_without_diarrhea" = "qval")

less_CC.vs.hc <- volcano_cross_CC_hc %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_CC.vs.control_without_diarrhea" = "coef",
         "qval_CC.vs.control_without_diarrhea" = "qval")

more_CC.vs.cd <- volcano_cross_CC_cd %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_CC.vs.cd" = "coef",
         "qval_CC.vs.cd" = "qval")

less_CC.vs.cd <- volcano_cross_CC_cd %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_CC.vs.cd" = "coef",
         "qval_CC.vs.cd" = "qval")

## Enriched species in CC vs both controls
more_both_CC <- inner_join(more_CC.vs.hc, more_CC.vs.cd)

## Depleted species in CC vs both controls
less_both_CC <- inner_join(less_CC.vs.hc, less_CC.vs.cd)


## Heatmap for LC and CC
altered_heat_LC <- volcano_cross_LC_hc %>% 
  select(1,3,6) %>% 
  rename("coef_LC.vs.Control_without_diarrhea" = "coef",
         "qval_LC.vs.Control_without_diarrhea" = "qval") %>% 
  left_join(select(volcano_cross_LC_cd, feature, coef, qval)) %>% 
  rename("coef_LC.vs.Chronic_diarrhea" = "coef",
         "qval_LC.vs.Chronic_diarrhea" = "qval")

altered_heat_CC <- volcano_cross_CC_hc %>% 
  select(1,3,6) %>% 
  rename("coef_CC.vs.Control_without_diarrhea" = "coef",
         "qval_CC.vs.Control_without_diarrhea" = "qval") %>% 
  left_join(select(volcano_cross_CC_cd, feature, coef, qval)) %>% 
  rename("coef_CC.vs.Chronic_diarrhea" = "coef",
         "qval_CC.vs.Chronic_diarrhea" = "qval")

altered_heat_LC_CC <- altered_heat_LC %>% left_join(altered_heat_CC) %>% as.data.frame()
altered_heat_LC_CC <- altered_heat_LC_CC %>% mutate(new_species = case_when(
  feature == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
  feature == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
  feature == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
  feature == "GGB80011_SGB15265" ~ "p_Firmicutes_GGB80011_SGB15265",
  feature == "GGB9524_SGB14923" ~ "p_Firmicutes_GGB9524_SGB14923",
  TRUE ~ feature
))

altered_heat$new_species <- rownames(altered_heat)
altered_heat_LC_CC <- altered_heat_LC_CC %>% right_join(altered_heat)


rownames(altered_heat_LC_CC) <- altered_heat_LC_CC$new_species
altered_heat_LC_CC <- altered_heat_LC_CC %>% select(-feature, -new_species)

b_altered_LC_CC <- altered_heat_LC_CC %>% select(9,11,1,3,5,7)
b_altered_LC_CC <- b_altered_LC_CC[custom_order,] %>% as.matrix()

b_altered_LC_CC_q <- altered_heat_LC_CC %>% select(10,12,2,4,6,8) 
b_altered_LC_CC_q <- b_altered_LC_CC_q[custom_order,] %>% as.matrix()

asterisks_matrix_altered_LC_CC <- matrix("", nrow = nrow(b_altered_LC_CC_q), ncol = ncol(b_altered_LC_CC_q))
asterisks_matrix_altered_LC_CC[b_altered_LC_CC_q < 0.05] <- "***"
asterisks_matrix_altered_LC_CC[b_altered_LC_CC_q >= 0.05 & b_altered_LC_CC_q < 0.1] <- "**"
asterisks_matrix_altered_LC_CC[b_altered_LC_CC_q >= 0.1 & b_altered_LC_CC_q < 0.25] <- "*"


# Create breaks with 0 in the middle
breaks <- seq(-3, 3, length.out = 51)

# Ensure the breaks include 0
mid <- which.min(abs(breaks - 0))

# Create a custom color palette
my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))

# Plot the heatmap with reordered rows
heatmap_altered_LC_CC <- pheatmap(b_altered_LC_CC,
                                  cluster_rows = FALSE,
                                  cluster_cols = FALSE,
                                  display_numbers = asterisks_matrix_altered_LC_CC,
                                  fontsize_number = 10,
                                  color = my_color_palette,
                                  breaks = breaks,
                                  annotation_legend = TRUE,
                                  angle_col = 90)

# ggsave("heatmap_altered_subtypes.pdf", heatmap_altered_LC_CC, width = 6.5, height = 11.5, dpi=1200, limitsize = FALSE)




####---------------  Species abundance: MC vs controls in women (Cross-sectional)----------------#####
# Preparing dataframes
taxo.rel_women <- taxo.rel %>% inner_join(select(alpha.cross.sectional_active_LC_CC, ID, mc_all_label)) %>% 
  filter(sex == "Female")
taxo.rel_women <- taxo.rel_women %>% as.data.frame()
rownames(taxo.rel_women) <- taxo.rel_women$ID

taxo.rel_maaslin_cross_women <- taxo.rel_women %>% select(1:467)
taxo.rel_maaslin_cross_women <- taxo.rel_maaslin_cross_women[, more_or_less]
metadata_maaslin_cross_women <- taxo.rel_women %>% select(468:477)

# maaslin_cross_women <- Maaslin2(input_data = taxo.rel_maaslin_cross_women,
#                                 input_metadata = metadata_maaslin_cross_women,
#                                 output = "maaslin_cross_women",
#                                 min_abundance = 0,
#                                 min_prevalence = 0,
#                                 normalization = "NONE", #already relative abundance
#                                 transform = "LOG",
#                                 analysis_method = "LM",
#                                 max_significance = 0.25,
#                                 fixed_effects = c("mc_all_label","age","bmi","disease_status"),
#                                 reference = c("mc_all_label,MC","disease_status,remission"))


volcano_cross_women <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross_women/all_results.tsv")
volcano_cross_women <- volcano_cross_women %>% filter(metadata=="mc_all_label") %>% select(1,3,4,5,8,9) 

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/volcano_cross_women.xlsx"
# write.xlsx(volcano_cross_CC, file_path, rowNames = FALSE)

### Volcano plot for CC vs HC
volcano_cross_women_hc <- volcano_cross_women %>% filter(value=="Healthy control")
volcano_cross_women_hc$coef <- -(volcano_cross_women_hc$coef)

volcano_cross_women_cd <- volcano_cross_women %>% filter(value=="Chronic diarrhea")
volcano_cross_women_cd$coef <- -(volcano_cross_women_cd$coef)


### Find species that are different in MC compared to both controls
more_women.vs.hc <- volcano_cross_women_hc %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.control_without_diarrhea_women" = "coef",
         "qval_MC.vs.control_without_diarrhea_women" = "qval")

less_women.vs.hc <- volcano_cross_women_hc %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.control_without_diarrhea_women" = "coef",
         "qval_MC.vs.control_without_diarrhea_women" = "qval")

more_women.vs.cd <- volcano_cross_women_cd %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.cd_women" = "coef",
         "qval_MC.vs.cd_women" = "qval")

less_women.vs.cd <- volcano_cross_women_cd %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.cd_women" = "coef",
         "qval_MC.vs.cd_women" = "qval")

## Enriched species in CC vs both controls
more_both_MC_women <- inner_join(more_women.vs.hc, more_women.vs.cd)

## Depleted species in CC vs both controls
less_both_MC_women <- inner_join(less_women.vs.hc, less_women.vs.cd)





####---------------  Species abundance: MC vs controls in men (Cross-sectional)----------------#####
# Preparing dataframes
taxo.rel_men <- taxo.rel %>% inner_join(select(alpha.cross.sectional_active_LC_CC, ID, mc_all_label)) %>% 
  filter(sex == "Male")
taxo.rel_men <- taxo.rel_men %>% as.data.frame()
rownames(taxo.rel_men) <- taxo.rel_men$ID

taxo.rel_maaslin_cross_men <- taxo.rel_men %>% select(1:467)
taxo.rel_maaslin_cross_men <- taxo.rel_maaslin_cross_men[, more_or_less]
metadata_maaslin_cross_men <- taxo.rel_men %>% select(468:477)

# maaslin_cross_men <- Maaslin2(input_data = taxo.rel_maaslin_cross_men,
#                               input_metadata = metadata_maaslin_cross_men,
#                               output = "maaslin_cross_men",
#                               min_abundance = 0,
#                               min_prevalence = 0,
#                               normalization = "NONE", #already relative abundance
#                               transform = "LOG",
#                               analysis_method = "LM",
#                               max_significance = 0.25,
#                               fixed_effects = c("mc_all_label","age","bmi","disease_status"),
#                               reference = c("mc_all_label,MC","disease_status,remission"))


volcano_cross_men <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_cross_men/all_results.tsv")
volcano_cross_men <- volcano_cross_men %>% filter(metadata=="mc_all_label") %>% select(1,3,4,5,8,9) 

# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/volcano_cross_men.xlsx"
# write.xlsx(volcano_cross_CC, file_path, rowNames = FALSE)

### Volcano plot for CC vs HC
volcano_cross_men_hc <- volcano_cross_men %>% filter(value=="Healthy control")
volcano_cross_men_hc$coef <- -(volcano_cross_men_hc$coef)

volcano_cross_men_cd <- volcano_cross_men %>% filter(value=="Chronic diarrhea")
volcano_cross_men_cd$coef <- -(volcano_cross_men_cd$coef)


### Find species that are different in MC compared to both controls
more_men.vs.hc <- volcano_cross_men_hc %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.control_without_diarrhea_men" = "coef",
         "qval_MC.vs.control_without_diarrhea_men" = "qval")

less_men.vs.hc <- volcano_cross_men_hc %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.control_without_diarrhea_men" = "coef",
         "qval_MC.vs.control_without_diarrhea_men" = "qval")

more_men.vs.cd <- volcano_cross_men_cd %>% 
  filter(coef > 0.5 & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.cd_men" = "coef",
         "qval_MC.vs.cd_men" = "qval")

less_men.vs.cd <- volcano_cross_men_cd %>% 
  filter(coef <(-0.5) & qval < 0.25) %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.cd_men" = "coef",
         "qval_MC.vs.cd_men" = "qval")

## Enriched species in CC vs both controls
more_both_MC_men <- inner_join(more_men.vs.hc, more_men.vs.cd)

## Depleted species in CC vs both controls
less_both_MC_men <- inner_join(less_men.vs.hc, less_men.vs.cd)


## Heatmap for all subgroups
altered_heat_LC <- volcano_cross_LC_hc %>% 
  select(1,3,6) %>% 
  rename("coef_LC.vs.Control_without_diarrhea" = "coef",
         "qval_LC.vs.Control_without_diarrhea" = "qval") %>% 
  left_join(select(volcano_cross_LC_cd, feature, coef, qval)) %>% 
  rename("coef_LC.vs.Chronic_diarrhea" = "coef",
         "qval_LC.vs.Chronic_diarrhea" = "qval")

altered_heat_CC <- volcano_cross_CC_hc %>% 
  select(1,3,6) %>% 
  rename("coef_CC.vs.Control_without_diarrhea" = "coef",
         "qval_CC.vs.Control_without_diarrhea" = "qval") %>% 
  left_join(select(volcano_cross_CC_cd, feature, coef, qval)) %>% 
  rename("coef_CC.vs.Chronic_diarrhea" = "coef",
         "qval_CC.vs.Chronic_diarrhea" = "qval")

altered_heat_women <- volcano_cross_women_hc %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Control_without_diarrhea_women" = "coef",
         "qval_MC.vs.Control_without_diarrhea_women" = "qval") %>% 
  left_join(select(volcano_cross_women_cd, feature, coef, qval)) %>% 
  rename("coef_MC.vs.Chronic_diarrhea_women" = "coef",
         "qval_MC.vs.Chronic_diarrhea_women" = "qval")

altered_heat_men <- volcano_cross_men_hc %>% 
  select(1,3,6) %>% 
  rename("coef_MC.vs.Control_without_diarrhea_men" = "coef",
         "qval_MC.vs.Control_without_diarrhea_men" = "qval") %>% 
  left_join(select(volcano_cross_men_cd, feature, coef, qval)) %>% 
  rename("coef_MC.vs.Chronic_diarrhea_men" = "coef",
         "qval_MC.vs.Chronic_diarrhea_men" = "qval")



altered_heat_LC_CC_sex <- altered_heat_LC %>% left_join(altered_heat_CC) %>% left_join(altered_heat_women) %>% left_join(altered_heat_men) %>% as.data.frame()
altered_heat_LC_CC_sex <- altered_heat_LC_CC_sex %>% mutate(new_species = case_when(
  feature == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
  feature == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
  feature == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
  feature == "GGB80011_SGB15265" ~ "p_Firmicutes_GGB80011_SGB15265",
  feature == "GGB9524_SGB14923" ~ "p_Firmicutes_GGB9524_SGB14923",
  TRUE ~ feature
))

altered_heat$new_species <- rownames(altered_heat)
altered_heat_LC_CC_sex <- altered_heat_LC_CC_sex %>% right_join(altered_heat)


rownames(altered_heat_LC_CC_sex) <- altered_heat_LC_CC_sex$new_species
altered_heat_LC_CC_sex <- altered_heat_LC_CC_sex %>% select(-feature, -new_species)

b_altered_LC_CC_sex <- altered_heat_LC_CC_sex %>% select(17,19,1,3,5,7,9,11,13,15)
b_altered_LC_CC_sex <- b_altered_LC_CC_sex[custom_order,] %>% as.matrix()

b_altered_LC_CC_sex_q <- altered_heat_LC_CC_sex %>% select(18,20,2,4,6,8,10,12,14,16) 
b_altered_LC_CC_sex_q <- b_altered_LC_CC_sex_q[custom_order,] %>% as.matrix()

asterisks_matrix_altered_LC_CC_sex <- matrix("", nrow = nrow(b_altered_LC_CC_sex_q), ncol = ncol(b_altered_LC_CC_sex_q))
asterisks_matrix_altered_LC_CC_sex[b_altered_LC_CC_sex_q < 0.05] <- "***"
asterisks_matrix_altered_LC_CC_sex[b_altered_LC_CC_sex_q >= 0.05 & b_altered_LC_CC_sex_q < 0.1] <- "**"
asterisks_matrix_altered_LC_CC_sex[b_altered_LC_CC_sex_q >= 0.1 & b_altered_LC_CC_sex_q < 0.25] <- "*"

# Create breaks with 0 in the middle
breaks <- seq(-3.5, 3.5, length.out = 51)

# Ensure the breaks include 0
mid <- which.min(abs(breaks - 0))

# Create a custom color palette
my_color_palette <- c(colorRampPalette(c("navy", "white"))(mid - 1), "white", colorRampPalette(c("white", "firebrick3"))(51 - mid))

# Define column annotations
col_annotation <- data.frame(Subgroup = c("all","all","subtype","subtype","subtype","subtype","sex","sex","sex","sex"))

# Ensure the column annotation has the right row names (matching the heatmap columns)
rownames(col_annotation) <- colnames(b_altered_LC_CC_sex)

# Define gaps_col based on the transitions in the Subgroup annotation
gaps_col <- which(diff(as.numeric(factor(col_annotation$Subgroup))) != 0)

# Plot the heatmap with gaps between the specified columns
heatmap_altered_LC_CC_sex <- pheatmap(b_altered_LC_CC_sex,
                                      cluster_rows = FALSE,
                                      cluster_cols = FALSE,
                                      display_numbers = asterisks_matrix_altered_LC_CC_sex,
                                      fontsize_number = 10,
                                      color = my_color_palette,
                                      breaks = breaks,
                                      gaps_col = gaps_col,  # Add gaps between "subtype" and "sex" columns
                                      annotation_legend = FALSE,
                                      angle_col = 90)
heatmap_altered_LC_CC_sex
# ggsave("heatmap_altered_LC_CC_sex.pdf", heatmap_altered_LC_CC_sex, width = 9, height = 11.5, dpi=1200, limitsize = FALSE)




####---------------  Paired Wilcoxon test for species of interest (from cross-sectional study) in Longitudinal study ----------------#####
## Use information from cross-sectional to inform longitudinal maaslin
# Select rows with qval ≤ 0.25
list_cross_for_long <- volcano_cross %>% slice(1:70) %>% select(feature)
list_cross_for_long <- list_cross_for_long$feature

taxo_long <- df.long %>% left_join(all_taxo, by = c("ID" = "id"))
rownames(taxo_long) <- taxo_long$ID
disease_status <- taxo_long$disease_status

taxo_cross_for_long <- taxo_long[, list_cross_for_long, drop = FALSE]
taxo_cross_for_long$disease_status <- disease_status



## For multi-omics association
mgx_halla_long <- taxo_long %>% select(14:480)
# 
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/mgx_halla_long.xlsx"
# write.xlsx(mgx_halla_long, file_path, rowNames = TRUE)





# Create an empty list to store the results
p_values_cross_for_long <- list()

# Iterate over each species column in taxo_cross_for_long
for (species_col in 1:70) {
  # Extract the species name
  species_name <- colnames(taxo_cross_for_long)[species_col]
  # Perform the comparison and store the result
  comparison_result <- compare_means(
    formula = as.formula(paste(species_name, "~ disease_status")),
    data = taxo_cross_for_long,
    method = "wilcox.test",
    paired = TRUE,
    p.adjust.method = "none"
  )
  # Extract p-value from the comparison result
  p_value_cross_for_long <- comparison_result$p
  # Store the p-value along with the species name
  p_values_cross_for_long[[species_name]] <- p_value_cross_for_long
}

# Convert the list of p-values to a dataframe
Wilcoxon_cross_for_long <- data.frame(species = names(p_values_cross_for_long), p_value = unlist(p_values_cross_for_long))
Wilcoxon_cross_for_long <- Wilcoxon_cross_for_long[order(Wilcoxon_cross_for_long$p_value), ]
rownames(Wilcoxon_cross_for_long) <- 1:nrow(Wilcoxon_cross_for_long)

Wilcoxon_cross_for_long$q_value <- p.adjust(Wilcoxon_cross_for_long$p_value, method = "BH")

# Box plot for significant species (longitudinal)
for_box_plot <- Wilcoxon_cross_for_long[Wilcoxon_cross_for_long$q_value <= 0.25, ]
for_box_plot <- for_box_plot$species
species_for_box_plot <- taxo_cross_for_long[, for_box_plot, drop = FALSE]
species_for_box_plot$id <- rownames(species_for_box_plot)
species_for_box_plot <- as.data.frame(species_for_box_plot)
species_for_box_plot <- species_for_box_plot %>%
  left_join(select(alpha_with_disease_status, ID, disease_status), by = c("id" = "ID")) %>%
  select(id, disease_status, 1:7)
rownames(species_for_box_plot) <- species_for_box_plot$id

species_for_box_plot_long <- pivot_longer(species_for_box_plot,
                                          cols = 3:9,
                                          names_to = "species",
                                          values_to = "relative_abundance")


species_for_box_plot_long <- species_for_box_plot_long %>%
  mutate(new_species = case_when(
    species == "GGB9534_SGB14937" ~ "Clostridia_GGB9534_SGB14937",
    species == "GGB33586_SGB53517" ~ "Clostridia_GGB33586_SGB53517",
    species == "GGB9699_SGB15216" ~ "Oscillospiraceae_GGB9699_SGB15216",
    species == "GGB3510_SGB4687" ~ "Clostridiaceae_GGB3510_SGB4687",
    TRUE ~ species
  ))

species_order <- c("Clostridia_GGB9534_SGB14937", "Clostridia_GGB33586_SGB53517", "Clostridiales_bacterium_NSJ_32",
                          "Oscillospiraceae_GGB9699_SGB15216", "Clostridiaceae_GGB3510_SGB4687", "Nitrosopumilus_SGB14899", "Collinsella_SGB4121")

# Arrange the data frame by the custom order
species_for_box_plot_long <- species_for_box_plot_long %>%
  arrange(factor(new_species, levels = species_order))

# Box plot
qval.to.plot <- Wilcoxon_cross_for_long %>%
  slice(1:7) %>%
  mutate(qval.format = round(q_value, 3)) %>%
  mutate(group1 = "0") %>%
  mutate(group2 = "1")

qval.to.plot <- qval.to.plot %>%
  mutate(new_species = case_when(
    species == "GGB9534_SGB14937" ~ "Clostridia_GGB9534_SGB14937",
    species == "GGB33586_SGB53517" ~ "Clostridia_GGB33586_SGB53517",
    species == "GGB9699_SGB15216" ~ "Oscillospiraceae_GGB9699_SGB15216",
    species == "GGB3510_SGB4687" ~ "Clostridiaceae_GGB3510_SGB4687",
    TRUE ~ species))


# Calculate half of the minimum non-zero value
half_min_value_long <- min(species_for_box_plot_long$relative_abundance[species_for_box_plot_long$relative_abundance > 0]) / 2

# Replace zeros with the calculated half minimum value
species_for_box_plot_long <- species_for_box_plot_long %>%
  mutate(relative_abundance = ifelse(relative_abundance == 0, half_min_value_long, relative_abundance))

# Take log for relative_abundance
species_for_box_plot_long$log_relative_abundance <- log(species_for_box_plot_long$relative_abundance)

# Box plot for all altered species
box.long <- ggplot(species_for_box_plot_long, aes(factor(new_species, levels = species_order), y = log_relative_abundance)) +
  geom_boxplot(aes(fill = disease_status), position = position_dodge(width = 0.8), width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(fill = disease_status),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1, alpha = 0.4)+
  geom_text(data = qval.to.plot, aes(x = factor(new_species, levels = species_order), y = 2.5, label = paste("q =", qval.format)), vjust = 0, hjust = 0.5, size = 5, color = "black") +
  labs(x = "Species", y = "log(relative abundance(%))", fill = "MC acitivity") +
  ylim(-10,2.6)+
  scale_fill_manual(values = c("active" = "#900C3F", "remission" = "antiquewhite4")) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 8))

box.long
# ggsave(filename = "box.long.pdf",
#        plot = box.long, units = "in", width=15, height=7, dpi = 1200)





####---------------  Microbial pathway: MaAsLin2 regressions (Cross-sectional) ----------------#####
# Preparing dataframes
load("mc_input_2024-07-08.RData") 
all_path <- all_path %>% select(-id)

# Change col_names to path_1~path_482
path_col_names <- colnames(all_path)

num_cols <- ncol(all_path)
new_col_names <- paste("path", 1:num_cols, sep = "_")
colnames(all_path) <- new_col_names


# maaslin_path_cross <- Maaslin2(input_data = all_path,
#                                input_metadata = metadata_maaslin_cross,
#                                output = "maaslin_path_cross",
#                                min_abundance = 0.001,
#                                min_prevalence = 0.1,
#                                normalization = "NONE", #already relative abundance
#                                transform = "LOG",
#                                analysis_method = "LM",
#                                max_significance = 0.25,
#                                fixed_effects = c("mc_all","age","bmi","sex","disease_status"),
#                                reference = c("mc_all,1","sex,Female","disease_status,remission"))

heat_path_cross <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_path_cross/all_results.tsv")
heat_path_cross <- heat_path_cross %>% filter(metadata=="mc_all") %>% select(1,3,4,9)
heat_path_cross$coef <- -(heat_path_cross$coef)

# Change colnames back to the original path_col_names
name_mapping <- setNames(path_col_names, new_col_names)
colnames(all_path) <- name_mapping[colnames(all_path)]

heat_path_cross$feature <- name_mapping[heat_path_cross$feature]

# remove rows with feature that has "unclassified" or "g__" or "s__" 
heat_path_cross_filtered <- heat_path_cross %>%
  filter(!grepl("unclassified", feature)) %>%
  filter(!grepl("g__|s__", feature))

heat_path_cross_filtered$pathname <- heat_path_cross_filtered$feature

### heatmap: Cross-sectional
heat_path_cross_filtered$feature <- heat_path_cross_filtered$pathname
heat_path_cross_filtered$pathname <- NULL
heat_path_cross_filtered <- as.data.frame(heat_path_cross_filtered)

# Create new columns for coefficients for CD and HC
heat_path_cross_filtered$coef_HC <- NA
heat_path_cross_filtered$coef_CD <- NA

# Assign values based on conditions
heat_path_cross_filtered$coef_HC[heat_path_cross_filtered$value == 2] <- heat_path_cross_filtered$coef[heat_path_cross_filtered$value == 2]
heat_path_cross_filtered$coef_CD[heat_path_cross_filtered$value == 0] <- heat_path_cross_filtered$coef[heat_path_cross_filtered$value == 0]


# Create new columns for q-values for CD and HC
heat_path_cross_filtered$qval_HC <- NA
heat_path_cross_filtered$qval_CD <- NA

# Assign values based on conditions
heat_path_cross_filtered$qval_HC[heat_path_cross_filtered$value == 2] <- heat_path_cross_filtered$qval[heat_path_cross_filtered$value == 2]
heat_path_cross_filtered$qval_CD[heat_path_cross_filtered$value == 0] <- heat_path_cross_filtered$qval[heat_path_cross_filtered$value == 0]

heat_path_cross_filtered <-  heat_path_cross_filtered %>% select(-value, -coef, -qval)

# Group by "feature" and summarize
heat_path_cross_filtered_combined <- heat_path_cross_filtered %>%
  group_by(feature) %>%
  summarize(
    coef_HC = first(coef_HC[!is.na(coef_HC)]),
    coef_CD = first(coef_CD[!is.na(coef_CD)]),
    qval_HC = first(qval_HC[!is.na(qval_HC)]),
    qval_CD = first(qval_CD[!is.na(qval_CD)])
  )

heat_path_cross_filtered_combined <- as.data.frame(heat_path_cross_filtered_combined)
rownames(heat_path_cross_filtered_combined) <- heat_path_cross_filtered_combined$feature
heat_path_cross_filtered_combined$feature <- NULL


####---------------  Microbial pathway: Heatmap (Cross-sectional, q<0.25) ----------------#####
heat_path_cross_0.25 <- heat_path_cross_filtered_combined %>% filter(qval_HC<0.25 | qval_CD<0.25)
path_cross <- heat_path_cross_0.25 %>% select(1:2) %>% rename("MC_vs_HC" = "coef_HC", "MC_vs_CD" = "coef_CD")

path_cross_mat <- heat_path_cross_0.25 %>% select(1:2) %>% rename("MC_vs_HC" = "coef_HC", "MC_vs_CD" = "coef_CD") %>% as.matrix()

q_val_cross <- heat_path_cross_0.25 %>%
  select(3:4) %>%
  mutate(
    MC_vs_HC = ifelse(qval_HC >= 0.25, "",
                      ifelse(qval_HC >= 0.10, "*",
                             ifelse(qval_HC >= 0.05, "**", "***"))),
    MC_vs_CD = ifelse(qval_CD >= 0.25, "",
                      ifelse(qval_CD >= 0.10, "*",
                             ifelse(qval_CD >= 0.05, "**", "***")))
  ) %>%
  select(-qval_HC, -qval_CD) %>%
  as.matrix()

heat_path_cross_0.25_fig <- pheatmap(path_cross_mat,
                                     display_numbers = q_val_cross,
                                     fontsize_number = 10,
                                     fontsize_row = 7,
                                     cluster_rows = FALSE, 
                                     cluster_cols = FALSE,
                                     angle_col = 0)

heat_path_cross_0.25_fig



## Pathways enriched in MC as compared to both HC and CD
path_cross_more_both <- path_cross %>% filter(MC_vs_HC>0 & MC_vs_CD>0) %>% as.matrix()

q_val_cross_more_both <- heat_path_cross_0.25 %>%
  filter(coef_HC>0 & coef_CD>0) %>%
  select(3:4) %>%
  mutate(
    MC_vs_HC = ifelse(qval_HC >= 0.25, "",
                      ifelse(qval_HC >= 0.10, "*",
                             ifelse(qval_HC >= 0.05, "**", "***"))),
    MC_vs_CD = ifelse(qval_CD >= 0.25, "",
                      ifelse(qval_CD >= 0.10, "*",
                             ifelse(qval_CD >= 0.05, "**", "***")))
  ) %>%
  select(-qval_HC, -qval_CD) %>%
  as.matrix()

heat_path_cross_more_both_0.25 <- pheatmap(path_cross_more_both,
                                           display_numbers = q_val_cross_more_both,
                                           fontsize_number = 10,
                                           fontsize_row = 7,
                                           cluster_rows = FALSE, 
                                           cluster_cols = FALSE,
                                           angle_col = 0,
                                           color = colorRampPalette(c("#FFEDA0", "#FEB24C", "#F03B20"))(50))

heat_path_cross_more_both_0.25



## Pathways depleted in MC as compared to both HC and CD
path_cross_less_both <- path_cross %>% filter(MC_vs_HC<0 & MC_vs_CD<0) %>% as.matrix()

q_val_cross_less_both <- heat_path_cross_0.25 %>%
  filter(coef_HC<0 & coef_CD<0) %>%
  select(3:4) %>%
  mutate(
    MC_vs_HC = ifelse(qval_HC >= 0.25, "",
                      ifelse(qval_HC >= 0.10, "*",
                             ifelse(qval_HC >= 0.05, "**", "***"))),
    MC_vs_CD = ifelse(qval_CD >= 0.25, "",
                      ifelse(qval_CD >= 0.10, "*",
                             ifelse(qval_CD >= 0.05, "**", "***")))
  ) %>%
  select(-qval_HC, -qval_CD) %>%
  as.matrix()

heat_path_cross_less_both_0.25 <- pheatmap(path_cross_less_both,
                                           display_numbers = q_val_cross_less_both,
                                           fontsize_number = 10,
                                           fontsize_row = 7,
                                           cluster_rows = FALSE, 
                                           cluster_cols = FALSE,
                                           angle_col = 0,
                                           color = colorRampPalette(c("#D0E1F9", "#74C2E1", "#377EB8"))(50))

heat_path_cross_less_both_0.25




####---------------  HAllA for microbial pathways ----------------#####
## Altered microbial pathways compared to both controls
path_halla <- heat_path_cross_filtered_combined %>% 
  filter(qval_HC<0.25 & qval_CD<0.25) %>% 
  rename("MC_vs_HC" = "coef_HC", "MC_vs_CD" = "coef_CD", "qval_MC_vs_HC" = "qval_HC", "qval_MC_vs_CD" = "qval_CD") %>% 
  filter((MC_vs_HC<0 & MC_vs_CD<0) | (MC_vs_HC>0 & MC_vs_CD>0))


path_more_or_less <- rownames(path_halla)
mc_cd_hc <- full_join(mc_cd, mc_hc)

all_path$ID <- rownames(all_path)
path_halla_residual <- readRDS("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_path_cross/fits/residuals.rds") %>%
  t() %>% as.data.frame()

# Change colnames back to the original path_col_names
colnames(path_halla_residual) <- name_mapping[colnames(path_halla_residual)]
path_halla_residual$ID <- rownames(path_halla_residual)

taxo.rel_maaslin_cross$ID <- rownames(taxo.rel_maaslin_cross)

path_halla_residual_1 <- path_halla_residual %>% inner_join(select(taxo.rel_maaslin_cross, ID))
rownames(path_halla_residual_1) <- path_halla_residual_1$ID


path_halla_residual_1$id <- substr(path_halla_residual_1$ID, 1, 5)

path_halla_residual <- path_halla_residual_1 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(path_halla_residual) <- path_halla_residual$ID
path_halla_residual$ID <- NULL
path_halla_residual$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/path_halla_residual.xlsx"
# write.xlsx(path_halla_residual, file_path, rowNames = TRUE)


## All microbial pathways
# path_halla_full <- all_path 
# matching_ids <- rownames(path_halla)
# path_halla_full <- path_halla_full[path_halla_full$ID %in% matching_ids, ]
# path_halla_full$ID <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/path_halla_full.xlsx"
# write.xlsx(path_halla_full, file_path, rowNames = TRUE)



## Microbial pathways correlated with sphingolipids (rho > 0.2 or < -0.2)
path_halla_sphingolipids_0.2 <- read_excel("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/sphingolipids_0.2_paths.xlsx") %>%
  select(2) %>% rename("microbial_pathway" = "X_features")

path_halla_sphingolipids_0.2 <- path_halla_sphingolipids_0.2$microbial_pathway

all_path_2 <- all_path %>% inner_join(select(taxo.rel_maaslin_cross, ID)) %>% select(all_of(path_halla_sphingolipids_0.2), ID)
rownames(all_path_2) <- all_path_2$ID

all_path_2$id <- substr(all_path_2$ID, 1, 5)

path_halla_sphingolipids_0.2 <- all_path_2 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(path_halla_sphingolipids_0.2) <- path_halla_sphingolipids_0.2$ID
path_halla_sphingolipids_0.2$ID <- NULL
path_halla_sphingolipids_0.2$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/path_halla_sphingolipids_0.2.xlsx"
# write.xlsx(path_halla_sphingolipids_0.2, file_path, rowNames = TRUE)




## Microbial pathways correlated with lysophospholipids (rho > 0.2 or < -0.2)
path_halla_lysophospholipids_0.2 <- read_excel("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysophospholipids_0.2_paths.xlsx") %>%
  select(2) %>% rename("microbial_pathway" = "X_features")

path_halla_lysophospholipids_0.2 <- path_halla_lysophospholipids_0.2$microbial_pathway

all_path_3 <- all_path %>% inner_join(select(taxo.rel_maaslin_cross, ID)) %>% select(all_of(path_halla_lysophospholipids_0.2), ID)
rownames(all_path_3) <- all_path_3$ID

all_path_3$id <- substr(all_path_3$ID, 1, 5)

path_halla_lysophospholipids_0.2 <- all_path_3 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(path_halla_lysophospholipids_0.2) <- path_halla_lysophospholipids_0.2$ID
path_halla_lysophospholipids_0.2$ID <- NULL
path_halla_lysophospholipids_0.2$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/path_halla_lysophospholipids_0.2.xlsx"
# write.xlsx(path_halla_lysophospholipids_0.2, file_path, rowNames = TRUE)





## Microbial pathways correlated with lysoplasmalogens (rho > 0.2 or < -0.2)
path_halla_lysoplasmalogens_0.2 <- read_excel("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/lysoplasmalogens_0.2_paths.xlsx") %>%
  select(2) %>% rename("microbial_pathway" = "X_features")

path_halla_lysoplasmalogens_0.2 <- path_halla_lysoplasmalogens_0.2$microbial_pathway

all_path_4 <- all_path %>% inner_join(select(taxo.rel_maaslin_cross, ID)) %>% select(all_of(path_halla_lysoplasmalogens_0.2), ID)
rownames(all_path_4) <- all_path_4$ID

all_path_4$id <- substr(all_path_4$ID, 1, 5)

path_halla_lysoplasmalogens_0.2 <- all_path_4 %>% inner_join(select(mc_cd_hc, id), by=c("id" = "id"))
rownames(path_halla_lysoplasmalogens_0.2) <- path_halla_lysoplasmalogens_0.2$ID
path_halla_lysoplasmalogens_0.2$ID <- NULL
path_halla_lysoplasmalogens_0.2$id <- NULL
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/path_halla_lysoplasmalogens_0.2.xlsx"
# write.xlsx(path_halla_lysoplasmalogens_0.2, file_path, rowNames = TRUE)








####-------------- Microbial EC: MaAsLin2 regressions (Cross-sectional) --------------####
EC <- EC %>% select(-id)

# Change col_names to path_1~path_472
EC_col_names <- colnames(EC)

num_cols_EC <- ncol(EC)
new_col_names_EC <- paste("EC", 1:num_cols_EC, sep = "_")
colnames(EC) <- new_col_names_EC


# maaslin_EC_cross <- Maaslin2(input_data = EC,
#                              input_metadata = metadata_maaslin_cross,
#                              output = "maaslin_EC_cross",
#                              min_abundance = 0.001,
#                              min_prevalence = 0.1,
#                              normalization = "NONE", #already relative abundance
#                              transform = "LOG",
#                              analysis_method = "LM",
#                              max_significance = 0.25,
#                              fixed_effects = c("mc_all","age","bmi","sex","disease_status"),
#                              reference = c("mc_all,1","sex,Female","disease_status,remission"))

heat_EC_cross <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_EC_cross/all_results.tsv")
heat_EC_cross <- heat_EC_cross %>% filter(metadata=="mc_all") %>% select(1,3,4,9)
heat_EC_cross$coef <- -(heat_EC_cross$coef)

# Change colnames back to the original EC_col_names
name_mapping_EC <- setNames(EC_col_names, new_col_names_EC)
colnames(EC) <- name_mapping_EC[colnames(EC)]

heat_EC_cross$feature <- name_mapping_EC[heat_EC_cross$feature]

# remove rows with feature that has "unclassified" or "g__" or "s__" 
heat_EC_cross_filtered <- heat_EC_cross %>%
  filter(!grepl("unclassified", feature)) %>%
  filter(!grepl("g__|s__", feature))



### heatmap: Cross-sectional
heat_EC_cross_filtered <- as.data.frame(heat_EC_cross_filtered)

# Create new columns for coefficients for CD and HC
heat_EC_cross_filtered$coef_HC <- NA
heat_EC_cross_filtered$coef_CD <- NA

# Assign values based on conditions
heat_EC_cross_filtered$coef_HC[heat_EC_cross_filtered$value == 2] <- heat_EC_cross_filtered$coef[heat_EC_cross_filtered$value == 2]
heat_EC_cross_filtered$coef_CD[heat_EC_cross_filtered$value == 0] <- heat_EC_cross_filtered$coef[heat_EC_cross_filtered$value == 0]


# Create new columns for q-values for CD and HC
heat_EC_cross_filtered$qval_HC <- NA
heat_EC_cross_filtered$qval_CD <- NA

# Assign values based on conditions
heat_EC_cross_filtered$qval_HC[heat_EC_cross_filtered$value == 2] <- heat_EC_cross_filtered$qval[heat_EC_cross_filtered$value == 2]
heat_EC_cross_filtered$qval_CD[heat_EC_cross_filtered$value == 0] <- heat_EC_cross_filtered$qval[heat_EC_cross_filtered$value == 0]

heat_EC_cross_filtered <-  heat_EC_cross_filtered %>% select(-value, -coef, -qval)

# Group by "feature" and summarize
heat_EC_cross_filtered_combined <- heat_EC_cross_filtered %>%
  group_by(feature) %>%
  summarize(
    coef_HC = first(coef_HC[!is.na(coef_HC)]),
    coef_CD = first(coef_CD[!is.na(coef_CD)]),
    qval_HC = first(qval_HC[!is.na(qval_HC)]),
    qval_CD = first(qval_CD[!is.na(qval_CD)])
  )

heat_EC_cross_filtered_combined <- as.data.frame(heat_EC_cross_filtered_combined)
rownames(heat_EC_cross_filtered_combined) <- heat_EC_cross_filtered_combined$feature
heat_EC_cross_filtered_combined$feature <- NULL


EC_halla <- heat_EC_cross_filtered_combined %>% 
  filter(qval_HC<0.25 & qval_CD<0.25) %>% 
  rename("MC_vs_HC" = "coef_HC", "MC_vs_CD" = "coef_CD", "qval_MC_vs_HC" = "qval_HC", "qval_MC_vs_CD" = "qval_CD") %>% 
  filter((MC_vs_HC<0 & MC_vs_CD<0) | (MC_vs_HC>0 & MC_vs_CD>0))






####---------------  Microbial pathway: MaAsLin2 regressions (Longitudinal) ----------------#####
# Preparing dataframes
long_id <- rownames(pcoa_taxo_long)
metadata_maaslin_long <- alpha_with_disease_status %>% filter(ID %in% long_id)
rownames(metadata_maaslin_long) <- metadata_maaslin_long$ID
metadata_maaslin_long <- as.data.frame(metadata_maaslin_long)

# maaslin_path_long <- Maaslin2(input_data = all_path,
#                               input_metadata = metadata_maaslin_long,
#                               output = "maaslin_path_long",
#                               min_abundance = 0.001,
#                               min_prevalence = 0.1,
#                               normalization = "NONE", #already relative abundance
#                               transform = "LOG",
#                               analysis_method = "LM",
#                               max_significance = 0.25,
#                               fixed_effects = c("age","bmi","sex","disease_status"),
#                               reference = c("sex,Female","disease_status,remission"))

heat_path_long <- read_tsv("/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/mc_R/maaslin_path_long/all_results.tsv")
heat_path_long <- heat_path_long %>% filter(metadata=="disease_status") %>% select(1,4,9) %>% slice(1:20)

# remove rows with feature that has "unclassified" or "g__" or "s__" 
heat_path_long_filtered <- heat_path_long %>%
  filter(!grepl("unclassified", feature)) %>%
  filter(!grepl("g__|s__", feature))

heat_path_long_filtered$pathname <- gsub("\\."," ",
                                         gsub("\\.\\."," ",heat_path_long_filtered$feature,),)

### heatmap: Longitudinal
heat_path_long_filtered$feature <- heat_path_long_filtered$pathname
heat_path_long_filtered$pathname <- NULL
heat_path_long_filtered <- as.data.frame(heat_path_long_filtered)

rownames(heat_path_long_filtered) <- heat_path_long_filtered$feature
heat_path_long_filtered$feature <- NULL


####---------------  Heatmap (Longitudinal, q<0.25) ----------------#####
path_long <- heat_path_long_filtered %>% select(1) %>% rename("Active_vs_Remission" = "coef") %>% as.matrix()
q_val_long <- heat_path_long_filtered %>%
  select(2) %>%
  mutate(
    Active_vs_Remission = ifelse(qval >= 0.25, "",
                                 ifelse(qval >= 0.10, "*",
                                        ifelse(qval >= 0.05, "**", "***")))
  ) %>%
  select(-qval) %>%
  as.matrix()

heat_path_long_0.25 <- pheatmap(path_long,
                                display_numbers = q_val_long,
                                fontsize_number = 20,
                                fontsize_row = 10,
                                cluster_rows = FALSE, 
                                cluster_cols = FALSE,
                                angle_col = 0)

heat_path_long_0.25

# ggsave(filename = "heat_path_long_0.25.png",
#        plot = heat_path_long_0.25, units = "in", width=8, height=7, dpi = 600)




####---------------  Scatter plot: microbe and metabolite ----------------#####
veillonella_parvula <- mgx_halla %>% select(Veillonella_parvula)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/veillonella_parvula_mgx.xlsx"
# write.xlsx(veillonella_parvula, file_path, rowNames = TRUE)

veillonella_rogosae <- mgx_halla %>% select(Veillonella_rogosae)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/veillonella_rogosae_mgx.xlsx"
# write.xlsx(veillonella_rogosae, file_path, rowNames = TRUE)

intestinibacter_bartlettii <- mgx_halla %>% select(Intestinibacter_bartlettii)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/intestinibacter_bartlettii_mgx.xlsx"
# write.xlsx(intestinibacter_bartlettii, file_path, rowNames = TRUE)

methylobacterium_SGB15164 <- mgx_halla %>% select(Methylobacterium_SGB15164)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/methylobacterium_SGB15164_mgx.xlsx"
# write.xlsx(methylobacterium_SGB15164, file_path, rowNames = TRUE)

collinsella_SGB4121 <- mgx_halla %>% select(Collinsella_SGB4121)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/collinsella_SGB4121_mgx.xlsx"
# write.xlsx(collinsella_SGB4121, file_path, rowNames = TRUE)

mediterraneibacter_butyricigenes <- mgx_halla %>% select(Mediterraneibacter_butyricigenes)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/mediterraneibacter_butyricigenes_mgx.xlsx"
# write.xlsx(mediterraneibacter_butyricigenes, file_path, rowNames = TRUE)

haemophilus_parainfluenzae <- mgx_halla %>% select(Haemophilus_parainfluenzae)
# file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_microbe/haemophilus_parainfluenzae_mgx.xlsx"
# write.xlsx(haemophilus_parainfluenzae, file_path, rowNames = TRUE)





# ####---------------  X Figure ALPHA DIVERSITY (MC vs HC) ----------------#####
# load("match_2024-05-06.RData") 
# ## Chao1 index
# # Cross-sectional: prepare patients list (At baseline, Active MC vs HC)
# df.cross.hc_mc_1 <- alpha.cross.sectional_active %>% filter(mc_all==1)
# df.cross.hc_mc_2 <- alpha_with_disease_status %>% filter(mc_all==2)
# df.cross.hc_mc <- rbind(df.cross.hc_mc_1, df.cross.hc_mc_2)
# df.cross.hc_mc$mc_all_label <- ifelse(startsWith(df.cross.hc_mc$ID, "G"), "Healthy_control", df.cross.hc_mc$mc_all_label)
# 
# # Compare Cross-sectional 
# stat.test_cross_chao1_hc <- compare_means(alpha.chao1 ~ mc_all_label, 
#                                           data=df.cross.hc_mc,
#                                           method = "wilcox.test", paired = F, p.adjust.method="none") 
# 
# bx.hc_mc <- ggboxplot(df.cross.hc_mc, y = "alpha.chao1", x = "mc_all_label",
#                       fill="mc_all_label", palette = c("#900C3F", "#189AF9"),
#                       xlab="Disease type", ylab="Alpha diversity (Chao1 index)", 
#                       ylim=c(0,400),font.x=12, font.y=12, font.tickslab=11, 
#                       legend = "none")
# 
# df.cross.hc_mc %>% filter(mc_all_label=="MC") # MC: n = 131
# df.cross.hc_mc %>% filter(mc_all_label=="Healthy_control") # Chronic diarrhea: n = 393
# 
# fig_chao1_cross_hc_mc <- bx.hc_mc +
#   stat_pvalue_manual(stat.test_cross_chao1_hc, label = "Wilcoxon p < {p.format}", y.position = 400)+
#   annotate("text", x=1, y=320, label = "n = 131") +
#   annotate("text", x=2, y=320, label = "n = 393")
# 
# fig_chao1_cross_hc_mc
# 
# # ggsave(filename = "fig_chao1_cross_hc_mc.png",
# #        plot=fig_chao1_cross_hc_mc, units = "in", height = 5, width = 5, dpi=600)



# ####---------------  X Figure PCoA (MC vs HC)  ----------------------#####
# species.abund.baseline.active <- alpha.cross.sectional_active[, 1] %>% left_join(species.abund, by = c("ID" = "id"))
# species.abund.baseline.active <- as.data.frame(species.abund.baseline.active)
# rownames(species.abund.baseline.active) <- species.abund.baseline.active$ID
# species.abund.baseline.active <- species.abund.baseline.active %>% rename(id = ID)
# species.abund_hc_mc <- full_join(species.abund.baseline.active, taxo_hc)
# rownames(species.abund_hc_mc) <- species.abund_hc_mc$id
# 
# # Replace all NA values in species.abund_hc_mc with 0
# species.abund_hc_mc[is.na(species.abund_hc_mc)] <- 0
# 
# 
# subset.species.abund_hc_mc <- species.abund_hc_mc[, !colnames(species.abund_hc_mc) %in% c("id")]
# 
# bray_hc_mc <- vegdist(subset.species.abund_hc_mc, "bray")  # pairwise dissmilarity: 0=same; 1=maximally dissimilar
# pca_hc_mc <- cmdscale(bray_hc_mc, eig = T)
# pcoap_hc_mc <- data.frame(pca_hc_mc$points)
# pcoap_hc_mc$ID <- rownames(pcoap_hc_mc)
# pcoap_hc_mc <- pcoap_hc_mc[order(pcoap_hc_mc$ID),]
# 
# # percentages of variation explained by PCO1 & 2
# eigs_hc_mc <- pca_hc_mc$eig
# pc1.pct_hc_mc <- eigs_hc_mc[1]/sum(eigs_hc_mc) 
# pc1.pct_hc_mc # pc1.pct = 10.6%
# pc2.pct_hc_mc <- eigs_hc_mc[2]/sum(eigs_hc_mc) 
# pc2.pct_hc_mc # pc2.pct = 6.4%
# 
# alpha.pcoa_hc_mc <- left_join(alpha_with_disease_status, pcoap, by = "ID")
# alpha.pcoa <- alpha.pcoa %>% 
#   mutate(mc_all_label = case_when(
#     mc_all == 1 ~ "MC",
#     mc_all == 0 ~ "Chronic diarrhea",
#     TRUE ~ NA_character_))
# 
# # Modify levels for disease type
# alpha.pcoa$mc_all_label <- fct_recode(alpha.pcoa$mc_all_label,
#                                       "Chronic diarrhea" = "Chronic diarrhea",
#                                       "Microscopic colitis" = "MC")
# 
# # Modify levels for disease status
# alpha.pcoa$disease_status <- fct_recode(alpha.pcoa$disease_status,
#                                         "Active" = "active",
#                                         "Remission" = "remission")
# 
# 
# 
# fig_PCoA_all <- ggscatter(alpha.pcoa, y = "X2", x = "X1",
#                           ylab = paste0('PCo2 (', round(pc2.pct * 100, 1), '%)'),
#                           xlab = paste0('PCo1 (', round(pc1.pct * 100, 1), '%)'),
#                           font.x = 12, font.y = 12, font.tickslab = 11,
#                           color = "mc_all_label", 
#                           palette = c("#900C3F", "#F9D937"), size = 2,
#                           legend = "bottom", legend.title = "Legend Title") +
#   stat_ellipse(aes(color = mc_all_label), linetype=2) +
#   guides(color = guide_legend(title = "Disease"))
# 
# fig_PCoA_all
# 
# ggsave(filename = "fig_PCoA_all.png",
#        plot=fig_PCoA_all, units = "in", height = 5, width = 5, dpi=600)


# ####--------------- X Enzymes: targeted pathways ----------------#####
# enzyme_ceramide_degradation <- all_path_full %>% select(contains("ceramide")) %>% select(1) %>% rename(ceramide_degradation = 1)
# enzyme_ceramide_degradation_subset <- enzyme_ceramide_degradation[rownames(enzyme_ceramide_degradation) %in% rownames(path_halla), , drop = FALSE]
# enzyme_ceramide_degradation <- as.data.frame(enzyme_ceramide_degradation_subset)
# 
# # file_path <- "/Users/sheng-yinchen/Library/CloudStorage/OneDrive-MassGeneralBrigham/CTEU/Projects/MC/metabolomics_MC/MC_metabolite_R/metabolite_enzyme/enzyme_ceramide_degradation_mgx.xlsx"
# # write.xlsx(enzyme_ceramide_degradation, file_path, rowNames = TRUE)
# ####---------------  X Heatmap (Cross-sectional, q<0.05) ----------------#####
# heat_path_cross_0.05 <- heat_path_cross_filtered_combined %>% filter(qval_HC<0.05 | qval_CD<0.05)
# path_cross_0.05 <- heat_path_cross_0.05 %>% select(1:2) %>% rename("MC_vs_HC" = "coef_HC", "MC_vs_CD" = "coef_CD")
# 
# path_cross_0.05_mat <- heat_path_cross_0.05 %>% select(1:2) %>% rename("MC_vs_HC" = "coef_HC", "MC_vs_CD" = "coef_CD") %>% as.matrix()
# q_val_cross_0.05 <- heat_path_cross_0.05 %>%
#   select(3:4) %>%
#   mutate(
#     MC_vs_HC = ifelse(qval_HC >= 0.05, "",
#                 ifelse(qval_HC >= 0.01, "*",
#                        ifelse(qval_HC >= 0.001, "**", "***"))),
#     MC_vs_CD = ifelse(qval_CD >= 0.05, "",
#                 ifelse(qval_CD >= 0.01, "*",
#                        ifelse(qval_CD >= 0.001, "**", "***")))
#   ) %>%
#   select(-qval_HC, -qval_CD) %>%
#   as.matrix()
# 
# heat_path_cross_0.05 <- pheatmap(path_cross_0.05,
#                                  display_numbers = q_val_cross_0.05,
#                                  fontsize_number = 10,
#                                  fontsize_row = 7,
#                                  cluster_rows = FALSE, 
#                                  cluster_cols = FALSE,
#                                  angle_col = 0)
# heat_path_cross_0.05
# 
# # ggsave(filename = "heat_path_cross_0.05.png",
# #        plot = heat_path_cross_0.05, units = "in", width=10, height=10, dpi = 600)




### Box plot for MC vs CD (only select the top 30 species)
box_cd <- volcano_cross_cd[volcano_cross_cd$qval <= 0.25, ] %>% slice(1:30)
box_cd <- box_cd$feature
species_box_cd <- taxo.rel_maaslin_cross[, box_cd, drop = FALSE]
species_box_cd$id <- rownames(species_box_cd)
species_box_cd <- as.data.frame(species_box_cd)
species_box_cd <- species_box_cd %>%
  left_join(select(metadata_maaslin_cross, ID, mc_all), by = c("id" = "ID")) %>%
  select(id, mc_all, 1:30)
rownames(species_box_cd) <- species_box_cd$id
species_box_cd <- species_box_cd %>% filter(!mc_all==2)


species_box_cd_long <- pivot_longer(species_box_cd,
                                    cols = 3:32,
                                    names_to = "species",
                                    values_to = "relative_abundance")


species_box_cd_long <- species_box_cd_long %>%
  mutate(new_species = case_when(
    species == "GGB3175_SGB4191" ~ "f_Clostridiaceae_GGB3175_SGB4191",
    species == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
    species == "GGB3602_SGB4574" ~ "f_Lachnospiraceae_GGB3602_SGB4574",
    species == "GGB9713_SGB15249" ~ "f_Oscillospiraceae_GGB9713_SGB15249",
    TRUE ~ species
  ))

cd_order <- c("Haemophilus_parainfluenzae", "Intestinibacter_bartlettii", "Anaerostipes_hadrus",
              "Clostridium_sp_AF20_17LB", "Clostridium_spiroforme", "f_Clostridiaceae_GGB3175_SGB4191",
              "Veillonella_parvula", "Clostridium_sp_AF36_4", "Lachnospira_eligens",
              "Veillonella_rogosae", "Faecalibacterium_prausnitzii", "Veillonella_dispar",
              "Anaerostipes_caccae", "Bacteroides_nordii", "p_Firmicutes_GGB3612_SGB4882",
              "Alistipes_timonensis", "Roseburia_sp_BX1005", "Faecalimonas_umbilicata",
              "Clostridiales_bacterium", "Pseudoflavonifractor_capillosus", "Dorea_formicigenerans",
              "f_Lachnospiraceae_GGB3602_SGB4574", "Parasutterella_SGB9260", "Bacteroides_xylanisolvens",
              "Clostridium_sp_AF34_10BH", "Collinsella_SGB4121", "f_Oscillospiraceae_GGB9713_SGB15249",
              "Roseburia_inulinivorans", "Blautia_glucerasea", "Senegalimassilia_anaerobia")

# Arrange the data frame by the custom order
species_box_cd_long <- species_box_cd_long %>%
  arrange(factor(new_species, levels = cd_order))

# Box plot
cd_qval.to.plot <- volcano_cross_cd %>%
  slice(1:30) %>%
  mutate(qval.format = round(qval, 3)) %>%
  mutate(group1 = "0") %>%
  mutate(group2 = "1") %>%
  mutate(qval_text = ifelse(qval < 0.001, "q < 0.001", paste0("q =", qval.format)))

cd_qval.to.plot <- cd_qval.to.plot %>%
  mutate(new_species = case_when(
    feature == "GGB3175_SGB4191" ~ "f_Clostridiaceae_GGB3175_SGB4191",
    feature == "GGB3612_SGB4882" ~ "p_Firmicutes_GGB3612_SGB4882",
    feature == "GGB3602_SGB4574" ~ "f_Lachnospiraceae_GGB3602_SGB4574",
    feature == "GGB9713_SGB15249" ~ "f_Oscillospiraceae_GGB9713_SGB15249",
    TRUE ~ feature
  ))

species_box_cd_long$mc_all_label <- "Chronic_diarrhea"
species_box_cd_long$mc_all_label[species_box_cd_long$mc_all == 1] <- "MC"
species_box_cd_long$mc_all_label <- factor(species_box_cd_long$mc_all_label, levels = c("MC", "Chronic_diarrhea"))

bx.species_cd <- ggplot(species_box_cd_long, aes(x = factor(new_species, levels = cd_order), y = relative_abundance)) +
  geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 0.8), width = 0.5) +
  geom_text(data = cd_qval.to.plot, aes(x = factor(new_species, levels = cd_order), y = 1.0, label = qval_text), vjust = 0, hjust = 0.5, size = 3, color = "black") +
  labs(x = "Species", y = "Relative abundance (%)", fill = "MC or Chronic diarrhea") +
  scale_fill_manual(values = c("MC" = "#900C3F", "Chronic_diarrhea" = "#F9D937")) +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1), axis.text.y = element_text(size = 8)) +
  ylim(0, 1.0)

bx.species_cd

# ggsave(filename = "bx.species_cd.png",
#        plot = bx.species_cd, units = "in", width=18, height=10, dpi = 600)



# ### Box plot for MC vs HC (only select the top 30 species)
# box_hc <- volcano_cross_hc[volcano_cross_hc$qval <= 0.25, ] %>% slice(1:30)
# box_hc <- box_hc$feature
# species_box_hc <- taxo.rel_maaslin_cross[, box_hc, drop = FALSE]
# species_box_hc$id <- rownames(species_box_hc)
# species_box_hc <- as.data.frame(species_box_hc)
# species_box_hc <- species_box_hc %>%
#   left_join(select(metadata_maaslin_cross, ID, mc_all), by = c("id" = "ID")) %>%
#   select(id, mc_all, 1:30)
# rownames(species_box_hc) <- species_box_hc$id
# species_box_hc <- species_box_hc %>% filter(!mc_all==0)
# 
# 
# species_box_hc_long <- pivot_longer(species_box_hc,
#                                     cols = 3:32,
#                                     names_to = "species",
#                                     values_to = "relative_abundance")
# 
# 
# species_box_hc_long <- species_box_hc_long %>%
#   mutate(new_species = case_when(
#     species == "GGB3653_SGB4964" ~ "f_Lachnospiraceae_GGB3653_SGB4964",
#     species == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
#     species == "GGB9775_SGB15395" ~ "p_Firmicutes_GGB9775_SGB15395",
#     species == "GGB2653_SGB3574" ~ "p_Firmicutes_GGB2653_SGB3574",
#     species == "GGB3061_SGB4063" ~ "p_Firmicutes_GGB3061_SGB4063",
#     species == "GGB9770_SGB15390" ~ "p_Firmicutes_GGB9770_SGB15390",
#     species == "GGB32900_SGB53446" ~ "f_Oscillospiraceae_GGB32900_SGB53446",
#     species == "GGB2980_SGB3962" ~ "f_Eubacteriales_Family_XIII_Incertae_Sedis_GGB2980_SGB3962",
#     species == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
#     species == "GGB9707_SGB15229" ~ "f_Oscillospiraceae_GGB9707_SGB15229",
#     species == "GGB33928_SGB15225" ~ "f_Oscillospiraceae_GGB33928_SGB15225",
#     TRUE ~ species
#   ))
# 
# hc_order <- c(
#   "Clostridiales_Family_XIII_bacterium_BX16",
#   "Dysosmobacter_sp_NSJ_60",
#   "Methylobacterium_SGB15164",
#   "f_Lachnospiraceae_GGB3653_SGB4964",
#   "Lawsonibacter_asaccharolyticus",
#   "f_Clostridia_unclassified_GGB33586_SGB53517",
#   "Oscillospiraceae_bacterium_Marseille_Q3528",
#   "p_Firmicutes_GGB9775_SGB15395",
#   "Anaerotruncus_rubiinfantis",
#   "Rothia_mucilaginosa",
#   "p_Firmicutes_GGB2653_SGB3574",
#   "p_Firmicutes_GGB3061_SGB4063",
#   "Clostridiales_bacterium",
#   "p_Firmicutes_GGB9770_SGB15390",
#   "Collinsella_SGB4121",
#   "Alistipes_putredinis",
#   "Mediterraneibacter_butyricigenes",
#   "Eggerthella_lenta",
#   "f_Oscillospiraceae_GGB32900_SGB53446",
#   "Haemophilus_parainfluenzae",
#   "Veillonella_dispar",
#   "f_Eubacteriales_Family_XIII_Incertae_Sedis_GGB2980_SGB3962",
#   "Intestinibacter_bartlettii",
#   "Clostridiaceae_bacterium_Marseille_Q4143",
#   "Dysosmobacter_welbionis",
#   "Clostridium_symbiosum",
#   "f_Clostridia_unclassified_GGB9534_SGB14937",
#   "Veillonella_parvula",
#   "f_Oscillospiraceae_GGB9707_SGB15229",
#   "f_Oscillospiraceae_GGB33928_SGB15225"
# )
# 
# # Arrange the data frame by the custom order
# species_box_hc_long <- species_box_hc_long %>%
#   arrange(factor(new_species, levels = hc_order))
# 
# # Box plot
# hc_qval.to.plot <- volcano_cross_hc %>%
#   slice(1:30) %>%
#   mutate(qval.format = round(qval, 3)) %>%
#   mutate(group1 = "0") %>%
#   mutate(group2 = "1") %>%
#   mutate(qval_text = ifelse(qval < 0.001, "q < 0.001", paste0("q =", qval.format)))
# 
# 
# hc_qval.to.plot <- hc_qval.to.plot %>%
#   mutate(new_species = case_when(
#     feature == "GGB3653_SGB4964" ~ "f_Lachnospiraceae_GGB3653_SGB4964",
#     feature == "GGB33586_SGB53517" ~ "f_Clostridia_unclassified_GGB33586_SGB53517",
#     feature == "GGB9775_SGB15395" ~ "p_Firmicutes_GGB9775_SGB15395",
#     feature == "GGB2653_SGB3574" ~ "p_Firmicutes_GGB2653_SGB3574",
#     feature == "GGB3061_SGB4063" ~ "p_Firmicutes_GGB3061_SGB4063",
#     feature == "GGB9770_SGB15390" ~ "p_Firmicutes_GGB9770_SGB15390",
#     feature == "GGB32900_SGB53446" ~ "f_Oscillospiraceae_GGB32900_SGB53446",
#     feature == "GGB2980_SGB3962" ~ "f_Eubacteriales_Family_XIII_Incertae_Sedis_GGB2980_SGB3962",
#     feature == "GGB9534_SGB14937" ~ "f_Clostridia_unclassified_GGB9534_SGB14937",
#     feature == "GGB9707_SGB15229" ~ "f_Oscillospiraceae_GGB9707_SGB15229",
#     feature == "GGB33928_SGB15225" ~ "f_Oscillospiraceae_GGB33928_SGB15225",
#     TRUE ~ feature
#   ))
# 
# species_box_hc_long$mc_all_label <- "Healthy_control"
# species_box_hc_long$mc_all_label[species_box_hc_long$mc_all == 1] <- "MC"
# species_box_hc_long$mc_all_label <- factor(species_box_hc_long$mc_all_label, levels = c("MC", "Healthy_control"))
# 
# bx.species_hc <- ggplot(species_box_hc_long, aes(x = factor(new_species, levels = hc_order), y = relative_abundance)) +
#   geom_boxplot(aes(fill = mc_all_label), position = position_dodge(width = 0.8), width = 0.5) +
#   geom_text(data = hc_qval.to.plot, aes(x = factor(new_species, levels = hc_order), y = 1.0, label = qval_text), vjust = 0, hjust = 0.5, size = 3, color = "black") +
#   labs(x = "Species", y = "Relative abundance (%)", fill = "MC or Healthy control") +
#   scale_fill_manual(values = c("MC" = "#900C3F", "Healthy_control" = "#189AF9")) +
#   theme_minimal() +
#   theme(legend.position = "top", axis.text.x = element_text(angle = 70, hjust = 1), axis.text.y = element_text(size = 8)) +
#   ylim(0, 1.0)
# 
# bx.species_hc
# 
# # ggsave(filename = "bx.species_hc.png",
# #        plot = bx.species_hc, units = "in", width=18, height=10, dpi = 600)
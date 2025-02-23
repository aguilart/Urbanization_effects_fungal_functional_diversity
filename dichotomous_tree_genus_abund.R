

######################################################################################################
##########  DATA ANALYSIS FUNCTIONAL TURNOVER (GENUS ABUNDANCE)   ####################
######################################################################################################



# This will be close to 1. It will not be one because fungi without genera assinged were removed in
# update analysis
sum(distinct(genus_site_funct_long[genus_site_funct_long$sites == "TAM-U3-200819-A",
                      c("genus", "sites", "abundance")])$abundance)

sum(genus_site_funct_long$abundance[protu_df1$sites == "TAM-U3-200819-A"])#This is also 1


# Dichotomous tree -> Genera abundance

free_vs_symb_abund <-
  genus_site_funct_long %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  filter(!is.na(functional_group))

free_vs_symb_abund$functional_group[-grep("sapro",
                                         free_vs_symb_abund$functional_group,
                                         ignore.case = T)] <- "host_associated"
free_vs_symb_abund$functional_group[grep("sapro",
                                        free_vs_symb_abund$functional_group,
                                        ignore.case = T)] <- "free_living"
free_vs_symb_abund <- 
  distinct(free_vs_symb_abund) %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(abundance)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq)

free_vs_symb_abund <- left_join(free_vs_symb_abund,
                               comm_wt_means_spore[, c(1:6)])

# 2. Wood sapro vs Sapro

sapro_vs_wsapro_abund <-
  genus_site_funct_long %>% 
  filter(!is.na(functional_group)) %>% 
  filter(grepl("sapro", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = recode(functional_group, Dung_Saprotroph = "Saprotroph")) %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(abundance)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

sapro_vs_wsapro_abund$Wood_Saprotroph[which(is.na(sapro_vs_wsapro_abund$Wood_Saprotroph))] <- 0

# 3. Comensalists vs symbionts
com_vs_symb_abund <-
  genus_site_funct_long %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  mutate(functional_group=case_when(
    !grepl("Epiphyte", functional_group)~"Symbionts",
    grepl("Epiphyte", functional_group)~"Commensalists")) %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(abundance)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

com_vs_symb_abund$Commensalists[which(is.na(com_vs_symb_abund$Commensalists))] <- 0
com_vs_symb_abund$Commensalists[which(is.na(com_vs_symb_abund$Symbionts))] <- 0


# 4. Mutualists vs pathogens

mut_vs_path_abund <-
  genus_site_funct_long %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  filter(!grepl("Epiphyte", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = case_when(
    !grepl("pathogen|parasite|trapping", functional_group, ignore.case = T) ~ "Mutualists",
    grepl("pathogen|parasite|trapping", functional_group, ignore.case = T)~"Pathogens")) %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(abundance)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

# Y <- cbind(mut_vs_path_abund$Mutualists,
#            mut_vs_path_abund$Pathogens)

# 5- Pathogens big (plants & animals) vs small (fungi & protists)
pa_vs_fp_abund <-
  genus_site_funct_long %>% 
  #mutate(abundance = as.numeric(abundance)) %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  filter(!grepl("Epiphyte", functional_group, ignore.case = T)) %>% 
  filter(grepl("pathogen|parasite|trapping", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = case_when(
    !grepl("lichen|fungal|algae", functional_group, ignore.case = T) ~ "Path_Plant_Animals",
    grepl("lichen|fungal|algae", functional_group, ignore.case = T)~"Path_Fungi_Protists")) %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(abundance)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

pa_vs_fp_abund$Path_Fungi_Protists[which(is.na(pa_vs_fp_abund$Path_Fungi_Protists))] <- 0


# 6. Plant pathogens vs all animal pathogens
plant_vs_animal_abund <-
  genus_site_funct_long %>% 
  #mutate(abundance = as.numeric(abundance)) %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  filter(!grepl("Epiphyte", functional_group, ignore.case = T)) %>% 
  filter(grepl("pathogen|parasite|trapping", functional_group, ignore.case = T)) %>% 
  filter(!grepl("lichen|fungal|algae", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = case_when(
    !grepl("Plant", functional_group, ignore.case = T) ~ "Animal_pathogen",
    grepl("Plant", functional_group, ignore.case = T)~"Plant_pathogen")) %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(abundance)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

plant_vs_animal_abund$Animal_pathogen[which(is.na(plant_vs_animal_abund$Animal_pathogen))] <- 0



# Calculating effect sizes

Y <- cbind(free_vs_symb_abund$free_living,
           free_vs_symb_abund$host_associated)

free_vs_symb_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                          family = binomial,
                          data = free_vs_symb_abund)

#

Y <- cbind(sapro_vs_wsapro_abund$Wood_Saprotroph,
           sapro_vs_wsapro_abund$Saprotroph)

sapro_vs_wsapro_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                             family = binomial,
                             data = sapro_vs_wsapro_abund)

#

Y <- cbind(com_vs_symb_abund$Symbionts,
           com_vs_symb_abund$Commensalists)

com_vs_symb_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                         family = binomial,
                         data = com_vs_symb_abund)

#
Y <- cbind(mut_vs_path_abund$Mutualists,
           mut_vs_path_abund$Pathogens)

mut_vs_path_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                         family = binomial,
                         data = mut_vs_path_abund)

#
Y <- cbind(pa_vs_fp_abund$Path_Plant_Animals,
           pa_vs_fp_abund$Path_Fungi_Protists)

pa_vs_fp_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                      family = binomial,
                      data = pa_vs_fp_abund)

#
Y <- cbind(plant_vs_animal_abund$Animal_pathogen,
           plant_vs_animal_abund$Plant_pathogen)

plant_vs_animal_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                             family = binomial,
                             data = plant_vs_animal_abund)
####

library(ggeffects)

CI_nu <- ggpredict(free_vs_symb_m_nu, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_free_living"
CI_nu_1 <- CI_nu

Y <- cbind(free_vs_symb_abund$host_associated, free_vs_symb_abund$free_living)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = free_vs_symb_abund)

CI_nu <- ggpredict(m, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_host_associated"
CI_nu_1 <- rbind(CI_nu_1, CI_nu)

#
CI_nu <- ggpredict(sapro_vs_wsapro_m_nu, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Wood_Saprotroph"
CI_nu_2 <- CI_nu

Y <- cbind(sapro_vs_wsapro_abund$Saprotroph, sapro_vs_wsapro_abund$Wood_Saprotroph)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = sapro_vs_wsapro_abund)

CI_nu <- ggpredict(m, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Saprotroph"
CI_nu_2 <- rbind(CI_nu_2, CI_nu)

#
CI_nu <- ggpredict(com_vs_symb_m_nu, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Symbionts"
CI_nu_3 <- CI_nu


# I get an error here
Y <- cbind(com_vs_symb_abund$Commensalists, com_vs_symb_abund$Symbionts)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = com_vs_symb_abund)

CI_nu <- ggpredict(m, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Commensalists"
CI_nu_3 <- rbind(CI_nu_3, CI_nu)

#
CI_nu <- ggpredict(mut_vs_path_m_nu, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Mutualists"
CI_nu_4 <- CI_nu

Y <- cbind(mut_vs_path_abund$Pathogens, mut_vs_path_abund$Mutualists)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = mut_vs_path_abund)

CI_nu <- ggpredict(m, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Pathogens"
CI_nu_4 <- rbind(CI_nu_4, CI_nu)

#
CI_nu <- ggpredict(pa_vs_fp_m_nu, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Path_Plant_Animals"
CI_nu_5 <- CI_nu

Y <- cbind(pa_vs_fp_abund$Path_Fungi_Protists, pa_vs_fp_abund$Path_Plant_Animals)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = pa_vs_fp_abund)

CI_nu <- ggpredict(m, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Path_Fungi_Protists"
CI_nu_5 <- rbind(CI_nu_5, CI_nu)

#
CI_nu <- ggpredict(plant_vs_animal_m_nu, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Animal_pathogen"
CI_nu_6 <- CI_nu

Y <- cbind(plant_vs_animal_abund$Plant_pathogen, plant_vs_animal_abund$Animal_pathogen)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = plant_vs_animal_abund)

CI_nu <- ggpredict(m, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_Plant_pathogen"
CI_nu_6 <- rbind(CI_nu_6, CI_nu)

rm(CI_nu)

effect_size_nu <-
  bind_rows(
    CI_nu_1 %>% 
      select(Type, predicted, group, comparison) %>% 
      pivot_wider(names_from = c(Type, group),
                  values_from = predicted),
    CI_nu_2 %>% 
      select(Type, predicted, group, comparison) %>% 
      pivot_wider(names_from = c(Type, group),
                  values_from = predicted),
    CI_nu_3 %>% 
      select(Type, predicted, group, comparison) %>% 
      pivot_wider(names_from = c(Type, group),
                  values_from = predicted),
    CI_nu_4 %>% 
      select(Type, predicted, group, comparison) %>% 
      pivot_wider(names_from = c(Type, group),
                  values_from = predicted),
    CI_nu_5 %>% 
      select(Type, predicted, group, comparison) %>% 
      pivot_wider(names_from = c(Type, group),
                  values_from = predicted),
    CI_nu_6 %>% 
      select(Type, predicted, group, comparison) %>% 
      pivot_wider(names_from = c(Type, group),
                  values_from = predicted)
  )

effect_size_nu$effect_air <- effect_size_nu$Air_Natural - effect_size_nu$Air_Urban 
effect_size_nu$effect_soil <- effect_size_nu$Soil_Natural - effect_size_nu$Soil_Urban


#effect_size_nu$x <- c("B", "H", "C", "A", "F", "I", "E", "G", "F", "H", "G", "E")
effect_size_nu$xend <- c(3, 27, 6, 0, 21, 30, 18, 24, 21, 27, 24, 18)
effect_size_nu$yend <- c(18, 18, 15, 15, 15, 15, 12, 12, 6, 6, 0, 0)

effect_size_nu$x <- c(15,15,3,3,27,27,21,21,24,24,21,21)
effect_size_nu$y <- c(21,21,18,18,18,18,15,15,12,12,6,6)

effect_size_nu$comparison <- gsub("prob_", "", effect_size_nu$comparison)

effect_size_nu$group <- rep(LETTERS[1:6], each = 2)

effect_size_nu$dir_air <- "#9ECAE1"# Blue, positive to cities
effect_size_nu$dir_air[which(effect_size_nu$effect_air>0)] <- "#FC9272"# Red, positive to natural

effect_size_nu$dir_soil <- "#3182BD"
effect_size_nu$dir_soil[which(effect_size_nu$effect_soil>0)] <- "#DE2D26"

m <- 0.4

effect_size_nu$x_air <- effect_size_nu$x - m
effect_size_nu$x_air[which(effect_size_nu$effect_air<0)] <- 
  effect_size_nu$x[which(effect_size_nu$effect_air<0)] + m


effect_size_nu$xend_air <- effect_size_nu$xend - m
effect_size_nu$xend_air[which(effect_size_nu$effect_air<0)] <- 
  effect_size_nu$xend[which(effect_size_nu$effect_air<0)] + m

#
effect_size_nu$x_soil <- effect_size_nu$x + m
effect_size_nu$x_soil[which(effect_size_nu$effect_soil<0)] <- 
  effect_size_nu$x[which(effect_size_nu$effect_soil<0)] - m

effect_size_nu$xend_soil <- effect_size_nu$xend + m
effect_size_nu$xend_soil[which(effect_size_nu$effect_soil<0)] <- 
  effect_size_nu$xend[which(effect_size_nu$effect_soil<0)] - m

#
effect_size_nu$x_air[c(1,2)] <- c(14,16)
effect_size_nu$xend_air[c(1,2)] <- c(2,28)
effect_size_nu$x_soil[c(1,2)] <- c(16,14)
effect_size_nu$xend_soil[c(1,2)] <- c(4,26)

effect_size_nu$comparison[6] <- "Commensalists"
effect_size_nu$comparison[9] <- "Pathogens of\nPlants & Animals"
effect_size_nu$comparison[10] <- "Pathogens of\nFungi & Protists"


# Finally plotting

pdf("Figures/decision_tree_genus_abundance.pdf",
    width = 8.515,
    height = 11.015)

effect_size_nu %>% 
  bind_rows(data.frame(comparison = "Functional\ndiversity",
                       xend = 15,#"E"
                       yend = 21,
                       size = 0.1)) %>% 
  ggplot() +
  
  geom_segment(aes(x=x_air, xend=xend_air, y=y, yend=yend, size = abs(effect_air),
                   color = dir_air)) +
  
  geom_segment(aes(x=x_soil, xend=xend_soil, y=y, yend=yend, size = abs(effect_soil),
                   color = dir_soil))+
  
  geom_point(aes(xend,yend), size=20, shape=16, col= "gray") +
  
  geom_text(aes(xend, yend, label = comparison, fontface = "bold")) +
  
  scale_color_identity() +
  guides(size = guide_legend(title="Effect size")) +
  
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank())

dev.off()
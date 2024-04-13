# Dichotomous tree -> PROTU presence 

# I wanted to stop using this one and instead getting a tree based on the abundance of genera
# that turned out to be tricky from the glmer side and as MArch 2024 I got an error that stops me
# from at least getting the dichotomous tree.

# This version is based on the number of PROTU. It is an "as good as it gets" proxy for the abundance
# of genera, as genera that have been assigned to many PROTU, will increase the counts to a guild

free_vs_symb_data <-
  site_funct_group_count %>% 
  filter(!is.na(functional_group))

free_vs_symb_data$functional_group[-grep("sapro",
                                         free_vs_symb_data$functional_group,
                                         ignore.case = T)] <- "host_associated"
free_vs_symb_data$functional_group[grep("sapro",
                                        free_vs_symb_data$functional_group,
                                        ignore.case = T)] <- "free_living"
free_vs_symb_data <- 
  free_vs_symb_data %>% 
  mutate(value = as.numeric(value)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(value)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq)

free_vs_symb_data <- left_join(free_vs_symb_data,
                               comm_wt_means_spore[, c(1:6)])

# 593
Y <- cbind(free_vs_symb_data$free_living,
           free_vs_symb_data$host_associated)

free_vs_symb_m_all <-glmer(Y ~ Type * Ecosystem + (1|Location/plot),
                           family = binomial,
                           data = free_vs_symb_data)

free_vs_symb_m_t_rt <-glmer(Y ~ Type + Ecosystem + (1|Location/plot),
                            family = binomial,
                            data = free_vs_symb_data)

free_vs_symb_m_rt <-glmer(Y ~ Ecosystem + (1|Location/plot),
                          family = binomial,
                          data = free_vs_symb_data)

free_vs_symb_m_t <-glmer(Y ~ Type + (1|Location/plot),
                         family = binomial,
                         data = free_vs_symb_data)

free_vs_symb_0 <-glmer(Y ~ 1 + (1|Location/plot),
                       family = binomial,
                       data = free_vs_symb_data)

comp_mod_free_sym <- rbind(
  anova(free_vs_symb_0,
        free_vs_symb_m_rt, method = "ML"),
  anova(free_vs_symb_0,
        free_vs_symb_m_t, method = "ML"),
  anova(free_vs_symb_m_all,
        free_vs_symb_m_t_rt,
        free_vs_symb_0, method = "ML"))

comp_mod_free_sym <- rapply(comp_mod_free_sym, classes = "numeric", how = "replace",
                            function(x){round(x, 2)})
comp_mod_free_sym$model <- row.names(comp_mod_free_sym)
comp_mod_free_sym <- comp_mod_free_sym[c(1,4,2,6,7) , c(9, 1:8)]
write.table(comp_mod_free_sym, "model_comparisons/1_comp_mod_free_sym.txt",
            row.names = F, sep = ";")

# 2. Wood sapro vs Sapro

sapro_vs_wsapro_data <-
  site_funct_group_count %>% 
  #mutate(value = as.numeric(value)) %>% 
  filter(!is.na(functional_group)) %>% 
  filter(grepl("sapro", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = recode(functional_group, Dung_Saprotroph = "Saprotroph")) %>% 
  mutate(value = as.numeric(value)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(value)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

Y <- cbind(sapro_vs_wsapro_data$Wood_Saprotroph,
           sapro_vs_wsapro_data$Saprotroph)

sapro_vs_wsapro_m_all <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                              family = binomial,
                              data = sapro_vs_wsapro_data)

sapro_vs_wsapro_m_t_rt <-glmer(Y ~ Type + Ecosystem + (1|Location/plot),
                               family = binomial,
                               data = sapro_vs_wsapro_data)

sapro_vs_wsapro_m_rt <-glmer(Y ~ Ecosystem + (1|Location/plot),
                             family = binomial,
                             data = sapro_vs_wsapro_data)

sapro_vs_wsapro_m_t <-glmer(Y ~ Type + (1|Location/plot),
                            family = binomial,
                            data = sapro_vs_wsapro_data)

sapro_vs_wsapro_0 <-glmer(Y ~ 1 + (1|Location/plot),
                          family = binomial,
                          data = sapro_vs_wsapro_data)

comp_mod_sapro_wsapro <- rbind(
  anova(sapro_vs_wsapro_0,
        sapro_vs_wsapro_m_rt, method = "ML"),
  anova(sapro_vs_wsapro_0,
        sapro_vs_wsapro_m_t, method = "ML"),
  anova(sapro_vs_wsapro_m_all,
        sapro_vs_wsapro_m_t_rt,
        sapro_vs_wsapro_0, method = "ML"))

comp_mod_sapro_wsapro <- rapply(comp_mod_sapro_wsapro, classes = "numeric", how = "replace",
                                function(x){round(x, 2)})
comp_mod_sapro_wsapro$model <- row.names(comp_mod_sapro_wsapro)
comp_mod_sapro_wsapro <- comp_mod_sapro_wsapro[c(1,4,2,6,7) , c(9, 1:8)]
write.table(comp_mod_sapro_wsapro, "model_comparisons/2_comp_mod_sapro_wsapro.txt",
            row.names = F, sep = ";")

# 3. Comensalists vs symbionts
com_vs_symb_data <-
  site_funct_group_count %>% 
 # mutate(value = as.numeric(value)) %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  mutate(functional_group=case_when(
    !grepl("Epiphyte", functional_group)~"Symbionts",
    grepl("Epiphyte", functional_group)~"Commensalists")) %>% 
  mutate(value = as.numeric(value)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(value)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

com_vs_symb_data$Commensalists[which(is.na(com_vs_symb_data$Commensalists))] <- 0
com_vs_symb_data$Commensalists[which(is.na(com_vs_symb_data$Symbionts))] <- 0

Y <- cbind(com_vs_symb_data$Symbionts,
           com_vs_symb_data$Commensalists)

com_vs_symb_m_all <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                          family = binomial,
                          data = com_vs_symb_data)

com_vs_symb_m_t_rt <-glmer(Y ~ Type + Ecosystem + (1|Location/plot),
                           family = binomial,
                           data = com_vs_symb_data)

com_vs_symb_m_rt <-glmer(Y ~ Ecosystem + (1|Location/plot),
                         family = binomial,
                         data = com_vs_symb_data)

com_vs_symb_m_t <-glmer(Y ~ Type + (1|Location/plot),
                        family = binomial,
                        data = com_vs_symb_data)

com_vs_symb_0 <-glmer(Y ~ 1 + (1|Location/plot),
                      family = binomial,
                      data = com_vs_symb_data)


comp_mod_com_symb <- rbind(
  anova(com_vs_symb_0,
        com_vs_symb_m_rt, method = "ML"),
  anova(com_vs_symb_0,
        com_vs_symb_m_t, method = "ML"),
  anova(com_vs_symb_m_all,
        com_vs_symb_m_t_rt,
        com_vs_symb_0, method = "ML"))

comp_mod_com_symb <- rapply(comp_mod_com_symb, classes = "numeric", how = "replace",
                            function(x){round(x, 2)})
comp_mod_com_symb$model <- row.names(comp_mod_com_symb)
comp_mod_com_symb <- comp_mod_com_symb[c(1,4,2,6,7) , c(9, 1:8)]
write.table(comp_mod_com_symb, "model_comparisons/3_comp_mod_com_symb.txt",
            row.names = F, sep = ";")

# 4. Mutualists vs pathogens

mut_vs_path_data <-
  site_funct_group_count %>% 
  #mutate(value = as.numeric(value)) %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  filter(!grepl("Epiphyte", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = case_when(
    !grepl("pathogen|parasite|trapping", functional_group, ignore.case = T) ~ "Mutualists",
    grepl("pathogen|parasite|trapping", functional_group, ignore.case = T)~"Pathogens")) %>% 
  mutate(value = as.numeric(value)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(value)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

Y <- cbind(mut_vs_path_data$Mutualists,
           mut_vs_path_data$Pathogens)

mut_vs_path_m_all <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                          family = binomial,
                          data = mut_vs_path_data)

mut_vs_path_m_t_rt <-glmer(Y ~ Type + Ecosystem + (1|Location/plot),
                           family = binomial,
                           data = mut_vs_path_data)

mut_vs_path_m_rt <-glmer(Y ~ Ecosystem + (1|Location/plot),
                         family = binomial,
                         data = mut_vs_path_data)

mut_vs_path_m_t <-glmer(Y ~ Type + (1|Location/plot),
                        family = binomial,
                        data = mut_vs_path_data)

mut_vs_path_0 <-glmer(Y ~ 1 + (1|Location/plot),
                      family = binomial,
                      data = mut_vs_path_data)

comp_mod_mut_path <- rbind(
  anova(mut_vs_path_0,
        mut_vs_path_m_rt, method = "ML"),
  anova(mut_vs_path_0,
        mut_vs_path_m_t, method = "ML"),
  anova(mut_vs_path_m_all,
        mut_vs_path_m_t_rt,
        mut_vs_path_0, method = "ML"))

comp_mod_mut_path <- rapply(comp_mod_mut_path, classes = "numeric", how = "replace",
                            function(x){round(x, 2)})
comp_mod_mut_path$model <- row.names(comp_mod_mut_path)
comp_mod_mut_path <- comp_mod_mut_path[c(1,4,2,6,7) , c(9, 1:8)]
write.table(comp_mod_mut_path, "model_comparisons/4_comp_mod_mut_path.txt",
            row.names = F, sep = ";")

# 5- Pathogens big (plants & animals) vs small (fungi & protists)
pa_vs_fp_data <-
  site_funct_group_count %>% 
  #mutate(value = as.numeric(value)) %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  filter(!grepl("Epiphyte", functional_group, ignore.case = T)) %>% 
  filter(grepl("pathogen|parasite|trapping", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = case_when(
    !grepl("lichen|fungal|algae", functional_group, ignore.case = T) ~ "Path_Plant_Animals",
    grepl("lichen|fungal|algae", functional_group, ignore.case = T)~"Path_Fungi_Protists")) %>% 
  mutate(value = as.numeric(value)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(value)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

pa_vs_fp_data$Path_Fungi_Protists[which(is.na(pa_vs_fp_data$Path_Fungi_Protists))] <- 0

Y <- cbind(pa_vs_fp_data$Path_Plant_Animals,
           pa_vs_fp_data$Path_Fungi_Protists)

pa_vs_fp_m_all <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                       family = binomial,
                       data = pa_vs_fp_data)

pa_vs_fp_m_t_rt <-glmer(Y ~ Type + Ecosystem + (1|Location/plot),
                        family = binomial,
                        data = pa_vs_fp_data)

pa_vs_fp_m_rt <-glmer(Y ~ Ecosystem + (1|Location/plot),
                      family = binomial,
                      data = pa_vs_fp_data)

pa_vs_fp_m_t <-glmer(Y ~ Type + (1|Location/plot),
                     family = binomial,
                     data = pa_vs_fp_data)

pa_vs_fp_0 <-glmer(Y ~ 1 + (1|Location/plot),
                   family = binomial,
                   data = pa_vs_fp_data)

comp_mod_pa_fp <- rbind(
  anova(pa_vs_fp_0,
        pa_vs_fp_m_rt, method = "ML"),
  anova(pa_vs_fp_0,
        pa_vs_fp_m_t, method = "ML"),
  anova(pa_vs_fp_m_all,
        pa_vs_fp_m_t_rt,
        pa_vs_fp_0, method = "ML"))

comp_mod_pa_fp <- rapply(comp_mod_pa_fp, classes = "numeric", how = "replace",
                         function(x){round(x, 2)})
comp_mod_pa_fp$model <- row.names(comp_mod_pa_fp)
comp_mod_pa_fp <- comp_mod_pa_fp[c(1,4,2,6,7) , c(9, 1:8)]
write.table(comp_mod_pa_fp, "model_comparisons/5_comp_mod_pa_fp.txt",
            row.names = F, sep = ";")

# 6. Plant pathogens vs all animal pathogens
plant_vs_animal_data <-
  site_funct_group_count %>% 
  #mutate(value = as.numeric(value)) %>% 
  filter(!is.na(functional_group)) %>% 
  filter(!grepl("sapro", functional_group, ignore.case = T)) %>% 
  filter(!grepl("Epiphyte", functional_group, ignore.case = T)) %>% 
  filter(grepl("pathogen|parasite|trapping", functional_group, ignore.case = T)) %>% 
  filter(!grepl("lichen|fungal|algae", functional_group, ignore.case = T)) %>% 
  mutate(functional_group = case_when(
    !grepl("Plant", functional_group, ignore.case = T) ~ "Animal_pathogen",
    grepl("Plant", functional_group, ignore.case = T)~"Plant_pathogen")) %>% 
  mutate(value = as.numeric(value)) %>% 
  distinct() %>% 
  group_by(sites, functional_group) %>%
  summarise(freq = sum(value)) %>% 
  pivot_wider(names_from = functional_group, values_from = freq) %>% 
  left_join(comm_wt_means_spore[, c(1:6)])

plant_vs_animal_data$Animal_pathogen[which(is.na(plant_vs_animal_data$Animal_pathogen))] <- 0

Y <- cbind(plant_vs_animal_data$Animal_pathogen,
           plant_vs_animal_data$Plant_pathogen)

plant_vs_animal_m_all <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                              family = binomial,
                              data = plant_vs_animal_data)

plant_vs_animal_m_t_rt <-glmer(Y ~ Type + Ecosystem + (1|Location/plot),
                               family = binomial,
                               data = plant_vs_animal_data)

plant_vs_animal_m_rt <-glmer(Y ~ Ecosystem + (1|Location/plot),
                             family = binomial,
                             data = plant_vs_animal_data)

plant_vs_animal_m_t <-glmer(Y ~ Type + (1|Location/plot),
                            family = binomial,
                            data = plant_vs_animal_data)

plant_vs_animal_0 <-glmer(Y ~ 1 + (1|Location/plot),
                          family = binomial,
                          data = plant_vs_animal_data)

comp_mod_plant_animal <- rbind(
  anova(plant_vs_animal_0,
        plant_vs_animal_m_rt, method = "ML"),
  anova(plant_vs_animal_0,
        plant_vs_animal_m_t, method = "ML"),
  anova(plant_vs_animal_m_all,
        plant_vs_animal_m_t_rt,
        plant_vs_animal_0, method = "ML"))

comp_mod_plant_animal <- rapply(comp_mod_plant_animal, classes = "numeric", how = "replace",
                                function(x){round(x, 2)})
comp_mod_plant_animal$model <- row.names(comp_mod_plant_animal)
comp_mod_plant_animal <- comp_mod_plant_animal[c(1,4,2,6,7) , c(9, 1:8)]
write.table(comp_mod_plant_animal, "model_comparisons/6_comp_mod_plant_animal.txt",
            row.names = F, sep = ";")

write.table(rbind(comp_mod_free_sym, comp_mod_sapro_wsapro,
                  comp_mod_com_symb, comp_mod_mut_path,
                  comp_mod_pa_fp, comp_mod_plant_animal), "model_comparisons/all_comp.txt",
            row.names = F, sep = ";")



# For the tree
Y <- cbind(free_vs_symb_data$free_living,
           free_vs_symb_data$host_associated)

free_vs_symb_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                          family = binomial,
                          data = free_vs_symb_data)

#

Y <- cbind(sapro_vs_wsapro_data$Wood_Saprotroph,
           sapro_vs_wsapro_data$Saprotroph)

sapro_vs_wsapro_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                             family = binomial,
                             data = sapro_vs_wsapro_data)

#

Y <- cbind(com_vs_symb_data$Symbionts,
           com_vs_symb_data$Commensalists)

com_vs_symb_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                         family = binomial,
                         data = com_vs_symb_data)

#
Y <- cbind(mut_vs_path_data$Mutualists,
           mut_vs_path_data$Pathogens)

mut_vs_path_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                         family = binomial,
                         data = mut_vs_path_data)

#
Y <- cbind(pa_vs_fp_data$Path_Plant_Animals,
           pa_vs_fp_data$Path_Fungi_Protists)

pa_vs_fp_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                      family = binomial,
                      data = pa_vs_fp_data)

#
Y <- cbind(plant_vs_animal_data$Animal_pathogen,
           plant_vs_animal_data$Plant_pathogen)

plant_vs_animal_m_nu <-glmer(Y ~ Type*Ecosystem + (1|Location/plot),
                             family = binomial,
                             data = plant_vs_animal_data)




####
library(ggeffects)

CI_nu <- ggpredict(free_vs_symb_m_nu, c("Type", "Ecosystem"))
CI_nu <- as.data.frame(CI_nu)
names(CI_nu)[1] <- "Type"
CI_nu$comparison <- "prob_free_living"
CI_nu_1 <- CI_nu

Y <- cbind(free_vs_symb_data$host_associated, free_vs_symb_data$free_living)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = free_vs_symb_data)

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

Y <- cbind(sapro_vs_wsapro_data$Saprotroph, sapro_vs_wsapro_data$Wood_Saprotroph)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = sapro_vs_wsapro_data)

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

Y <- cbind(com_vs_symb_data$Commensalists, com_vs_symb_data$Symbionts)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = com_vs_symb_data)

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

Y <- cbind(mut_vs_path_data$Pathogens, mut_vs_path_data$Mutualists)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = mut_vs_path_data)

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

Y <- cbind(pa_vs_fp_data$Path_Fungi_Protists, pa_vs_fp_data$Path_Plant_Animals)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = pa_vs_fp_data)

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

Y <- cbind(plant_vs_animal_data$Plant_pathogen, plant_vs_animal_data$Animal_pathogen)
m <-glmer(Y ~ Type*Ecosystem + (1|Location/plot), family = binomial, data = plant_vs_animal_data)

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
# effect_size_nu$xend <- c(3, 27, 6, 0, 21, 30, 18, 24, 21, 27, 24, 18)
# effect_size_nu$yend <- c(18, 18, 15, 15, 15, 15, 12, 12, 6, 6, 0, 0)


effect_size_nu$xend <- c(22.366, 152.93, 38.834, 5.895, 114.698, 191.162, 75.503, 152.269, 114.378, 190.84, 152.53, 76.069)
effect_size_nu$yend <- c(43.189, 43.189, 78.846, 78.846, 78.846, 78.846, 114.503, 114.503, 150.160, 150.160, 185.817, 185.817)


effect_size_nu$x <- c(88.240, 88.240, 22.366, 22.366, 152.930, 152.930, 114.378, 114.378, 152.53, 152.53, 114.378, 114.378)
effect_size_nu$y <- c(6.749, 6.749, 43.189, 43.189, 43.189, 43.189, 78.846, 78.846, 114.503, 114.503, 150.160, 150.160)

# effect_size_nu$x <- c(15,15,3,3,27,27,21,21,24,24,21,21)
# effect_size_nu$y <- c(21,21,18,18,18,18,15,15,12,12,6,6)

effect_size_nu$icon_air <-  "Figures/building_air.png"
effect_size_nu$icon_air[which(effect_size_nu$effect_air>0)] <-  "Figures/tree_air.png"
#effect_size_nu$icon_air[which(effect_size_nu$effect_air>0)] <-  NA

effect_size_nu$icon_soil <-  "Figures/building_soil.png"
effect_size_nu$icon_soil[which(effect_size_nu$effect_air>0)] <-  "Figures/tree_soil.png"
#effect_size_nu$icon_soil[which(effect_size_nu$effect_air>0)] <-  NA

effect_size_nu$comparison <- gsub("prob_", "", effect_size_nu$comparison)

# effect_size_nu$group <- rep(LETTERS[1:6], each = 2)
# 
effect_size_nu$dir_air <- "#9ECAE1"# Blue, positive to cities
effect_size_nu$dir_air[which(effect_size_nu$effect_air>0)] <- "#FC9272"# Red, positive to natural

effect_size_nu$dir_soil <- "#3182BD"
effect_size_nu$dir_soil[which(effect_size_nu$effect_soil>0)] <- "#DE2D26"



m <- 2

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
# effect_size_nu$x_air[c(1,2)] <- c(14,16)
# effect_size_nu$xend_air[c(1,2)] <- c(2,28)
# effect_size_nu$x_soil[c(1,2)] <- c(16,14)
# effect_size_nu$xend_soil[c(1,2)] <- c(4,26)

effect_size_nu$comparison[6] <- "Commensalists"
effect_size_nu$comparison[9] <- "Pathogens of\nPlants & Animals"
effect_size_nu$comparison[10] <- "Pathogens of\nFungi & Protists"



pdf("Figures/decision_tree_PROTU_freq.pdf",
    width = 8.515,
    height = 11.015)

effect_size_nu %>% 
  bind_rows(data.frame(comparison = "Functional\ndiversity",
                       xend = 88.240,#"E"
                       yend = 6.749,
                       size = 0.1)) %>% 
  ggplot() +
  
  geom_segment(aes(x=x_air, xend=xend_air, y=y*-1, yend=yend*-1, size = abs(effect_air),
                   color = dir_air)) +
  
  geom_segment(aes(x=x_soil, xend=xend_soil, y=y*-1, yend=yend*-1, size = abs(effect_soil),
                   color = dir_soil))+
  
  geom_point(aes(xend,yend*-1), size=20, shape=16, col= "gray") +
  
  geom_text(aes(xend, yend*-1, label = comparison, fontface = "bold")) +
  
  scale_color_identity() +
  guides(size = guide_legend(title="Effect size")) +
  
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank())

dev.off()




library(ggimage)
effect_size_nu %>% 
  bind_rows(data.frame(comparison = "Functional\ndiversity",
                       x = 15,#"E"
                       y = 21,
                       size = 0.1)) %>% 
  ggplot() +
  geom_point(aes(x,y), size=20, shape=16, col= "gray")+
  
  geom_image(aes(x-1.5, y-1, image = icon_air, size = I(abs(effect_air))), by = "height") +
  geom_image(aes(x+1.5, y-1, image = icon_soil, size = I(abs(effect_soil))), by = "height") +
  geom_text(aes(x, y, label = comparison, fontface = "bold")) +
  
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Another approach

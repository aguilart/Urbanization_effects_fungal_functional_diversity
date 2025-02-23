
######################################################################################################
##########  DATA ANALYSIS SEXUAL SPORE SIZE SHIFTS (ABUNDANCE WEIGHTED)   #############################################
######################################################################################################


library(lme4)
library(lmerTest)
library(performance)
library(tidyverse)
library(ggdist)


get_table <- 
  function(x){ z <- summary(x)
  y = z$coefficients
  
  colnames(y) <- gsub("\\.", "_", colnames(y))
  colnames(y) <- gsub(" ", "_", colnames(y))
  colnames(y) <- gsub("Pr\\(>\\|t\\|\\)", "p_value", colnames(y))
  colnames(y) <- tolower(colnames(y))
  
  y = as.data.frame(y)
  y
  
  }


# Load data in long format: each genera present per site and ecosystem  and their
# associated functional guild and spore size

spore_function_long <- readRDS("spore_function_long.RDS")



lichen_path_air_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Lichen_Parasite == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

lichen_path_soil_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Lichen_Parasite == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

# plant pathogens

plant_path_air_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Plant_Pathogen == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

plant_path_soil_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Plant_Pathogen == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

# human pathogens

human_path_air_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Human_Pathogen == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

human_path_soil_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Human_Pathogen == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

# Wood saprotrophs

wood_sapro_air_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Wood_Saprotroph == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

wood_sapro_soil_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Wood_Saprotroph == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

# Saprotrophs

sapro_air_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Saprotroph == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

sapro_soil_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Saprotroph == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

# Ectomycorrhizal

ecto_air_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Ectomycorrhizal == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

ecto_soil_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Ectomycorrhizal == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

# Endophyte

endo_air_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Endophyte == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

endo_soil_sex_spore_abund <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Endophyte == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T))

tabla_sex_spores_abund <- 
  bind_rows(
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Saprotrophs"),
    get_table(sapro_air_sex_spore_abund),
    get_table(sapro_soil_sex_spore_abund),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Wood Saprotrophs"),
    get_table(wood_sapro_air_sex_spore_abund),
    get_table(wood_sapro_soil_sex_spore_abund),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Endophytes"),
    get_table(endo_air_sex_spore_abund),
    get_table(endo_soil_sex_spore_abund),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Ectomycorrhizas"),
    get_table(ecto_air_sex_spore_abund),
    get_table(ecto_soil_sex_spore_abund),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Lichen pathogens"),
    get_table(lichen_path_air_sex_spore_abund),
    get_table(lichen_path_soil_sex_spore_abund),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Human pathogens"),
    get_table(human_path_air_sex_spore_abund),
    get_table(human_path_soil_sex_spore_abund),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Plant pathogens"),
    get_table(plant_path_air_sex_spore_abund),
    get_table(plant_path_soil_sex_spore_abund))


tabla_sex_spores_abund <- tabla_sex_spores_abund[-grep("Intercept", rownames(tabla_sex_spores_abund)), ]
tabla_sex_spores_abund$model <- rownames(tabla_sex_spores_abund)
tabla_sex_spores_abund <- tabla_sex_spores_abund[, c(6, 1:5)]
rownames(tabla_sex_spores_abund) <- NULL


m <- which(!grepl("Ecosystem", tabla_sex_spores_abund$model))
tabla_sex_spores_abund$model[m+1] <- "Urban air - Natural air"
tabla_sex_spores_abund$model[m+2] <- "Urban soil - Natural soil"

tabla_sex_spores_abund <-
  bind_cols(tabla_sex_spores_abund,
            data.frame(r2 = c(
              NA,
              r2_nakagawa(sapro_air_sex_spore_abund)$R2_marginal,
              r2_nakagawa(sapro_soil_sex_spore_abund)$R2_marginal,
              NA,
              r2_nakagawa(wood_sapro_air_sex_spore_abund)$R2_marginal,
              r2_nakagawa(wood_sapro_soil_sex_spore_abund)$R2_marginal,
              NA,
              r2_nakagawa(endo_air_sex_spore_abund)$R2_marginal,
              r2_nakagawa(endo_soil_sex_spore_abund)$R2_marginal,
              NA,
              r2_nakagawa(ecto_air_sex_spore_abund)$R2_marginal,
              r2_nakagawa(ecto_soil_sex_spore_abund)$R2_marginal,
              NA,
              r2_nakagawa(lichen_path_air_sex_spore_abund)$R2_marginal,
              r2_nakagawa(lichen_path_soil_sex_spore_abund)$R2_marginal,
              NA,
              r2_nakagawa(human_path_air_sex_spore_abund)$R2_marginal,
              r2_nakagawa(human_path_soil_sex_spore_abund)$R2_marginal,
              NA,
              r2_nakagawa(plant_path_air_sex_spore_abund)$R2_marginal,
              r2_nakagawa(plant_path_soil_sex_spore_abund)$R2_marginal)
            ))

tabla_sex_spores_abund <- rapply(tabla_sex_spores_abund, classes = "numeric", how = "replace",
                                 function(x){round(x, 2)})

tabla_sex_spores_abund <- tabla_sex_spores_abund[, c(1, 7, 2:6)]

write.table(tabla_sex_spores_abund,
            "results_sex_spores_abund.txt", sep =";",
            row.names = F)


pdf("Figures/boxplots_spores_mean_abund.pdf",
   width = 11.015 ,
   height = 17.03)#,   units = "in", res = 300)
cowplot::plot_grid(
bind_rows(
  spore_function_long %>% 
    filter(Saprotroph == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Saprotroph"),
  
  spore_function_long %>% 
    filter(Wood_Saprotroph == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Wood_Saprotroph"),
  
  spore_function_long %>% 
    filter(Endophyte == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Endophyte"),
  
  spore_function_long %>% 
    filter(Ectomycorrhizal == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Ectomycorrhizal"),
  
  spore_function_long %>% 
    filter(Lichen_Parasite == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Lichen_Parasite"),
  
  spore_function_long %>% 
    filter(Human_Pathogen == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Human_Pathogen"),
  
  spore_function_long %>% 
    filter(Plant_Pathogen == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Plant_Pathogen")) %>% 
  
  ggplot() +
  aes(x = Ecosystem, y = sexual_spore_volume_wt, fill = Ecosystem) +
  
  # geom_boxplot(alpha=0.9, outlier.shape =NA) + 
  # geom_jitter(position = position_jitter(seed = 1)) +
  
  geom_boxplot(alpha=0.9, outlier.shape =NA) + 
  geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  
  labs(y = "Community-weighted sexual spore mean",
       title = "Soil samples weighted mean (abundance)") +
  scale_fill_manual(values = c("#5CB85Cff","#EEA236FF" )) +
 
  facet_wrap(.~guild, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 20, face = "bold")),


### Air

bind_rows(
  spore_function_long %>% 
    filter(Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Saprotroph"),
  
  spore_function_long %>% 
    filter(Wood_Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Wood_Saprotroph"),
  
  spore_function_long %>% 
    filter(Endophyte == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Endophyte"),
  
  spore_function_long %>% 
    filter(Ectomycorrhizal == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Ectomycorrhizal"),
  
  spore_function_long %>% 
    filter(Lichen_Parasite == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Lichen_Parasite"),
  
  spore_function_long %>% 
    filter(Human_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Human_Pathogen"),
  
  spore_function_long %>% 
    filter(Plant_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), sum, na.rm = T) %>% 
    mutate(guild = "Plant_Pathogen")) %>% 
  
  ggplot() +
  aes(x = Ecosystem, y = sexual_spore_volume_wt, fill = Ecosystem) +
  
  # geom_boxplot(alpha=0.9, outlier.shape =NA) + 
  # geom_jitter(position = position_jitter(seed = 1)) +
  
  geom_boxplot(alpha=0.9, outlier.shape =NA) + 
  geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  scale_fill_manual(values = c("#5CB85Cff","#EEA236FF" )) +
  
  labs(y = "Community-weighted sexual spore mean",
       title = "Air samples weighted mean (abundance)") +
  facet_wrap(.~guild, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 20, face = "bold")),

nrow = 2, ncol = 1, align = "v")
dev.off()

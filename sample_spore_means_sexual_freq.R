
library(lme4)
library(lmerTest)
library(performance)
library(tidyverse)


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


#### Analysis with frequencies: sexual spores ####
# Lichens

lichen_path_air_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Lichen_Parasite == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))


lichen_path_soil_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Lichen_Parasite == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

# plant pathogens

plant_path_air_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Plant_Pathogen == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

plant_path_soil_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Plant_Pathogen == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

# human pathogens

human_path_air_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Human_Pathogen == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

human_path_soil_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Human_Pathogen == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))


# Wood saprotrophs

wood_sapro_air_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Wood_Saprotroph == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

wood_sapro_soil_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Wood_Saprotroph == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

# Saprotrophs

sapro_air_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Saprotroph == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

sapro_soil_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Saprotroph == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

# Ectomycorrhizal

ecto_air_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Ectomycorrhizal == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

ecto_soil_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Ectomycorrhizal == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

# Endophyte

endo_air_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Endophyte == 1) %>% 
         filter(Type == "Air") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

endo_soil_sex_spore_freq <-
  lmer(sexual_spore_volume_wt ~ Ecosystem + (1|Location/plot),
       data = spore_function_long %>% 
         filter(Endophyte == 1) %>% 
         filter(Type == "Soil") %>% 
         #filter(genus != "Calonectria") %>% 
         filter(!is.na(sexual_spore_volume)) %>% 
         mutate(asexual_spore_volume_wt = frequency*log10(asexual_spore_volume)) %>% 
         mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
         group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
         summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), mean, na.rm = T))

###


tabla_sex_spores_freq <- 
  bind_rows(
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Saprotrophs"),
    get_table(sapro_air_sex_spore_freq),
    get_table(sapro_soil_sex_spore_freq),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Wood Saprotrophs"),
    get_table(wood_sapro_air_sex_spore_freq),
    get_table(wood_sapro_soil_sex_spore_freq),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Endophytes"),
    get_table(endo_air_sex_spore_freq),
    get_table(endo_soil_sex_spore_freq),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Ectomycorrhizas"),
    get_table(ecto_air_sex_spore_freq),
    get_table(ecto_soil_sex_spore_freq),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Lichen pathogens"),
    get_table(lichen_path_air_sex_spore_freq),
    get_table(lichen_path_soil_sex_spore_freq),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Human pathogens"),
    get_table(human_path_air_sex_spore_freq),
    get_table(human_path_soil_sex_spore_freq),
    data.frame(estimate = NA, std__error = NA, df = NA, t_value= NA, p_value = NA, row.names = "Plant pathogens"),
    get_table(plant_path_air_sex_spore_freq),
    get_table(plant_path_soil_sex_spore_freq))


tabla_sex_spores_freq <- tabla_sex_spores_freq[-grep("Intercept", rownames(tabla_sex_spores_freq)), ]
tabla_sex_spores_freq$model <- rownames(tabla_sex_spores_freq)
tabla_sex_spores_freq <- tabla_sex_spores_freq[, c(6, 1:5)]
rownames(tabla_sex_spores_freq) <- NULL


m <- which(!grepl("Ecosystem", tabla_sex_spores_freq$model))
tabla_sex_spores_freq$model[m+1] <- "Urban air - Natural air"
tabla_sex_spores_freq$model[m+2] <- "Urban soil - Natural soil"

tabla_sex_spores_freq <-
  bind_cols(tabla_sex_spores_freq,
            data.frame(r2 = c(
              NA,
              r2_nakagawa(sapro_air_sex_spore_freq)$R2_marginal,
              r2_nakagawa(sapro_soil_sex_spore_freq)$R2_marginal,
              NA,
              r2_nakagawa(wood_sapro_air_sex_spore_freq)$R2_marginal,
              r2_nakagawa(wood_sapro_soil_sex_spore_freq)$R2_marginal,
              NA,
              r2_nakagawa(endo_air_sex_spore_freq)$R2_marginal,
              r2_nakagawa(endo_soil_sex_spore_freq)$R2_marginal,
              NA,
              r2_nakagawa(ecto_air_sex_spore_freq)$R2_marginal,
              r2_nakagawa(ecto_soil_sex_spore_freq)$R2_marginal,
              NA,
              r2_nakagawa(lichen_path_air_sex_spore_freq)$R2_marginal,
              r2_nakagawa(lichen_path_soil_sex_spore_freq)$R2_marginal,
              NA,
              r2_nakagawa(human_path_air_sex_spore_freq)$R2_marginal,
              r2_nakagawa(human_path_soil_sex_spore_freq)$R2_marginal,
              NA,
              r2_nakagawa(plant_path_air_sex_spore_freq)$R2_marginal,
              r2_nakagawa(plant_path_soil_sex_spore_freq)$R2_marginal)
            ))

tabla_sex_spores_freq <- rapply(tabla_sex_spores_freq, classes = "numeric", how = "replace",
                                function(x){round(x, 2)})

tabla_sex_spores_freq <- tabla_sex_spores_freq[, c(1, 7, 2:6)]

write.table(tabla_sex_spores_freq,
            "results_sex_spores_freq.txt", sep =";",
            row.names = F)

my_theme <-
theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"))


library(ggsci)

pdf("Figures/halfviolin_spores_mean_freq_soil.pdf",
    width = 12,
    height = 8)

cowplot::plot_grid(
bind_rows(
spore_function_long %>% 
  filter(Saprotroph == 1) %>% 
  filter(Type == "Soil") %>% 
  filter(!is.na(sexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
  summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
  mutate(guild = "Saprotroph"),

spore_function_long %>% 
  filter(Wood_Saprotroph == 1) %>% 
  filter(Type == "Soil") %>% 
  filter(!is.na(sexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
  summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
  mutate(guild = "Wood_Saprotroph")) %>% 
  
  ggplot() +
  aes(x = Ecosystem, y = sexual_spore_volume_wt, fill = Ecosystem) +
  #geom_violin(trim = FALSE, alpha=0.9) + theme_classic()+
  # geom_boxplot(alpha=0.9) + 
  stat_halfeye(aes(thickness = after_stat(pdf*n)), scale = 0.5) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA) +
  #geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  labs(y = "Community-wide sexual spore mean",
       title = "Free-living") +
  scale_fill_manual(values = c("#5CB85Cff","#EEA236FF" )) +
  facet_wrap(.~guild, scales = "free", nrow = 1) +
  my_theme,

bind_rows(
spore_function_long %>% 
  filter(Endophyte == 1) %>% 
  filter(Type == "Soil") %>% 
  filter(!is.na(sexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
  summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
  mutate(guild = "Endophyte"),

spore_function_long %>% 
  filter(Ectomycorrhizal == 1) %>% 
  filter(Type == "Soil") %>% 
  filter(!is.na(sexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
  summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
  mutate(guild = "Ectomycorrhizal")) %>% 
  
ggplot() +
  aes(x = Ecosystem, y = sexual_spore_volume_wt, fill = Ecosystem) +
  #geom_violin(trim = FALSE, alpha=0.9) + theme_classic()+
  # geom_boxplot(alpha=0.9) + 
  stat_halfeye(aes(thickness = after_stat(pdf*n)), scale = 0.5) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA) +
  #geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  labs(y = "Community-wide sexual spore mean",
       title = "Mutualists") +
  scale_fill_manual(values = c("#5CB85Cff","#EEA236FF" )) +
  facet_wrap(.~guild, scales = "free", nrow = 1) +
  my_theme,

bind_rows(
spore_function_long %>% 
  filter(Lichen_Parasite == 1) %>% 
  filter(Type == "Soil") %>% 
  filter(!is.na(sexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
  summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
  mutate(guild = "Lichen_Parasite"),

spore_function_long %>% 
  filter(Human_Pathogen == 1) %>% 
  filter(Type == "Soil") %>% 
  filter(!is.na(sexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
  summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
  mutate(guild = "Human_Pathogen"),

spore_function_long %>% 
  filter(Plant_Pathogen == 1) %>% 
  filter(Type == "Soil") %>% 
  filter(!is.na(sexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
  summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
  mutate(guild = "Plant_Pathogen")) %>% 

  ggplot() +
  aes(x = Ecosystem, y = sexual_spore_volume_wt, fill = Ecosystem) +
  #geom_violin(trim = FALSE, alpha=0.9) + theme_classic()+
 # geom_boxplot(alpha=0.9) + 
  stat_halfeye(aes(thickness = after_stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA) +
  #geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  
  #geom_text(aes(label = sites, color = Location), position = position_jitter(seed = 1)) +
  #scale_fill_manual(values=c("steelblue", "orange")) +
  labs(y = "Community-wide spore mean",
       title = "Pathogens") +
  scale_fill_manual(values = c("#5CB85Cff","#EEA236FF" )) +
  facet_wrap(.~guild, scales = "free", nrow = 1) +
  my_theme,
nrow = 3, ncol = 1, align = "v")

dev.off()

###

pdf("Figures/halfviolin_spores_mean_freq_air.pdf",
    width = 11.015 ,
    height = 8.515)
bind_rows(
  spore_function_long %>% 
    filter(Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Saprotroph"),
  
  spore_function_long %>% 
    filter(Wood_Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Wood_Saprotroph"),
  
  spore_function_long %>% 
    filter(Endophyte == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Endophyte"),
  
  spore_function_long %>% 
    filter(Ectomycorrhizal == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Ectomycorrhizal"),
  
  spore_function_long %>% 
    filter(Lichen_Parasite == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Lichen_Parasite"),
  
  spore_function_long %>% 
    filter(Human_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Human_Pathogen"),
  
  spore_function_long %>% 
    filter(Plant_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Plant_Pathogen")) %>% 
  
  ggplot() +
  aes(x = Ecosystem, y = sexual_spore_volume_wt, fill = Ecosystem) +
  #geom_violin(trim = FALSE, alpha=0.9) + theme_classic()+
  #geom_boxplot(alpha=0.9) + 
  
  stat_halfeye() +
  geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  
  #geom_text(aes(label = sites, color = Location), position = position_jitter(seed = 1)) +
  #scale_fill_manual(values=c("steelblue", "orange")) +
  labs(y = "Community-wide spore mean",
       title = "Air samples based on genus presence/absence") +
  scale_fill_manual(values = c("#5CB85Cff","#EEA236FF" )) +
  facet_wrap(.~guild, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"))
dev.off()

# As png

png("Figures/halfviolin_spores_mean_freq_air.png",
    width = 11.015 ,
    height = 8.515,
    units = "in", res = 300)
bind_rows(
  spore_function_long %>% 
    filter(Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Saprotroph"),
  
  spore_function_long %>% 
    filter(Wood_Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Wood_Saprotroph"),
  
  spore_function_long %>% 
    filter(Endophyte == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Endophyte"),
  
  spore_function_long %>% 
    filter(Ectomycorrhizal == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Ectomycorrhizal"),
  
  spore_function_long %>% 
    filter(Lichen_Parasite == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Lichen_Parasite"),
  
  spore_function_long %>% 
    filter(Human_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Human_Pathogen"),
  
  spore_function_long %>% 
    filter(Plant_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(sexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
    group_by(Type, Location, Ecosystem, plot, sample, sites) %>%
    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T) %>% 
    mutate(guild = "Plant_Pathogen")) %>% 
  
  ggplot() +
  aes(x = Ecosystem, y = sexual_spore_volume_wt, fill = Ecosystem) +
  #geom_violin(trim = FALSE, alpha=0.9) + theme_classic()+
  #geom_boxplot(alpha=0.9) + 
  
  stat_halfeye() +
  geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  
  #geom_text(aes(label = sites, color = Location), position = position_jitter(seed = 1)) +
  #scale_fill_manual(values=c("steelblue", "orange")) +
  labs(y = "Community-wide spore mean",
       title = "Air samples based on genus presence/absence") +
  scale_fill_manual(values = c("#5CB85Cff","#EEA236FF" )) +
  facet_wrap(.~guild, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        title = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 20, face = "bold"))
dev.off()



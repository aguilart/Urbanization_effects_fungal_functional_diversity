
######################################################################################################
#####################  DATA ANALYSIS BIOTIC HOMOGENIZATION   #############################################
######################################################################################################

# Loading OTU table in long format with the genus asociated and functional group

site_funct_group_count <- readRDS("site_funct_group_count.RDS")

# Loading libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(metafor)

functional_richness <-    
  left_join(
    site_funct_group_count[,-1] %>% 
      #site_funct_group_count_c %>% 
      filter(!is.na(functional_group)) %>% 
      distinct(),
    samples[,c("ID", "Location", "plot")] %>% rename(sites = ID)) %>% 
  select(genus, Type, Ecosystem, functional_group, Location, plot) %>% 
  distinct()

functional_richness$Location[-grep("Hel|Joen|Tamp|Laht", functional_richness$Location)] <- "Jyvaskyla"

functional_richness_wide <-
  functional_richness %>% 
  group_by(Location, Type, Ecosystem, functional_group, plot) %>% 
  tally() %>% 
  pivot_wider(names_from = functional_group, values_from = n, values_fill = 0)

soil <- which(functional_richness_wide$Type == "Soil")
air <- which(functional_richness_wide$Type == "Air")

x <- functional_richness_wide[soil,]


dis <- vegdist(x[, c(5:length(names(x)))], method = "euclidean")
mod <- betadisper(dis, x$Ecosystem)
#mod <-anova(mod)
mod

x <- functional_richness_wide[air,]

dis <- vegdist(x[, c(5:length(names(x)))], method = "euclidean")
mod <- betadisper(dis, x$Ecosystem)
#mod <-anova(mod)
mod


functional_betadisper_soil <-
  lapply(
    split(functional_richness_wide[soil,], functional_richness_wide$Location[soil]),
    function(x){
      dis <- vegdist(x[, c(5:length(names(x)))], method = "euclidean")
      mod <- betadisper(dis, x$Ecosystem)
      #mod <-anova(mod)
      mod
    })
#
results_beta_soil_within <- 
  bind_rows(lapply(functional_betadisper_soil, anova), .id = "id")
results_beta_soil_within <- as.data.frame(results_beta_soil_within[c(1,3,5,7,9), ])

results_beta_soil_within <- rapply(results_beta_soil_within,
                                   classes = "numeric", how = "replace",
                                   function(x){round(x, 2)})
write.table(results_beta_soil_within,
            "model_comparisons/results_beta_soil_within.txt",
            row.names = F, sep = ";")
#
results_beta_air_within <- 
  bind_rows(lapply(functional_betadisper_air, anova), .id = "id")
results_beta_air_within <- as.data.frame(results_beta_air_within[c(1,3,5,7,9), ])

results_beta_air_within <- rapply(results_beta_air_within,
                                  classes = "numeric", how = "replace",
                                  function(x){round(x, 2)})
write.table(results_beta_air_within,
            "model_comparisons/results_beta_air_within.txt",
            row.names = F, sep = ";")



functional_betadisper_air <-
  lapply(
    split(functional_richness_wide[air,], functional_richness_wide$Location[air]),
    function(x){
      dis <- vegdist(x[, c(5:length(names(x)))], method = "euclidean")
      mod <- betadisper(dis, x$Ecosystem)
      # mod <-anova(mod)
      mod
    })

#
results_beta_soil_across <- 
  bind_rows(lapply(functional_betadisper_soil_across, anova), .id = "id")
results_beta_soil_across <- as.data.frame(results_beta_soil_across[c(1,3), ])

results_beta_soil_across <- rapply(results_beta_soil_across,
                                   classes = "numeric", how = "replace",
                                   function(x){round(x, 2)})
write.table(results_beta_soil_across,
            "model_comparisons/results_beta_soil_across.txt",
            row.names = F, sep = ";")
#
results_beta_air_across <- 
  bind_rows(lapply(functional_betadisper_air_across, anova), .id = "id")
results_beta_air_across <- as.data.frame(results_beta_air_across[c(1,3), ])

results_beta_air_across <- rapply(results_beta_air_across,
                                  classes = "numeric", how = "replace",
                                  function(x){round(x, 2)})
write.table(results_beta_air_across,
            "model_comparisons/results_beta_air_across.txt",
            row.names = F, sep = ";")


pdf("Figures/Barplot_funct_freq_genus.pdf",
    width = 12, height = 6)
distinct(site_funct_group_count[, c(3:6)]) %>% 
  filter(!is.na(functional_group)) %>% 
  group_by(Type, Ecosystem, functional_group) %>% 
  tally() %>%
  mutate(freq = 100*(n / sum(n))) %>% 
  filter(freq > 0.1) %>% 
  ggplot() +
  aes(y = freq, fill = Ecosystem, reorder_within(x = functional_group, -freq, list(Type))) + 
  scale_fill_brewer(palette = "Dark2") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(.~Type, scales = "free",drop = T) +
  labs(y= "Percentage of functional groups",
       title = "Presence/absence of functional groups per Ecosystem based on number of genera")+
  scale_x_reordered() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = "bottom")
dev.off()

png("Figures/Barplot_funct_freq_genus.png",
    width = 16,
    height = 8,
    units = "in", res = 300)
distinct(site_funct_group_count[, c(3:6)]) %>% 
  filter(!is.na(functional_group)) %>% 
  group_by(Type, Ecosystem, functional_group) %>% 
  tally() %>%
  mutate(freq = 100*(n / sum(n))) %>% 
  filter(freq > 0.1) %>% 
  ggplot() +
  aes(y = freq, fill = Ecosystem, reorder_within(x = functional_group, -freq, list(Type))) + 
  scale_fill_brewer(palette = "Dark2") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(.~Type, scales = "free",drop = T) +
  labs(y= "Percentage of functional groups",
       title = "Presence/absence of functional groups per Ecosystem based on number of genera")+
  scale_x_reordered() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = "bottom")
dev.off()


#####################################################################################################
library(ape)
library(nlme)

spore_function_long$plot_simple <- 
  gsub("HEL-|JOE-|JYV-|LAH-|TAM-","", spore_function_long$plot)

spore_function_long$site_simple <- "A"
spore_function_long$site_simple[grepl("-N", spore_function_long$sites)&grepl("180819", spore_function_long$sites)] <- "B"
spore_function_long$site_simple[grepl("-N", spore_function_long$sites)&grepl("220819|230819", spore_function_long$sites)] <- "C"

spore_function_long$site_simple[grepl("-U", spore_function_long$sites)&grepl("200819", spore_function_long$sites)] <- "B"
spore_function_long$site_simple[grepl("-U", spore_function_long$sites)&grepl("210819", spore_function_long$sites)] <- "C"
samples$ID[grepl("HEL-N1", samples$ID)&grepl("S", samples$ID)]
samples$ID[grepl("HEL-N2", samples$ID)&grepl("S", samples$ID)]
samples$ID[grepl("HEL-N3", samples$ID)&grepl("S", samples$ID)]


unique(str_extract(samples$ID[grep("U", samples$ID)], "\\-[[:digit:]]+\\-"))

sex_spore_var_comp <- 
  spore_function_long[-which(is.na(spore_function_long$sexual_spore_volume)), ]

nat_s <- which(sex_spore_var_comp$Ecosystem == "Natural"&
                 sex_spore_var_comp$Type == "Soil")

urb_s <- which(sex_spore_var_comp$Ecosystem == "Urban"&
                 sex_spore_var_comp$Type == "Soil")

for_within_homogeneity <-
  left_join(
    sex_spore_var_comp[nat_s,] %>% 
      # filter(Saprotroph == 1) %>% 
      # filter(Wood_Saprotroph == 1) %>% 
      # filter(Plant_Pathogen == 1) %>% 
      # filter(Lichen_Parasite == 1) %>% 
      # filter(Endophyte == 1) %>% 
      # filter(Ectomycorrhizal == 1) %>% 
      filter(Human_Pathogen == 1) %>% 
      mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
      #group_by(Location, plot_simple) %>%  
      select(genus, Location, sexual_spore_volume_wt) %>% distinct %>% 
      group_by(Location) %>%
      summarise(mean_size_natural = mean(sexual_spore_volume_wt, na.rm = T),
                sd_size_natural = sd(sexual_spore_volume_wt, na.rm = T),
                n_natural = n()),# %>% 
    #mutate(plot_simple = gsub("N","P", plot_simple), ID = paste(Location, plot_simple, sep = "_")),
    
    sex_spore_var_comp[urb_s,] %>% 
      # filter(Saprotroph == 1) %>% 
      # filter(Wood_Saprotroph == 1) %>% 
      # filter(Plant_Pathogen == 1) %>% 
      # filter(Lichen_Parasite == 1) %>% 
      # filter(Endophyte == 1) %>% 
      # filter(Ectomycorrhizal == 1) %>% 
      filter(Human_Pathogen == 1) %>% 
      mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
      #group_by(Location, plot_simple) %>%  
      select(genus, Location, sexual_spore_volume_wt) %>% distinct %>% 
      group_by(Location) %>%
      summarise(mean_size_urban = mean(sexual_spore_volume_wt, na.rm = T),
                sd_size_urban = sd(sexual_spore_volume_wt, na.rm = T),
                n_urban = n()) %>% 
      #mutate(plot_simple = gsub("U","P", plot_simple), ID = paste(Location, plot_simple, sep = "_")) %>% ungroup() %>% select(ID, mean_size_urban, sd_size_urban, n_urban), by = "ID")
      ungroup() %>% select(Location, mean_size_urban, sd_size_urban, n_urban), by = "Location")

# Making the test
for_within_homogeneity$CVR <-
  Calc.lnCVR(CMean = for_within_homogeneity$mean_size_urban,
             CSD = for_within_homogeneity$sd_size_urban,
             CN = for_within_homogeneity$n_urban,
             EMean = for_within_homogeneity$mean_size_natural,
             ESD = for_within_homogeneity$sd_size_natural,
             EN = for_within_homogeneity$n_natural)

# Calculate the sampling error for lnCVR using 
# the function Calc.var.lnCVR.
for_within_homogeneity$V.CVR<-Calc.var.lnCVR(CMean = for_within_homogeneity$mean_size_urban,
                                             CSD = for_within_homogeneity$sd_size_urban,
                                             CN = for_within_homogeneity$n_urban,
                                             EMean = for_within_homogeneity$mean_size_natural,
                                             ESD = for_within_homogeneity$sd_size_natural,
                                             EN = for_within_homogeneity$n_natural,
                                             Equal.E.C.Corr = F)

response <- for_within_homogeneity$CVR
var <- for_within_homogeneity$V.CVR

# Run the meta-analysis.
model <- rma(response, vi=var)

# Observe the models results.
summary(model)

# model_sapro <- model # This one is not significant, still urban is more variable than natural
# model_ws <- model # This one is not significant, the only one where urban is less variable than urban
# model_pa <- model # This one is not sigificant, still urban is more variable than natural
# model_lp <- model # This is also significant, urban is more variable than natural # Only one that remain significant after only consider city wide comparisons (no plot- plot or site-site within city comparisons)
# model_endo <- model # This is also significant, urban is more variable than natural
#model_ecto <- model # This one is signficant, urban is more variable than natural

model_hp <- model # This one is signficant, urban is more variable than natural



getting_values <- function(x){
  y <- data.frame(estimate = x$beta,
                  upper = x$ci.ub,
                  lower = x$ci.lb)
  y}

variation_within <- bind_rows(lapply(
  list (model_hp,
        model_lp, model_pa,
        model_ecto, model_endo,
        model_sapro, model_ws), getting_values))

variation_within$group <- c("1_Human pathogens",
                            "2_Lichen parasites",
                            "3_Plant pathogens",
                            "4_Ectomycorrhizal",
                            "5_Endophytes",
                            "6_Saprotrophs",
                            "7_Wood saprotrophs")

variation_within <- 
  rapply(variation_within, classes = "numeric", how = "replace",
         function(x){round(x, 2)})


variation_within <- variation_within[,c(4, 1:3)]
names(variation_within)[2] <- "estimated lnCVR"

write.table(variation_within,
            "results_within_city_variation.txt", sep =";",
            row.names = F)

pdf("Figures/forest_plot_homogeneity.pdf",
    width = 11.015 ,
    height = 8.515)
variation_within %>% 
  ggplot() +
  aes(x = `estimated lnCVR`, y = group) +
  geom_point(size = 4) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 2) +
  labs(x = "ln(Coeff. Urban / Coeff. Natural)", y = "Functional groups") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))
dev.off()

#####################################################################################################

# Across city variation in urban relative to natural
library(car)

across_cities_variation <- 
  bind_rows(list(
    # Ectomycorrhizal
    leveneTest(lm(sexual_spore_volume_wt ~ Ecosystem,
                  data = spore_function_long %>% 
                    filter(Ectomycorrhizal == 1) %>% 
                    filter(Type == "Soil") %>% 
                    #filter(genus != "Calonectria") %>% 
                    filter(!is.na(sexual_spore_volume)) %>% 
                    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
                    group_by(Type, Location, Ecosystem) %>%
                    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T))),
    # Endophyte
    leveneTest(lm(sexual_spore_volume_wt ~ Ecosystem,
                  data = spore_function_long %>% 
                    filter(Endophyte == 1) %>% 
                    filter(Type == "Soil") %>% 
                    #filter(genus != "Calonectria") %>% 
                    filter(!is.na(sexual_spore_volume)) %>% 
                    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
                    group_by(Type, Location, Ecosystem) %>%
                    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T))),
    
    # Lichen parasite
    leveneTest(lm(sexual_spore_volume_wt ~ Ecosystem,
                  data = spore_function_long %>% 
                    filter(Lichen_Parasite == 1) %>% 
                    filter(Type == "Soil") %>% 
                    #filter(genus != "Calonectria") %>% 
                    filter(!is.na(sexual_spore_volume)) %>% 
                    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
                    group_by(Type, Location, Ecosystem) %>%
                    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T))),
    
    # Human pathogen
    leveneTest(lm(sexual_spore_volume_wt ~ Ecosystem,
                  data = spore_function_long %>% 
                    filter(Human_Pathogen == 1) %>% 
                    filter(Type == "Soil") %>% 
                    #filter(genus != "Calonectria") %>% 
                    filter(!is.na(sexual_spore_volume)) %>% 
                    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
                    group_by(Type, Location, Ecosystem) %>%
                    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T))),
    
    # Plant pathogen
    leveneTest(lm(sexual_spore_volume_wt ~ Ecosystem,
                  data = spore_function_long %>% 
                    filter(Plant_Pathogen == 1) %>% 
                    filter(Type == "Soil") %>% 
                    #filter(genus != "Calonectria") %>% 
                    filter(!is.na(sexual_spore_volume)) %>% 
                    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
                    group_by(Type, Location, Ecosystem) %>%
                    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T))),
    
    # Saprotroph
    leveneTest(lm(sexual_spore_volume_wt ~ Ecosystem,
                  data = spore_function_long %>% 
                    filter(Saprotroph == 1) %>% 
                    filter(Type == "Soil") %>% 
                    #filter(genus != "Calonectria") %>% 
                    filter(!is.na(sexual_spore_volume)) %>% 
                    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
                    group_by(Type, Location, Ecosystem) %>%
                    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T))),
    
    #Wood saprotroph
    leveneTest(lm(sexual_spore_volume_wt ~ Ecosystem,
                  data = spore_function_long %>% 
                    filter(Wood_Saprotroph == 1) %>% 
                    filter(Type == "Soil") %>% 
                    #filter(genus != "Calonectria") %>% 
                    filter(!is.na(sexual_spore_volume)) %>% 
                    mutate(sexual_spore_volume_wt = frequency*log10(sexual_spore_volume)) %>% 
                    group_by(Type, Location, Ecosystem) %>%
                    summarise_at(c("sexual_spore_volume_wt"), mean, na.rm = T))) )  )


#Printing the results of this tests:
rownames(across_cities_variation) <- NULL
across_cities_variation <- as.data.frame(across_cities_variation[c(1,3,5,7,9,11,13), ])
across_cities_variation$group <- c("Ectomycorrhizal", "Endophyte", "Lichen parasite", "Human pathogen",
                                   "Plant pathogen", "Saprotroph", "Wood saprotroph")

across_cities_variation <-
  rapply(across_cities_variation, classes = "numeric", how = "replace",
         function(x){round(x, 2)})

across_cities_variation <- across_cities_variation[, c(4,2,3)]

write.table(across_cities_variation,
            "results_across_city_variation.txt", sep =";",
            row.names = F)

####################################################################################################



#Sexual

# Saprotrophs
png("Figures/heatmap_genus_abund_asex_sapro.png",
    width = 11.015 ,
    height = 17.03,
    units = "in", res = 300)
cowplot::plot_grid(
  spore_function_long %>% 
    filter(Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
    geom_tile(aes(fill = asexual_spore_volume_wt)) +
    scale_size(range = c(0.5, 10)) +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Saprotroph abundance air") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  
  spore_function_long %>% 
    filter(Saprotroph == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
    geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Saprotroph abundance soil") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  nrow = 2, ncol = 1, align = "v")
dev.off()

#Wood saprotrophs
png("Figures/heatmap_genus_abund_asex_wood_sapro.png",
    width = 11.015 ,
    height = 17.03,
    units = "in", res = 300)
cowplot::plot_grid(
  spore_function_long %>% 
    filter(Wood_Saprotroph == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Wood saprotroph abundance air") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  
  spore_function_long %>% 
    filter(Wood_Saprotroph == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Wood saprotroph abundance soil") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  nrow = 2, ncol = 1, align = "v")
dev.off()

# Endophytes

png("Figures/heatmap_genus_abund_asex_endo.png",
    width = 11.015 ,
    height = 17.03,
    units = "in", res = 300)
cowplot::plot_grid(
  spore_function_long %>% 
    filter(Endophyte == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Endophyte abundance air") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  
  spore_function_long %>% 
    filter(Endophyte == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Endophyte abundance soil") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  nrow = 2, ncol = 1, align = "v")
dev.off()

## Lichen parasites
png("Figures/heatmap_genus_abund_asex_lichenP.png",
    width = 11.015 ,
    height = 17.03,
    units = "in", res = 300)
cowplot::plot_grid(
  spore_function_long %>% 
    filter(Lichen_Parasite == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Lichen parasite abundance air") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  
  spore_function_long %>% 
    filter(Lichen_Parasite == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Lichen parasite abundance soil") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  nrow = 2, ncol = 1, align = "v")
dev.off()

### Ectomycorrhizal

png("Figures/heatmap_genus_abund_asex_ecto.png",
    width = 11.015 ,
    height = 17.03,
    units = "in", res = 300)
cowplot::plot_grid(
  spore_function_long %>% 
    filter(Ectomycorrhizal == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Ectomycorrhizal abundance air") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  
  spore_function_long %>% 
    filter(Ectomycorrhizal == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Ectomycorrhizal abundance soil") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  nrow = 2, ncol = 1, align = "v")
dev.off()

### Human pathogens
png("Figures/heatmap_genus_abund_asex_human_path.png",
    width = 11.015 ,
    height = 17.03,
    units = "in", res = 300)
cowplot::plot_grid(
  spore_function_long %>% 
    filter(Human_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
    geom_tile(aes(fill =asexual_spore_volume_wt)) +
    
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Human_Pathogen abundance air") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  
  spore_function_long %>% 
    filter(Human_Pathogen == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance))+
    geom_tile(aes(fill =asexual_spore_volume_wt)) +
    
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Human_Pathogen abundance soil") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  nrow = 2, ncol = 1, align = "v")
dev.off()

### Plant pathogens
png("Figures/heatmap_genus_abund_asex_plant_path.png",
    width = 11.015 ,
    height = 17.03,
    units = "in", res = 300)
cowplot::plot_grid(
  spore_function_long %>% 
    filter(Plant_Pathogen == 1) %>% 
    filter(Type == "Air") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance)) +
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Plant_Pathogen abundance air") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  
  spore_function_long %>% 
    filter(Plant_Pathogen == 1) %>% 
    filter(Type == "Soil") %>% 
    filter(!is.na(asexual_spore_volume)) %>% 
    mutate(asexual_spore_volume_wt = log10(asexual_spore_volume)) %>% 
    mutate(sexual_spore_volume_wt = log10(sexual_spore_volume)) %>% 
    # group_by(Type, Location, Ecosystem, plot, sites, genus) %>%  
    # summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt", "abundance"),
    #              mean, na.rm = T) %>% 
    ggplot() +
    aes(x = paste(Location, sites), y = reorder(genus, abundance))+
        geom_tile(aes(fill =asexual_spore_volume_wt)) +
    #geom_point(aes(size = sexual_spore_volume_wt, color = sexual_spore_volume_wt)) +
    #geom_text(aes(label = genus)) +
    scale_size(range = c(0.5, 10)) +
    #facet_grid(. ~ Type, scales = "fixed") +
    facet_wrap(. ~ Ecosystem, scales = "free_x") +
    scale_fill_viridis_c(option = "turbo", direction = 1) +
    labs(title = "Plant_Pathogen abundance soil") +
    ylab(label = "Relative abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"),
  nrow = 2, ncol = 1, align = "v")
dev.off()


######################################################################################################
##################################      DATA ASSEMBLY    #############################################
######################################################################################################

rm(list = ls()); gc()
library(tidyverse)

#For the moment, this function appears twice in the code
organizando<-function(df){
  
  #Cleaning step
  m<-apply(df,1,paste,collapse="_")
  if(any(duplicated(m))){
    m<-which(duplicated(m))
    df<-df[-m,]}
  
  #Merging step
  s1<-split(df,df[,1])
  s2<-lapply(s1,function(df){
    data.frame(accepted_name=unique(df[,1]),
               guild=paste(df[,2],collapse = " "))
  })
  
  s3<-do.call("rbind",s2)
  
  #Eliminating duplicated words
  s3$guild<-
    #trial<-
    sapply(s3$guild,function(x){
      paste(
        unique(c(unlist(strsplit(x," ")))),
        collapse = " ")  
    })
  
  #Ordering alphabetically
  s3$guild<-unname(sapply(s3$guild, function(x) {
    paste(sort(strsplit(
      x[1], ' ')[[1]]), collapse=' ')} ))
  
  s3
} 

# ----------------------- Loading the Spore data
Spore_data<-readRDS("C:/Users/Carlos/Documents/Professional/SporeSizeAcrossFungKingdom/SporeSizeFungalKingdom/output/Spore_data_12Nov21.RDS")
#On Agugust 2023, I noticed some entries for conidia are too big and incorrect. In the meantime, the 
# easiest solution is to exclude all entries for conidia with very big width values which are equal to their lengths

m <- 
  which(
    (
      (Spore_data$SporeType == "Meiospores" & Spore_data$phylum == "Ascomycota") & 
        Spore_data$spore_width>99 & 
        (Spore_data$spore_width == Spore_data$spore_length)
    )|
      (
        (Spore_data$SporeType == "Mitospores" & Spore_data$phylum == "Ascomycota") & 
          Spore_data$spore_width>27 & 
          (Spore_data$spore_width == Spore_data$spore_length)
      ) 
  )


Spore_data <- Spore_data[-m, ]

#saveRDS(Spore_data, "Spore_data.RDS")
#rm(BirdEgg,FunGuildData,df_species_byPrimerSet,kewData,FungalTaxanomy_col,l)


a <- unique(paste(Spore_data$order,Spore_data$family, sep = "_"))
b <- unique(Spore_data$family)

a <- str_split(a,"_")
a <-sapply(a,function(x){x[2]})

a[duplicated(a)]

#Using the taxonomy that the (global) data from Otso and Nerea have
Spore_data$order[Spore_data$family=="Mycosphaerellaceae"] <- "Capnodiales"
Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
Spore_data$order[Spore_data$family=="Elaphomycetaceae"] <- "Eurotiales"
Spore_data$order[Spore_data$family=="Physciaceae"] <- "Teloschistales"
Spore_data$order[Spore_data$family=="Trapeliaceae"] <- "Baeomycetales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"

Spore_data$SporeType <- gsub(" ","_",Spore_data$SporeType)
spore_genus <- organizando(Spore_data[,c("genus","SporeType")])
spore_genus$guild <- gsub("\\) c\\(","",spore_genus$guild)
spore_genus$guild <- gsub("\\)","",spore_genus$guild)
spore_genus$guild <- gsub('"'," ",spore_genus$guild)
spore_genus$guild <- gsub("\\,","",spore_genus$guild)
spore_genus$guild <- trimws(spore_genus$guild)
spore_genus$guild <- gsub("c\\(","",spore_genus$guild)
spore_genus$guild <- gsub("\\s+"," ",spore_genus$guild)
names(spore_genus)[2] <- "SporeType"

Genus_Spore_data<-
  Spore_data %>% 
  group_by(phylum,
           class,order,family,genus,SporeType) %>% 
  summarize_at(c("SporeVolume","spore_length","spore_width",
                 "SporeArea","Q_ratio"), median)

#----------------------------------------- Load community data (Finland data) ---------------------------
source("Finland_data/S1_S_and_ABU_models_CA.R")
#detach("package:lme4", unload = TRUE)
detach("package:MASS", unload = TRUE)

#Original PROTU table with associated taxonomy
protuID_tax_org <- readRDS("finland_data.RDS")
protuID_tax_org$canonicalname <- trimws(protuID_tax_org$canonicalname)

#Corrected taxonomy according to Index Fungorum (see Taxonomy_Finland.R)
corr_tax <- readRDS("corr_tax_comb.RDS"); m <- apply(corr_tax,1,paste,collapse="_"); m <- which(duplicated(m)); corr_tax <- corr_tax[-m,]


#PROTU table + corrected taxonomy
protuID_tax <- left_join(protuID_tax_org[,c("PROTU","HR","canonicalname")],
                         corr_tax)#removing status and matchtype

rm(protuID_tax_org,corr_tax, tr)

#Site x OTU matrix (174 sites x 79155 PROTUs)
site_protu_matrix <- protu_df# the matrix object "protu" is the original one and contains the same info


#OTU x site matrix (the tranpose of the site + otu matrix) (79155 PROTU´s x 174 sites)
protu_df <- t(protu_df)
protu_df <- as.data.frame(protu_df)
protu_df$species <- rownames(protu_df)
rownames(protu_df) <- NULL
m <- length(names(protu_df))
protu_df <- protu_df[, c(m, 1:(m - 1))]#this is a protu x sites matrix
protu_site_matrix <- protu_df

#Long version of site x species matrix with only three columns: species, sites and abundance
protu_df1 <- gather(protu_df, key = sites,
                    value = abundance, names(protu_df)[2]:last(names(protu_df))); names(protu_df1)[1] <- "PROTU"

#adding site data (whether it comes from soil, air or natural or urban)
protu_df1 <- left_join(protu_df1,samples %>% 
                         rename(sites = ID))

#protu_df1 preserves the same info as the original protu, it was only transposed and new data
#has been added. Importantly, the total abundance of species for each site remains one as it is 
#happening in the original protu. Examples:

sum(protu[which(rownames(protu) == "TAM-U3-200819-A"), ])#This is 1
sum(protu_df1$abundance[protu_df1$sites == "TAM-U3-200819-A"])#This is also 1

sum(protu[which(rownames(protu) == "HEL-N1-170819-S"), ])#This is 1
sum(protu_df1$abundance[protu_df1$sites =="HEL-N1-170819-S"])#This is also 1
#

protu_site_matrix_long <- protu_df1
rm(protu_df, protu_df1, otu, otu.dna); gc()

#------------------------ Loading the guild data 
source("C:/Users/Carlos/Documents/Fungal_LifeStyle_Database/Loading_GuildData.R")#;rm(m,rev,rev2)
# saveRDS(Guild_Data, "Guild_Data.RDS")
# saveRDS(taxonomy_guilds, "taxonomy_guilds.RDS")
#taxonomy_guilds <- readRDS("taxonomy_guilds.RDS")
rm(Guild_Data_o)
#Guild_Data$trophicMode<-gsub("-"," ",Guild_Data$trophicMode)

#How many have assigned guild at species level
#Guild_Data$guild[which(Guild_Data$guild == "")] <- NA
#Guild_Data$accepted_name[which(Guild_Data$accepted_name == "NA")] <- NA
Guild_Data$accepted_name <- gsub("_", " ", Guild_Data$accepted_name)
Guild_Data$names_to_use <- gsub("_", " ", Guild_Data$names_to_use)

# I need to better define "Soft rots"
Guild_Data$guild[
  grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)] <- 
  gsub("Wood_Saprotroph", "Soft_Saprotroph",
       Guild_Data$guild[
         grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)])

#####  On april 17 I sent this version to Syklar and Brendan (Finland group)####
# source("C:/Users/Carlos/Documents/Fungal_LifeStyle_Database/Loading_GuildData.R")#;rm(m,rev,rev2)
# rm(Guild_Data_o)
# Guild_Data$accepted_name <- gsub("_", " ", Guild_Data$accepted_name)
# Guild_Data$names_to_use <- gsub("_", " ", Guild_Data$names_to_use)
# Guild_Data$guild[
#   grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)] <- 
#   gsub("Wood_Saprotroph", "Soft_wood_Saprotroph",
#        Guild_Data$guild[
#          grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)])
# write.csv(Guild_Data, "Guild_Data_17Ap23.csv", row.names = F)
#####

organizando<-function(df){
  
  #Cleaning step
  m<-apply(df,1,paste,collapse="_")
  if(any(duplicated(m))){
    m<-which(duplicated(m))
    df<-df[-m,]}
  
  #Merging step
  s1<-split(df,df[,1])
  s2<-lapply(s1,function(df){
    data.frame(accepted_name=unique(df[,1]),
               guild=paste(df[,2],collapse = " "))
  })
  
  s3<-do.call("rbind",s2)
  
  #Eliminating duplicated words
  s3$guild<-
    #trial<-
    sapply(s3$guild,function(x){
      paste(
        unique(c(unlist(strsplit(x," ")))),
        collapse = " ")  
    })
  
  #Ordering alphabetically
  s3$guild<-unname(sapply(s3$guild, function(x) {
    paste(sort(strsplit(
      x[1], ' ')[[1]]), collapse=' ')} ))
  
  s3
} 


df_tax <- Guild_Data[,c("names_to_use","guild","accepted_name")]
m <- which(duplicated(taxonomy_guilds$accepted_name))#For some reason, even the same Species Fungorum
df_tax <- left_join(df_tax,taxonomy_guilds[-m,c("kingdom","phylum","class",
                                                "order","family","genus",
                                                "accepted_name")] )
df_tax <- as.data.frame(df_tax)

guilds_genus <- organizando(df_tax[,c("genus","guild")])
guilds_genus$guild_count<-contando(guilds_genus)
guilds_genus$tax_level<-"B_genus"

m <- apply(taxonomy_guilds[,c("phylum","genus")],1,paste,collapse="_")
m <- which(duplicated(m))

guilds_genus <-
  left_join(guilds_genus,taxonomy_guilds[-m,c("phylum","genus")] %>% rename(accepted_name=genus))

#####

# On January 22 Nerea checked the assignation of guilds to the genera found in the 
# Finnish data (only). Check the script Manual check Nerea on the story behind this
# check
new_check_na <- read.csv2("new_guild_check_NA.csv")
new_check_na$guild_revised_NA[which(new_check_na$guild_revised_NA=="")] <- 
  new_check_na$guild[which(new_check_na$guild_revised_NA=="")]
new_check_na$guild_revised_NA <- gsub("Wood_sa", "Wood_Sa", new_check_na$guild_revised_NA)

guild_genus<-new_check_na[, c(1, 3)]; rm(new_check_na)
names(guild_genus)[2] <- "guild"

# However, Nerea still included some of the genera considered soft-rots as saprotrophs
# in this version I am considered soft-rots as "saprotroph" and not "wood-saprotrophs"

guild_genus$guild[
  which(guild_genus$accepted_name%in%guilds_genus$accepted_name[grep("Soft",guilds_genus$guild)])] <- 
  gsub("Wood_Saprotroph", "Saprotroph",
       guild_genus$guild[
         which(guild_genus$accepted_name%in%guilds_genus$accepted_name[grep("Soft",guilds_genus$guild)])])

# Replacing the ones coming from the Nerea check
guild_genus <- organizando(guild_genus)

# Replacing the ones coming from the Guild_database

guilds_genus$guild[grep("Soft",guilds_genus$guild)] <- 
  gsub("Soft_Saprotroph", "Saprotroph", guilds_genus$guild[grep("Soft",guilds_genus$guild)])
guilds_genus <- organizando(guilds_genus)

# Updating wood saprotrophs entries ("soft" saprotrophs should not be included as wood saprotrophs)

names(guilds_genus)[1] <- "genus"

more <-
  unique(protuID_tax$genus[
    which(!protuID_tax$genus%in%guild_genus$accepted_name)])

more2 <- 
  guilds_genus[which(guilds_genus$genus%in%more), c("genus", "guild")]

more2 <- unique(more2)
more2$guild[which(more2$genus=="Cylindrocarpon")] <- "Plant_Pathogen Saprotroph"
more2$guild[which(more2$genus=="Furcasterigmium")] <- "Saprotroph"
more2$guild[which(more2$genus=="Paraphoma")] <- "Plant_Pathogen Saprotroph"
more2$guild[which(more2$genus=="Wardomyces")] <- "Saprotroph"
names(more2)[1] <- "accepted_name"

#Cylindrocarpon not a wood_saprotroph
#Furcasterigmium not a wood_saprotroph
#Paraphoma not a wood_saprotroph
#Wardomyces not a wood_saprotroph

guild_genus <- bind_rows(guild_genus, more2)
rm(guilds_genus);gc()

#Transforming the guild table into wide format (each column is a function)
library(qdapTools)

dat <- guild_genus
m <- length(names(dat))
dat <- cbind(dat, str_split_fixed(dat$guild, " ", n = Inf))
n <- length(names(dat))
names(dat)[(m + 1):n] <- paste(rep("guild", (n - m)), 1:(n - m), sep = "_")

dat <- as.matrix(dat)
dat[which(dat == "")] <- NA
dat <- as.data.frame(dat)

probando <- dat
probando <- pivot_longer(cols = c(3:length(names(probando))), probando)#guild_1 to guild_7
probando$name <- NULL
probando <- probando[-which(is.na(probando$value)), ]
probando$guild <- NULL
#probando<-probando[-which(duplicated(paste(probando$accepted_name,probando$value))),]

probando <- pivot_wider(probando,names_from = value)
probando <- as.data.frame(probando)
#probando <- as.matrix(probando)

probando[, -1][!is.na(probando[, -1])] <- 1
probando[is.na(probando)] <- 0
functional_dat <- probando; rm(probando)
names(functional_dat)[1] <- "genus"

genus_funct_matrix <- functional_dat; rm(functional_dat, dat, more, more2)

#1 Add the Protu ID to the genus_funct_matrix
protuID_funct <- left_join(protuID_tax[,c("PROTU","genus")],
                           genus_funct_matrix, na_matches = "never")


#2 Add now the functions
protuID_site_matrix_funct <- left_join(protu_site_matrix %>% rename(PROTU = species),
                                       protuID_funct, na_matches = "never")

sum(protuID_site_matrix_funct$`HEL-N1-170819-A`)#all these ones sum 1

#------------------------------------ Merging community data with Functional Data -----------------------------------------------------------------------------

#Second method: this one is based on presence absence only. That is whether a Protu
#is present or not. This one makes a bit more sense (although info is lost about
#the abundance of the PROTU) as I can recalculate the abundance of each functional group
#per site, then use this info in a stalked bar

#For this I need a table with: 
#1, The info the relative frequency of each functional group (y axis)
#2. the names of the sites (x axis)
#2. Column with the name of the functional group (fill)
#3. A column with the type soil or air (x axis)
#4. A column with the location (natural or urban) (facet)

step1 <-
  left_join(protu_site_matrix_long, protuID_tax[,c("PROTU","genus")],
            na_matches = "never") %>% 
  filter(abundance != 0) %>% 
  select(PROTU, sites, genus, Type, Ecosystem)

m <- length(names(protuID_site_matrix_funct))

step2 <-
  left_join(step1,
            protuID_site_matrix_funct[-which(duplicated(protuID_site_matrix_funct$genus)),
                                      c(176:m)])#the genus and the functional groups

m <- length(names(step2))
step3 <- pivot_longer(cols = c(6:m), step2)#These are all the guilds
step3 <- step3[ - which(step3$value == 0), ]
step3$name[which(is.na(step3$value))] <- NA
names(step3)[6] <- "functional_group"
site_funct_group_count <- step3; rm(step1, step2, step3)

#-------------------------------------------------- Merging with Spore Data -----------------------------------------------------------------------------

ProtuID_spores<-
  left_join(
    protuID_tax[,c("PROTU","genus")],
    Genus_Spore_data %>% 
      ungroup() %>% 
      select(genus,SporeType,SporeVolume) %>%
      filter(SporeType=="Mitospores"))

names(ProtuID_spores)[4]<-"asexual_spore_volume"
ProtuID_spores$SporeType<-NULL

ProtuID_spores<-
  left_join(
    ProtuID_spores,
    Genus_Spore_data %>% 
      ungroup() %>% 
      select(genus,SporeType,SporeVolume) %>%
      filter(SporeType=="Meiospores") %>% 
      select(genus,SporeVolume))

names(ProtuID_spores)[4]<-"sexual_spore_volume"
ProtuID_spores$SporeType<-NULL
#rm(Genus_Spore_data)


#Adding spore data to the protu x site matrix (long version)
protu_site_matrix_long<-left_join(protu_site_matrix_long,ProtuID_spores[,c("PROTU","asexual_spore_volume")])
protu_site_matrix_long<-left_join(protu_site_matrix_long,ProtuID_spores[,c("PROTU","sexual_spore_volume")])
protu_site_matrix_long<-left_join(protu_site_matrix_long,ProtuID_spores[,c("PROTU","genus")])

#This is 1
sum(protu_site_matrix_long$abundance[protu_site_matrix_long$site=="HEL-N1-170819-A"])

comm_wt_means_spore <- protu_site_matrix_long %>% 
  filter(abundance != 0) %>%
  mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, RefinedType, plot, Ecosystem, sites) %>% 
  summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T)

comm_wt_means_spore$plot <- gsub("HEL-|JOE-|JYV-|LAH-|TAM-", "", comm_wt_means_spore$plot)
comm_wt_means_spore$plot <- gsub("HEL-|JOE-|JYV-|LAH-|TAM-", "", comm_wt_means_spore$plot)
comm_wt_means_spore$plot <- gsub("U|N", "P", comm_wt_means_spore$plot)

#comm_wt_means_spore$RefinedType <- gsub("Urban|Natural", "", comm_wt_means_spore$RefinedType)
comm_wt_means_spore$RefinedType[which(comm_wt_means_spore$RefinedType == "UrbanEdge")] <- "UrbanaEdge"

comm_wt_means_spore$Location[-grep("Hel|Joen|Lah|Tamp", comm_wt_means_spore$Location)] <- "Jyvaskyla"

rm(list = ls()); gc()
library(tidyverse)

#For the moment, this function appears twice in the code
organizando<-function(df){
  
  #Cleaning step
  m<-apply(df,1,paste,collapse="_")
  if(any(duplicated(m))){
    m<-which(duplicated(m))
    df<-df[-m,]}
  
  #Merging step
  s1<-split(df,df[,1])
  s2<-lapply(s1,function(df){
    data.frame(accepted_name=unique(df[,1]),
               guild=paste(df[,2],collapse = " "))
  })
  
  s3<-do.call("rbind",s2)
  
  #Eliminating duplicated words
  s3$guild<-
    #trial<-
    sapply(s3$guild,function(x){
      paste(
        unique(c(unlist(strsplit(x," ")))),
        collapse = " ")  
    })
  
  #Ordering alphabetically
  s3$guild<-unname(sapply(s3$guild, function(x) {
    paste(sort(strsplit(
      x[1], ' ')[[1]]), collapse=' ')} ))
  
  s3
} 

# ----------------------- Loading the Spore data
Spore_data<-readRDS("C:/Users/Carlos/Documents/Professional/SporeSizeAcrossFungKingdom/SporeSizeFungalKingdom/output/Spore_data_12Nov21.RDS")
#On Agugust 2023, I noticed some entries for conidia are too big and incorrect. In the meantime, the 
# easiest solution is to exclude all entries for conidia with very big width values which are equal to their lengths

m <- 
which(
  (
  (Spore_data$SporeType == "Meiospores" & Spore_data$phylum == "Ascomycota") & 
    Spore_data$spore_width>99 & 
    (Spore_data$spore_width == Spore_data$spore_length)
  )|
    (
      (Spore_data$SporeType == "Mitospores" & Spore_data$phylum == "Ascomycota") & 
        Spore_data$spore_width>27 & 
        (Spore_data$spore_width == Spore_data$spore_length)
    ) 
  )
  

Spore_data <- Spore_data[-m, ]

#saveRDS(Spore_data, "Spore_data.RDS")
#rm(BirdEgg,FunGuildData,df_species_byPrimerSet,kewData,FungalTaxanomy_col,l)


a <- unique(paste(Spore_data$order,Spore_data$family, sep = "_"))
b <- unique(Spore_data$family)

a <- str_split(a,"_")
a <-sapply(a,function(x){x[2]})

a[duplicated(a)]

#Using the taxonomy that the (global) data from Otso and Nerea have
Spore_data$order[Spore_data$family=="Mycosphaerellaceae"] <- "Capnodiales"
Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
Spore_data$order[Spore_data$family=="Elaphomycetaceae"] <- "Eurotiales"
Spore_data$order[Spore_data$family=="Physciaceae"] <- "Teloschistales"
Spore_data$order[Spore_data$family=="Trapeliaceae"] <- "Baeomycetales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"
# Spore_data$order[Spore_data$family=="Venturiaceae"] <- "Venturiales"

Spore_data$SporeType <- gsub(" ","_",Spore_data$SporeType)
spore_genus <- organizando(Spore_data[,c("genus","SporeType")])
spore_genus$guild <- gsub("\\) c\\(","",spore_genus$guild)
spore_genus$guild <- gsub("\\)","",spore_genus$guild)
spore_genus$guild <- gsub('"'," ",spore_genus$guild)
spore_genus$guild <- gsub("\\,","",spore_genus$guild)
spore_genus$guild <- trimws(spore_genus$guild)
spore_genus$guild <- gsub("c\\(","",spore_genus$guild)
spore_genus$guild <- gsub("\\s+"," ",spore_genus$guild)
names(spore_genus)[2] <- "SporeType"

Genus_Spore_data<-
  Spore_data %>% 
  group_by(phylum,
           class,order,family,genus,SporeType) %>% 
  summarize_at(c("SporeVolume","spore_length","spore_width",
                 "SporeArea","Q_ratio"), median)

#----------------------------------------- Load community data (Finland data) ---------------------------
source("Finland_data/S1_S_and_ABU_models_CA.R")
#detach("package:lme4", unload = TRUE)
detach("package:MASS", unload = TRUE)

#Original PROTU table with associated taxonomy
protuID_tax_org <- readRDS("finland_data.RDS")
protuID_tax_org$canonicalname <- trimws(protuID_tax_org$canonicalname)

#Corrected taxonomy according to Index Fungorum (see Taxonomy_Finland.R)
corr_tax <- readRDS("corr_tax_comb.RDS"); m <- apply(corr_tax,1,paste,collapse="_"); m <- which(duplicated(m)); corr_tax <- corr_tax[-m,]


#PROTU table + corrected taxonomy
protuID_tax <- left_join(protuID_tax_org[,c("PROTU","HR","canonicalname")],
                         corr_tax)#removing status and matchtype

rm(protuID_tax_org,corr_tax, tr)

#Site x OTU matrix (174 sites x 79155 PROTUs)
site_protu_matrix <- protu_df# the matrix object "protu" is the original one and contains the same info


#OTU x site matrix (the tranpose of the site + otu matrix) (79155 PROTU´s x 174 sites)
protu_df <- t(protu_df)
protu_df <- as.data.frame(protu_df)
protu_df$species <- rownames(protu_df)
rownames(protu_df) <- NULL
m <- length(names(protu_df))
protu_df <- protu_df[, c(m, 1:(m - 1))]#this is a protu x sites matrix
protu_site_matrix <- protu_df

#Long version of site x species matrix with only three columns: species, sites and abundance
protu_df1 <- gather(protu_df, key = sites,
                    value = abundance, names(protu_df)[2]:last(names(protu_df))); names(protu_df1)[1] <- "PROTU"

#adding site data (whether it comes from soil, air or natural or urban)
protu_df1 <- left_join(protu_df1,samples %>% 
                         rename(sites = ID))

#protu_df1 preserves the same info as the original protu, it was only transposed and new data
#has been added. Importantly, the total abundance of species for each site remains one as it is 
#happening in the original protu. Examples:

sum(protu[which(rownames(protu) == "TAM-U3-200819-A"), ])#This is 1
sum(protu_df1$abundance[protu_df1$sites == "TAM-U3-200819-A"])#This is also 1

sum(protu[which(rownames(protu) == "HEL-N1-170819-S"), ])#This is 1
sum(protu_df1$abundance[protu_df1$sites =="HEL-N1-170819-S"])#This is also 1
#

protu_site_matrix_long <- protu_df1
rm(protu_df, protu_df1, otu, otu.dna); gc()

#------------------------ Loading the guild data 
source("C:/Users/Carlos/Documents/Fungal_LifeStyle_Database/Loading_GuildData.R")#;rm(m,rev,rev2)
# saveRDS(Guild_Data, "Guild_Data.RDS")
# saveRDS(taxonomy_guilds, "taxonomy_guilds.RDS")
#taxonomy_guilds <- readRDS("taxonomy_guilds.RDS")
rm(Guild_Data_o)
#Guild_Data$trophicMode<-gsub("-"," ",Guild_Data$trophicMode)

#How many have assigned guild at species level
#Guild_Data$guild[which(Guild_Data$guild == "")] <- NA
#Guild_Data$accepted_name[which(Guild_Data$accepted_name == "NA")] <- NA
Guild_Data$accepted_name <- gsub("_", " ", Guild_Data$accepted_name)
Guild_Data$names_to_use <- gsub("_", " ", Guild_Data$names_to_use)

# I need to better define "Soft rots"
Guild_Data$guild[
grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)] <- 
gsub("Wood_Saprotroph", "Soft_Saprotroph",
  Guild_Data$guild[
    grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)])

#####  On april 17 I sent this version to Syklar and Brendan (Finland group)####
# source("C:/Users/Carlos/Documents/Fungal_LifeStyle_Database/Loading_GuildData.R")#;rm(m,rev,rev2)
# rm(Guild_Data_o)
# Guild_Data$accepted_name <- gsub("_", " ", Guild_Data$accepted_name)
# Guild_Data$names_to_use <- gsub("_", " ", Guild_Data$names_to_use)
# Guild_Data$guild[
#   grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)] <- 
#   gsub("Wood_Saprotroph", "Soft_wood_Saprotroph",
#        Guild_Data$guild[
#          grep("List_of_Fungi_in_Soft-Rot_Tissue", Guild_Data$citationSource)])
# write.csv(Guild_Data, "Guild_Data_17Ap23.csv", row.names = F)
#####

organizando<-function(df){
  
  #Cleaning step
  m<-apply(df,1,paste,collapse="_")
  if(any(duplicated(m))){
    m<-which(duplicated(m))
    df<-df[-m,]}
  
  #Merging step
  s1<-split(df,df[,1])
  s2<-lapply(s1,function(df){
    data.frame(accepted_name=unique(df[,1]),
               guild=paste(df[,2],collapse = " "))
  })
  
  s3<-do.call("rbind",s2)
  
  #Eliminating duplicated words
  s3$guild<-
    #trial<-
    sapply(s3$guild,function(x){
      paste(
        unique(c(unlist(strsplit(x," ")))),
        collapse = " ")  
    })
  
  #Ordering alphabetically
  s3$guild<-unname(sapply(s3$guild, function(x) {
    paste(sort(strsplit(
      x[1], ' ')[[1]]), collapse=' ')} ))
  
  s3
} 


df_tax <- Guild_Data[,c("names_to_use","guild","accepted_name")]
m <- which(duplicated(taxonomy_guilds$accepted_name))#For some reason, even the same Species Fungorum
df_tax <- left_join(df_tax,taxonomy_guilds[-m,c("kingdom","phylum","class",
                                              "order","family","genus",
                                              "accepted_name")] )
df_tax <- as.data.frame(df_tax)

guilds_genus <- organizando(df_tax[,c("genus","guild")])
guilds_genus$guild_count<-contando(guilds_genus)
guilds_genus$tax_level<-"B_genus"

m <- apply(taxonomy_guilds[,c("phylum","genus")],1,paste,collapse="_")
m <- which(duplicated(m))

guilds_genus <-
  left_join(guilds_genus,taxonomy_guilds[-m,c("phylum","genus")] %>% rename(accepted_name=genus))

#####

# On January 22 Nerea checked the assignation of guilds to the genera found in the 
# Finnish data (only). Check the script Manual check Nerea on the story behind this
# check
new_check_na <- read.csv2("new_guild_check_NA.csv")
new_check_na$guild_revised_NA[which(new_check_na$guild_revised_NA=="")] <- 
  new_check_na$guild[which(new_check_na$guild_revised_NA=="")]
new_check_na$guild_revised_NA <- gsub("Wood_sa", "Wood_Sa", new_check_na$guild_revised_NA)

guild_genus<-new_check_na[, c(1, 3)]; rm(new_check_na)
names(guild_genus)[2] <- "guild"

# However, Nerea still included some of the genera considered soft-rots as saprotrophs
# in this version I am considered soft-rots as "saprotroph" and not "wood-saprotrophs"

guild_genus$guild[
  which(guild_genus$accepted_name%in%guilds_genus$accepted_name[grep("Soft",guilds_genus$guild)])] <- 
gsub("Wood_Saprotroph", "Saprotroph",
guild_genus$guild[
  which(guild_genus$accepted_name%in%guilds_genus$accepted_name[grep("Soft",guilds_genus$guild)])])

# Replacing the ones coming from the Nerea check
guild_genus <- organizando(guild_genus)

# Replacing the ones coming from the Guild_database

guilds_genus$guild[grep("Soft",guilds_genus$guild)] <- 
gsub("Soft_Saprotroph", "Saprotroph", guilds_genus$guild[grep("Soft",guilds_genus$guild)])
guilds_genus <- organizando(guilds_genus)

# Updating wood saprotrophs entries ("soft" saprotrophs should not be included as wood saprotrophs)

names(guilds_genus)[1] <- "genus"

more <-
  unique(protuID_tax$genus[
    which(!protuID_tax$genus%in%guild_genus$accepted_name)])

more2 <- 
  guilds_genus[which(guilds_genus$genus%in%more), c("genus", "guild")]

more2 <- unique(more2)
more2$guild[which(more2$genus=="Cylindrocarpon")] <- "Plant_Pathogen Saprotroph"
more2$guild[which(more2$genus=="Furcasterigmium")] <- "Saprotroph"
more2$guild[which(more2$genus=="Paraphoma")] <- "Plant_Pathogen Saprotroph"
more2$guild[which(more2$genus=="Wardomyces")] <- "Saprotroph"
names(more2)[1] <- "accepted_name"

#Cylindrocarpon not a wood_saprotroph
#Furcasterigmium not a wood_saprotroph
#Paraphoma not a wood_saprotroph
#Wardomyces not a wood_saprotroph

guild_genus <- bind_rows(guild_genus, more2)
rm(guilds_genus);gc()

#Transforming the guild table into wide format (each column is a function)
library(qdapTools)

dat <- guild_genus
m <- length(names(dat))
dat <- cbind(dat, str_split_fixed(dat$guild, " ", n = Inf))
n <- length(names(dat))
names(dat)[(m + 1):n] <- paste(rep("guild", (n - m)), 1:(n - m), sep = "_")

dat <- as.matrix(dat)
dat[which(dat == "")] <- NA
dat <- as.data.frame(dat)

probando <- dat
probando <- pivot_longer(cols = c(3:length(names(probando))), probando)#guild_1 to guild_7
probando$name <- NULL
probando <- probando[-which(is.na(probando$value)), ]
probando$guild <- NULL
#probando<-probando[-which(duplicated(paste(probando$accepted_name,probando$value))),]

probando <- pivot_wider(probando,names_from = value)
probando <- as.data.frame(probando)
#probando <- as.matrix(probando)

probando[, -1][!is.na(probando[, -1])] <- 1
probando[is.na(probando)] <- 0
functional_dat <- probando; rm(probando)
names(functional_dat)[1] <- "genus"

genus_funct_matrix <- functional_dat; rm(functional_dat, dat, more, more2)

#1 Add the Protu ID to the genus_funct_matrix
protuID_funct <- left_join(protuID_tax[,c("PROTU","genus")],
                           genus_funct_matrix, na_matches = "never")


#2 Add now the functions
protuID_site_matrix_funct <- left_join(protu_site_matrix %>% rename(PROTU = species),
                                       protuID_funct, na_matches = "never")

sum(protuID_site_matrix_funct$`HEL-N1-170819-A`)#all these ones sum 1

#------------------------------------ Merging community data with Functional Data -----------------------------------------------------------------------------

#Second method: this one is based on presence absence only. That is whether a Protu
#is present or not. This one makes a bit more sense (although info is lost about
#the abundance of the PROTU) as I can recalculate the abundance of each functional group
#per site, then use this info in a stalked bar

#For this I need a table with: 
#1, The info the relative frequency of each functional group (y axis)
#2. the names of the sites (x axis)
#2. Column with the name of the functional group (fill)
#3. A column with the type soil or air (x axis)
#4. A column with the location (natural or urban) (facet)

step1 <-
  left_join(protu_site_matrix_long, protuID_tax[,c("PROTU","genus")],
            na_matches = "never") %>% 
  filter(abundance != 0) %>% 
  select(PROTU, sites, genus, Type, Ecosystem)

m <- length(names(protuID_site_matrix_funct))

step2 <-
  left_join(step1,
            protuID_site_matrix_funct[-which(duplicated(protuID_site_matrix_funct$genus)),
                                      c(176:m)])#the genus and the functional groups

m <- length(names(step2))
step3 <- pivot_longer(cols = c(6:m), step2)#These are all the guilds
step3 <- step3[ - which(step3$value == 0), ]
step3$name[which(is.na(step3$value))] <- NA
names(step3)[6] <- "functional_group"
site_funct_group_count <- step3; rm(step1, step2, step3)

#-------------------------------------------------- Merging with Spore Data -----------------------------------------------------------------------------

ProtuID_spores<-
  left_join(
    protuID_tax[,c("PROTU","genus")],
    Genus_Spore_data %>% 
      ungroup() %>% 
      select(genus,SporeType,SporeVolume) %>%
      filter(SporeType=="Mitospores"))

names(ProtuID_spores)[4]<-"asexual_spore_volume"
ProtuID_spores$SporeType<-NULL

ProtuID_spores<-
  left_join(
    ProtuID_spores,
    Genus_Spore_data %>% 
      ungroup() %>% 
      select(genus,SporeType,SporeVolume) %>%
      filter(SporeType=="Meiospores") %>% 
      select(genus,SporeVolume))

names(ProtuID_spores)[4]<-"sexual_spore_volume"
ProtuID_spores$SporeType<-NULL
#rm(Genus_Spore_data)


#Adding spore data to the protu x site matrix (long version)
protu_site_matrix_long<-left_join(protu_site_matrix_long,ProtuID_spores[,c("PROTU","asexual_spore_volume")])
protu_site_matrix_long<-left_join(protu_site_matrix_long,ProtuID_spores[,c("PROTU","sexual_spore_volume")])
protu_site_matrix_long<-left_join(protu_site_matrix_long,ProtuID_spores[,c("PROTU","genus")])

#This is 1
sum(protu_site_matrix_long$abundance[protu_site_matrix_long$site=="HEL-N1-170819-A"])

comm_wt_means_spore <- protu_site_matrix_long %>% 
  filter(abundance != 0) %>%
  mutate(asexual_spore_volume_wt = abundance*log10(asexual_spore_volume)) %>% 
  mutate(sexual_spore_volume_wt = abundance*log10(sexual_spore_volume)) %>% 
  group_by(Type, Location, RefinedType, plot, Ecosystem, sites) %>% 
  summarise_at(c("asexual_spore_volume_wt", "sexual_spore_volume_wt"), sum, na.rm = T)

comm_wt_means_spore$plot <- gsub("HEL-|JOE-|JYV-|LAH-|TAM-", "", comm_wt_means_spore$plot)
comm_wt_means_spore$plot <- gsub("HEL-|JOE-|JYV-|LAH-|TAM-", "", comm_wt_means_spore$plot)
comm_wt_means_spore$plot <- gsub("U|N", "P", comm_wt_means_spore$plot)

#comm_wt_means_spore$RefinedType <- gsub("Urban|Natural", "", comm_wt_means_spore$RefinedType)
comm_wt_means_spore$RefinedType[which(comm_wt_means_spore$RefinedType == "UrbanEdge")] <- "UrbanaEdge"

comm_wt_means_spore$Location[-grep("Hel|Joen|Lah|Tamp", comm_wt_means_spore$Location)] <- "Jyvaskyla"

#####

# Loading
library(lme4)
library(lmerTest)
library(performance)
library(metafor)

spore_function <- protuID_site_matrix_funct[, c(1:176)] %>% 
  group_by(genus) %>% 
  summarise_at((vars(`HEL-N1-170819-A`:`TAM-U3-210819-S`)), sum)

sum(spore_function$`HEL-N1-170819-A`) # This is still 1

spore_function_long <- gather(spore_function, key = sites,
                              value = abundance, names(spore_function)[2]:last(names(spore_function)))#; names(spore_function_long)[1] <- "genus"

#adding site data (whether it comes from soil, air or natural or urban)
spore_function_long <- left_join(spore_function_long, 
                                 samples %>% 
                                   rename(sites = ID))

spore_function_long <- 
  spore_function_long[which(spore_function_long$abundance > 0), ]

sum(spore_function_long$abundance[spore_function_long$sites == "HEL-N1-170819-A"])

spore_function_long$frequency <- spore_function_long$abundance
spore_function_long$frequency[which(spore_function_long$frequency > 0)] <- 1
spore_function_long <- spore_function_long[-which(is.na(spore_function_long$genus)), ] # here I lose some fungi
spore_function_long$Location[-grep("Hel|Joe|Lah|Tamp", spore_function_long$Location)] <- "Jyvaskyla"

# Adding spore data
spore_function_long <-
  left_join(
    spore_function_long,
    Genus_Spore_data %>% 
      ungroup() %>% 
      select(genus,SporeType,SporeVolume) %>%
      filter(SporeType=="Meiospores") %>% 
      select(genus,SporeVolume)); names(spore_function_long)[length(names(spore_function_long))] <- "sexual_spore_volume"

spore_function_long <-
  left_join(
    spore_function_long,
    Genus_Spore_data %>% 
      ungroup() %>% 
      select(genus,SporeType,SporeVolume) %>%
      filter(SporeType=="Mitospores") %>% 
      select(genus,SporeVolume)); names(spore_function_long)[length(names(spore_function_long))] <- "asexual_spore_volume"

# Adding the functional data (it contains the functions in wide format, useful for spore functions analysis)
spore_function_long_temp <-
  left_join(spore_function_long, genus_funct_matrix)

# Having the functions in long format (useful for analyisis of turnover of functional groups)
genus_site_funct_long <-
  left_join(spore_function_long,
            distinct(site_funct_group_count[,c("sites","genus","functional_group")]))

# Renaming just because downstream analysis requires it

spore_function_long <- spore_function_long_temp
rm(spore_function_long_temp); gc()



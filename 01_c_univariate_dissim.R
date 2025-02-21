

##***********************************##
##* Maria Angeles Pérez-Navarro
##* King's College University
##* AlienImpacts project
##* March 2022
##* *********************************##

library(dplyr)
library(readr)
library(tidyr)
library("methods")
# library(textreadr)#
library(stringr)
#library(radiant.data)
library(tidyverse)
library(taxize)
#library(Taxonstand)
library(ade4)


###******************************************###
###******1. PLOT LEVEL DISSIMILARITIES******####
###******************************************###

#** check complete cases for main variables ####

abun_traits <- read_delim("../../output/tables/abundances/ab_traits_complete_df1.csv", 
                          delim=";", col_names=T)

abun_traits <- abun_traits%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D",
    Transect_old=="1"~"A",
    Transect_old=="2"~"B",
    Transect_old=="3"~"C",
    Transect_old=="4"~"D",
    Transect_old=="5"~"E",
    Transect_old=="6"~"F",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ),
  plot_code=paste0(Year,"_",
                   OldField, "_",
                   Transect, "_",
                   Plot))

misce <- abun_traits%>%
  filter(genus%in%c("Miscellaneous","miscellaneous",
                    "Unknown", "unknown", "woody",
                    "Woody", "Mosses", "mosses",
                    "animal", "ant", "anthill",
                    "anthills", "Bare", "dung",
                    "fern", "Forb", "Fungi",
                    "fungi", "gooher", "goopher",
                    "gopher", "mosses/lichens",
                    "gophermound","Grass", "hairy",
                    "Hairy", "Leaves", "Lichens",
                    "lichens","Man")|
           species%in%c("Pine cones",
                        "Pine needles",
                        "Quercus elipsoidalis grub"))%>%
  dplyr::select(species)%>%
  distinct()


abun_traits <- abun_traits%>%
  mutate(no_plant = case_when(
    species%in%c(misce$species)~"remove",
    !species%in%c(misce$species)~"keep"
  ))


abun_traits <- abun_traits%>%
  dplyr::select(-c(genus, sp, sub_sp))%>%
  mutate(genus = word(species, 1, sep = fixed(" ")),
         sp = word(species, 2, sep = fixed(" ")),
         sub_sp= word(species, 3, sep = fixed(" ")))# correct genus and sp names


abun_traits%>%
  group_by(no_plant)%>%
  summarise(mean_cover=mean(cover_percent))

abun_traits <- abun_traits%>%
  mutate(transect_code=paste0(Year,"_", OldField,"_", Transect),
         species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)),
         sp_code=paste0(species_code,"_",
                        plot_code))

nrow(abun_traits)
abun_traits%>%
  dplyr::select(-cover_percent)%>%
  distinct()%>%
  nrow()

#* summarise cover data at species level within plot 
#* in case it may have some 
#* error in the original database. 
#* Actually there are 46 repeated species within plot

abun_traits1 <- abun_traits%>%
  dplyr::select(-cover_percent)%>%
  distinct()#table at species x plot level

abun_cover <- abun_traits%>%
  rename(cover_percent1= cover_percent)%>%
  group_by(plot_code, species)%>%
  summarise(cover_percent=sum(cover_percent1))%>%
  ungroup()

abun_traits_new <- abun_traits1%>%
  left_join(abun_cover, by=c("species","plot_code"))

nrow(abun_traits1)
nrow(abun_cover)
nrow(abun_traits_new)

abun_traits <- abun_traits_new
rm(abun_traits1, abun_cover, abun_traits_new)

names(abun_traits)

#* Cover percent data are not scaled by 100% create column of total
#* measured percent per plot and total after discounting non plant species
#* Then add new column of cover proportion -don't use that to compare between plots-

old_plot_cover <- abun_traits%>%
  group_by(plot_code)%>%
  summarise(old_total_cover=sum(cover_percent))

new_plot_cover <- abun_traits%>%
  filter(no_plant=="keep")%>%
  group_by(plot_code)%>%
  summarise(new_total_cover=sum(cover_percent))

abun_traits1 <- list(abun_traits, old_plot_cover, new_plot_cover)%>%
  purrr::reduce(left_join, by="plot_code")


nrow(abun_traits1)
nrow(abun_traits)

abun_traits <- abun_traits1
rm(abun_traits1)

abun_traits <- abun_traits%>%
  mutate(real_cover_percent=
           (cover_percent*100)/old_total_cover,
         new_cover_percent=
           (cover_percent*100)/new_total_cover)

names(abun_traits)

z <- abun_traits%>%
  dplyr::select(Year, OldField,
         Transect, Plot,
         plot_code, species,
         cover_percent,
         old_total_cover,
         new_total_cover,
         real_cover_percent,
         new_cover_percent)%>%
  arrange(plot_code)

abun_traits%>%
  dplyr::select(real_cover_percent)%>%
  pull()%>%hist()

abun_traits%>%
  filter(no_plant=="keep")%>%
  dplyr::select(new_cover_percent)%>%
  pull()%>%hist()

## split complete cases dataset

basic_info <- abun_traits%>%
  dplyr::select(Exp, Year, OldField, 
         Transect, Plot, YearAb,
         plot_code, transect_code,
         species, genus, sp, sub_sp, 
         catford_species, 
         old_total_cover,
         new_total_cover,
         cover_percent,
         real_cover_percent,
         new_cover_percent,
         Functional_group, Duration, 
         Lifeform, Pathway, Taxon, 
         Specid, Origin, Family,
         sp_code, no_plant)

# df_traits <- abun_traits%>%
#   dplyr::select(species, SLA, seed_mass,
#          leaf_area,LDMC,
#          height,leaf_fresh_mass,
#          leaf_dry_mass)

complete_traits <- abun_traits%>%
  dplyr::select(Exp, Year, OldField, 
         Transect, Plot, YearAb,
         plot_code, transect_code,
         species, genus, sp, sub_sp, 
         catford_species,
         species_code, sp_code,
         old_total_cover,
         new_total_cover,
         cover_percent,
         real_cover_percent,
         new_cover_percent,
         no_plant,
         Functional_group, Duration, 
         Lifeform, Pathway, Taxon, 
         Specid, Origin, Family,
         SLA, seed_mass,
         leaf_area,LDMC,
         height,leaf_fresh_mass,
         leaf_dry_mass)%>%
  filter(complete.cases(SLA, seed_mass,
                        leaf_area,LDMC,
                        height,leaf_fresh_mass,
                        leaf_dry_mass))

nrow(abun_traits)
nrow(basic_info)
abun_traits%>%filter(no_plant=="keep")%>%
  nrow()
nrow(complete_traits)

# data summary of representativeness at table, plot and transect level

basic_info%>%
  dplyr::select(-c(Exp, Year, OldField, 
            Transect, Plot, YearAb, 
            plot_code, transect_code,
            cover_percent))%>%
  distinct()%>%
  nrow()

complete_traits%>%
  dplyr::select(species, SLA, seed_mass,
         leaf_area,LDMC,
         height,leaf_fresh_mass,
         leaf_dry_mass)%>%
  distinct()%>%
  nrow()

sort(unique(complete_traits$species))
nrow(complete_traits)
nrow(abun_traits)

#************************************************##
#** prepare df with the dplyr::selected set of traits ####
#************************************************##

# abun_traits <- abun_traits%>%
#   filter(no_plant=="keep")

# prepare species data for pca

complete_traits <- complete_traits%>%
  filter(genus!="Quercus")

complete_traits_sp <- complete_traits%>%
  dplyr::select(species, SLA, seed_mass,
         leaf_area,LDMC,
         height,leaf_fresh_mass,
         leaf_dry_mass)%>%
  distinct()%>%
  tibble::column_to_rownames('species')#only complete cases


#******************************##
#** get distances univariant ####
#******************************##

#**** using community trait centroid ####

source('niche_required_functions.R')
complete_traits

# table preparation

table <- complete_traits%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))%>%
  mutate(sp_code=paste0(species_code,"_",
                        plot_code))

centroid_df <-  table%>%
  dplyr::select(names(basic_info),
                species_code,
                sp_code,
                SLA, seed_mass,
                leaf_area,
                LDMC,
                height, 
                leaf_fresh_mass,
                leaf_dry_mass)%>%
  distinct()

plot_list <- unique(centroid_df$plot_code)
fun_traits <- c("SLA", "seed_mass",
                "leaf_area",
                "LDMC",
                "height", 
                "leaf_fresh_mass",
                "leaf_dry_mass")

table_plot_list <- list()

for(i in 1:length(plot_list)) {
  
  table_plot <- table%>%
    filter(plot_code== plot_list[i])
  
  species_list <- unique(table_plot$sp_code)
  table_plot_sp_list <- list()
  
  for(j in 1:length(species_list)){
    
    table_plot_sp <- table_plot%>%
      filter(sp_code== species_list[j])
    
    table_plot_rest <- table_plot%>%
      filter(species!= species_list[j])%>%
      dplyr::select(c(sp_code,species, 
               cover_percent,
               all_of(fun_traits)))
    
    species_id <- table_plot%>%
      filter(sp_code== species_list[j])%>%
      dplyr::select(sp_code)%>%as.character()
    
    trait_centroid_df <- tibble(
      centroid_SLA=COGravity(x=table_plot_rest$SLA, 
                             wt=table_plot_rest$cover_percent)[1],#same to real_cover_percent
      centroid_seed_mass=COGravity(x=table_plot_rest$seed_mass, 
                                   wt=table_plot_rest$cover_percent)[1],
      centroid_leaf_area=COGravity(x=table_plot_rest$leaf_area, 
                                   wt=table_plot_rest$cover_percent)[1],
      centroid_LDMC=COGravity(x=table_plot_rest$LDMC, 
                              wt=table_plot_rest$cover_percent)[1],
      centroid_height=COGravity(x=table_plot_rest$height, 
                                wt=table_plot_rest$cover_percent)[1],
      centroid_leaf_fresh_mass=COGravity(x=table_plot_rest$leaf_fresh_mass, 
                             wt=table_plot_rest$cover_percent)[1],
      centroid_leaf_dry_mass=COGravity(x=table_plot_rest$leaf_dry_mass, 
                             wt=table_plot_rest$cover_percent)[1])%>%
      as.data.frame()
    
    table_plot_sp <- cbind(table_plot_sp,
                           trait_centroid_df)
    
    table_plot_sp <- table_plot_sp%>%
      mutate(dissim_SLA=SLA-centroid_SLA, # important to keep first the target species first so it is easier to discuss that the target species is hierarchically superior
             dissim_SLA_abs=abs(SLA - centroid_SLA),
             dissim_seed_mass=seed_mass -centroid_seed_mass,
             dissim_seed_mass_abs=abs(seed_mass-centroid_seed_mass),
             dissim_leaf_area=leaf_area-centroid_leaf_area,
             dissim_leaf_area_abs=abs(leaf_area-centroid_leaf_area),
             dissim_LDMC=LDMC-centroid_LDMC,
             dissim_LDMC_abs=abs(LDMC-centroid_LDMC),
             dissim_height=height-centroid_height,
             dissim_height_abs=abs(height-centroid_height),
             dissim_leaf_fresh_mass=leaf_fresh_mass-centroid_leaf_fresh_mass,
             dissim_leaf_fresh_mass_abs=abs(leaf_fresh_mass-centroid_leaf_fresh_mass),
             dissim_leaf_dry_mass=leaf_dry_mass-centroid_leaf_dry_mass,
             dissim_leaf_dry_mass_abs=abs(leaf_dry_mass-centroid_leaf_dry_mass))
    
    
    # table_plot$sp_code <- factor(table_plot$sp_code)
    # require(RColorBrewer)
    # cols <- colorRampPalette(brewer.pal(11, "Spectral"))
    # display.brewer.pal(n=11, name="Spectral")
    # cols1 <- brewer.pal(n=11, name="Spectral")
    # scales::show_col(cols1)
    # mypal <- cols1[c(4, 8, 10)]
    # 
    # gg_i <- ggplot(data=table_plot_rest,
    #                aes(x=SLA)) +
    #   geom_bar(aes(y=cover_percent,
    #                fill=species),
    #            stat="identity",
    #            position="dodge", 
    #            width = 0.35)+
    #   scale_fill_manual(values = cols1)+
    #   geom_vline(xintercept=table_plot_sp$SLA,
    #              color="red")+
    #   geom_vline(xintercept=table_plot_sp$centroid_SLA,
    #              color="black")+
    #   labs(color = "Species")+
    #   guides(color=guide_legend(ncol=2))+
    #   xlab("SLA")+
    #   ggtitle(plot_list[i])+
    #   theme_linedraw()+
    #   theme(plot.title = element_text(hjust = 0.5),
    #         legend.position = "bottom",
    #         legend.title = element_text(size = 9),
    #         legend.text  = element_text(size = 8),
    #         legend.key.size = unit(5, "mm"),
    #         panel.grid.major=element_blank(),
    #         panel.grid.minor=element_blank(),
    #         axis.text=element_text(color="black"),
    #         text=element_text(family="serif"))
    # 
    # 
    # ggsave(paste0("../../output/community_diagrams/", plot_list[i], "univariant.png"), 
    #        plot=gg_i, width = 9, height = 9, dpi = 300, units="cm")
    # ggsave(paste0("../../output/community_diagrams/", plot_list[i], "univariant.pdf"),
    #        plot=gg_i, width = 9, height = 9, units="cm")
    
    
    table_plot_sp_list[[j]] <- table_plot_sp
    
    
  }
  
  table_plot_new <- table_plot_sp_list%>%
    bind_rows()
  table_plot_list[[i]] <- table_plot_new
  
  print(i/length(plot_list)*100) 
  
  
}


centroid_univ_df <- table_plot_list%>%
  bind_rows()

nrow(table)
nrow(centroid_univ_df)

setdiff(table$sp_code, centroid_univ_df$sp_code)

data.table::fwrite(centroid_univ_df,
                   "../../output/tables/dissim/dissimilarity_centroid_univariate_plot_2024.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)



#**** using weighted pairwaise distances ####


source('niche_required_functions.R')
complete_traits

# table preparation

table <- complete_traits%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))%>%
  mutate(sp_code=paste0(species_code,"_",
                        plot_code))

centroid_df <-  table%>%
  dplyr::select(names(basic_info),
                species_code,
                sp_code,
                SLA, seed_mass,
                leaf_area,
                LDMC,
                height, 
                leaf_fresh_mass,
                leaf_dry_mass)%>%
  distinct()

plot_list <- unique(centroid_df$plot_code)
fun_traits <- c("SLA", "seed_mass",
                "leaf_area",
                "LDMC",
                "height", 
                "leaf_fresh_mass",
                "leaf_dry_mass")

table_plot_list <- list()

for(i in 1:length(plot_list)) {
  
  table_plot <- table%>%
    filter(plot_code== plot_list[i])
  
  species_list <- unique(table_plot$sp_code)
  table_plot_sp_list <- list()
  
  for(j in 1:length(species_list)){
    
    table_plot_sp <- table_plot%>%
      filter(sp_code== species_list[j])
    
    table_plot_rest <- table_plot%>%
      filter(sp_code!= species_list[j])
    
    if(nrow(table_plot_rest)>=1){
      
      trait_sp <- table_plot_sp%>%
        dplyr::select(all_of(fun_traits))%>%
        distinct()
      
      table_plot_rest <- table_plot_rest%>%
        mutate(dissim_SLA=abs(SLA-trait_sp$SLA),
               dissim_seed_mass=abs(seed_mass-trait_sp$seed_mass),
               dissim_leaf_area=abs(leaf_area-trait_sp$leaf_area),
               dissim_LDMC= abs(LDMC-trait_sp$LDMC),
               dissim_height=abs(height-trait_sp$height),
               dissim_leaf_fresh_mass=abs(leaf_fresh_mass-trait_sp$leaf_fresh_mass),
               dissim_leaf_dry_mass=abs(leaf_dry_mass-trait_sp$leaf_dry_mass))
      
      av_dissim <- table_plot_rest%>%
        dplyr::select(c(sp_code,cover_percent,
                 starts_with("dissim")))%>%
        summarise_at(vars(-c(sp_code,cover_percent)),
                     funs(weighted.mean(., w=cover_percent)))#also real_cover_percent
      
    }else{
      
     av_dissim <- tibble(
       dissim_SLA=NA_real_,
       dissim_seed_mass=NA_real_,
       dissim_leaf_area=NA_real_,
       dissim_height=NA_real_,
       dissim_LDMC=NA_real_,
       dissim_leaf_fresh_mass=NA_real_,
       dissim_leaf_dry_mass=NA_real_
     )
      
      
    }
   
    table_plot_sp <- cbind(table_plot_sp,
                           av_dissim)
    
    table_plot_sp_list[[j]] <- table_plot_sp
    
    
  }
  
  table_plot_new <- table_plot_sp_list%>%
    bind_rows()
  table_plot_list[[i]] <- table_plot_new
 
  print(i/length(plot_list)*100) 
  
}

pairwise_df <- table_plot_list%>%
  bind_rows()

nrow(table)
nrow(pairwise_df)

setdiff(table$sp_code, pairwise_df$sp_code)
  
data.table::fwrite(pairwise_df,
                   "../../output/tables/dissim/dissimilarity_univariant_pairwise_abs_plot.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)


###**********************************************###
###******2. TRANSECT LEVEL DISSIMILARITIES******####
###**********************************************###

#** check complete cases for main variables ####

abun_traits <- read_delim("../../output/tables/abundances/ab_traits_complete_df1.csv", 
                          delim=";", col_names=T)

abun_traits <- abun_traits%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D",
    Transect_old=="1"~"A",
    Transect_old=="2"~"B",
    Transect_old=="3"~"C",
    Transect_old=="4"~"D",
    Transect_old=="5"~"E",
    Transect_old=="6"~"F",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ),
  plot_code=paste0(Year,"_",
                   OldField, "_",
                   Transect, "_",
                   Plot))

misce <- abun_traits%>%
  filter(genus%in%c("Miscellaneous","miscellaneous",
                    "Unknown", "unknown", "woody",
                    "Woody", "Mosses", "mosses",
                    "animal", "ant", "anthill",
                    "anthills", "Bare", "dung",
                    "fern", "Forb", "Fungi",
                    "fungi", "gooher", "goopher",
                    "gopher", "mosses/lichens",
                    "gophermound","Grass", "hairy",
                    "Hairy", "Leaves", "Lichens",
                    "lichens","Man")|
           species%in%c("Pine cones",
                        "Pine needles",
                        "Quercus elipsoidalis grub"))%>%
  dplyr::select(species)%>%
  distinct()


abun_traits <- abun_traits%>%
  mutate(no_plant = case_when(
    species%in%c(misce$species)~"remove",
    !species%in%c(misce$species)~"keep"
  ))


abun_traits <- abun_traits%>%
  dplyr::select(-c(genus, sp, sub_sp))%>%
  mutate(genus = word(species, 1, sep = fixed(" ")),
         sp = word(species, 2, sep = fixed(" ")),
         sub_sp= word(species, 3, sep = fixed(" ")))# correct genus and sp names


abun_traits%>%
  group_by(no_plant)%>%
  summarise(mean_cover=mean(cover_percent))

abun_traits <- abun_traits%>%
  mutate(transect_code=paste0(Year,"_", OldField,"_", Transect),
         species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)),
         sp_code=paste0(species_code,"_",
                        plot_code))

nrow(abun_traits)
abun_traits%>%
  dplyr::select(-cover_percent)%>%
  distinct()%>%
  nrow()

#* summarise cover data at species level within plot 
#* in case it may have some 
#* error in the original database. 
#* Actually there are 46 repeated species within plot

abun_traits1 <- abun_traits%>%
  dplyr::select(-cover_percent)%>%
  distinct()#table at species x plot level

abun_cover <- abun_traits%>%
  rename(cover_percent1= cover_percent)%>%
  group_by(plot_code, species)%>%
  summarise(cover_percent=sum(cover_percent1))%>%
  ungroup()

abun_traits_new <- abun_traits1%>%
  left_join(abun_cover, by=c("species","plot_code"))

nrow(abun_traits1)
nrow(abun_cover)
nrow(abun_traits_new)

abun_traits <- abun_traits_new
rm(abun_traits1, abun_cover, abun_traits_new)

names(abun_traits)
nrow(abun_traits)

#* Cover percent data are not scaled by 100% create column of total
#* measured percent per plot and total after discounting non plant species
#* Then add new column of cover proportion -don't use that to compare between plots-

old_plot_cover <- abun_traits%>%
  group_by(plot_code)%>%
  summarise(old_total_cover=sum(cover_percent))

new_plot_cover <- abun_traits%>%
  filter(no_plant=="keep")%>%
  group_by(plot_code)%>%
  summarise(new_total_cover=sum(cover_percent))

abun_traits1 <- list(abun_traits, old_plot_cover, new_plot_cover)%>%
  purrr::reduce(left_join, by="plot_code")


nrow(abun_traits1)
nrow(abun_traits)

abun_traits <- abun_traits1
rm(abun_traits1)

abun_traits <- abun_traits%>%
  mutate(real_cover_percent=
           (cover_percent*100)/old_total_cover,#scale by 100% of the sampled surface per plot
         new_cover_percent=
           (cover_percent*100)/new_total_cover)#scale by surface of study plants

names(abun_traits)

z <- abun_traits%>%
  dplyr::select(Year, OldField,
         Transect, Plot,
         plot_code, species,
         cover_percent,
         old_total_cover,
         new_total_cover,
         real_cover_percent,
         new_cover_percent)%>%
  arrange(plot_code)

abun_traits%>%
  dplyr::select(real_cover_percent)%>%
  pull()%>%hist()

abun_traits%>%
  filter(no_plant=="keep")%>%
  dplyr::select(new_cover_percent)%>%
  pull()%>%hist()

## split complete cases dataset

basic_info <- abun_traits%>%
  dplyr::select(Exp, Year, OldField, 
         Transect, Plot, YearAb,
         plot_code, transect_code,
         species, genus, sp, sub_sp, 
         catford_species, 
         old_total_cover,
         new_total_cover,
         cover_percent,
         real_cover_percent,
         new_cover_percent,
         Functional_group, Duration, 
         Lifeform, Pathway, Taxon, 
         Specid, Origin, Family,
         sp_code, no_plant)

# df_traits <- abun_traits%>%
#   dplyr::select(species, SLA, seed_mass,
#          leaf_area,LDMC,
#          height,leaf_fresh_mass,
#          leaf_dry_mass)

complete_traits <- abun_traits%>%
  dplyr::select(Exp, Year, OldField, 
         Transect, Plot, YearAb,
         plot_code, transect_code,
         species, genus, sp, sub_sp, 
         catford_species,
         species_code, sp_code,
         old_total_cover,
         new_total_cover,
         cover_percent,
         real_cover_percent,
         new_cover_percent,
         no_plant,
         Functional_group, Duration, 
         Lifeform, Pathway, Taxon, 
         Specid, Origin, Family,
         SLA, seed_mass,
         leaf_area,LDMC,
         height,leaf_fresh_mass,
         leaf_dry_mass)%>%
  filter(complete.cases(SLA, seed_mass,
                        leaf_area,LDMC,
                        height,leaf_fresh_mass,
                        leaf_dry_mass))

nrow(abun_traits)
nrow(basic_info)
abun_traits%>%filter(no_plant=="keep")%>%
  nrow()
nrow(complete_traits)

# data summary of representativeness at table, plot and transect level

basic_info%>%
  dplyr::select(-c(Exp, Year, OldField, 
            Transect, Plot, YearAb, 
            plot_code, transect_code,
            cover_percent))%>%
  distinct()%>%
  nrow()

complete_traits%>%
  dplyr::select(species, SLA, seed_mass,
         leaf_area,LDMC,
         height,leaf_fresh_mass,
         leaf_dry_mass)%>%
  distinct()%>%
  nrow()

sort(unique(complete_traits$species))
nrow(complete_traits)
nrow(abun_traits)

complete_traits <- complete_traits%>%
  filter(genus!="Quercus")

#******************************##
#** get distances univariant ####
#******************************##

#**** using community trait centroid ####

source('niche_required_functions.R')
complete_traits

complete_traits <- complete_traits%>%
  filter(genus!="Quercus")


# table preparation

table_cat <- complete_traits%>%
  dplyr::select(-c(Plot, plot_code, 
            sp_code,
            old_total_cover,
            new_total_cover, 
            cover_percent,
            real_cover_percent,
            new_cover_percent))%>%
  distinct()%>%
  mutate(sp_code=paste0(species_code,"_",
                        transect_code))

nrow(table_cat)

nplot_trans <- complete_traits%>%
  dplyr::select(transect_code, Plot)%>%
  distinct()%>%
  group_by(transect_code)%>%
  summarise(n=n())%>%
  ungroup()%>%
  mutate(transect_cover=n*100)

hist((nplot_trans$transect_cover))
unique(nplot_trans$n)
nplot_trans%>%
  filter(n<25)

table_cov <- complete_traits%>%
  dplyr::select(transect_code,
         species_code,
         old_total_cover,
         new_total_cover, 
         cover_percent,
         real_cover_percent,
         new_cover_percent)%>%
  group_by(transect_code, species_code)%>%
  summarise_all(sum, na.rm=T)%>%
  ungroup()%>%
  mutate(sp_code=paste0(species_code,"_",
                        transect_code))

table_cov1 <- table_cov%>%
  left_join(nplot_trans, by="transect_code")%>%
  mutate(mean_cover_percent=cover_percent/transect_cover)#cover percent is sum at transect level, so devide by total transect cover
#do not use the rest of cover values

nrow(table_cov)
nrow(table_cov1)
table_cov <- table_cov1

hist(complete_traits$cover_percent)
hist(table_cov$mean_cover_percent)
hist(table_cov$cover_percent)
hist(table_cov$transect_cover)

table <- table_cat%>%
  left_join(table_cov%>%
              dplyr::select(setdiff(names(table_cov),
                             names(table_cat)),
                     sp_code), by="sp_code")

nrow(table)
names(table)

basic_vars <- names(basic_info)[!names(basic_info)%in%c("Plot",
                                                        "plot_code")]

centroid_df <-  table%>%
  dplyr::select(basic_vars,
                species_code,
                sp_code,
                SLA, seed_mass,
                leaf_area,
                LDMC,
                height, 
                leaf_fresh_mass,
                leaf_dry_mass)%>%
  distinct()

nrow(centroid_df)

transect_list <- unique(centroid_df$transect_code)
fun_traits <- c("SLA", "seed_mass",
                "leaf_area",
                "LDMC",
                "height", 
                "leaf_fresh_mass",
                "leaf_dry_mass")

table_transect_list <- list()



for(i in 1:length(transect_list)) {
  
  table_transect <- table%>%
    filter(transect_code== transect_list[i])
  
  species_list <- unique(table_transect$sp_code)
  table_transect_sp_list <- list()
  
  for(j in 1:length(species_list)){
    
    table_transect_sp <- table_transect%>%
      filter(sp_code== species_list[j])
    
    table_transect_rest <- table_transect%>%
      filter(species!= species_list[j])%>%
      dplyr::select(c(sp_code,species, 
               mean_cover_percent,
               all_of(fun_traits)))
    
    species_id <- table_transect%>%
      filter(sp_code== species_list[j])%>%
      dplyr::select(sp_code)%>%as.character()
    
    trait_centroid_df <- tibble(
      centroid_SLA=COGravity(x=table_transect_rest$SLA, 
                             wt=table_transect_rest$mean_cover_percent)[1],#same to real_cover_percent
      centroid_seed_mass=COGravity(x=table_transect_rest$seed_mass, 
                                   wt=table_transect_rest$mean_cover_percent)[1],
      centroid_leaf_area=COGravity(x=table_transect_rest$leaf_area, 
                                   wt=table_transect_rest$mean_cover_percent)[1],
      centroid_LDMC=COGravity(x=table_transect_rest$LDMC, 
                              wt=table_transect_rest$mean_cover_percent)[1],
      centroid_height=COGravity(x=table_transect_rest$height, 
                                wt=table_transect_rest$mean_cover_percent)[1],
      centroid_leaf_fresh_mass=COGravity(x=table_transect_rest$leaf_fresh_mass, 
                                         wt=table_transect_rest$mean_cover_percent)[1],
      centroid_leaf_dry_mass=COGravity(x=table_transect_rest$leaf_dry_mass, 
                                       wt=table_transect_rest$mean_cover_percent)[1])%>%
      as.data.frame()
    
    table_transect_sp <- cbind(table_transect_sp,
                               trait_centroid_df)
    
    table_transect_sp <- table_transect_sp%>%
      mutate(dissim_SLA=SLA-centroid_SLA,#important to keep target species value before to simplify interpretation (the target species is higher or lower than the rest of the community)
             dissim_SLA_abs=abs(SLA-centroid_SLA),#abs is comparable to multivariate distance, but it worth to keep the value of the difference 
             dissim_seed_mass=seed_mass-centroid_seed_mass,
             dissim_seed_mass_abs=abs(seed_mass-centroid_seed_mass),
             dissim_leaf_area=leaf_area-centroid_leaf_area,
             dissim_leaf_area_abs=abs(leaf_area-centroid_leaf_area),
             dissim_LDMC=LDMC-centroid_LDMC,
             dissim_LDMC_abs=abs(LDMC-centroid_LDMC),
             dissim_height=height-centroid_height,
             dissim_height_abs=abs(height-centroid_height),
             dissim_leaf_fresh_mass=leaf_fresh_mass-centroid_leaf_fresh_mass,
             dissim_leaf_fresh_mass_abs=abs(leaf_fresh_mass-centroid_leaf_fresh_mass),
             dissim_leaf_dry_mass=leaf_dry_mass-centroid_leaf_dry_mass,
             dissim_leaf_dry_mass_abs=abs(leaf_dry_mass-centroid_leaf_dry_mass))
    
    
    # table_transect$sp_code <- factor(table_transect$sp_code)
    # require(RColorBrewer)
    # cols <- colorRampPalette(brewer.pal(11, "Spectral"))
    # display.brewer.pal(n=11, name="Spectral")
    # cols1 <- brewer.pal(n=11, name="Spectral")
    # scales::show_col(cols1)
    # mypal <- cols1[c(4, 8, 10)]
    # 
    # gg_i <- ggplot(data=table_transect_rest,
    #                aes(x=SLA)) +
    #   geom_bar(aes(y=cover_percent,
    #                fill=species),
    #            stat="identity",
    #            position="dodge", 
    #            width = 0.35)+
    #   scale_fill_manual(values = cols1)+
    #   geom_vline(xintercept=table_transect_sp$SLA,
    #              color="red")+
    #   geom_vline(xintercept=table_transect_sp$centroid_SLA,
    #              color="black")+
    #   labs(color = "Species")+
    #   guides(color=guide_legend(ncol=2))+
    #   xlab("SLA")+
    #   ggtitle(transect_list[i])+
    #   theme_linedraw()+
    #   theme(plot.title = element_text(hjust = 0.5),
    #         legend.position = "bottom",
    #         legend.title = element_text(size = 9),
    #         legend.text  = element_text(size = 8),
    #         legend.key.size = unit(5, "mm"),
    #         panel.grid.major=element_blank(),
    #         panel.grid.minor=element_blank(),
    #         axis.text=element_text(color="black"),
    #         text=element_text(family="serif"))
    # 
    # 
    # ggsave(paste0("../../output/community_diagrams/", transect_list[i], "univariant.png"), 
    #        plot=gg_i, width = 9, height = 9, dpi = 300, units="cm")
    # ggsave(paste0("../../output/community_diagrams/", transect_list[i], "univariant.pdf"),
    #        plot=gg_i, width = 9, height = 9, units="cm")
    
    
    table_transect_sp_list[[j]] <- table_transect_sp
    
    
  }
  
  table_transect_new <- table_transect_sp_list%>%
    bind_rows()
  table_transect_list[[i]] <- table_transect_new
  
  print(i/length(transect_list)*100) 
  
}


centroid_univ_df <- table_transect_list%>%
  bind_rows()

nrow(table)
nrow(centroid_univ_df)

setdiff(table$sp_code, centroid_univ_df$sp_code)

data.table::fwrite(centroid_univ_df,
                   "../../output/tables/dissim/dissimilarity_centroid_univariate_transect_2024.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)



#**** using weighted pairwaise distances ####


source('niche_required_functions.R')
complete_traits

complete_traits%>%
  filter(genus!="Quercus")

# table preparation

table_cat <- complete_traits%>%
  dplyr::select(-c(Plot, plot_code, 
            sp_code,
            old_total_cover,
            new_total_cover, 
            cover_percent,
            real_cover_percent,
            new_cover_percent))%>%
  distinct()%>%
  mutate(sp_code=paste0(species_code,"_",
                        transect_code))

nrow(table_cat)

nplot_trans <- complete_traits%>%
  dplyr::select(transect_code, Plot)%>%
  distinct()%>%
  group_by(transect_code)%>%
  summarise(n=n())%>%
  ungroup()%>%
  mutate(transect_cover=n*100)

hist((nplot_trans$transect_cover))
unique(nplot_trans$n)
nplot_trans%>%
  filter(n<25)

table_cov <- complete_traits%>%
  dplyr::select(transect_code,
         species_code,
         old_total_cover,
         new_total_cover, 
         cover_percent,
         real_cover_percent,
         new_cover_percent)%>%
  group_by(transect_code, species_code)%>%
  summarise_all(sum, na.rm=T)%>%
  ungroup()%>%
  mutate(sp_code=paste0(species_code,"_",
                        transect_code))

table_cov1 <- table_cov%>%
  left_join(nplot_trans, by="transect_code")%>%
  mutate(mean_cover_percent=cover_percent/transect_cover)#cover percent is sum at transect level, so devide by total transect cover
#do not use the rest of cover values

nrow(table_cov)
nrow(table_cov1)
table_cov <- table_cov1

hist(complete_traits$cover_percent)
hist(table_cov$mean_cover_percent)
hist(table_cov$cover_percent)
hist(table_cov$transect_cover)

table <- table_cat%>%
  left_join(table_cov%>%
              dplyr::select(setdiff(names(table_cov),
                             names(table_cat)),
                     sp_code), by="sp_code")

nrow(table)
names(table)

basic_vars <- names(basic_info)[!names(basic_info)%in%c("Plot",
                                                        "plot_code")]

centroid_df <-  table%>%
  dplyr::select(basic_vars,
                species_code,
                sp_code,
                SLA, seed_mass,
                leaf_area,
                LDMC,
                height, 
                leaf_fresh_mass,
                leaf_dry_mass)%>%
  distinct()

nrow(centroid_df)

transect_list <- unique(centroid_df$transect_code)
fun_traits <- c("SLA", "seed_mass",
                "leaf_area",
                "LDMC",
                "height", 
                "leaf_fresh_mass",
                "leaf_dry_mass")

table_transect_list <- list()


for(i in 1:length(transect_list)) {
  
  table_transect <- table%>%
    filter(transect_code== transect_list[i])
  
  species_list <- unique(table_transect$sp_code)
  table_transect_sp_list <- list()
  
  for(j in 1:length(species_list)){
    
    table_transect_sp <- table_transect%>%
      filter(sp_code== species_list[j])
    
    table_transect_rest <- table_transect%>%
      filter(sp_code!= species_list[j])
    
    if(nrow(table_transect_rest)>=1){
      
      trait_sp <- table_transect_sp%>%
        dplyr::select(all_of(fun_traits))%>%
        distinct()
      
      table_transect_rest <- table_transect_rest%>%
        mutate(dissim_SLA=abs(SLA-trait_sp$SLA),
               dissim_seed_mass=abs(seed_mass-trait_sp$seed_mass),
               dissim_leaf_area=abs(leaf_area-trait_sp$leaf_area),
               dissim_LDMC= abs(LDMC-trait_sp$LDMC),
               dissim_height=abs(height-trait_sp$height),
               dissim_leaf_fresh_mass=abs(leaf_fresh_mass-trait_sp$leaf_fresh_mass),
               dissim_leaf_dry_mass=abs(leaf_dry_mass-trait_sp$leaf_dry_mass))
      
      av_dissim <- table_transect_rest%>%
        dplyr::select(c(sp_code,mean_cover_percent,
                 starts_with("dissim")))%>%
        summarise_at(vars(-c(sp_code,mean_cover_percent)),
                     funs(weighted.mean(., w=mean_cover_percent)))#also real_cover_percent
      
    }else{
      
      av_dissim <- tibble(
        dissim_SLA=NA_real_,
        dissim_seed_mass=NA_real_,
        dissim_leaf_area=NA_real_,
        dissim_LDMC= NA_real_,
        dissim_height=NA_real_,
        dissim_leaf_fresh_mass=NA_real_,
        dissim_leaf_dry_mass=NA_real_
      )
      
      
    }
    
    table_transect_sp <- cbind(table_transect_sp,
                               av_dissim)
    
    table_transect_sp_list[[j]] <- table_transect_sp
    
    
  }
  
  table_transect_new <- table_transect_sp_list%>%
    bind_rows()
  table_transect_list[[i]] <- table_transect_new
  
  print(i/length(transect_list)*100) 
  
  
}

pairwise_df <- table_transect_list%>%
  bind_rows()

nrow(table)
nrow(pairwise_df)

setdiff(table$sp_code, pairwise_df$sp_code)

data.table::fwrite(pairwise_df,
                   "../../output/dissimilarity_univariant_pairwise_abs_transect.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)




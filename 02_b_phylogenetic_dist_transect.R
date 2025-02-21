



##***********************************##
##* Maria Angeles Pérez-Navarro
##* King's College University
##* AlienImpacts project
##* Phylogenetic distances transect
##* April 2022
##* *********************************##

library(dplyr)
library(readr)
library(tidyr)
library(XML)
library(EML)#
library("methods")
library(textreadr)#
library(stringr)
library(radiant.data)
library(tidyverse)
library(taxize)
library(Taxonstand)
library(ade4)
library(ape)
library(phytools)
library(phangorn)
library(geiger)
library(V.PhyloMaker)
#devtools::install_github("jinyizju/V.PhyloMaker")
library(dendextend)
library(phylogram)
library(circlize)
#devtools::install_github("BlakeRMills/MetBrewer") 
library(MetBrewer)
library(RColorBrewer)


#** check complete cases for main variables ####

abun_traits <- read_delim("../../output/tables/abundances/ab_traits_complete_df1.csv", 
                          delim=";", col_names=T)



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
  select(species)%>%
  distinct()


abun_traits <- abun_traits%>%
  mutate(no_plant = case_when(
    species%in%c(misce$species)~"remove",
    !species%in%c(misce$species)~"keep"
  ))


abun_traits <- abun_traits%>%
  select(-c(genus, sp, sub_sp))%>%
  mutate(genus = word(species, 1, sep = fixed(" ")),
         sp = word(species, 2, sep = fixed(" ")),
         sub_sp= word(species, 3, sep = fixed(" ")))# correct genus and sp names

abun_traits <- abun_traits%>%
  mutate(species1=case_when(
    sp== "sp."~ genus,
    sp!= "sp."~ species
  ))


# Species names in the plant list

#* check that for searching accepted names
#* https://cran.r-project.org/web/packages/Taxonstand/Taxonstand.pdf
#* when the species is at genus level it has problem to
#* find the species. Try to remove sp from species name

# 
# species <- unique(abun_traits$species1)
# new_names_list <- list()
# 
# for (i in 1:35){
# 
#   start <- ((i-1)*10)+1
#   end <- i*10
#   sp_list <- species[c(start:end)]
# 
#   r1 <- TPL(sp_list, corr = TRUE)
#   str(r1)
#   head(r1)
#   if(is.numeric(r1$Date)){
#     r1$Date <- as.Date(r1$Date, origin="2022-02-15")
#   }else{ }
# 
# 
#   new_names_list[[i]] <- r1
# 
# }
# 
# new_names <- bind_rows(new_names_list)
# names(new_names)
# head(new_names)
# nrow(new_names)
# 
# 
# tpl_names <- new_names%>%
#   filter(!is.na(Taxon))%>%
#   mutate(tpl_species=case_when(
#     is.na(New.Genus)~ paste0(Taxon, " sp."),
#     !is.na(New.Genus)&is.na(New.Species)~paste0(New.Genus,
#                                         " sp."),
#     !is.na(New.Genus)&!is.na(New.Species)~paste0(New.Genus,
#                                                " ", New.Species)
#   ),
#   tpl_genus=case_when(
#     is.na(New.Genus)~Taxon,
#     !is.na(New.Genus)~New.Genus
#   ))%>%
#   select(Taxon,
#          Family,
#          tpl_genus,
#          New.Species,
#          tpl_species)%>%
#   rename(taxon=Taxon,
#          tpl_family=Family,
#          tpl_sp=New.Species)
# 
# 
# data.table::fwrite(tpl_names,
#                    "../../data/tpl_names.csv",
#                    sep=";", dec=".", col.names=T, row.names=F, quote=F)


tpl_names <- read_delim("../../data/tpl_names.csv", 
                        delim=";", col_names=T)


abun_traits1 <- abun_traits%>%
  left_join(tpl_names, by=c("species1"="taxon"))%>%
  distinct()

nrow(abun_traits1)
nrow(abun_traits)

abun_traits1%>%
  select(species, tpl_species,
         genus, tpl_genus,
         sp, tpl_sp)

abun_traits <- abun_traits1
rm(abun_traits1)

abun_traits%>%group_by(Year, OldField, Transect)

abun_traits <- abun_traits%>%
  mutate(transect_code=paste0(Year,"_", OldField,"_", Transect))

nrow(abun_traits)
abun_traits%>%
  select(-cover_percent)%>%
  distinct()%>%
  nrow()

#* summarise cover data at species level within plot 
#* in case it may have some 
#* error in the original database. 
#* Actually there are 46 repeated species within plot

abun_traits1 <- abun_traits%>%
  select(-cover_percent)%>%
  distinct()

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

#** add new cover values up to 100% and up to 
#* total cover of known plants

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


#* set basic info

basic_info <- abun_traits%>%
  select(Exp, Year, OldField, 
         Transect, Plot, YearAb,
         plot_code, transect_code,
         species, genus, sp, sub_sp, 
         catford_species, 
         no_plant,
         tpl_family, tpl_genus,
         tpl_sp, tpl_species,
         cover_percent, real_cover_percent,
         new_cover_percent,
         old_total_cover,
         new_total_cover,
         Functional_group, Duration, 
         Lifeform, Pathway, Taxon, 
         Specid, Origin, Family,
         species1)

#*get phylogeny using phylo.maker from V.Phylomaker using Smith and Brown 2018

sp_names <- basic_info%>%
  select(Specid, Family,
         species, genus,
         sp, sub_sp,
         tpl_family,
         tpl_genus,
         tpl_sp,
         tpl_species,
         no_plant,
         species1)%>%
  distinct()%>%
  filter(!is.na(Family))%>%
  filter(Family!="Unknown")%>%
  mutate(tpl_family=na_if(tpl_family, ""))%>%
  mutate(family_def=case_when(
    is.na(tpl_family)~Family,
    !is.na(tpl_family)~tpl_family # some tpl families have na, replace by old family name
  ),
  species_def=case_when(
    sp=="sp."~ substring(tpl_species, 1,
                         nchar(tpl_species)-4),
    sp!="sp."~ tpl_species
  ))%>%
  arrange(tpl_species)

nrow(sp_names)
names(sp_names)
head(sp_names)

sp_file <- sp_names%>%
  filter(no_plant=="keep")%>%
  select(species_def,
         tpl_genus,
         family_def)%>%
  rename(species=species_def,
         genus=tpl_genus,
         family=family_def)%>%
  distinct()

#* adapt family names to phylo.maker format. 
#* Instead of Compositae ~ Asteraceae and 
#* Instead of Leguminosae ~ Fabaceae

sp_file1 <- sp_file%>%
  mutate(family_new=case_when(
    family=="Aceraceae"~ "Sapindaceae",
    family=="Compositae" ~ "Asteraceae",
    family=="Leguminosae" ~ "Fabaceae",
    family!="Leguminosae"|
      family!="Compositae"|
      family!="Aceraceae"~ family
  ))%>%
  select(-family)%>%
  rename(family=family_new)%>%
  distinct()

sp_file <- sp_file1
nrow(sp_file)
rm(sp_file1)

unique(is.na(sp_file$species))

#* calculate phylogenetic tree using V.phyloMaker package based on 
#* smith and brown 2018 plant phylogeny
#* code for plotting https://r-charts.com/part-whole/circular-dendrogram/

whole_tree <- phylo.maker(sp_file, scenarios=c("S1","S2","S3"))

sce1_tree <- whole_tree$scenario.1

par(mfrow = c(1, 1))
plot.phylo(as.phylo(sce1_tree), cex = 0.15, main = "scenario.1",
           type="fan")

plot.phylo(whole_tree$scenario.3, cex = 0.15, main = "scenario.3",
           type="fan")

mycols <- colorRampPalette(brewer.pal(11, "Spectral"))(15)
mycols1 <- c(met.brewer("Lakota",10),met.brewer("Cross",5))
scales::show_col(c(mycols, mycols1))

dend <- as.dendrogram(sce1_tree)

dend <- dend%>%
  color_branches(k=15)%>%
  color_labels(k=15)

dend1 <- dend%>%  #k= number of clusters
  set("branches_k_color", value = rev(mycols1), k=15) %>%
  set("labels_colors", value = rev(mycols1), k=15)%>%
  set("labels_cex", 0.45) 

phylo_plot <- circlize_dendrogram(dend1,
                                  dend_track_height = 0.85) 


ggsave("./phylo_circle.png",
       plot=phylo_plot, 
       width = 15, height = 15, units="cm",
       dpi=600)# do not work for dendograms¿


# check number of species in the tree

length(sce1_tree$tip.label)
sort(sce1_tree$tip.label)

#* test subset tree

sel_species <- c("foeniculum",
                 "llatum",
                 "munis",
                 "noidea")# in case of having uncomplete names

sel_label<-sapply(sel_species,grep,sce1_tree$tip.label)%>%unlist()
pr_species <- sce1_tree$tip.label[sel_label]

pr_tree<-drop.tip(sce1_tree,
                  setdiff(sce1_tree$tip.label,pr_species))
plotTree(pr_tree)

cophenetic.phylo(pr_tree)

col_bin <- ifelse(grepl(sel_species, 
                        sce1_tree$tip.label),
                  "blue","orange")
col_bin2 <- ifelse(sce1_tree$tip.label%in%pr_species, 
                   "red", "black")

plot.phylo(sce1_tree, cex = 0.25,
           main = "scenario.1",
           type="fan", tip.color=col_bin2)

#********************************##
#** get phylogenetic distances ####
#********************************##

# table preparation

#* we should apply this approach to the
#* subset of species with functional data for the selected
#* traits as phylogenetic distance per species is 
#* obtained as the average of pairwise distances between the 
#* target species and the rest of species of the community 
#* weighted by the species cover in the community. So,
#* including a different subset of species will give a 
#* different value of phylogenetic distances between whole species
#* table dataset and complete cases dataset


basic_info2 <- basic_info%>%
  distinct()%>%
  filter(!is.na(Family))%>%
  filter(Family!="Unknown")%>%
  mutate(tpl_family=na_if(tpl_family, ""))%>%
  mutate(family_def=case_when(
    is.na(tpl_family)~Family,
    !is.na(tpl_family)~tpl_family # some tpl families have na, replace by old family name
  ),
  species_def=case_when(
    sp=="sp."~ substring(tpl_species, 1,
                         nchar(tpl_species)-4),
    sp!="sp."~ tpl_species
  ))%>%
  mutate(family_def = as.factor(family_def),
         family_def = forcats::fct_recode(family_def,
                                          `Asteraceae`="Compositae",
                                          `Fabaceae`= "Leguminosae",
                                          `Sapindaceae`="Aceraceae"))%>%
  select(-tpl_family,
         -tpl_species)%>%
  rename(tpl_family=family_def,
         tpl_species=species_def)%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))%>%
  mutate(sp_code=paste0(species_code,"_",
                        plot_code))



complete_traits <- abun_traits%>%
  select(Exp, Year, OldField, 
         Transect, Plot, YearAb,
         plot_code, transect_code,
         species, genus, sp, sub_sp, 
         catford_species, 
         cover_percent, 
         real_cover_percent,
         new_cover_percent,
         old_total_cover,
         new_total_cover,
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
                        leaf_dry_mass))%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))%>%
  mutate(sp_code=paste0(species_code,"_",
                        plot_code))

complete_traits <- complete_traits%>%
  filter(genus!="Quercus")

com_df <- complete_traits%>% 
  left_join(basic_info2%>%
              select(sp_code,
                     setdiff(names(basic_info2),
                             names(complete_traits))),
            by="sp_code")

nrow(com_df)
nrow(complete_traits)
names(com_df)
unique(is.na(com_df$tpl_species))

#####
table_cat <- com_df%>%
  select(-c(Plot, plot_code, 
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

nplot_trans <- com_df%>%
  select(transect_code, Plot)%>%
  distinct()%>%
  group_by(transect_code)%>%
  summarise(n=n())%>%
  ungroup()%>%
  mutate(transect_cover=n*100)

hist((nplot_trans$transect_cover))
unique(nplot_trans$n)
nplot_trans%>%
  filter(n<25)

table_cov <- com_df%>%
  select(transect_code,
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

hist(com_df$cover_percent)
hist(table_cov$mean_cover_percent)
hist(table_cov$cover_percent)
hist(table_cov$transect_cover)

table <- table_cat%>%
  left_join(table_cov%>%
              select(setdiff(names(table_cov),
                             names(table_cat)),
                     sp_code), by="sp_code")

nrow(table)
sort(names(table))

transect_df_list <- list()
transect_list <- unique(table$transect_code)


for(i in 1:length(transect_list)) {
  
  table_plot <- table%>%
    filter(transect_code== transect_list[i])
  
  sp_list <- unique(table_plot$sp_code)
  species_list0 <- unique(table_plot$tpl_species)
  species_list <- gsub(" ", "_", species_list0)
  
  sp_df_list <- list()
  
  if(nrow(table_plot)<2){
    
    for(j in 1:length(sp_list)){
      
      table_plot_sp <- table_plot%>%
        filter(sp_code== sp_list[j])
      
      table_plot_sp <- table_plot_sp%>%
        mutate(phylo_w_dist=NA_real_)
      
      sp_df_list[[j]] <- table_plot_sp
      
    }
    
  }else{
    
    cut_tree<-drop.tip(sce1_tree,
                       setdiff(sce1_tree$tip.label,species_list))
    plotTree(cut_tree)
    
    phylo_dist <- cophenetic.phylo(cut_tree)%>%
      as.data.frame()%>%
      rownames_to_column()
    
    for(j in 1:length(sp_list)){
      
      table_plot_sp <- table_plot%>%
        filter(sp_code== sp_list[j])
      
      table_plot_rest <- table_plot%>%
        filter(sp_code!= sp_list[j])
      
      dist_x <- phylo_dist%>%
        filter(rowname==species_list[j])%>%
        select(-rowname)%>%
        t()%>%
        data.frame()%>%
        rownames_to_column()%>%
        filter(rowname!=species_list[j])%>%
        mutate(tpl_species=gsub("_", " ", rowname))
      
      names(dist_x)[2] <- "dist"
      
      dist_x2 <- table_plot_rest%>%
        left_join(dist_x, by="tpl_species")
      
      dist_x2%>%
        select(transect_code,
               tpl_species,
               dist,
               mean_cover_percent)
      
      phylo_w_dist <- weighted.mean(dist_x2$dist,
                                    dist_x2$mean_cover_percent)# also real_cover_percent
      
      table_plot_sp <- table_plot_sp%>%
        mutate(phylo_w_dist=phylo_w_dist)
      
      sp_df_list[[j]] <- table_plot_sp
      
    } #loop end
    
  }# else end
  
  
  table_plot2 <- bind_rows(sp_df_list)
  transect_df_list[[i]] <- table_plot2
  
  #*add plot for community species
  
  # native_sp <- table_plot2%>%
  #   filter(Origin=="Native")%>%
  #   select(tpl_species)%>%
  #   mutate(tpl_species1=gsub(" ", "_", tpl_species))
  #   
  # invasive_sp <- table_plot2%>%
  #   filter(Origin!="Native")%>%
  #   select(tpl_species)%>%
  #   mutate(tpl_species1=gsub(" ", "_", tpl_species))
  # 
  # 
  # col_bin <- ifelse(grepl(sel_species, sce1_tree$tip.label),"blue","orange")
  # col_bin2 <- ifelse(sce1_tree$tip.label%in%pr_species, "red", "black")
  # 
  # sp_tree <- tibble(
  #   species=sce1_tree$tip.label
  # )%>%
  #   mutate(color=case_when(
  #     species%in%native_sp$tpl_species1~"green",
  #     species%in%invasive_sp$tpl_species1~"red",
  #     species!=native_sp$tpl_species1&species!=invasive_sp$tpl_species1~"black"
  #   ),
  #   cex=case_when(
  #     species%in%native_sp$tpl_species1~0.4,
  #     species%in%invasive_sp$tpl_species1~0.4,
  #     species!=native_sp$tpl_species1&species!=invasive_sp$tpl_species1~0.25
  #   ))
  # 
  # sp_tree%>%
  #   filter(color!="black")
  # 
  # col_bin <- sp_tree$color
  # cex_bin <- sp_tree$cex
  # 
  # 
  # plot.phylo(sce1_tree, cex = cex_bin, main = "scenario.1",
  #            type="fan", tip.color=col_bin)
  # 
  
  
}


phylo_dist_all <- bind_rows(transect_df_list)
nrow(phylo_dist_all)
nrow(table)

data.table::fwrite(phylo_dist_all,
                   "../../output/tables/dissim/phylo_distances_transect_new.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)



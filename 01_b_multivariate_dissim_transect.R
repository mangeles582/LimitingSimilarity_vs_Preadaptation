



##***********************************##
##* Maria Angeles Pérez-Navarro
##* King's College University
##* AlienImpacts project
##* Jan 2022
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
  select(-cover_percent)%>%
  distinct()%>%
  nrow()

#* summarise cover data at species level within plot 
#* in case it may have some 
#* error in the original database. 
#* Actually there are 46 repeated species within plot

abun_traits1 <- abun_traits%>%
  select(-cover_percent)%>%
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
  select(Year, OldField,
         Transect, Plot,
         plot_code, species,
         cover_percent,
         old_total_cover,
         new_total_cover,
         real_cover_percent,
         new_cover_percent)%>%
  arrange(plot_code)

abun_traits%>%
  select(real_cover_percent)%>%
  pull()%>%hist()

abun_traits%>%
  filter(no_plant=="keep")%>%
  select(new_cover_percent)%>%
  pull()%>%hist()

## split complete cases dataset

basic_info <- abun_traits%>%
  select(Exp, Year, OldField, 
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
#   select(species, SLA, seed_mass,
#          leaf_area,LDMC,
#          height,leaf_fresh_mass,
#          leaf_dry_mass)

complete_traits <- abun_traits%>%
  select(Exp, Year, OldField, 
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
  select(-c(Exp, Year, OldField, 
            Transect, Plot, YearAb, 
            plot_code, transect_code,
            cover_percent))%>%
  distinct()%>%
  nrow()

complete_traits%>%
  select(species, SLA, seed_mass,
         leaf_area,LDMC,
         height,leaf_fresh_mass,
         leaf_dry_mass)%>%
  distinct()%>%
  nrow()

sort(unique(complete_traits$species))
nrow(complete_traits)
nrow(abun_traits)

##* 110 from 310 species have data for all the selected traits 30% aprox
##* 101678 rows from 152203 rows of abundances dataset have trait data
##* ie 0.67%

# plot level

res_plot <- abun_traits%>%
  select(plot_code)%>%
  unique()%>%
  nrow()

row_complete_plot <- complete_traits%>%
  group_by(plot_code)%>%
  summarise(complete = n())%>%
  ungroup()

row_total_plot <- abun_traits%>%
  group_by(plot_code)%>%
  summarise(total = n())%>%
  ungroup()

row_total_plant_plot <- abun_traits%>%
  filter(no_plant=="keep")%>%
  group_by(plot_code)%>%
  summarise(total_plant = n())%>%
  ungroup()

nrow_sum_plot0 <- row_total_plot%>%
  left_join(row_total_plant_plot, by="plot_code")%>%
  replace_na(list(complete=0))

nrow_sum_plot <- nrow_sum_plot0%>%
  left_join(row_complete_plot, by="plot_code")%>%
  replace_na(list(complete=0))

nrow_sum_plot <- nrow_sum_plot%>%
  mutate(repres_percent= complete/total_plant,
         level="plot",
         variable="rows")

cover_complete_plot <- complete_traits%>%
  group_by(plot_code)%>%
  summarise(complete = sum(cover_percent))%>%
  replace_na(list(complete=0))

cover_total_plot <- abun_traits%>%
  group_by(plot_code)%>%
  summarise(total = sum(cover_percent))

cover_total_plant_plot <- abun_traits%>%
  filter(no_plant=="keep")%>%
  group_by(plot_code)%>%
  summarise(total_plant = sum(cover_percent))

cover_sum_plot0 <- cover_total_plot%>%
  left_join(cover_total_plant_plot, by="plot_code")%>%
  replace_na(list(complete=0))

cover_sum_plot <- cover_sum_plot0%>%
  left_join(cover_complete_plot, by="plot_code")

cover_sum_plot <- cover_sum_plot%>%
  mutate(repres_percent= complete/total_plant,
         level="plot",
         variable="cover")## check how much abundant are not species items

## transect level


res_transect <- abun_traits%>%
  select(transect_code)%>%
  unique()%>%
  nrow()

row_complete_transect <- complete_traits%>%
  group_by(transect_code, species)%>%
  summarise(cover_percent_sum= sum(cover_percent))%>%
  ungroup()%>%
  group_by(transect_code)%>%
  summarise(complete = n())%>%
  ungroup()

row_total_transect <- abun_traits%>%
  group_by(transect_code, species)%>%
  summarise(cover_percent_sum= sum(cover_percent))%>%
  ungroup()%>%
  group_by(transect_code)%>%
  summarise(total = n())%>%
  ungroup()

row_total_plant_transect <- abun_traits%>%
  filter(no_plant=="keep")%>%
  group_by(transect_code, species)%>%
  summarise(cover_percent_sum= sum(cover_percent))%>%
  ungroup()%>%
  group_by(transect_code)%>%
  summarise(total_plant = n())%>%
  ungroup()

nrow_sum_transect0 <- row_total_transect%>%
  left_join(row_total_plant_transect, by="transect_code")%>%
  replace_na(list(complete=0))

nrow_sum_transect <- nrow_sum_transect0%>%
  left_join(row_complete_transect, by="transect_code")%>%
  replace_na(list(complete=0))

nrow_sum_transect <- nrow_sum_transect%>%
  mutate(repres_percent= complete/total_plant,
         level="transect",
         variable="rows")

cover_complete_transect <- complete_traits%>%
  group_by(transect_code)%>%
  summarise(complete = sum(cover_percent))%>%
  replace_na(list(complete=0))

cover_total_transect <- abun_traits%>%
  group_by(transect_code)%>%
  summarise(total = sum(cover_percent))

cover_total_plant_transect <- abun_traits%>%
  filter(no_plant=="keep")%>%
  group_by(transect_code)%>%
  summarise(total_plant = sum(cover_percent))

cover_sum_transect0 <- cover_total_transect%>%
  left_join(cover_total_plant_transect, by="transect_code")

cover_sum_transect <- cover_sum_transect0%>%
  left_join(cover_complete_transect, by="transect_code")

cover_sum_transect <- cover_sum_transect%>%
  mutate(repres_percent= complete/total_plant,
         level="transect",
         variable="cover")

## add hist plots

list_sum <- list(nrow_sum_plot, cover_sum_plot,
                 nrow_sum_transect, cover_sum_transect)

sum_all <- bind_rows(list_sum)

hist(nrow_sum_plot$total_plant)
hist(nrow_sum_transect$total_plant)

row_sum_all <- sum_all%>%
  filter(variable=="rows")%>%
  select(plot_code, total_plant, complete, level)%>%
  pivot_longer(!c(plot_code, level),
               names_to="dataset", 
               values_to="n")

(total_rows <- ggplot(data=row_sum_all,
                      aes(x=n, fill=dataset))+
    geom_histogram(alpha=0.75, position = 'identity')+
    scale_fill_manual(values=c("#009966",
                               "#CCFF99"))+
    facet_wrap(.~level, scales="free")+
    ylab ( "count") +
    xlab ( "number of rows") +
    theme_bw()+
    theme(#legend.position="none",
      legend.background = element_blank(),
      legend.title = element_text(size = 11), 
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      axis.text=element_text(color="black",
                             size=10),
      axis.title=element_text(color="black",
                              size=13),
      strip.text = element_text(color="black",
                                size=15,  
                                family="serif",
                                hjust=0.5),
      text=element_text(family="serif"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())
)

(total_all <- ggplot(data=sum_all,
                     aes(x=repres_percent))+
    geom_histogram(fill="#009966")+
    facet_grid(level~variable, scales="free")+
    geom_vline(xintercept=0.5, linetype="dashed", 
               color = "red", size=0.4)+
    geom_vline(xintercept=0.75, linetype="dashed", 
               color = "orange", size=0.4)+
    ylab ( "count") +
    xlab ( "percentage of species with trait data") +
    theme_bw()+
    theme(#legend.position="none",
      legend.background = element_blank(),
      legend.title = element_text(size = 11), 
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.5, "lines"),
      axis.text=element_text(color="black",
                             size=10),
      axis.title=element_text(color="black",
                              size=13),
      strip.text = element_text(color="black",
                                size=15,  
                                family="serif",
                                hjust=0.5),
      text=element_text(family="serif"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank())
)


#* estimate number of plots with less than 2 species
#* these data will be lost as it is imposible to calculate any distance

row_sum_all%>%
  filter(dataset=="complete")%>%
  filter(n<=1)%>%
  group_by(level)%>%
  summarise(nrow=n())


#* save plot cover and nrow represntativeness

#* new cover percent shouldn't be use as it is not accurate in terms of
#* cover neither in terms of frequency -as there more unseampled species-
#* cover percent is the percentage in relation to the total sampled surface
#* so, it helps to compare in terms of plant sizes
#* and real cover percent is the percentage in relation to the total items surveyed
#* it give information about the frequencies of species regarding other
#* items in the plot. It will be useful if we are interested on relative 
#* abundances and want to scale up to other transect or field level keeping the 
#* focus on relative abundances


nrow_sum_plot
cover_sum_plot

# data.table::fwrite(nrow_sum_plot,
#                    "../../output/tables/abundances/nrow_representativeness.csv",
#                    sep=";", dec=".", col.names=T, row.names=F, quote=F)
# 
# 
# data.table::fwrite(cover_sum_plot,
#                    "../../output/tables/abundances/cover_representativeness.csv",
#                    sep=";", dec=".", col.names=T, row.names=F, quote=F)
# 


#* Regarding to the percentage of cover with trait data per plot
#* we see that there's a wide variability and that a huge amount of plot
#* have trait data for 50% of plant cover or more
#* In case of considering species representation per plot the percentage of 
#* representativeness is considerably higher. 
#* In addition it is not necessary to work at transect level as the percentage
#* cover or species with data do not increase.

#* add accumulation density curve of species lacking complete
#* trait data for the selected functional traits

# abun_traits <- abun_traits%>%
#   filter(no_plant=="keep")

repres_species <- unique(complete_traits$species_code)

unrepres_species <- abun_traits%>%
  filter(no_plant=="keep")%>%
  filter(!species_code%in%repres_species)%>%
  filter()

unique(unrepres_species$species)

unrepres_list <- unrepres_species$species

sp_df <- unrepres_species%>%
  filter(species==unrepres_list[104])%>%
  arrange(new_cover_percent)

hist(sp_df$new_cover_percent)
plot(sort(sp_df$new_cover_percent))

acum_df <- sp_df%>%
  select(plot_code, new_cover_percent)%>%
  mutate(new_cover_percent1=c(new_cover_percent[-1],NA)
  )%>%
  mutate(increase=new_cover_percent1-new_cover_percent)

ggplot(data=acum_df, aes(x=reorder(plot_code, new_cover_percent),
                         y=new_cover_percent, 
                         group=1))+# if not geom_line does not work
  geom_line()+
  geom_point()+
  ggtitle(unrepres_list[104])+
  theme_classic()

#* consider to run a loop with this code to look for more abundant
#* species still lacking trait information to increase the representativeness

#*******************************************************##
#** prepare trait pca with the selected set of traits ####
#*******************************************************##

# abun_traits <- abun_traits%>%
#   filter(no_plant=="keep")

# prepare species data for pca

complete_traits <- complete_traits%>%
  filter(genus!="Quercus")

complete_traits_sp <- complete_traits%>%
  select(species, SLA, seed_mass,
         leaf_area,LDMC,
         height,leaf_fresh_mass,
         leaf_dry_mass)%>%
  distinct()%>%
  tibble::column_to_rownames('species')#only complete cases


# df_traits_sp <- df_traits%>%
#   distinct()%>%
#   tibble::column_to_rownames('species')%>%
#   filter(!(is.na(SLA)&is.na(seed_mass)&
#              is.na(leaf_area)&is.na(LDMC)&
#              is.na(height)&is.na(leaf_fresh_mass)&
#              is.na(leaf_dry_mass)))#it has na values for some traits


# pca cannot contain na entries

dudi_trait <- dudi.pca(complete_traits_sp, 
                       scale=T, scannf = F, nf=5)

summary(dudi_trait)
screeplot(dudi_trait)
dudi_trait$li
dudi_trait$co
s.corcircle(dudi_trait$co[, c(1,2)])# if you don't extract specific axis it shows the two first
par(mfrow=c(2,2))
factoextra::fviz_pca_var(dudi_trait, col.var = "contrib", repel=T,
                         gradient.cols=c("#00AFBB", "#E7B800", "#FF0000"))

factoextra::fviz_pca_var(dudi_trait, axes=c(1,3),
                         col.var = "contrib", repel=T,
                         gradient.cols=c("#00AFBB", "#E7B800", "#FF0000"))

factoextra::fviz_pca_var(dudi_trait, axes=c(2,3),
                         col.var = "contrib", repel=T,
                         gradient.cols=c("#00AFBB", "#E7B800", "#FF0000"))




#save(dudi_trait, file = "../../output/dudi_trait.RData")
load("../../output/dudi_trait.RData")

## join abundance table and trait axis values

dudi_df <- dudi_trait$li%>%
  rownames_to_column(var="species")%>%
  as_tibble()

join_pca_traits <- complete_traits%>%
  left_join(dudi_df, by="species")

nrow(join_pca_traits)
nrow(complete_traits)

#*************************##
#** get distances in 2d ####
#*************************##

#**** using community trait centroid ####

source('niche_required_functions.R')
join_pca_traits

# table preparation

join_pca_traits <- join_pca_traits%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))

## create table at transect level

names(join_pca_traits)

pca_traits_expand <- join_pca_traits%>%
  group_by(transect_code)%>%
  expand(species_code)

join_pca_traits2 <- pca_traits_expand%>%
  left_join(join_pca_traits%>%
              select(transect_code,
                     Plot, plot_code),
            by=c("transect_code"))%>%
  distinct()%>%
  mutate(sp_code=paste0(species_code,"_", plot_code))

cover_mean <- join_pca_traits2%>%
  left_join(join_pca_traits%>%
              select(sp_code,
                     cover_percent,
                     species),
            by="sp_code")%>%
  mutate(cover_percent = replace_na(cover_percent, 0))%>%
  group_by(transect_code, species_code)%>%
  summarise(mean_cover=mean(cover_percent))#dis can be used alterntively to table_cov

table_cat <- join_pca_traits%>%
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

nplot_trans <- join_pca_traits%>%
  select(transect_code, Plot)%>%
  distinct()%>%
  group_by(transect_code)%>%
  summarise(n=n())%>%
  ungroup()%>%
  mutate(transect_cover=n*100)#do not multiply by 100

hist((nplot_trans$transect_cover))
unique(nplot_trans$n)
nplot_trans%>%
  filter(n<25)

table_cov <- join_pca_traits%>%
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

hist(cover_mean$mean_cover)
hist(table_cov1$mean_cover_percent)
cor(cover_mean$mean_cover, table_cov1$mean_cover_percent)

cover_mean1 <- cover_mean%>%
  left_join(table_cov1%>%select(transect_code, 
                                species_code,
                                mean_cover_percent,
                                cover_percent),
            by=c("transect_code", "species_code"))

nrow(table_cov)
nrow(table_cov1)
table_cov <- table_cov1

hist(join_pca_traits$cover_percent)
hist(table_cov$mean_cover_percent)
hist(table_cov$cover_percent)
hist(table_cov$transect_cover)
 
table <- table_cat%>%
  left_join(table_cov%>%
              select(setdiff(names(table_cov),
                             names(table_cat)),
                     sp_code), by="sp_code")

nrow(table)

basic_vars <- names(basic_info)[!names(basic_info)%in%c("Plot",
                                                        "plot_code")]

centroid_df <-  table%>%
  dplyr::select(basic_vars,
                species_code,
                sp_code,
                Axis1, Axis2,
                Axis3, Axis4, 
                Axis5)%>%
  distinct()


transect_list <- unique(centroid_df$transect_code)

for(i in 1:length(transect_list)) {
  
  table_transect <- table%>%
    filter(transect_code== transect_list[i])
  
  species_list <- unique(table_transect$species)
  
  for(j in 1:length(species_list)){
    
    table_transect_rest <- table_transect%>%
      filter(species!= species_list[j])
    
    # species_id <- unique(table_transect_sp$species_code)
    
    species_id <- table_transect%>%
      filter(species== species_list[j])%>%
      select(sp_code)%>%as.character()
    
    trait_centroid <- COGravity(x=table_transect_rest$Axis1, 
                                y=table_transect_rest$Axis2, 
                                wt=table_transect_rest$mean_cover_percent)#same to real_cover_percent
    
    centroid_df[centroid_df$sp_code==species_id, "trait_centroid_x"] <- trait_centroid[1]
    centroid_df[centroid_df$sp_code==species_id, "trait_centroid_y"] <- trait_centroid[3]
    
    centroid_df_transect <- centroid_df%>%
      filter(sp_code== species_id)
    dissimilarity <- GMINdistance(centroid_df_transect[, c("trait_centroid_x",
                                                           "trait_centroid_y")],
                                  centroid_df_transect[, c("Axis1",
                                                           "Axis2")])[1]
    centroid_df[centroid_df$sp_code==species_id, "dissim"]<- dissimilarity
    
    
    # table_transect$sp_code <- factor(table_transect$sp_code)
    # require(RColorBrewer)
    # cols <- colorRampPalette(brewer.pal(11, "Spectral"))
    # display.brewer.pal(n=11, name="Spectral")
    # cols1 <- brewer.pal(n=11, name="Spectral")
    # scales::show_col(cols1)
    # mypal <- cols1[c(4, 8, 10)]
    
    # gg_i <- ggplot(table_transect %>%
    #                  dplyr::mutate(cover_percent = case_when(
    #                    cover_percent == 0~NA_real_,
    #                    cover_percent != 0~cover_percent))) +
    #   geom_point( aes(x=Axis1, y=Axis2,
    #                   color=Origin, size=cover_percent*10
    #   ), alpha =.95)+
    #   scale_color_manual(values = mypal)+
    #  #scale_color_viridis_d(option="plasma")+
    #   labs(color = "Species")+
    #   guides(color=guide_legend(ncol=2)) +
    #   xlim(min(table_transect$Axis1)-2,
    #        max(table_transect$Axis1)+2)+
    #   ylim(min(table_transect$Axis2)-2,
    #        max(table_transect$Axis2)+2)+
    #   scale_size(range = c(2, 7))+
    #   scale_alpha_continuous("Cover_percentage",
    #                          range=c(0.2, 0.9))+
    #   guides(size = FALSE,
    #          alpha = FALSE
    #   )+
    #   geom_point(data=centroid_df_transect,
    #              aes(x=trait_centroid_x,
    #                  y=trait_centroid_y),
    #              color="black", size=2, shape=19)+
    #   # geom_point(data=centroid_df_transect,
    #   #            aes(x=Axis1, y=Axis2),
    #   #            color=col_sp, size=2, shape=17)+
    #   xlab("Functional Axis 1 (42.6%)")+
    #   ylab("Functional Axis 2 (17.1%)")+
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
    
    
    # ggsave(paste0("../../output/community_diagrams/", transect_list[i], ".png"), 
    #        plot=gg_i, width = 9, height = 9, dpi = 300, units="cm")
    # ggsave(paste0("../../output/community_diagrams/", transect_list[i], ".pdf"),
    #        plot=gg_i, width = 9, height = 9, units="cm")
    
  }
  
  
  
}

nrow(centroid_df)
data.table::fwrite(centroid_df,
                   "../../output/tables/dissim/dissimilarity_centroid2d_transect_new.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)



#**** using weighted pairwaise distances ####


source('niche_required_functions.R')
join_pca_traits

# table preparation

join_pca_traits <- join_pca_traits%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))

## create table at transect level

names(join_pca_traits)

table_cat <- join_pca_traits%>%
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

nplot_trans <- join_pca_traits%>%
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

table_cov <- join_pca_traits%>%
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

hist(join_pca_traits$cover_percent)
hist(table_cov$mean_cover_percent)
hist(table_cov$cover_percent)
hist(table_cov$transect_cover)

table <- table_cat%>%
  left_join(table_cov%>%
              select(setdiff(names(table_cov),
                             names(table_cat)),
                     sp_code), by="sp_code")

nrow(table)

basic_vars <- names(basic_info)[!names(basic_info)%in%c("Plot",
                                                        "plot_code")]

centroid_df <-  table%>%
  dplyr::select(basic_vars,
                species_code,
                sp_code,
                Axis1, Axis2,
                Axis3, Axis4, 
                Axis5)%>%
  distinct()


transect_list <- unique(centroid_df$transect_code)

for(i in 1:length(transect_list)) {
  
  table_transect <- table%>%
    filter(transect_code== transect_list[i])
  
  species_list <- unique(table_transect$sp_code)
  
  for(j in 1:length(species_list)){
    
    table_transect_sp <- table_transect%>%
      filter(sp_code== species_list[j])
    
    table_transect_rest <- table_transect%>%
      filter(sp_code!= species_list[j])
    
    if(nrow(table_transect_rest)<1){
      centroid_df[centroid_df$sp_code==species_list[j], "dissim"] <- NA
    }else{
      table_transect_rest <- table_transect_rest%>%
        mutate(dissim=GMINdistance(table_transect_sp[,c("Axis1", 
                                                        "Axis2")],
                                   table_transect_rest[,c("Axis1",
                                                          "Axis2")])$min.dist)
      
      av_dissim <- weighted.mean(table_transect_rest$dissim,
                                 table_transect_rest$mean_cover_percent)#also real_cover_percent
      
      centroid_df[centroid_df$sp_code==species_list[j], "dissim"] <- av_dissim
    }
    
  }
 
  print(i) 
  
}


nrow(centroid_df)
data.table::fwrite(centroid_df,
                   "../../output/tables/dissim/dissimilarity_weight_distances2d_transect_new.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)


#*************************##
#** get distances in 3d ####
#*************************##

#**** using community trait centroid ####

source('niche_required_functions.R')
join_pca_traits

# table preparation

join_pca_traits <- join_pca_traits%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))

## create table at transect level

names(join_pca_traits)

table_cat <- join_pca_traits%>%
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

nplot_trans <- join_pca_traits%>%
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

table_cov <- join_pca_traits%>%
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

hist(join_pca_traits$cover_percent)
hist(table_cov$mean_cover_percent)
hist(table_cov$cover_percent)
hist(table_cov$transect_cover)

table <- table_cat%>%
  left_join(table_cov%>%
              select(setdiff(names(table_cov),
                             names(table_cat)),
                     sp_code), by="sp_code")

nrow(table)

basic_vars <- names(basic_info)[!names(basic_info)%in%c("Plot",
                                                        "plot_code")]

centroid_df <-  table%>%
  dplyr::select(basic_vars,
                species_code,
                sp_code,
                Axis1, Axis2,
                Axis3, Axis4, 
                Axis5)%>%
  distinct()


transect_list <- unique(centroid_df$transect_code)

for(i in 1:length(transect_list)) {
  
  table_transect <- table%>%
    filter(transect_code== transect_list[i])
  
  species_list <- unique(table_transect$species)
  
  for(j in 1:length(species_list)){
    
    table_transect_sp <- table_transect%>%
      filter(species!= species_list[j])
    
    # species_id <- unique(table_transect_sp$species_code)
    
    species_id <- table_transect%>%
      filter(species== species_list[j])%>%
      select(sp_code)%>%as.character()
    
    trait_centroid <- COGravity(x=table_transect_sp$Axis1, 
                                y=table_transect_sp$Axis2, 
                                z=table_transect_sp$Axis3,
                                wt=table_transect_sp$mean_cover_percent)#cover_percent
    
    centroid_df[centroid_df$sp_code==species_id, "trait_centroid_x"] <- trait_centroid[1]
    centroid_df[centroid_df$sp_code==species_id, "trait_centroid_y"] <- trait_centroid[3]
    centroid_df[centroid_df$sp_code==species_id, "trait_centroid_z"] <- trait_centroid[5]
    
    centroid_df_transect <- centroid_df%>%
      filter(sp_code== species_id)
    dissimilarity <- GMINdistance(centroid_df_transect[,c("trait_centroid_x",
                                                          "trait_centroid_y",
                                                          "trait_centroid_z")],
                                  centroid_df_transect[,c("Axis1",
                                                          "Axis2",
                                                          "Axis3")]
    )[1]
    
    centroid_df[centroid_df$sp_code==species_id, "dissim"]<- dissimilarity
    
    # require("scatterplot3d")
    # require(plotly)
    # 
    # table_transect$sp_code <- factor(table_transect$sp_code)
    # require(RColorBrewer)
    # cols <- colorRampPalette(brewer.pal(11, "Spectral"))
    # display.brewer.pal(n=11, name="Spectral")
    # cols1 <- brewer.pal(n=11, name="Spectral")
    # mypal <- cols1[c(4, 8, 10)]
    # 
    # transect3d <- table_transect %>%
    #                    dplyr::mutate(
    #                      cover_percent = case_when(
    #                      cover_percent == 0~NA_real_,
    #                      cover_percent != 0~cover_percent),
    #                      origin2=case_when(
    #                        Origin=="Native"~"Native",
    #                        Origin=="Introduced"~"Introduced",
    #                        Origin!="Native"&
    #                          Origin!="Introduced"~"Introduced")
    #                      )
    # 
    # unique(transect3d$origin2)
    # transect3d$origin2 <- as.factor(transect3d$origin2)
    # transect3d <- transect3d%>%
    #   dplyr::select(species,
    #                 origin2,
    #                 Axis1,Axis2,Axis3,
    #                 cover_percent)%>%
    #   mutate(type=case_when(
    #     species==species_list[j]~"target sp",
    #     species!=species_list[j]~"rest sp"
    #   ))
    # 
    # # rbind with trait_centroid
    # 
    # centroid_plot <- trait_centroid%>%
    #   as.data.frame()%>%t()%>%
    #   as.data.frame()%>%
    #   dplyr::select(COGx, COGy, COGz)%>%
    #   rename(Axis1=COGx,
    #          Axis2=COGy,
    #          Axis3=COGz)%>%
    #   mutate(species=paste0("centroid_", species_list[j]),
    #          origin2=NA_real_,
    #          cover_percent=150,
    #          type="centroid")%>%
    #   dplyr::select(species,
    #                 origin2,
    #                 Axis1,Axis2,Axis3,
    #                 cover_percent,
    #                 type)%>%
    #   remove_rownames()
    #   
    # 
    # transect3d <- rbind(transect3d, centroid_plot)
    # transect3d <- transect3d%>%
    #   mutate(scale=scales::rescale(transect3d$cover_percent, c(3, 30)))
    # 
    # plot_ly(data=transect3d,
    #         x=~Axis1,
    #         y=~Axis2,
    #         z=~Axis3,
    #         type="scatter3d",
    #         mode="markers",
    #         color=~type,
    #         # size=~scale,
    #         marker=list(size=~scale))
    # 
    # 
    # colors <- c("#FDAE61", "#3288BD")
    # colors <- colors[as.numeric(transect3d$origin2)]
    # 
    # png(paste0("../../output/community_diagrams/", transect_list[i], 
    #            res=2000, pointsize= 0.5, units="mm", "_3d.png"))
    # 
    # 
    # scatterplot3d(transect3d[,c("Axis1", "Axis2", "Axis3")],
    #               main="3D multivariate functional distances",
    #               xlab = "Funct Axis1 (45.4%)",
    #               ylab = "",
    #               zlab = "Funct Axis2 (19.8%)",
    #               cex.symbol=scales::rescale(transect3d$cover_percent, c(1, 10)),
    #               pch=16,
    #               color=colors,
    #               alpha=0.7)
    # 
    # text(x = 4.7, y = 0.7, "Funct Axis3 (14.4%)", srt = 29)
    # 
    # legend("bottom", legend = levels(transect3d$origin2),
    #        col =  c("#FDAE61", "#3288BD"), 
    #        pch = 16, 
    #        inset = -0.3, xpd = TRUE, horiz = TRUE)
    # 
    # dev.off()
    # 
    
    
  }
  
  print(i)
  
}

nrow(centroid_df)

data.table::fwrite(centroid_df,
                   "../../output/tables/dissim/dissimilarity_centroid3d_transect_new.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)



#**** using weighted pairwaise distances ####


source('niche_required_functions.R')
join_pca_traits

# table preparation

join_pca_traits <- join_pca_traits%>%
  mutate(species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)))

## create table at transect level

names(join_pca_traits)

table_cat <- join_pca_traits%>%
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

nplot_trans <- join_pca_traits%>%
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

table_cov <- join_pca_traits%>%
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

hist(join_pca_traits$cover_percent)
hist(table_cov$mean_cover_percent)
hist(table_cov$cover_percent)
hist(table_cov$transect_cover)

table <- table_cat%>%
  left_join(table_cov%>%
              select(setdiff(names(table_cov),
                             names(table_cat)),
                     sp_code), by="sp_code")

nrow(table)

basic_vars <- names(basic_info)[!names(basic_info)%in%c("Plot",
                                                        "plot_code")]

centroid_df <-  table%>%
  dplyr::select(basic_vars,
                species_code,
                sp_code,
                Axis1, Axis2,
                Axis3, Axis4, 
                Axis5)%>%
  distinct()


transect_list <- unique(centroid_df$transect_code)


for(i in 1:length(transect_list)) {
  
  table_transect <- table%>%
    filter(transect_code== transect_list[i])
  
  species_list <- unique(table_transect$sp_code)
  
  for(j in 1:length(species_list)){
    
    table_transect_sp <- table_transect%>%
      filter(sp_code== species_list[j])
    
    table_transect_rest <- table_transect%>%
      filter(sp_code!= species_list[j])
    
    if(nrow(table_transect_rest)<1){
      centroid_df[centroid_df$sp_code==species_list[j], "dissim"] <- NA
    }else{
      table_transect_rest <- table_transect_rest%>%
        mutate(dissim=GMINdistance(table_transect_sp[,c("Axis1", 
                                                        "Axis2",
                                                        "Axis3")],
                                   table_transect_rest[,c("Axis1",
                                                          "Axis2",
                                                          "Axis3")])$min.dist)
      
      av_dissim <- weighted.mean(table_transect_rest$dissim,
                                 table_transect_rest$mean_cover_percent)#real_cover_percent
      
      centroid_df[centroid_df$sp_code==species_list[j], "dissim"] <- av_dissim
    }
    
  }

  print(i)  
  
}

nrow(centroid_df)

data.table::fwrite(centroid_df,
                   "../../output/tables/dissim/dissimilarity_weight_distances3d_transect.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)




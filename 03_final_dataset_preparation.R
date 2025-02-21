


##***********************************##
##* Maria Angeles Pérez-Navarro
##* King's College University
##* AlienImpacts project
##* Statistic Analyses
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
library(lme4)
library(lmerTest)
library(DHARMa)
#devtools::install_github("EnquistLab/RNSR", ref="master", subdir="NSR") 
library(NSR)
library(ggpointdensity)
library(HH)
library(performance)


#**********************************###
#***** 1 prepare data tables******####
#**********************************###


#**********************************###
#******1.1.plot level tables******####
#**********************************###


coord <- read_delim( "../../data/abundances/location_matrix_byplot.csv", 
                    delim=",", col_names=T)

unique(coord$transect)

coord <- coord%>%
  dplyr::select(-c(xloc, yloc))%>%
  rename(Field=field,
         Transect_old=transect,
         Plot=plot)%>% 
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D"))

field_info <- read_delim( "../../data/abundances/field_information.csv", 
                     delim=",", col_names=T)

unique(field_info$Transect)

field_info <- field_info%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D"))

field_info <- coord%>%
  left_join(field_info, by=c("Field", "Transect"))

field_info <- field_info%>%
  rename(OldField=Field)

phylo_dis <- read_delim( "../../output/tables/dissim/phylo_distances_plot.csv", 
                          delim=";", col_names=T)

# fun_dis <- read_delim( "../../output/tables/dissim/dissimilarity_centroid3d_plot.csv", 
#                           delim=";", col_names=T)

fun_dis <- read_delim( "../../output/tables/dissim/dissimilarity_weight_distances2d_plot.csv",
                          delim=";", col_names=T)

setdiff(fun_dis$OldField, field_info$OldField)

fun_dis <- fun_dis%>% 
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
  ))%>%
  dplyr::select(-c(species_code,
          sp_code, genus, 
          sp, sub_sp,
          plot_code))%>%
  mutate(genus = word(species, 1, sep = fixed(" ")),
         sp = word(species, 2, sep = fixed(" ")),
         sub_sp= word(species, 3, sep = fixed(" ")))%>%#temporal
  mutate(plot_code=paste0(Year,"_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year,"_",
                          OldField, "_",
                          Transect),
         transect_code2=paste0(OldField, "_",
                               Transect),
         species_code= paste0(substring(genus, 1, 3), "_",
                                        substring(sp, 1,3)),
         sp_code=paste0(species_code,"_",
                                 plot_code))

fun_dis <- fun_dis%>%
  left_join(field_info, 
            by=c("OldField", "Transect", "Plot"))


phylo_dis <- phylo_dis%>%
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
  ))%>%
  mutate(plot_code=paste0(Year,"_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year,"_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                           Transect),
         sp_code=paste0(species_code,"_",
                        plot_code))

nrow(fun_dis)
nrow(phylo_dis)
unique(is.na(phylo_dis$OldField))
unique(is.na(fun_dis$OldField))

sort(unique(phylo_dis$OldField))
sort(unique(fun_dis$OldField))


#*consider substitute by dissimilarity_centroid3d
#*dissimilarity_weight_distances2d
#*dissimilarity_weight_distances3d

repres <- read_delim( "../../output/tables/abundances/cover_representativeness.csv", 
                      delim=";", col_names=T)

repres <- repres%>%
  filter(level=="plot")

#*consider substitute by nrow_representativeness

yearab_df <- read_delim( "../../data/abundances/e014_Year of abandon.txt")

burn_df <- read_delim( "../../data/abundances/e014_Burned info.txt")

burn_df <- burn_df%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~ "A",
    Transect_old=="R"~ "B",
    Transect_old=="W"~ "C",
    Transect_old=="Y"~ "D",
    Transect_old=="G "~"A",
    Transect_old=="R "~"B",
    Transect_old=="W "~"C",
    Transect_old=="Y "~"D",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ))

#***** add other field environmental data ********

carbon <-read_delim(file="../../data/abundances/e014_Soil carbon data.txt")
nitro <-read_delim(file="../../data/abundances/e014_Soil nitrogen data.txt")
omatter <-read_delim(file="../../data/abundances/e014_Soil organic matter data.txt")
light <-read_delim(file="../../data/abundances/e014_Percent light penetration data.txt")

names(carbon)
names(nitro)
names(omatter)
names(light)

carbon <- carbon%>%
  rename(OldField=field,
         Transect_old=transect,
         depth_c_cm="Depth (cm)",
         Plot= quadrat, 
         YearAb_carbon=CultivatedHistorical,
         YearAb_photo_carbon=CultivatedPhoto)%>%
  mutate(Transect_old=as.character(Transect_old))%>%
  mutate(Transect = case_when(
    Transect_old=="1"~"A",
    Transect_old=="2"~"B",
    Transect_old=="3"~"C",
    Transect_old=="4"~"D",
    Transect_old=="5"~"E",
    Transect_old=="6"~"F"
  ))%>%
  mutate(plot_code=paste0(Year, "_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year, "_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect))

nitro <- nitro%>%
  rename(OldField=field,
         Transect_old=transect,
         depth_n_cm="Depth (cm)",
         Plot= quadrat)%>%
  mutate(Transect_old=as.character(Transect_old))%>%
  mutate(Transect = case_when(
    Transect_old=="1"~"A",
    Transect_old=="2"~"B",
    Transect_old=="3"~"C",
    Transect_old=="4"~"D",
    Transect_old=="5"~"E",
    Transect_old=="6"~"F"
  ))%>%
  mutate(plot_code=paste0(Year, "_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year, "_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect))%>%
  dplyr::select(-c(CultivatedPhoto,
            CultivatedHistorical,
            LastCrop,
            Area_ha))

omatter <- omatter%>%
  rename(Experiment="Experiment number",
         OldField="Old field designator",
         Sampling_date_om="Sampling date",
         Plot=Quadrat)%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ))%>%
  mutate(Year=paste0("19", 
                     substr(as.character(Sampling_date_om), start=1,stop=2)),
    plot_code=paste0(Year, "_",
                     OldField, "_",
                     Transect, "_",
                     Plot),
    plot_code2=paste0(OldField, "_",
                      Transect, "_",
                      Plot),
    transect_code=paste0(Year, "_",
                         OldField, "_",
                         Transect),
    transect_code2=paste0(OldField, "_",
                          Transect))%>%
  mutate(Year=as.numeric(Year))%>%
  dplyr::select(-c(CultivatedHistorical,
            LastCrop,
            Area_ha,
            Experiment))

light <- light%>%
  rename(Experiment="Experiment number",
         OldField="Old field code",
         Transect_old="Transect code",
         Sampling_date_light="Sampling date",
         Year="Sampling year",
         Plot=Quadrat,
         YearAb_light=CultivatedHistorical)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ))%>%
  mutate(plot_code=paste0(Year, "_",
                             OldField, "_",
                             Transect, "_",
                             Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year, "_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect))%>%
  dplyr::select(-c(LastCrop,
            Area_ha,
            Experiment))

head(carbon)
head(nitro)
head(omatter)
head(light)

unique(fun_dis$Year)
unique(carbon$Year)
unique(nitro$Year)
unique(omatter$Year)
unique(light$Year)

sort(unique(fun_dis$OldField))
sort(unique(carbon$OldField))
sort(unique(nitro$OldField))
sort(unique(omatter$OldField))
sort(unique(light$OldField))

sort(unique(fun_dis$Transect))
sort(unique(carbon$Transect))
sort(unique(nitro$Transect))
sort(unique(omatter$Transect))
sort(unique(light$Transect))

sort(unique(fun_dis$transect_code2))
sort(unique(carbon$transect_code2))
sort(unique(nitro$transect_code2))
sort(unique(omatter$transect_code2))
sort(unique(light$transect_code2))

setdiff(unique(carbon$transect_code2), unique(fun_dis$transect_code2))
setdiff(unique(nitro$transect_code2), unique(fun_dis$transect_code2))
setdiff(unique(omatter$transect_code2), unique(fun_dis$transect_code2))
setdiff(unique(light$transect_code2), unique(fun_dis$transect_code2))

sort(unique(fun_dis$Plot))
sort(unique(carbon$Plot))
sort(unique(nitro$Plot))
sort(unique(omatter$Plot))
sort(unique(light$Plot))

names(carbon)
names(nitro)
names(omatter)
names(light)

local_env <- list(carbon, nitro, omatter, light)%>%
  purrr::reduce(full_join, by=c("Year","OldField",
                                "Transect","Plot",
                                "plot_code", "plot_code2",
                                "transect_code",
                                "transect_code2"))%>%
  distinct()

z <- local_env%>%
  dplyr::select(plot_code2, YearAb_carbon,
         YearAb_light)


#* as there is no data for all plot and year combinations we need
#* to average them. Check if it is possible at plot level or 
#* should be at transect level

fun_dis%>%
  dplyr::select(OldField,
         Transect,
         Plot)%>%
  distinct()

nitro%>%
  dplyr::select(OldField,
         Transect,
         Plot)%>%
  distinct()


setdiff(fun_dis$plot_code2, nitro$plot_code2)
setdiff(fun_dis$OldField, nitro$OldField)

sort(unique(fun_dis$OldField))
sort(unique(fun_dis$Transect))

sort(unique(nitro$OldField))
sort(unique(nitro$Transect))

sort(unique(omatter$OldField))
sort(unique(omatter$Transect))

sort(unique(light$OldField))
sort(unique(light$Transect))

## add carbon content

carbon_plot_av <- carbon%>%
  group_by(plot_code2)%>%##get plot average across years
  summarise(avg_plot_c_perc=mean(ProportionC, na.rm=T),
            sd_plot_c_perc=sd(ProportionC, na.rm=T))

carbon_plot_av <- carbon%>%
  dplyr::select(OldField,
         Transect,
         transect_code,
         plot_code,
         transect_code2, 
         plot_code2)%>%
  distinct()%>%
  left_join(carbon_plot_av, by="plot_code2")

carbon_transect_av <- carbon_plot_av%>%#avoid to give more weight to plots sampled more years
  group_by(transect_code2)%>%
  summarise(avg_transect_c_perc=mean(avg_plot_c_perc, na.rm=T),
            sd_transect_c_perc=sd(avg_plot_c_perc, na.rm=T))

carbon_field_av <- carbon_plot_av%>%
  group_by(OldField)%>%
  summarise(avg_field_c_perc=mean(avg_plot_c_perc, na.rm=T),
            sd_field_c_perc=sd(avg_plot_c_perc, na.rm=T))#sd is 0.3 of the mean

carbon_plot_av%>%filter(grepl("_NA", plot_code2))

fun_dis_env0 <- fun_dis%>%
  left_join(carbon_plot_av%>%
              dplyr::select(plot_code2, 
                     setdiff(names(carbon_plot_av),
                             names(fun_dis))),
            by="plot_code2")

fun_dis_env1 <- fun_dis_env0%>%
  left_join(carbon_transect_av, by="transect_code2")

fun_dis_env <- fun_dis_env1%>%
  left_join(carbon_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(carbon_perc=case_when(
    !is.na(avg_plot_c_perc)~avg_plot_c_perc,
    (is.na(avg_plot_c_perc)&
       !is.na(avg_transect_c_perc))~avg_transect_c_perc,
    (is.na(avg_plot_c_perc)&
       is.na(avg_transect_c_perc))~avg_field_c_perc
  ))

rm(fun_dis_env0, fun_dis_env1)

nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)
fun_dis_env%>%filter(grepl("_NA", plot_code2))


## add nitrogen content

nitro_plot_av <- nitro%>%
  group_by(plot_code2)%>%
  summarise(avg_plot_n_perc=mean(ProportionN, na.rm=T),
            sd_plot_n_perc=sd(ProportionN, na.rm=T))

nitro_plot_av <- nitro%>%
  dplyr::select(OldField,
         Transect,
         transect_code,
         plot_code,
         transect_code2, 
         plot_code2)%>%
  distinct()%>%
  left_join(nitro_plot_av, by="plot_code2")

nitro_transect_av <- nitro_plot_av%>%
  group_by(transect_code2)%>%
  summarise(avg_transect_n_perc=mean(avg_plot_n_perc, na.rm=T),
            sd_transect_n_perc=sd(avg_plot_n_perc, na.rm=T))
  
nitro_field_av <- nitro_plot_av%>%
  group_by(OldField)%>%
  summarise(avg_field_n_perc=mean(avg_plot_n_perc, na.rm=T),
            sd_field_n_perc=sd(avg_plot_n_perc, na.rm=T))#sd is 0.14 of the mean


fun_dis_env0 <- fun_dis_env%>%
  left_join(nitro_plot_av%>%
              dplyr::select(plot_code2,
                     setdiff(names(nitro_plot_av),
                             names(fun_dis_env))), 
            by="plot_code2")

fun_dis_env1 <- fun_dis_env0%>%
  left_join(nitro_transect_av, by="transect_code2")

fun_dis_env <- fun_dis_env1%>%
  left_join(nitro_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(nitro_perc=case_when(
    !is.na(avg_plot_n_perc)~avg_plot_n_perc,
    (is.na(avg_plot_n_perc)&
      !is.na(avg_transect_n_perc))~avg_transect_n_perc,
    (is.na(avg_plot_n_perc)&
      is.na(avg_transect_n_perc))~avg_field_n_perc
  ))

rm(fun_dis_env0, fun_dis_env1)


nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)
fun_dis_env%>%filter(grepl("_NA", plot_code2))


## add organic matter content

omatter_plot_av <- omatter%>%
  group_by(plot_code2)%>%
  summarise(avg_plot_om_perc=mean(`Percent organic matter`, na.rm=T),
            sd_plot_om_perc=sd(`Percent organic matter`, na.rm=T))

omatter_transect_av <- omatter%>%# it is not needed to average plot level data at it has only one year sampled
  group_by(transect_code2)%>%
  summarise(avg_transect_om_perc=mean(`Percent organic matter`, na.rm=T),
            sd_transect_om_perc=sd(`Percent organic matter`, na.rm=T))

omatter_field_av <- omatter%>%
  group_by(OldField)%>%
  summarise(avg_field_om_perc=mean(`Percent organic matter`, na.rm=T),
            sd_field_om_perc=sd(`Percent organic matter`, na.rm=T))#sd is quite high


fun_dis_env0 <- fun_dis_env%>%
  left_join(omatter_plot_av, by="plot_code2")

fun_dis_env1 <- fun_dis_env0%>%
  left_join(omatter_transect_av, by="transect_code2")

fun_dis_env <- fun_dis_env1%>%
  left_join(omatter_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(omatter_perc=case_when(
    !is.na(avg_plot_om_perc)~avg_plot_om_perc,
    (is.na(avg_plot_om_perc)&
       !is.na(avg_transect_om_perc))~avg_transect_om_perc,
    (is.na(avg_plot_om_perc)&
       is.na(avg_transect_om_perc))~avg_field_om_perc
  ))

rm(fun_dis_env0, fun_dis_env1)


nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)
fun_dis_env%>%filter(grepl("_NA", plot_code2))


## add light penetration

light_plot_av <- light%>%
  group_by(plot_code2)%>%
  summarise(avg_plot_light_perc=mean(`Percent light penetration`, na.rm=T),
            sd_plot_light_perc=sd(`Percent light penetration`, na.rm=T))

light_plot_av <- light%>%
  dplyr::select(OldField,
         Transect,
         transect_code,
         plot_code,
         transect_code2, 
         plot_code2)%>%
  distinct()%>%
  left_join(light_plot_av, by="plot_code2")


light_transect_av <- light_plot_av%>%
  group_by(transect_code2)%>%
  summarise(avg_transect_light_perc=mean(avg_plot_light_perc, na.rm=T),
            sd_transect_light_perc=sd(avg_plot_light_perc, na.rm=T))

light_field_av <- light_plot_av%>%
  group_by(OldField)%>%
  summarise(avg_field_light_perc=mean(avg_plot_light_perc, na.rm=T),
            sd_field_light_perc=sd(avg_plot_light_perc, na.rm=T))#sd is half of the mean


fun_dis_env0 <- fun_dis_env%>%
  left_join(light_plot_av%>%
              dplyr::select(plot_code2,
                     setdiff(names(light_plot_av),
                             names(fun_dis_env))), 
            by="plot_code2")

fun_dis_env1 <- fun_dis_env0%>%
  left_join(light_transect_av, by="transect_code2")

fun_dis_env <- fun_dis_env1%>%
  left_join(light_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(light_perc=case_when(
    !is.na(avg_plot_light_perc)~avg_plot_light_perc,
    (is.na(avg_plot_light_perc)&
       !is.na(avg_transect_light_perc))~avg_transect_light_perc,
    (is.na(avg_plot_light_perc)&
       is.na(avg_transect_light_perc))~avg_field_light_perc
  ))

rm(fun_dis_env0, fun_dis_env1)

nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)
fun_dis_env%>%filter(grepl("_NA", plot_code2))


#* remove all plot and transect level variables
#* keep field level data and across levels env data

names(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  dplyr::select(-c(avg_plot_c_perc,
            avg_plot_n_perc,
            avg_plot_om_perc,
            avg_plot_light_perc,
            sd_plot_c_perc,
            sd_plot_n_perc,
            sd_plot_om_perc,
            sd_plot_light_perc,
            avg_transect_c_perc,
            avg_transect_n_perc,
            avg_transect_om_perc,
            avg_transect_light_perc,
            sd_transect_c_perc,
            sd_transect_n_perc,
            sd_transect_om_perc,
            sd_transect_light_perc,
            sd_field_c_perc,
            sd_field_n_perc,
            sd_field_om_perc,
            sd_field_light_perc))

fun_dis_env
nrow(fun_dis_env)
nrow(fun_dis)

fun_dis <- fun_dis_env

#**** add plot-trait representativeness *****

#* check plot with hihgh representativeness, considering that
#* as representation percent higher than 75% and add this info as a
#* new column to subset datatable 

high_repres <- repres%>%
  filter(repres_percent>=0.5)

nrow(repres)
nrow(high_repres)

#* 15122/17063 = 88.6% of plots have trait data for more than 50% 
#* of plot cover

fun_dis <- fun_dis%>%
  mutate(repres=case_when(
    plot_code%in%high_repres$sample_code~"Y",
    !plot_code%in%high_repres$sample_code~"N"
  ))

a <- fun_dis%>%
  filter(repres=="Y")%>%
  dplyr::select(plot_code)%>%
  distinct()%>%
  nrow()


b <- fun_dis%>%
  filter(repres=="N")%>%
  dplyr::select(plot_code)%>%
  distinct()%>%
  nrow()

sum(a+b)

phylo_dis%>%
  dplyr::select(plot_code)%>%
  distinct()%>%
  nrow()

unique(is.na(fun_dis$repres))
nrow(repres)
rm(a,b)
#* the lack of coincidence between repres rows and complete_traits,
#* or dissimilarities data frames are due to the lack of some plots.
#* While repres contain all the plots present in the study
#* the other dataframes just contains those plots with at least one
#* species with complete trait data for the selected functional traits
#* that means that 17063-16944 only 119 plots have no data for any 
#* of their species

##******** join all datasets ********##

names(fun_dis)
names(phylo_dis)

nrow(fun_dis)
nrow(phylo_dis)

unique(setdiff(unique(phylo_dis$sp_code), unique(fun_dis$sp_code)))
unique(setdiff(unique(fun_dis$sp_code), unique(phylo_dis$sp_code)))

dis_df <- fun_dis%>%
  dplyr::select(sp_code,Plot,Transect,
         OldField, Year,
         setdiff(names(fun_dis), names(phylo_dis)))%>%
  left_join(phylo_dis, by=c("sp_code", "Plot",
                            "Transect", "OldField",
                            "Year"))%>%
  rename(fun_w_dist=dissim)%>%
  mutate(sp_code2=paste0(
    OldField,"_",
    Transect, "_",
    Plot, "_",
    tpl_species
  ))


nrow(dis_df)
names(dis_df)

unique(is.na(phylo_dis$tpl_species))
unique(is.na(phylo_dis$OldField))
unique(is.na(fun_dis$OldField))
unique(is.na(dis_df$OldField))
unique(is.na(dis_df$sp_code2))


#**** add time since of abandonment **********

#* correct error in YearAb variable which come from the
#* original e014 abundances datatable. The database contained species
#* within Year Ab variable, plant cover in species column and 
#* na y cover column. Nevertheless, we can obtain the lacking data from 
#* yearab_df dataset and knowing that year of abandonment is information
#* at field level

unique(is.na(dis_df$YearAbandoned))

dis_df%>%
  filter(is.na(YearAbandoned))%>%
  dplyr::select(OldField, Transect, 
         YearAb,YearAbandoned)%>%
  distinct()

dis_df <- dis_df%>%
  mutate(YearAbandoned=as.character(YearAbandoned))%>%
  mutate(yearab2=case_when(
    is.na(YearAbandoned) ~ YearAb,
    !is.na(YearAbandoned) ~ YearAbandoned
  ))

unique(is.na(dis_df$yearab2))

# test correspondence between Adam Clark and Cedar creek database

yab_filter <- dis_df%>%
  filter(grepl("[[:digit:]]", YearAb))%>%
  dplyr::select(OldField,YearAb, YearAbandoned)%>%
  distinct()%>%
  mutate(YearAb=as.numeric(YearAb),
         YearAbandoned=as.numeric(YearAbandoned))%>%
  drop_na()

cor(yab_filter$YearAb, yab_filter$YearAbandoned)

plot(yab_filter$YearAb~yab_filter$YearAbandoned)# plot 28 was abandon in 1991

dis_df <- dis_df%>%
  dplyr::select(-c(YearAb, YearAbandoned))%>%
  mutate(YearAb=as.numeric(yearab2))

#add time since abandon

dis_df <- dis_df%>%
  mutate(time_since_ab=Year-YearAb)

unique(is.na(dis_df$time_since_ab))
nrow(dis_df)

#**** add year of colonization and minimum residence time ******

#* return all the species that have appear at some point in each plot
#* check year of introduction in US https://www.invasivespeciesinfo.gov/terrestrial/plants

# inv_pres <- spxplot%>%
#   filter(Origin!="Native")

spxplot <- dis_df%>%
  dplyr::select(sp_code2,
         OldField,
         Transect, Plot,
         plot_code2,
         species, genus,
         sp, sub_sp,
         tpl_species,
         tpl_family,
         Origin)%>%
  distinct()

spxplot_list <- spxplot$sp_code2
colon_year_list <- list()

for(i in 1:length(spxplot_list)){
  
  dis_df_sp <- dis_df%>%
    filter(sp_code2==spxplot_list[i])
  
  year_pres <- unique(dis_df_sp$Year)
  colon_year <- min(year_pres)
  
  sel_plot_code <- unique(dis_df_sp$plot_code2)
  
  year_sampling_plot <- dis_df%>%
    filter(plot_code2==sel_plot_code)%>%
    dplyr::select(Year)%>%
    arrange(Year)%>%
    distinct()%>%pull()## check that as far as I'm not sure if all plots were sampled all the years
  
  mean_sampling_gap <- diff(year_sampling_plot)%>%
    mean()%>%
    round(.,0)
  
  first_plot_sample <- min(year_sampling_plot)
  colon_year_order <- which(year_sampling_plot==colon_year)
  pre_colon_sample <- if(colon_year_order==1){
    NA_real_}else{
      year_sampling_plot[colon_year_order-1]}
  
  min_year_df <- dis_df_sp%>%
    mutate(colon_year= colon_year,
    time_since_colon= Year-colon_year,
    first_plot_sample= first_plot_sample, 
    pre_colon_sample=pre_colon_sample,
    mean_sampling_gap=mean_sampling_gap)# case when do not work when pre_colon_sample si num(0), ie. colon_year_order ==1
  
  colon_year_list[[i]] <- min_year_df
  
  progress_item <- i/length(spxplot_list)*100
  print(progress_item)
  
}
  
colon_year_df <- colon_year_list%>%
  bind_rows()

nrow(colon_year_df)
nrow(dis_df)

#* add case_when for minimum residence time -mrt-.
#* We will estimate it as time since colon plus
#* half of the period between colon year and the sampling
#* year previous to the colon year, meaning that the species
#* arrived at plot in the middle of the period between 
#* field samplings. In case the species appeared in the plot
#* the first sampling year we will estimate colonization year
#* as year of abandonment plus average period between abandoment
#* and colonization for the given species. So, mrt is the 
#* difference between estimated colonization and current sample.
#* Lag between abandonment and colonization will be estimated
#* for species appearing after the first sample year and species
#* appearing the first sample year but close to year of 
#* abandonment -precise colon time -

plot(density(colon_year_df$mean_sampling_gap, na.rm=T))
gap_average <- modeest::mlv(colon_year_df$mean_sampling_gap, na.rm=T)

colon_year_df <- colon_year_df%>%
  mutate(colon_less_yearab=colon_year-YearAb,
         colon_estim_precision=case_when(
           (colon_year==first_plot_sample&#species appear the first sample year
           colon_less_yearab>gap_average)~"low",# first sample was long after the abandonment
           (colon_year>first_plot_sample|#species appear after the first sample year
            colon_less_yearab<=gap_average)~"high"# or appear the first sample year but soon after abandonment
         ))
         

imprecise_colon_df <- colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  dplyr::select(species, YearAb, Year,
         first_plot_sample,
         pre_colon_sample,
         colon_year, colon_less_yearab,
         colon_estim_precision,
         OldField, Transect, Plot,
         sp_code2)%>%distinct()# simplify

a <- imprecise_colon_df%>%
  nrow()

b <- colon_year_df%>%
  filter(colon_estim_precision=="high")%>%
  nrow()

c <- colon_year_df%>%
  filter(is.na(YearAb))%>%
  nrow()

nrow(colon_year_df)
a+b+c ##a>>b! this is important

sp_imprecise_sample <- unique(imprecise_colon_df$sp_code2)

sp_imprecise_target <- colon_year_df%>%
  filter(sp_code2%in%sp_imprecise_sample)%>%
  dplyr::select(species)%>%
  distinct()%>%
  arrange(species)%>%
  pull()

spcolon_lag_list <- list()

for(i in 1:length(sp_imprecise_target)){
  
  first_year_old <- colon_year_df%>%
    filter(species==sp_imprecise_target[i])%>%
    filter(colon_estim_precision=="low")%>%
    dplyr::select(species, YearAb, Year,
           first_plot_sample,
           pre_colon_sample,
           colon_year, colon_less_yearab,
           colon_estim_precision,
           OldField, Transect, Plot,
           sp_code2)#yearxsp
  
  imprecise <- nrow(first_year_old)
  
  nofirst_df <- colon_year_df%>%
    filter(species==sp_imprecise_target[i]#&
           # OldField==unique(sp_field_target$OldField)&
           # Transect==unique(sp_field_target$Transect)
           )%>%#if possible add transect if else
    filter(colon_estim_precision=="high")%>%
    dplyr::select(species, YearAb, Year,
           first_plot_sample,
           pre_colon_sample,
           colon_year, colon_less_yearab,
           colon_estim_precision,
           OldField, Transect, Plot,
           sp_code2)
  
  precise <- nrow(nofirst_df)
  
  na_yearab <- colon_year_df%>%
    filter(species==sp_imprecise_target[i])%>%
    filter(is.na(YearAb))%>%
    nrow()
  
  colon_year_df%>%
    filter(species==sp_imprecise_target[i])%>%
    nrow()
  
  sum(imprecise, precise, na_yearab)
  
  if(precise>1){
  
  plot(density(nofirst_df$colon_less_yearab),col="red",
       main=sp_imprecise_target[i])
  
  max_dens <- max(density(nofirst_df$colon_less_yearab)$y)
  more_likely_gap <- modeest::mlv(nofirst_df$colon_less_yearab, na.rm=T)%>%
    round(., 0)
  
  ggplot(data=nofirst_df, aes(x=colon_less_yearab))+
    geom_density(color="#33CC00", lwd=1)+
    xlim(c(0, max(nofirst_df$colon_less_yearab)))+
    xlab("Time lag between field abandon and colonization")+
    ylab("Density")+
    ggtitle(sp_imprecise_target[i])+
    annotate("text", x = max(nofirst_df$colon_less_yearab)-15,
             y = max_dens-max_dens*0.1, 
             label = paste0("More likely time lag= ", more_likely_gap),
             family="serif")+
    annotate("text", x = max(nofirst_df$colon_less_yearab)-15,
             y = max_dens-max_dens*0.2, 
             label = paste0("Nº imprecise colon= ", imprecise),
             family="serif")+
    annotate("text", x = max(nofirst_df$colon_less_yearab)-15,
             y = max_dens-max_dens*0.3, 
             label = paste0("Nº precise colon= ", precise),
             family="serif")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(family="serif"),
          axis.text = element_text(color="black"))} else { }
  
  spcolon_lag_df <- tibble(
    species=sp_imprecise_target[i])
  spcolon_lag_df <- spcolon_lag_df%>%
    mutate(colon_delay=case_when(
      precise>2 ~ more_likely_gap,
      precise==1~ nofirst_df$colon_less_yearab[1],
      precise==2~ round(mean(nofirst_df$colon_less_yearab),0),
      precise==0 ~ min(first_year_old$colon_less_yearab-(gap_average/2))
  ),
  colon_delay_mean=round(mean(nofirst_df$colon_less_yearab),0),
  min_delay_imprecise=min(first_year_old$colon_less_yearab-(gap_average/2)),
  n_imprecise=imprecise,
  n_precise=precise)#try not to use case when with a object not included in the table
  
  spcolon_lag_list[[i]] <- spcolon_lag_df
  # create a datatable with two columns species and more likely time lapse between abandon and colonization (mode)
  
}

sp_delay_df <- spcolon_lag_list%>%
  bind_rows()

length(sp_imprecise_target)
nrow(sp_delay_df)

hist(sp_delay_df$colon_delay)
sp_delay_df

#* add colon delay to year of abandon. In case of adding 
#* colon delay to YearAb lead to a estimate colonization year higher 
#* than the actual first observed year, we will use as estimated colon
#* date the colonization year minus half of the sample gap, as in rest
#* of situations, which is condidered as the most recent colon
#* date as possible in the rest of cases. 

colon_year_df1 <- colon_year_df%>%
  left_join(sp_delay_df%>%
              dplyr::select(species, colon_delay,
                     colon_delay_mean), by="species")%>%
  distinct()# spcolon_lag_df has more than one row

nrow(colon_year_df)
nrow(colon_year_df1)
names(colon_year_df1)

colon_year_df <- colon_year_df1
rm(colon_year_df1)

colon_year_df <- colon_year_df%>%
  mutate(YearAb_plus_delay1= YearAb+colon_delay,
         YearAb_plus_delaymean= YearAb+colon_delay_mean)

colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  filter(YearAb_plus_delay1>=colon_year)%>%
  nrow()# third part of total imprecise cases will be supposed to have arrived recently to the plot in case of using the mode of delay

colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  filter(YearAb_plus_delaymean>=colon_year)%>%
  nrow()# almost half of the imprecise cases will be supposed to have arrived recently to the plot in case of using average delay

colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  nrow()

mrt_df <- colon_year_df%>%
    mutate(estim_colon_year1=case_when(#1 is mode based
    colon_year> first_plot_sample~# also !is.na(pre_colon_sample)
      round(pre_colon_sample+(colon_year-pre_colon_sample)/2, 0),
    (colon_estim_precision=="high"&
      colon_year==first_plot_sample)~
      round(YearAb+colon_less_yearab/2, 0),
    (colon_estim_precision=="low"&
      YearAb_plus_delay1>=colon_year)~
      round(colon_year-gap_average/2, 0),# as recent colon as average recent colon
    (colon_estim_precision=="low"&
       YearAb_plus_delay1<colon_year)~
      round(YearAb+colon_delay,0)),
    estim_colon_year2=case_when(# 2 is average based
      colon_year> first_plot_sample~# also !is.na(pre_colon_sample)
      round(pre_colon_sample+(colon_year-pre_colon_sample)/2, 0),
    (colon_estim_precision=="high"&
      colon_year==first_plot_sample)~
      round(YearAb+colon_less_yearab/2, 0),
    (colon_estim_precision=="low"&
      YearAb_plus_delaymean>=colon_year)~
      round(colon_year-gap_average/2, 0),# as recent colon as average recent colon
    (colon_estim_precision=="low"&
      YearAb_plus_delaymean<colon_year)~
      round(YearAb+colon_delay_mean,0))
    )

hist(mrt_df$estim_colon_year1)
max(mrt_df$estim_colon_year1, na.rm=T)

hist(mrt_df$estim_colon_year2)
max(mrt_df$estim_colon_year2, na.rm=T)

mrt_df <- mrt_df%>%
  mutate(mrt1=Year-estim_colon_year1,
         mrt2=Year-estim_colon_year2)
  
hist(mrt_df$mrt1)
hist(log(mrt_df$mrt1))
min(mrt_df$mrt1, na.rm=T)

hist(mrt_df$mrt2)
hist(log(mrt_df$mrt2))
min(mrt_df$mrt2, na.rm=T)


zd <- mrt_df%>%
  filter(mrt1==0)%>%
  dplyr::select(OldField,
         plot_code, tpl_species,
         Origin, YearAb,
         first_plot_sample,
         Year, 
         colon_year,
         colon_delay,
         estim_colon_year1,
         estim_colon_year2,
         time_since_colon,
         pre_colon_sample,
         mrt1,
         mrt2)

zd%>%
  dplyr::select(OldField)%>%
  distinct()#601 field has mrt=0 since it was inmediately samplet after abandonment

mrt_df <- mrt_df%>%
  rename(first_sp_record=colon_year)

nrow(mrt_df)
nrow(dis_df)
  
## join with dissimilarities df

dis_df1 <- dis_df%>%
  left_join(mrt_df%>%
              dplyr::select(sp_code, first_sp_record,
                     colon_delay, colon_delay_mean,
                     estim_colon_year1, 
                     estim_colon_year2, mrt1,mrt2), 
            by="sp_code")%>%
  distinct()

nrow(dis_df1)
nrow(dis_df)
names(dis_df1)  
dis_df <- dis_df1

rm(dis_df1)


###****add also burned category********##

sort(names(field_info))
dis_df <- dis_df%>%
  rename(Burned=Burned.)
unique(dis_df$Burned)

dis_df%>%
  filter(is.na(Burned))%>%
  dplyr::select(OldField, Transect)%>%
  distinct()

## it lacks OldFields 22,29,69

burn_df
burn_df1 <- burn_df%>%
  dplyr::select(-ExperimentNumber)%>%
  mutate(OldField=as.numeric(OldField),
         Transect=str_remove_all(Transect, "[ ]"))%>%
  mutate(transect_code2=paste0(OldField, "_",
                               Transect))

sort(unique(burn_df1$OldField))
sort(unique(dis_df$OldField))
sort(unique(burn_df1$Transect))

burn_df <- burn_df1
rm(burn_df1)

dis_df <- dis_df%>%
  left_join(burn_df%>%
              dplyr::select(transect_code2, 
                     BurnTreatment), by="transect_code2")

dis_df <- dis_df%>%
  mutate(burn=case_when(
    is.na(Burned)~BurnTreatment,
    !is.na(Burned)~Burned
  ))

unique(dis_df$burn)
dis_df <- dis_df%>%
  dplyr::select(-Burned, BurnTreatment)%>%
  mutate(BurnTreatment=burn)%>%
  dplyr::select(-burn)

## add time since first burning

dis_df <- dis_df%>%
  mutate(time_since_burn=Year-FirstYearAfterBurning)


###**** add melting time (time since first invasive arrives) ********##


#* check the origin of introduced/or native species using gbif distribution
#* and https://plants.usda.gov/home/plantProfile?symbol=POPR 
#* just making two vectors may work and then case_when %in%

## substitute the unknown status for the next species

dis_df <- dis_df%>%
  mutate(origin=case_when(
    species=="Aristida sp."~"Native",
    species=="Helianthus sp."~"Native",
    species=="Oenothera sp." ~ "Native",
    species=="Rhus sp." ~"Native",
    species=="Setaria sp." ~"Introduced",
    species=="Stachys sp." ~ "Native", 
    species=="Tradescantia sp." ~ "Native",
    !species%in%c("Aristida sp.",
                  "Helianthus sp.",
                  "Oenonthera sp.",
                  "Rhus sp.",
                  "Setaria sp.",
                  "Stachys sp.",
                  "Tradescantia sp.")~Origin
  ))

zo <- dis_df%>%
  filter(sp=="sp."&
         !Origin%in%c("Native", "Introduced"))%>%
  dplyr::select(species, Origin, origin)%>%
  distinct()

dis_df <- dis_df%>%
  mutate(origin2=case_when(
    origin=="Introduced"~"Introduced",
    origin!="Introduced"~"Rest"
  ),
    origin3=case_when(
      origin=="Native"~"Native",
      origin!="Native"~"Rest"
  ))# consider all unknown species as natives


#* get time since invasion -melting time - THis should be
#* a plot level value. But in case of invasive species it is not
#* clear to me whether a plot value could be used -ie. what means 
#* a melting time higher than mrt for invasive species, does the arrival
#* of other previous invasive species affect the ulterior arrival of
#* other introduced species?- 
#* I will use a plot level data with maximum time since first
#* invasive species appeared in the plot. Most cases have only 
#* one invasive species other have no invasive species,
#* so value will be NA or negative value as 0 
#* is value for recent invasions. 


plot_inv_time <- dis_df%>%
  filter(origin2=="Introduced")%>%
  group_by(plot_code)%>%
  summarise(melting_time1=max(mrt1, na.rm=T),
            melting_time2=max(mrt2, na.rm=T))

dis_df1 <- dis_df%>%
  left_join(plot_inv_time%>%
              dplyr::select(plot_code, 
                     melting_time1,
                     melting_time2), 
            by="plot_code")

nrow(dis_df1)
nrow(dis_df)
names(dis_df1)  
dis_df <- dis_df1

rm(dis_df1)


unique(is.na(dis_df$melting_time1))
unique(is.na(dis_df$mrt1))

ze <- dis_df%>%filter(is.na(melting_time1))%>%
  dplyr::select(OldField, plot_code,
         tpl_species, 
         first_sp_record,
         estim_colon_year1,
         origin2,
         mrt1, 
         melting_time1)## some plots have no introduced species and have NA

#* replace NA out of field 600 and 601 - which has NA since 
#* they lack YearAb and were only sampled in 2016 -  by
#* - Inf has they have never been colonized by introduced species 

dis_df1 <- dis_df%>%
  mutate(melting_time1_inf=replace_na(
    melting_time1, -Inf
  ),
  melting_time2_inf=replace_na(
    melting_time2, -Inf))

dis_df1%>%
  filter(is.na(melting_time1_inf))%>%
  dplyr::select(OldField, 
         melting_time1_inf)%>%
  distinct()


nrow(dis_df1)
nrow(dis_df)

dis_df <- dis_df1
rm(dis_df1)

#*** add number of invasive species per plot ****

dis_df_ninv <- dis_df%>%
  group_by(origin2, plot_code)%>%
  summarise(n=n())

dis_df_ninv <- dis_df_ninv%>%
  pivot_wider(names_from=origin2,
              values_from=n)%>%
  rename(n_intro=Introduced,
         n_native=Rest)%>%
  replace_na(list(n_intro = 0,
                  n_native = 0))%>%
  mutate(intr_ratio=n_intro/(n_intro+n_native))

hist(dis_df_ninv$intr_ratio)

dis_df1 <- dis_df%>%
  left_join(dis_df_ninv, 
            by="plot_code")

nrow(dis_df1)
nrow(dis_df)
names(dis_df1)  
dis_df <- dis_df1

rm(dis_df1, dis_df_ninv)


#**** add phylogenetic similarity as variable which  has better ***
#**** distribution than dissimilarity ******

dis_df <- dis_df%>%
  mutate(phylo_w_sim=1000-phylo_w_dist)

nrow(dis_df)#old version was 106325
dis_df%>%
  filter(genus=="Quercus")%>%
  nrow()

##**** save def analytical dataframe ****##

# data.table::fwrite(dis_df,
#                    "../../output/tables/analyses/analysis_df_plot_ctr3d.csv",
#                    sep=";", dec=".", col.names=T, row.names=F, quote=F)

data.table::fwrite(dis_df,
                   "../../output/tables/analyses/analysis_df_plot_pair2d.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)



#**************************************###
#******1.2.transect level tables******####
#**************************************###

coord <- read_delim( "../../data/abundances/location_matrix_byplot.csv", 
                     delim=",", col_names=T)

unique(coord$transect)

coord <- coord%>%
  dplyr::select(-c(xloc, yloc))%>%
  rename(Field=field,
         Transect_old=transect,
         Plot=plot)%>% 
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D"))

coord_transect <- coord%>%
  filter(Plot==1)## choose coordinate of plot1 as representative of the transect

nrow(coord_transect)

field_info <- read_delim( "../../data/abundances/field_information.csv", 
                          delim=",", col_names=T)

unique(field_info$Transect)

field_info <- field_info%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D"))

nrow(field_info)

field_info <- coord_transect%>%
  dplyr::select(-c(Plot, Transect_old))%>%
  left_join(field_info, by=c("Field", "Transect"))

field_info <- field_info%>%
  rename(OldField=Field)
nrow(field_info)

phylo_dis <- read_delim( "../../output/tables/dissim/phylo_distances_transect.csv", 
                         delim=";", col_names=T)

# fun_dis <- read_delim( "../../output/tables/dissim/dissimilarity_centroid3d_transect.csv",
#                        delim=";", col_names=T)

fun_dis <- read_delim( "../../output/tables/dissim/dissimilarity_weight_distances3d_transect.csv",
                       delim=";", col_names=T)


setdiff(fun_dis$OldField, field_info$OldField)
unique(fun_dis$Transect)
unique(phylo_dis$Transect)

fun_dis <- fun_dis%>% 
  dplyr::select(-c(species_code,
                   sp_code, genus, 
                   sp, sub_sp))%>%
  mutate(genus = word(species, 1, sep = fixed(" ")),
         sp = word(species, 2, sep = fixed(" ")),
         sub_sp= word(species, 3, sep = fixed(" ")))%>%#temporal
  mutate(transect_code=paste0(Year,"_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect),
         species_code= paste0(substring(genus, 1, 3), "_",
                              substring(sp, 1,3)),
         sp_code=paste0(species_code,"_",
                        transect_code))

fun_dis <- fun_dis%>%
  left_join(field_info, 
            by=c("OldField", "Transect"))


phylo_dis <- phylo_dis%>%
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
  ))%>%
  mutate(transect_code=paste0(Year,"_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect),
         sp_code=paste0(species_code,"_",
                        transect_code))


nrow(fun_dis)
nrow(phylo_dis)
unique(is.na(phylo_dis$OldField))
unique(is.na(fun_dis$OldField))

all(sort(unique(phylo_dis$OldField)) == sort(unique(fun_dis$OldField)))


fun_dis <- fun_dis%>%#mean_cover_percent is missing, don't know why
  left_join(phylo_dis%>%
              dplyr::select(sp_code, mean_cover_percent),
            by="sp_code"
  )

unique(fun_dis$mean_cover_percent)

#*consider substitute by dissimilarity_centroid3d
#*dissimilarity_weight_distances2d
#*dissimilarity_weight_distances3d

repres <- read_delim( "../../output/tables/abundances/cover_representativeness.csv", 
                      delim=";", col_names=T)

repres <- repres%>%
  filter(level=="transect")

#*consider substitute by nrow_representativeness

yearab_df <- read_delim( "../../data/abundances/e014_Year of abandon.txt")

burn_df <- read_delim( "../../data/abundances/e014_Burned info.txt")

burn_df <- burn_df%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~ "A",
    Transect_old=="R"~ "B",
    Transect_old=="W"~ "C",
    Transect_old=="Y"~ "D",
    Transect_old=="G "~"A",
    Transect_old=="R "~"B",
    Transect_old=="W "~"C",
    Transect_old=="Y "~"D",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ))

#***** add other field environmental data ********

carbon <-read_delim(file="../../data/abundances/e014_Soil carbon data.txt")
nitro <-read_delim(file="../../data/abundances/e014_Soil nitrogen data.txt")
omatter <-read_delim(file="../../data/abundances/e014_Soil organic matter data.txt")
light <-read_delim(file="../../data/abundances/e014_Percent light penetration data.txt")

names(carbon)
names(nitro)
names(omatter)
names(light)

carbon <- carbon%>%
  rename(OldField=field,
         Transect_old=transect,
         depth_c_cm="Depth (cm)",
         Plot= quadrat, 
         YearAb_carbon=CultivatedHistorical,
         YearAb_photo_carbon=CultivatedPhoto)%>%
  mutate(Transect_old=as.character(Transect_old))%>%
  mutate(Transect = case_when(
    Transect_old=="1"~"A",
    Transect_old=="2"~"B",
    Transect_old=="3"~"C",
    Transect_old=="4"~"D",
    Transect_old=="5"~"E",
    Transect_old=="6"~"F"
  ))%>%
  mutate(plot_code=paste0(Year, "_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year, "_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect))

nitro <- nitro%>%
  rename(OldField=field,
         Transect_old=transect,
         depth_n_cm="Depth (cm)",
         Plot= quadrat)%>%
  mutate(Transect_old=as.character(Transect_old))%>%
  mutate(Transect = case_when(
    Transect_old=="1"~"A",
    Transect_old=="2"~"B",
    Transect_old=="3"~"C",
    Transect_old=="4"~"D",
    Transect_old=="5"~"E",
    Transect_old=="6"~"F"
  ))%>%
  mutate(plot_code=paste0(Year, "_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year, "_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect))%>%
  dplyr::select(-c(CultivatedPhoto,
            CultivatedHistorical,
            LastCrop,
            Area_ha))

omatter <- omatter%>%
  rename(Experiment="Experiment number",
         OldField="Old field designator",
         Sampling_date_om="Sampling date",
         Plot=Quadrat)%>%
  rename(Transect_old=Transect)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ))%>%
  mutate(Year=paste0("19", 
                     substr(as.character(Sampling_date_om), start=1,stop=2)),
         plot_code=paste0(Year, "_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year, "_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect))%>%
  mutate(Year=as.numeric(Year))%>%
  dplyr::select(-c(CultivatedHistorical,
            LastCrop,
            Area_ha,
            Experiment))

light <- light%>%
  rename(Experiment="Experiment number",
         OldField="Old field code",
         Transect_old="Transect code",
         Sampling_date_light="Sampling date",
         Year="Sampling year",
         Plot=Quadrat,
         YearAb_light=CultivatedHistorical)%>%
  mutate(Transect = case_when(
    Transect_old=="G"~"A",
    Transect_old=="R"~"B",
    Transect_old=="W"~"C",
    Transect_old=="Y"~"D",
    Transect_old%in%c("A","B","C",
                      "D","E","F")~Transect_old
  ))%>%
  mutate(plot_code=paste0(Year, "_",
                          OldField, "_",
                          Transect, "_",
                          Plot),
         plot_code2=paste0(OldField, "_",
                           Transect, "_",
                           Plot),
         transect_code=paste0(Year, "_",
                              OldField, "_",
                              Transect),
         transect_code2=paste0(OldField, "_",
                               Transect))%>%
  dplyr::select(-c(LastCrop,
            Area_ha,
            Experiment))

head(carbon)
head(nitro)
head(omatter)
head(light)

unique(fun_dis$Year)
unique(carbon$Year)
unique(nitro$Year)
unique(omatter$Year)
unique(light$Year)

sort(unique(fun_dis$OldField))
sort(unique(carbon$OldField))
sort(unique(nitro$OldField))
sort(unique(omatter$OldField))
sort(unique(light$OldField))

sort(unique(fun_dis$Transect))
sort(unique(carbon$Transect))
sort(unique(nitro$Transect))
sort(unique(omatter$Transect))
sort(unique(light$Transect))

sort(unique(fun_dis$transect_code2))
sort(unique(carbon$transect_code2))
sort(unique(nitro$transect_code2))
sort(unique(omatter$transect_code2))
sort(unique(light$transect_code2))

setdiff(unique(carbon$transect_code2), unique(fun_dis$transect_code2))
setdiff(unique(nitro$transect_code2), unique(fun_dis$transect_code2))
setdiff(unique(omatter$transect_code2), unique(fun_dis$transect_code2))
setdiff(unique(light$transect_code2), unique(fun_dis$transect_code2))

sort(unique(fun_dis$Plot))
sort(unique(carbon$Plot))
sort(unique(nitro$Plot))
sort(unique(omatter$Plot))
sort(unique(light$Plot))

names(carbon)
names(nitro)
names(omatter)
names(light)

local_env <- list(carbon, nitro, omatter, light)%>%
  purrr::reduce(full_join, by=c("Year","OldField",
                                "Transect","Plot",
                                "plot_code", "plot_code2",
                                "transect_code",
                                "transect_code2"))%>%
  distinct()

z <- local_env%>%
  dplyr::select(plot_code2, YearAb_carbon,
         YearAb_light)

#* as there is no data for all plot and year combinations we need
#* to average them. Check if it is possible at plot level or 
#* should be at transect level

fun_dis%>%
  dplyr::select(OldField,
         Transect)%>%
  distinct()

nitro%>%
  dplyr::select(OldField,
         Transect)%>%
  distinct()

sort(unique(fun_dis$OldField))
sort(unique(fun_dis$Transect))

sort(unique(nitro$OldField))
sort(unique(nitro$Transect))

sort(unique(omatter$OldField))
sort(unique(omatter$Transect))

sort(unique(light$OldField))
sort(unique(light$Transect))

## add carbon content

carbon_plot_av <- carbon%>%
  group_by(plot_code2)%>%##get plot average across years
  summarise(avg_plot_c_perc=mean(ProportionC, na.rm=T),
            sd_plot_c_perc=sd(ProportionC, na.rm=T))

carbon_plot_av <- carbon%>%
  dplyr::select(OldField,
         Transect,
         transect_code,
         plot_code,
         transect_code2,
         plot_code2)%>%
  distinct()%>%
  left_join(carbon_plot_av, by="plot_code2")

carbon_transect_av <- carbon_plot_av%>%#avoid to give more weight to plots sampled more years
  group_by(transect_code2)%>%
  summarise(avg_transect_c_perc=mean(avg_plot_c_perc, na.rm=T),
            sd_transect_c_perc=sd(avg_plot_c_perc, na.rm=T))

carbon_field_av <- carbon_plot_av%>%
  group_by(OldField)%>%
  summarise(avg_field_c_perc=mean(avg_plot_c_perc, na.rm=T),
            sd_field_c_perc=sd(avg_plot_c_perc, na.rm=T))#sd is 0.3 of the mean


fun_dis_env0 <- fun_dis%>%
  left_join(carbon_transect_av%>%
              dplyr::select(setdiff(names(carbon_transect_av),
                             names(fun_dis)),
                     transect_code2),
            by="transect_code2")

fun_dis_env <- fun_dis_env0%>%
  left_join(carbon_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(carbon_perc=case_when(
    !is.na(avg_transect_c_perc)~avg_transect_c_perc,
    is.na(avg_transect_c_perc)~avg_field_c_perc
  ))

rm(fun_dis_env0)

nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)

## add nitrogen content

nitro_plot_av <- nitro%>%
  group_by(plot_code2)%>%
  summarise(avg_plot_n_perc=mean(ProportionN, na.rm=T),
            sd_plot_n_perc=sd(ProportionN, na.rm=T))

nitro_plot_av <- nitro%>%
  dplyr::select(OldField,
         Transect,
         transect_code,
         plot_code,
         transect_code2, 
         plot_code2)%>%
  distinct()%>%
  left_join(nitro_plot_av, by="plot_code2")

nitro_transect_av <- nitro_plot_av%>%
  group_by(transect_code2)%>%
  summarise(avg_transect_n_perc=mean(avg_plot_n_perc, na.rm=T),
            sd_transect_n_perc=sd(avg_plot_n_perc, na.rm=T))

nitro_field_av <- nitro_plot_av%>%
  group_by(OldField)%>%
  summarise(avg_field_n_perc=mean(avg_plot_n_perc, na.rm=T),
            sd_field_n_perc=sd(avg_plot_n_perc, na.rm=T))#sd is 0.14 of the mean


fun_dis_env0 <- fun_dis_env%>%
  left_join(nitro_transect_av%>%
              dplyr::select(transect_code2,
                     setdiff(names(nitro_transect_av),
                             names(fun_dis_env))), 
            by="transect_code2")

fun_dis_env <- fun_dis_env0%>%
  left_join(nitro_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(nitro_perc=case_when(
    !is.na(avg_transect_n_perc)~avg_transect_n_perc,
     is.na(avg_transect_n_perc)~avg_field_n_perc
  ))

rm(fun_dis_env0)
nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)

## add organic matter content

omatter_plot_av <- omatter%>%
  group_by(plot_code2)%>%
  summarise(avg_plot_om_perc=mean(`Percent organic matter`, na.rm=T),
            sd_plot_om_perc=sd(`Percent organic matter`, na.rm=T))

omatter_transect_av <- omatter%>%# it is not needed to average plot level data at it has only one year sampled
  group_by(transect_code2)%>%
  summarise(avg_transect_om_perc=mean(`Percent organic matter`, na.rm=T),
            sd_transect_om_perc=sd(`Percent organic matter`, na.rm=T))

omatter_field_av <- omatter%>%
  group_by(OldField)%>%
  summarise(avg_field_om_perc=mean(`Percent organic matter`, na.rm=T),
            sd_field_om_perc=sd(`Percent organic matter`, na.rm=T))#sd is quite high


fun_dis_env0 <- fun_dis_env%>%
  left_join(omatter_transect_av, by="transect_code2")

fun_dis_env <- fun_dis_env0%>%
  left_join(omatter_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(omatter_perc=case_when(
    !is.na(avg_transect_om_perc)~avg_transect_om_perc,
     is.na(avg_transect_om_perc)~avg_field_om_perc
  ))

rm(fun_dis_env0)
nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)

## add light penetration

light_plot_av <- light%>%
  group_by(plot_code2)%>%
  summarise(avg_plot_light_perc=mean(`Percent light penetration`, na.rm=T),
            sd_plot_light_perc=sd(`Percent light penetration`, na.rm=T))

light_plot_av <- light%>%
  dplyr::select(OldField,
         Transect,
         transect_code,
         plot_code,
         transect_code2, 
         plot_code2)%>%
  distinct()%>%
  left_join(light_plot_av, by="plot_code2")


light_transect_av <- light_plot_av%>%
  group_by(transect_code2)%>%
  summarise(avg_transect_light_perc=mean(avg_plot_light_perc, na.rm=T),
            sd_transect_light_perc=sd(avg_plot_light_perc, na.rm=T))

light_field_av <- light_plot_av%>%
  group_by(OldField)%>%
  summarise(avg_field_light_perc=mean(avg_plot_light_perc, na.rm=T),
            sd_field_light_perc=sd(avg_plot_light_perc, na.rm=T))#sd is half of the mean

fun_dis_env0 <- fun_dis_env%>%
  left_join(light_transect_av, by="transect_code2")

fun_dis_env <- fun_dis_env0%>%
  left_join(light_field_av, by="OldField")

fun_dis_env <- fun_dis_env%>%
  mutate(light_perc=case_when(
         !is.na(avg_transect_light_perc)~avg_transect_light_perc,
          is.na(avg_transect_light_perc)~avg_field_light_perc
  ))

rm(fun_dis_env0)

nrow(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  distinct()

nrow(fun_dis_env)

#* remove all plot and transect level variables
#* keep field level data and across levels env data

names(fun_dis_env)

fun_dis_env <- fun_dis_env%>%
  dplyr::select(-c(avg_transect_c_perc,
            avg_transect_n_perc,
            avg_transect_om_perc,
            avg_transect_light_perc,
            sd_transect_c_perc,
            sd_transect_n_perc,
            sd_transect_om_perc,
            sd_transect_light_perc,
            sd_field_c_perc,
            sd_field_n_perc,
            sd_field_om_perc,
            sd_field_light_perc))

fun_dis_env
nrow(fun_dis_env)
nrow(fun_dis)

fun_dis <- fun_dis_env
unique(fun_dis$transect_code2)
unique(is.na(fun_dis$transect_code2))
unique(is.na(fun_dis$omatter_perc))

fun_dis%>%
  filter(is.na(light_perc))%>%
  dplyr::select(OldField, transect_code2,
                nitro_perc,
                carbon_perc, light_perc,
                omatter_perc)%>%
  distinct()# these are the fields lacking in envrionmental datasets

#**** add transect-trait representativeness *****

#* check transect with hihgh representativeness, considering that
#* as representation percent higher than 75% and add this info as a
#* new column to subset datatable 

unique(repres$level)

high_repres <- repres%>%
  filter(repres_percent>=0.5)

nrow(repres)
nrow(high_repres)
unique(repres$sample_code)

#* 678/684 = 99.1% of plots have trait data for more than 50% 
#* of plot cover

fun_dis <- fun_dis%>%
  mutate(repres=case_when(
    transect_code%in%high_repres$sample_code~"Y",
    !transect_code%in%high_repres$sample_code~"N"
  ))

a <- fun_dis%>%
  filter(repres=="Y")%>%
  dplyr::select(transect_code)%>%
  distinct()%>%
  nrow()


b <- fun_dis%>%
  filter(repres=="N")%>%
  dplyr::select(transect_code)%>%
  distinct()%>%
  nrow()

sum(a+b)
rm(a,b)

##******** join all datasets ********##

names(fun_dis)
names(phylo_dis)

nrow(fun_dis)
nrow(phylo_dis)

unique(setdiff(unique(phylo_dis$sp_code), unique(fun_dis$sp_code)))
unique(setdiff(unique(fun_dis$sp_code), unique(phylo_dis$sp_code)))
unique(setdiff(unique(fun_dis$transect_code2), unique(phylo_dis$transect_code2)))
all(sort(unique(phylo_dis$transect_code2)) == sort(unique(fun_dis$transect_code2)))


dis_df <- fun_dis%>%
  dplyr::select(sp_code,Transect,
         OldField, Year,
         setdiff(names(fun_dis), names(phylo_dis)))%>%
  left_join(phylo_dis, by=c("sp_code",
                            "Transect", 
                            "OldField",
                            "Year"))%>%
  rename(fun_w_dist=dissim)%>%
  mutate(sp_code2=paste0(
    OldField,"_",
    Transect, "_",
    tpl_species
  ))


nrow(dis_df)
sort(names(dis_df))

unique(is.na(phylo_dis$tpl_species))
unique(is.na(phylo_dis$OldField))
unique(is.na(fun_dis$OldField))
unique(is.na(dis_df$OldField))
unique(is.na(dis_df$sp_code2))


#**** add time since of abandonment **********

#* correct error in YearAb variable which come from the
#* original e014 abundances datatable. The database contained species
#* within Year Ab variable, plant cover in species column and 
#* na y cover column. Nevertheless, we can obtain the lacking data from 
#* yearab_df dataset and knowing that year of abandonment is information
#* at field level

unique(is.na(dis_df$YearAbandoned))

dis_df%>%
  filter(is.na(YearAbandoned))%>%
  dplyr::select(OldField, Transect, 
         YearAb,YearAbandoned)%>%
  distinct()

dis_df <- dis_df%>%
  mutate(YearAbandoned=as.character(YearAbandoned))%>%
  mutate(yearab2=case_when(
    is.na(YearAbandoned) ~ YearAb,
    !is.na(YearAbandoned) ~ YearAbandoned
  ))

unique(is.na(dis_df$yearab2))

# test correspondence between Adam Clark and Cedar creek database

yab_filter <- dis_df%>%
  filter(grepl("[[:digit:]]", YearAb))%>%
  dplyr::select(OldField,YearAb, YearAbandoned)%>%
  distinct()%>%
  mutate(YearAb=as.numeric(YearAb),
         YearAbandoned=as.numeric(YearAbandoned))%>%
  drop_na()

cor(yab_filter$YearAb, yab_filter$YearAbandoned)

plot(yab_filter$YearAb~yab_filter$YearAbandoned)# plot 28 was abandon in 1991

dis_df <- dis_df%>%
  dplyr::select(-c(YearAb, YearAbandoned))%>%
  mutate(YearAb=as.numeric(yearab2))

#add time since abandon

dis_df <- dis_df%>%
  mutate(time_since_ab=Year-YearAb)

unique(is.na(dis_df$time_since_ab))
nrow(dis_df)


#**** add year of colonization and minimum residence time ******

#* return all the species that have appear at some point in each plot
#* check year of introduction in US https://www.invasivespeciesinfo.gov/terrestrial/plants

# inv_pres <- spxplot%>%
#   filter(Origin!="Native")

spxplot <- dis_df%>%
  dplyr::select(sp_code2,
         OldField,
         Transect, 
         species, genus,
         sp, sub_sp,
         tpl_species,
         tpl_family,
         Origin)%>%
  distinct()

spxplot_list <- spxplot$sp_code2
colon_year_list <- list()

for(i in 1:length(spxplot_list)){
  
  dis_df_sp <- dis_df%>%
    filter(sp_code2==spxplot_list[i])
  
  year_pres <- unique(dis_df_sp$Year)
  colon_year <- min(year_pres)
  
  sel_transect_code <- unique(dis_df_sp$transect_code2)
  
  year_sampling_transect <- dis_df%>%
    filter(transect_code2==sel_transect_code)%>%
    dplyr::select(Year)%>%
    arrange(Year)%>%
    distinct()%>%pull()## check that as far as I'm not sure if all transects were sampled all the years
  
  mean_sampling_gap <- diff(year_sampling_transect)%>%
    mean()%>%
    round(.,0)
  
  first_transect_sample <- min(year_sampling_transect)
  colon_year_order <- which(year_sampling_transect==colon_year)
  pre_colon_sample <- if(colon_year_order==1){
    NA_real_}else{
      year_sampling_transect[colon_year_order-1]}
  
    min_year_df <- dis_df_sp%>%
    mutate(colon_year= colon_year,
           time_since_colon= Year-colon_year,
           first_transect_sample= first_transect_sample, 
           pre_colon_sample=pre_colon_sample,
           mean_sampling_gap=mean_sampling_gap)# case when do not work when pre_colon_sample si num(0), ie. colon_year_order ==1
  
  colon_year_list[[i]] <- min_year_df
  
  progress_item <- i/length(spxplot_list)*100
  print(progress_item)
  
}

colon_year_df <- colon_year_list%>%
  bind_rows()

nrow(colon_year_df)
nrow(dis_df)

#* add case_when for minimum residence time -mrt-.
#* We will estimate it as time since colon plus
#* half of the period between colon year and the sampling
#* year previous to the colon year, meaning that the species
#* arrived at plot in the middle of the period between 
#* field samplings. In case the species appeared in the plot
#* the first sampling year we will estimate colonization year
#* as year of abandonment plus average period between abandoment
#* and colonization for the given species. So, mrt is the 
#* difference between estimated colonization and current sample.
#* Lag between abandonment and colonization will be estimated
#* for species appearing after the first sample year and species
#* appearing the first sample year but close to year of 
#* abandonment -precise colon time -

plot(density(colon_year_df$mean_sampling_gap, na.rm=T))
gap_average <- modeest::mlv(colon_year_df$mean_sampling_gap, na.rm=T)

colon_year_df <- colon_year_df%>%
  mutate(colon_less_yearab=colon_year-YearAb,
         colon_estim_precision=case_when(
           (colon_year==first_transect_sample&#species appear the first sample year
              colon_less_yearab>gap_average)~"low",# first sample was long after the abandonment
           (colon_year>first_transect_sample|#species appear after the first sample year
              colon_less_yearab<=gap_average)~"high"# or appear the first sample year but soon after abandonment
         ))

unique(colon_year_df$colon_year)

### continue working here!

imprecise_colon_df <- colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  dplyr::select(species, YearAb, Year,
         first_transect_sample,
         pre_colon_sample,
         colon_year, colon_less_yearab,
         colon_estim_precision,
         OldField, Transect,
         sp_code2)%>%distinct()# simplify

a <- imprecise_colon_df%>%
  nrow()

b <- colon_year_df%>%
  filter(colon_estim_precision=="high")%>%
  nrow()

c <- colon_year_df%>%
  filter(is.na(YearAb))%>%
  nrow()

nrow(colon_year_df)
a+b+c ##a>>b! this is important

sp_imprecise_sample <- unique(imprecise_colon_df$sp_code2)

sp_imprecise_target <- colon_year_df%>%
  filter(sp_code2%in%sp_imprecise_sample)%>%
  dplyr::select(species)%>%
  distinct()%>%
  arrange(species)%>%
  pull()

spcolon_lag_list <- list()

for(i in 1:length(sp_imprecise_target)){
  
  first_year_old <- colon_year_df%>%
    filter(species==sp_imprecise_target[i])%>%
    filter(colon_estim_precision=="low")%>%
    dplyr::select(species, YearAb, Year,
           first_transect_sample,
           pre_colon_sample,
           colon_year, colon_less_yearab,
           colon_estim_precision,
           OldField, Transect,
           sp_code2)#yearxsp
  
  imprecise <- nrow(first_year_old)
  
  nofirst_df <- colon_year_df%>%
    filter(species==sp_imprecise_target[i]#&
           # OldField==unique(sp_field_target$OldField)&
           # Transect==unique(sp_field_target$Transect)
    )%>%#if possible add transect if else
    filter(colon_estim_precision=="high")%>%
    dplyr::select(species, YearAb, Year,
           first_transect_sample,
           pre_colon_sample,
           colon_year, colon_less_yearab,
           colon_estim_precision,
           OldField, Transect,
           sp_code2)
  
  precise <- nrow(nofirst_df)
  
  na_yearab <- colon_year_df%>%
    filter(species==sp_imprecise_target[i])%>%
    filter(is.na(YearAb))%>%
    nrow()
  
  colon_year_df%>%
    filter(species==sp_imprecise_target[i])%>%
    nrow()
  
  sum(imprecise, precise, na_yearab)
  
  if(precise>1){
    
    plot(density(nofirst_df$colon_less_yearab),col="red",
         main=sp_imprecise_target[i])
    
    max_dens <- max(density(nofirst_df$colon_less_yearab)$y)
    more_likely_gap <- modeest::mlv(nofirst_df$colon_less_yearab, na.rm=T)%>%
      round(., 0)
    
    ggplot(data=nofirst_df, aes(x=colon_less_yearab))+
      geom_density(color="#33CC00", lwd=1)+
      xlim(c(0, max(nofirst_df$colon_less_yearab)))+
      xlab("Time lag between field abandon and colonization")+
      ylab("Density")+
      ggtitle(sp_imprecise_target[i])+
      annotate("text", x = max(nofirst_df$colon_less_yearab)-15,
               y = max_dens-max_dens*0.1, 
               label = paste0("More likely time lag= ", more_likely_gap),
               family="serif")+
      annotate("text", x = max(nofirst_df$colon_less_yearab)-15,
               y = max_dens-max_dens*0.2, 
               label = paste0("Nº imprecise colon= ", imprecise),
               family="serif")+
      annotate("text", x = max(nofirst_df$colon_less_yearab)-15,
               y = max_dens-max_dens*0.3, 
               label = paste0("Nº precise colon= ", precise),
               family="serif")+
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5),
            text=element_text(family="serif"),
            axis.text = element_text(color="black"))} else { }
  
  spcolon_lag_df <- tibble(
    species=sp_imprecise_target[i])
  spcolon_lag_df <- spcolon_lag_df%>%
    mutate(colon_delay=case_when(
      precise>2 ~ more_likely_gap,
      precise==1~ nofirst_df$colon_less_yearab[1],
      precise==2~ round(mean(nofirst_df$colon_less_yearab),0),
      precise==0 ~ min(first_year_old$colon_less_yearab-(gap_average/2))
    ),
    colon_delay_mean=round(mean(nofirst_df$colon_less_yearab),0),
    min_delay_imprecise=min(first_year_old$colon_less_yearab-(gap_average/2)),
    n_imprecise=imprecise,
    n_precise=precise)#try not to use case when with a object not included in the table
  
  spcolon_lag_list[[i]] <- spcolon_lag_df
  # create a datatable with two columns species and more likely time lapse between abandon and colonization (mode)
  
}

sp_delay_df <- spcolon_lag_list%>%
  bind_rows()

length(sp_imprecise_target)
nrow(sp_delay_df)

hist(sp_delay_df$colon_delay)
sp_delay_df

#* add colon delay to year of abandon. In case of adding 
#* colon delay to YearAb lead to a estimate colonization year higher 
#* than the actual first observed year, we will use as estimated colon
#* date the colonization year minus half of the sample gap, as in rest
#* of situations, which is condidered as the most recent colon
#* date as possible in the rest of cases. 

colon_year_df1 <- colon_year_df%>%
  left_join(sp_delay_df%>%
              dplyr::select(species, colon_delay,
                     colon_delay_mean), by="species")%>%
  distinct()# spcolon_lag_df has more than one row

nrow(colon_year_df)
nrow(colon_year_df1)
names(colon_year_df1)

colon_year_df <- colon_year_df1
rm(colon_year_df1)

colon_year_df <- colon_year_df%>%
  mutate(YearAb_plus_delay1= YearAb+colon_delay,#mode
         YearAb_plus_delaymean= YearAb+colon_delay_mean)

colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  filter(YearAb_plus_delay1>=colon_year)%>%
  nrow()# almost third part of total imprecise cases will be supposed to have arrived recently to the plot in case of using the mode of delay

colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  filter(YearAb_plus_delaymean>=colon_year)%>%
  nrow()# almost half of the imprecise cases will be supposed to have arrived recently to the plot in case of using average delay

colon_year_df%>%
  filter(colon_estim_precision=="low")%>%
  nrow()

mrt_df <- colon_year_df%>%
  mutate(estim_colon_year1=case_when(#1 is mode based
    colon_year> first_transect_sample~# also !is.na(pre_colon_sample)
      round(pre_colon_sample+(colon_year-pre_colon_sample)/2, 0),
    (colon_estim_precision=="high"&
       colon_year==first_transect_sample)~
      round(YearAb+colon_less_yearab/2, 0),
    (colon_estim_precision=="low"&
       YearAb_plus_delay1>=colon_year)~
      round(colon_year-gap_average/2, 0),# as recent colon as average recent colon
    (colon_estim_precision=="low"&
       YearAb_plus_delay1<colon_year)~
      round(YearAb+colon_delay,0)),
    estim_colon_year2=case_when(# 2 is average based
      colon_year> first_transect_sample~# also !is.na(pre_colon_sample)
        round(pre_colon_sample+(colon_year-pre_colon_sample)/2, 0),
      (colon_estim_precision=="high"&
         colon_year==first_transect_sample)~
        round(YearAb+colon_less_yearab/2, 0),
      (colon_estim_precision=="low"&
         YearAb_plus_delaymean>=colon_year)~
        round(colon_year-gap_average/2, 0),# as recent colon as average recent colon
      (colon_estim_precision=="low"&
         YearAb_plus_delaymean<colon_year)~
        round(YearAb+colon_delay_mean,0))
  )

hist(mrt_df$estim_colon_year1)
max(mrt_df$estim_colon_year1, na.rm=T)

hist(mrt_df$estim_colon_year2)
max(mrt_df$estim_colon_year2, na.rm=T)

mrt_df <- mrt_df%>%
  mutate(mrt1=Year-estim_colon_year1,#mode
         mrt2=Year-estim_colon_year2)#mean

hist(mrt_df$mrt1)
hist(log(mrt_df$mrt1))
min(mrt_df$mrt1, na.rm=T)

hist(mrt_df$mrt2)
hist(log(mrt_df$mrt2))
min(mrt_df$mrt2, na.rm=T)


zd <- mrt_df%>%
  filter(mrt1==0)%>%
  dplyr::select(OldField,
         tpl_species,
         Origin, YearAb,
         first_transect_sample,
         Year, 
         colon_year,
         colon_delay,
         estim_colon_year1,
         estim_colon_year2,
         time_since_colon,
         pre_colon_sample,
         mrt1,
         mrt2)

zd%>%
  dplyr::select(OldField)%>%
  distinct()#601 field has mrt=0 since it was immediately sampled after abandonment

mrt_df <- mrt_df%>%
  rename(first_sp_record=colon_year)

nrow(mrt_df)
nrow(dis_df)

## join with dissimilarities df

dis_df1 <- dis_df%>%
  left_join(mrt_df%>%
              dplyr::select(sp_code, first_sp_record,
                     colon_delay, colon_delay_mean,
                     estim_colon_year1, 
                     estim_colon_year2, mrt1,mrt2), 
            by="sp_code")%>%
  distinct()

nrow(dis_df1)
nrow(dis_df)
names(dis_df1)  
dis_df <- dis_df1

rm(dis_df1)


###****add also burned category********##

sort(names(field_info))
dis_df <- dis_df%>%
  rename(Burned=Burned.)
unique(dis_df$Burned)

dis_df%>%
  filter(is.na(Burned))%>%
  dplyr::select(OldField, Transect)%>%
  distinct()

## it lacks OldFields 22,29,69

burn_df
burn_df1 <- burn_df%>%
  dplyr::select(-ExperimentNumber)%>%
  mutate(OldField=as.numeric(OldField),
         Transect=str_remove_all(Transect, "[ ]"))%>%
  mutate(transect_code2=paste0(OldField, "_",
                               Transect))

sort(unique(burn_df1$OldField))
sort(unique(dis_df$OldField))
sort(unique(burn_df1$Transect))

burn_df <- burn_df1
rm(burn_df1)

dis_df <- dis_df%>%
  left_join(burn_df%>%
              dplyr::select(transect_code2, 
                     BurnTreatment), by="transect_code2")

dis_df <- dis_df%>%
  mutate(burn=case_when(
    is.na(Burned)~BurnTreatment,
    !is.na(Burned)~Burned
  ))

unique(dis_df$burn)
dis_df <- dis_df%>%
  dplyr::select(-Burned, BurnTreatment)%>%
  mutate(BurnTreatment=burn)%>%
  dplyr::select(-burn)

## add time since first burning

dis_df <- dis_df%>%
  mutate(time_since_burn=Year-FirstYearAfterBurning)


###**** add melting time (time since first invasive arrives) ********##


#* check the origin of introduced/or native species using gbif distribution
#* and https://plants.usda.gov/home/plantProfile?symbol=POPR 
#* just making two vectors may work and then case_when %in%

## substitute the unknown status for the next species

dis_df <- dis_df%>%
  mutate(origin=case_when(
    species=="Aristida sp."~"Native",
    species=="Helianthus sp."~"Native",
    species=="Oenothera sp." ~ "Native",
    species=="Rhus sp." ~"Native",
    species=="Setaria sp." ~"Introduced",
    species=="Stachys sp." ~ "Native", 
    species=="Tradescantia sp." ~ "Native",
    !species%in%c("Aristida sp.",
                  "Helianthus sp.",
                  "Oenonthera sp.",
                  "Rhus sp.",
                  "Setaria sp.",
                  "Stachys sp.",
                  "Tradescantia sp.")~Origin
  ))

zo <- dis_df%>%
  filter(sp=="sp."&
           !Origin%in%c("Native", "Introduced"))%>%
  dplyr::select(species, Origin, origin)%>%
  distinct()

dis_df <- dis_df%>%
  mutate(origin2=case_when(
    origin=="Introduced"~"Introduced",
    origin!="Introduced"~"Rest"
  ),
  origin3=case_when(
    origin=="Native"~"Native",
    origin!="Native"~"Rest"
  ))# consider all unknown species as natives


#* get time since invasion -melting time - THis should be
#* a plot level value. But in case of invasive species it is not
#* clear to me whether a plot value could be used -ie. what means 
#* a melting time higher than mrt for invasive species, does the arrival
#* of other previous invasive species affect the ulterior arrival of
#* other introduced species?- 
#* I will use a plot level data with maximum time since first
#* invasive species appeared in the plot. Most cases have only 
#* one invasive species other have no invasive species,
#* so value will be NA or negative value as 0 
#* is value for recent invasions. 


transect_inv_time <- dis_df%>%
  filter(origin2=="Introduced")%>%
  group_by(transect_code)%>%
  summarise(melting_time1=max(mrt1, na.rm=T),
            melting_time2=max(mrt2, na.rm=T))

dis_df1 <- dis_df%>%
  left_join(transect_inv_time%>%
              dplyr::select(transect_code, 
                     melting_time1,
                     melting_time2), 
            by="transect_code")

nrow(dis_df1)
nrow(dis_df)
names(dis_df1)  
dis_df <- dis_df1

rm(dis_df1)


unique(is.na(dis_df$melting_time1))
unique(is.na(dis_df$mrt1))

ze <- dis_df%>%filter(is.na(melting_time1))%>%
  dplyr::select(OldField, transect_code,
         tpl_species, 
         first_sp_record,
         estim_colon_year1,
         origin,
         mrt1, 
         melting_time1)## some plots have no introduced species and have NA

#* replace NA out of field 600 and 601 - which has NA since 
#* they lack YearAb and were only sampled in 2016 -  by
#* - Inf has they have never been colonized by introduced species 

dis_df1 <- dis_df%>%
  mutate(melting_time1_inf=replace_na(
    melting_time1, -Inf
  ),
  melting_time2_inf=replace_na(
    melting_time2, -Inf))

dis_df1%>%
  filter(is.na(melting_time1_inf))%>%
  dplyr::select(OldField, 
         melting_time1_inf)%>%
  distinct()


nrow(dis_df1)
nrow(dis_df)

dis_df <- dis_df1
rm(dis_df1)

#*** add number of invasive species per plot ****

dis_df_ninv <- dis_df%>%
  group_by(origin2, transect_code)%>%
  summarise(n=n())

dis_df_ninv <- dis_df_ninv%>%
  pivot_wider(names_from=origin2,
              values_from=n)%>%
  rename(n_intro=Introduced,
         n_native=Rest)%>%
  replace_na(list(n_intro = 0,
                  n_native = 0))%>%
  mutate(intr_ratio=n_intro/(n_intro+n_native))

hist(dis_df_ninv$intr_ratio)

dis_df1 <- dis_df%>%
  left_join(dis_df_ninv, 
            by="transect_code")

nrow(dis_df1)
nrow(dis_df)
names(dis_df1)  
dis_df <- dis_df1

rm(dis_df1, dis_df_ninv)


#**** add phylogenetic similarity as variable which  has better ***
#**** distribution than dissimilarity ******

dis_df <- dis_df%>%
  mutate(phylo_w_sim=1000-phylo_w_dist)

#** scale mean cover from "tanto por 1" to percentage

dis_df <- dis_df%>%
  mutate(mean_cover_percent=mean_cover_percent*100)
  

##**** save def analytical dataframe ****##

# data.table::fwrite(dis_df,
#                    "../../output/tables/analyses/analysis_df_dissim_transect_ctr3d.csv",
#                    sep=";", dec=".", col.names=T, row.names=F, quote=F)

data.table::fwrite(dis_df,
                   "../../output/tables/analyses/analysis_df_dissim_transect_pair3d.csv",
                   sep=";", dec=".", col.names=T, row.names=F, quote=F)



####** end ***######
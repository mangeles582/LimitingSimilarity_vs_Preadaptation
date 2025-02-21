

### just for representation



# figure s1 ---------------------------------------------------------------


join_traits <- read_delim("../../output/tables/functional/clean_catford_other_databases1.csv", 
                          delim=";", col_names=T)# also tratis_df_list[[2]], but it may change if I add new tables


creek_df_join <- read_delim("../../output/tables/abundances/clean_creek_abundances.csv", 
                            delim=";", col_names=T)


## correct spelling errors and differences in each database

sort(setdiff(creek_df_join$species, join_traits$species))
sort(setdiff(join_traits$species, creek_df_join$species))
sort(unique(creek_df_join$species))
sort(unique(join_traits$species))

creek_df_join <- creek_df_join%>%
  mutate(species = as.factor(species),
         species = forcats::fct_recode(species,
                                       `Achillea millefolium lanulosa`="Achillea millefolium(lanulosa)",
                                       `Andropogon gerardii` =  "Andropogon gerardi",
                                       `Alnus sp.`= "alnus sp.",
                                       `Ambrosia artemisiifolia` ="ambrosia artemisifolia",
                                       `Artemisia campestris caudata`= "Artemisia caudata campestris",
                                       `Artemisia campestris caudata`= "Artemisia (caudata) campestris",
                                       `Botrychium sp.`="Botrichium sp.",
                                       `Carex sp.`= "carex sp.",
                                       `Carex sp.`= "Carex  sp.",
                                       `Andropogon gerardii`="Dead but still standing Andropogon",
                                       `Erigeron annuus` = "Erigeron anuus",
                                       `Euphorbia maculata`= "Euphorbia (supina) maculata",
                                       `Fragaria vesca`="Fragaria vesca (americana)",
                                       `Lactuca serriola scariola`= "Lactuca serriola (scariola)",
                                       `Polygonatum biflorum` = "Polygonatum canaliculatum",
                                       `Quercus ellipsoidalis` = "Quercus borealis-ellipsoidalis",
                                       `Rosa sp.` = "rosa sp.",
                                       `Rosa sp.` = "rosa sp",
                                       `Rubus sp.`= "rubus sp.",
                                       `Rudbeckia hirta` = "rudbekia hirta",
                                       `Plantago purshii patagonica`="Plantago (purshii) patagonica",
                                       `Tragopogon dubius major` = "Tragopogon dubius (major)",
                                       `Setaria lutescens glauca`="Setaria lutescens (glauca)",
                                       `Stipa sp.`="stipa sp.",
                                       `Solidago sp.`="solidago sp."))

join_traits <- join_traits%>%
  filter(!species%in%c("Carex pensylvanica",
                       "Carex spp."))%>%# are wrong writtend and only have na
  mutate(species = as.factor(species),
         species = forcats::fct_recode(species,
                                       `Achillea millefolium lanulosa`="Achillea millefolium(lanulosa)",
                                       `Ambrosia artemisiifolia` ="Ambrosia artemisiaefolia",
                                       `Andropogon gerardii` =  "Andropogon gerardi",
                                       `Artemisia campestris caudata`= "Artemisia (caudata) campestris",
                                       `Cyperus schweinitzii` = "Cyperus Schweinitzii",
                                       `Plantago purshii patagonica` = "Plantago (purshii) patagonica",
                                       `Quercus ellipsoidalis`= "Quercus elipsoidalis",
                                       `Plantago purshii patagonica`="Plantago (purshii) patagonica",
                                       `Setaria lutescens glauca`="Setaria lutescens (glauca)",
                                       `Tragopogon dubius major` = "Tragopogon dubius (major)",
         ))%>%
  distinct()# check if there are not repeated species due to NAs for some repeated species

unique(duplicated(join_traits$species))

dup_traits_sp <- join_traits%>%
  filter(duplicated(species))%>%
  select(species)

zz <- join_traits%>%
  filter(species%in%dup_traits_sp$species)%>%
  select(species, lifeform, SLA, LDMC,
         height, leaf_area, leaf_fresh_mass,
         leaf_dry_mass, seed_mass,
         trait_db)

## remove ducplicates with less data (averaging would mix different databases)

join_traits1 <- join_traits%>%
  filter(!((species=="Ambrosia artemisiifolia"&is.na(LDMC)|
              species=="Andropogon gerardii"&trait_db=="Other"|
              species=="Cyperus schweinitzii"& is.na(LDMC)|
              species=="Quercus ellipsoidalis" & is.na(LDMC))))

zzz <- join_traits1%>%
  filter(species%in%dup_traits_sp$species)%>%
  select(species, lifeform, SLA, LDMC,
         height, leaf_area, leaf_fresh_mass,
         leaf_dry_mass, seed_mass,
         trait_db)

nrow(join_traits)
nrow(join_traits1)# should be 4 species less


join_traits <- join_traits1
rm(join_traits1)

## add genus and epithet column to each database

creek_df_join <- creek_df_join%>%
  mutate(genus = word(catford_species, 1, sep = fixed(" ")),
         sp = word(catford_species, 2, sep = fixed(" ")),
         sub_sp= word(catford_species, 3, sep = fixed(" ")))

join_traits <- join_traits%>%
  mutate(genus = word(species, 1, sep = fixed(" ")),
         sp = word(species, 2, sep = fixed(" ")),
         sub_sp= word(species, 3, sep = fixed(" ")))


## join both databases


ab_traits_join <- creek_df_join%>%
  left_join(join_traits%>%select(-c(genus, sp, sub_sp, func_grp,
                                    duration, lifeform, pathway,
                                    origin)), 
            by="species")%>%
  distinct()

nrow(creek_df_join)
nrow(ab_traits_join)

sort(unique(ab_traits_join$species))



ab_na_traits <- creek_df_join%>%
  select(catford_species, species, genus, sp, sub_sp, Origin)%>%
  distinct()%>%
  left_join(join_traits%>%select(-c(genus, sp, sub_sp, func_grp,
                                    duration, lifeform, pathway,
                                    origin)),
            by=c("species"="species"))%>%
  select(c(catford_species, genus, sp, sub_sp,
           !contains("sd")))%>%
  #filter(rowSums(is.na(select(., contains("mean")))) != ncol(select(.,contains("mean"))))%>%
  arrange(catford_species)


ab_na_traits <- ab_na_traits%>%
  dplyr::filter(!c(genus%in%c("animal", "ant","Bare", "dung", 
                              "fern","Forb","fungi",
                              "Fungi", "hairy", 
                              "gopher","Grass","Leaves",
                              "Lichens", "Man", "miscellaneous",
                              "Miscellaneous", "Mosses", 
                              "mosses/lichens", "pinach",
                              "unknown", "Unknown",
                              "Woody", "woody")))

nrow(ab_na_traits)
table(ab_na_traits$Origin)

ab_na_traits2 <- creek_df_join%>%
  select(catford_species, genus, sp, sub_sp, Origin)%>%
  distinct()%>%
  left_join(join_traits%>%select(-c(genus, sp, sub_sp)),
            by=c("catford_species"="species"))%>%
  select(c(catford_species, genus, sp, sub_sp,
           !contains("sd")))%>%
  arrange(catford_species)


ab_na_traits2 <- ab_na_traits2%>%
  dplyr::filter(!c(genus%in%c("animal", "ant", "anthill",
                              "anthills", "Bare", "dung", 
                              "fern","Forb","fungi",
                              "Fungi", "hairy", "gooher",
                              "goopher", "gopher","Grass",
                              "Leaves","Lichens", "lichens",
                              "Man", "miscellaneous",
                              "Miscellaneous", "miscellaneous",
                              "Mosses", "mosses",
                              "mosses/lichens", "pinach", "Pine",
                              "unknown", "Unknown",
                              "Woody", "woody")))

trait_repres0 <- ab_na_traits2%>%
  #select(c(Origin, contains("mean")))%>%
  group_by(Origin)%>%
  summarise_all(funs(sum(!is.na(.))))%>%
  select(-c(catford_species, genus, sp, sub_sp, origin))


#trait_repres <- colSums(is.na(ab_na_traits[names(ab_na_traits)[grepl("mean", names(ab_na_traits))]]))#also

trait_repres1 <- trait_repres0%>%
  pivot_longer(!Origin, names_to = "trait_name", values_to = "sp_sum")%>%
  mutate(sp_percentage=sp_sum/nrow(ab_na_traits2))#%>%
# mutate(trait=case_when(
#   grepl("_wmean", trait_name) ~ substr(trait_name,1,  nchar(trait_name)-6),
#   grepl("_mean", trait_name) ~ substr(trait_name,1, nchar(trait_name)-5)
# ))


sp_repres0 <- ab_na_traits2%>%
  select(-c(genus, sp, sub_sp, Origin,
            func_grp, duration, lifeform,
            pathway, origin))%>%
  #distinct()%>%
  group_by(catford_species)%>%
  summarise_all(funs(sum(is.na(.))))%>%
  ungroup()%>%
  mutate(na_sum = rowSums(across(where(is.numeric))))%>%
  select(catford_species, na_sum)%>%
  mutate(trait_sum=(ncol(ab_na_traits)-5)-na_sum,
         trait_percentage= trait_sum/(ncol(ab_na_traits)-5))%>%
  rename(species=catford_species)%>%
  arrange(species)

sp_repres00 <- sp_repres0%>%
  left_join(ab_na_traits2%>%select(catford_species, Origin), 
            by=c("species"="catford_species"))

sp_repres0 <- sp_repres00

sp_repres <- table(sp_repres0$trait_sum)%>%
  data.frame()%>%
  rename(trait_sum=Var1,
         freq=Freq)

trait_repres1$trait <- factor(trait_repres1$trait_name, 
                              levels= c(as.vector(trait_repres1$trait_name))%>%unique())

trait_repres1 <- trait_repres1%>%
  mutate(sp_origin=case_when(
    Origin=="Introduced"~"Introduced",
    Origin=="Native"~"Native",
    !Origin%in%c("Introduced", "Native")~"Unknown"
  ))

(trait_level <- ggplot(trait_repres1%>%
                         filter(!trait_name%in%c(
                           "catford_species", "genus", "sp", "trait_db"
                         )),
                       aes(x=factor(trait), 
                           y=sp_percentage*100,
                           fill=sp_origin))+
    geom_bar(stat="identity", width=0.7)+
    coord_flip()+
    scale_fill_manual(values=c("#FDE725FF",
                               "#481567FF",
                               "#2D708EFF",
                               "#3CBB75FF"),
                      na.value = "transparent")+
    xlab("number of traits with data")+
    ylab("percentage of species")+
    theme_minimal()
)


# figure s2 --------------------------------------------------------


abun_df <- read_delim("../../output/tables/abundances/clean_creek_abundances.csv", 
                      delim=";", col_names=T)
abun_df <- abun_df%>%
  mutate(genus = word(species, 1, sep = fixed(" ")),
         sp1 = word(species, 2, sep = fixed(" ")),
         sub_sp= word(species, 3, sep = fixed(" ")))

misce <- abun_df%>%
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


abun_df<- abun_df%>%
  mutate(no_plant = case_when(
    species%in%c(misce$species)~"remove",
    !species%in%c(misce$species)~"keep"
  ))

abun_df <- abun_df%>%
  filter(no_plant=="keep")


cadotte_df <- read_delim("../../data/trait_pilot/other_data/Cadotte_et_al_2009_traits_E120.csv", 
                         delim=",", col_names=T)

cadotte_df <- cadotte_df%>%
  mutate(species=gsub("_", " ", Species))

cadotte_df <- cadotte_df%>%
  filter(species%in%abun_df$species)

nrow(cadotte_df)
unique(cadotte_df$species)# check if database is already species averages


e133_df <- read_delim("../../data/trait_pilot/other_data/E133_leaf data.csv", 
                      delim=",", col_names=T)

e133_df$species_comillas <- e133_df$species

e133_df <- e133_df%>%
  mutate(species=substr(species_comillas, 2, nchar(species_comillas)-1))%>%
  dplyr::select(-species_comillas)


nrow(e133_df)
sort(e133_df$species)# check if database is already species averages
e133_df[duplicated(e133_df$species),]# there is one repeated species

e133_df_sum <- e133_df%>%
  group_by(species)%>%
  select(-sp)%>%
  summarise_all(list(mean), na.rm=T)

nrow(e133_df_sum)
sort(e133_df_sum$species)# check if database is already species averages

e133_df <- e133_df_sum

e133_df <- e133_df%>%
  filter(species%in%abun_df$species)

craine_df <- read_delim("../../data/trait_pilot/other_data/Joe Craine plant trait data.csv", 
                        delim=",", col_names=T)# we do not need to join with this as there are no extra species

nrow(craine_df)
sort(unique(craine_df$Species))# check if database is already species averages

craine_df_sum <- craine_df%>%
  filter(Treatment=="LoN")%>%
  group_by(Species)%>%
  dplyr::select(-c(Tube, Functional,
                   Treatment, Harvest))%>%#there is no control treatment so we'll keep the low Nitrogen one
  summarise_all(list(mean), na.rm=T)

nrow(craine_df_sum)
sort(craine_df_sum$Species)# check if database is already species averages

craine_df <- craine_df_sum
craine_df <- craine_df%>%
  rename(species=Species)

craine_df <- craine_df%>%
  filter(species%in%abun_df$species)

master_df <- read_delim("../../data/trait_pilot/other_data/Mastertraits_2008.csv", 
                        delim=",", col_names=T)

master_df <- master_df%>%
  mutate(species=gsub("_", " ", name))

nrow(master_df)
sort(master_df$species)# check if database is already species averages

master_df <- master_df%>%
  filter(species%in%abun_df$species)

willis_df <- read_delim("../../data/trait_pilot/other_data/Willis_et_al_2010_traits_E133.csv", 
                        delim=",", col_names=T)

willis_df <- willis_df%>%
  mutate(species=gsub("_", " ", Genus_species))

nrow(willis_df)
sort(willis_df$species)# check if database is already species averages

willis_df <- willis_df%>%
  filter(species%in%abun_df$species)


names(cadotte_df)
names(e133_df)
names(master_df)
names(willis_df)

## homogeneize names and units with Catford

names(cadotte_df)

cadotte_df <- cadotte_df%>%
  mutate(height=`height_(m)`*100,#catford is in cm
         SLA= `SLA_cm2/g`*0.1,#catford is in mm2/mg
         leaf_area=Area_of_leaf_blade_cm2*100)%>%#catford is in mm2
  rename(PA=`P/A`,# catford has no perimeter
         PAL=`P_A*L`,
         seed_mass= `Seed_weight_(g)`,# catford is also in g
  )%>%
  mutate(species=gsub("_", " ", Species))%>%
  select(-c(`height_(m)`,
            `SLA_cm2/g`,
            Area_of_leaf_blade_cm2))

names(cadotte_df)


e133_df <- e133_df%>%
  rename(leaf_area=area)%>%
  select(-SLA_cm2_g)

names(e133_df)


names(willis_df)[2:12] <- substring(names(willis_df)[2:12], 5)

willis_df <- willis_df%>%
  mutate(SLA=`SLA (cm^2/g)`*0.1,# Catford is in mm2/mg
         LDMC=(`DLW (g)`*1000)/`FLW (g)`# estimate it from dry mass and fresh mass mg/g
  )%>%
  rename(fresh_mass_nopet = `FLW w/o petiole (g)`,
         fresh_mass =  `FLW (g)`,
         dry_mass_nopet =  `DLW w/o Petiole (g)`,
         dry_mass = `DLW (g)`,
         LWF = `LWF`,
         leaf_length = `Leaf Length (mm)`,
         leaf_width = `Leaf Width (mm)`,
         petiol_length = `Petiol Length (mm)`,
         leaf_area = `Area (mm^2)`,# Catford is in mm2 as well
         leaf_perim = `Perimeter (mm)`)%>%
  select(-`SLA (cm^2/g)`)

names(master_df)

master_df <- master_df%>%
  mutate(area=Area, #catford is in mm2, MAster metadata say it is in cm2 but actually it should be in mm2 (Area/Area_n)*100
         SLA = SLA_pet*0.1 , #catford is in mm2/mg 
         SLA_nopet = SLA_nopet*0.1, #catford is in mm2/mg 
         SLA_greater3 = SLA_greater3*0.1, #catford is in mm2/mg 
         height=height*100)%>%
  rename(log_height=logheight,
         root_depth=rootdpth,
         log_root_depth= logrootdpth,
         leaf_area= area, 
         log_leaf_area= logArea,
         leaf_area_n= Area_n,
         leaf_length= Length,
         leaf_length_n= Length_n,
         leaf_width =  Width,      
         leaf_width_n= Width_n,
         log_SLA_pet = logSLA_pet,
         log_PA =  logPA,
         log_PAL = logPLA,
         leaf_perim = perim,
         leaf_perim_n = perim_n,
         PAL_greater3 = PLA_greater3,
         PAL=PLA,
         dry_mass=drymass_pet,
         dry_mass_n= drymass_pet_n,
         dry_mass_nopet = drymass_nopet,
         dry_mass_nopet_n = drymass_nopet_n)%>%
  select(-c(Area, SLA_pet))


names(cadotte_df) <- paste0(names(cadotte_df), "_cadotte")
names(e133_df) <- paste0(names(e133_df), "_e133")
names(willis_df) <- paste0(names(willis_df), "_willis")
names(master_df) <- paste0(names(master_df), "_master")

cadotte_df <- cadotte_df%>%
  rename(species= species_cadotte)

e133_df <- e133_df%>%
  rename(species= species_e133)

willis_df <- willis_df%>%
  rename(species= species_willis)

master_df <- master_df%>%
  rename(species= species_master)

other_traits_list <- list(cadotte_df, e133_df, 
                          master_df, willis_df)

other_traits_df <- other_traits_list %>% 
  purrr::reduce(full_join, by = c("species"))

other_traits_df %>%
  select(where(is.character))

other_traits_df$petiol_length_willis <- as.numeric(other_traits_df$petiol_length_willis)

chr_vars <- other_traits_df %>%
  select(where(is.character))%>%
  select(-species)

other_traits_longer <- other_traits_df %>%
  select(-names(chr_vars))%>%
  pivot_longer(!species, names_to = "trait_db", values_to = "value")

study_table0 <- other_traits_longer%>%
  group_by(trait_db)%>%
  summarise(no_na_sp=sum(!is.na(value)),
            na_sp=sum(is.na(value)),
            n=n())%>%
  mutate(
    trait=case_when(
      grepl("cadotte", trait_db) ~ substring(trait_db,1, nchar(trait_db)-8),
      grepl("e133", trait_db) ~ substring(trait_db,1, nchar(trait_db)-5),
      grepl("willis", trait_db) ~ substring(trait_db, 1, nchar(trait_db)-7 ),
      grepl("master", trait_db) ~ substring(trait_db, 1, nchar(trait_db)-7 )
    ),
    dataset=case_when(
      grepl("cadotte", trait_db) ~ "cadotte",
      grepl("e133", trait_db) ~ "e133",
      grepl("willis", trait_db) ~ "willis" ,
      grepl("master", trait_db) ~ "master" 
    ))


other_traits_longersum <- other_traits_longer%>%
  group_by(trait_db)%>%
  summarise(mean_value=mean(value, na.rm=T))%>%
  mutate(
    trait=case_when(
      grepl("cadotte", trait_db) ~ substring(trait_db,1, nchar(trait_db)-8),
      grepl("e133", trait_db) ~ substring(trait_db,1, nchar(trait_db)-5),
      grepl("willis", trait_db) ~ substring(trait_db, 1, nchar(trait_db)-7 ),
      grepl("master", trait_db) ~ substring(trait_db, 1, nchar(trait_db)-7 )
    ),
    dataset=case_when(
      grepl("cadotte", trait_db) ~ "cadotte",
      grepl("e133", trait_db) ~ "e133",
      grepl("willis", trait_db) ~ "willis" ,
      grepl("master", trait_db) ~ "master" 
    ))%>%
  arrange(trait)

## compare databases value distribution for main traits
## compare also with catford database


catford_traits <- read_delim("../../output/tables/functional/e275_e93_trait_compilation.csv", 
                             delim=";", col_names=T)

sort(unique(other_traits_longer$trait_db))

catford_traits <- catford_traits%>%
  filter(species%in%abun_df$species)

other_traits_longer2 <- other_traits_longer%>%
  mutate( dataset=case_when(
    grepl("cadotte", trait_db) ~ "cadotte",
    grepl("e133", trait_db) ~ "e133",
    grepl("willis", trait_db) ~ "willis" ,
    grepl("master", trait_db) ~ "master" 
  ),
  trait=case_when(
    grepl("cadotte", trait_db)~ substring(trait_db,1, nchar(trait_db)-8),
    grepl("e133", trait_db)~ substring(trait_db,1, nchar(trait_db)-5),
    grepl("willis", trait_db)~ substring(trait_db,1, nchar(trait_db)-7),
    grepl("master", trait_db)~ substring(trait_db,1, nchar(trait_db)-7)
  ))%>%
  select(species, trait, dataset, value)%>%
  filter(trait%in%c("SLA", "seed_mass",
                    "leaf_area","LDMC",
                    "height", "fresh_mass",
                    "dry_mass"))%>%
  mutate(trait=case_when(
    trait=="fresh_mass" ~ "leaf_fresh_mass",
    trait=="dry_mass" ~ "leaf_dry_mass",
    !trait%in%c("fresh_mass", "dry_mass") ~ trait
  ))

head(catford_traits)
sort(names(catford_traits))

catford_traits2 <- catford_traits%>%
  # select(species, SLA_mean, seed_mass_wmean,
  #        LDMC_mean, leaf_dry_mass_mean,
  #        leaf_area_mean, leaf_fresh_mass_mean,
  #        height_mean)%>%
  select(-c(func_grp, duration, lifeform, pathway, origin,
            disp_length_mean, disp_length_sd, root_mass2_mean,
            root_depth2_mean, infl_only_mean))%>%
  select(!contains(c("_sd", "_std_")))%>%
  rename(root_mass_mean=root_mass1_mean,
         root_depth_mean=root_depth1_mean)%>%
  pivot_longer(!species, names_to = "trait_db", 
               values_to = "value")%>%
  mutate(trait=case_when(
    grepl("_mean", trait_db)~ substring(trait_db,1, nchar(trait_db)-5),
    grepl("_wmean", trait_db)~ substring(trait_db,1, nchar(trait_db)-6),
  ), 
  dataset="catford")%>%
  select(species, trait, dataset, value)

cat_study_table <- catford_traits2%>%
  group_by(trait)%>%
  summarise(no_na_sp=sum(!is.na(value)),
            na_sp=sum(is.na(value)),
            n=n())%>%
  mutate(dataset="Catford")

study_table <- study_table0%>%
  select(trait, no_na_sp,
         na_sp, n,
         dataset)%>%
  rbind(cat_study_table)

(study_n_sp <- ggplot(study_table%>%
                        filter(!is.na(trait))%>%
                        filter(!grepl("_n", trait))%>%
                        filter(!grepl("log_", trait))%>%
                        filter(!trait%in%c("uCC",
                                           "uSM",
                                           "uFF",
                                           "uNM"))%>%
                        mutate(trait=fct_reorder(trait, 
                                                 no_na_sp)),
                      aes(x=factor(trait), 
                          y=no_na_sp,
                          fill=dataset))+
    geom_bar(stat="identity", width=0.7)+
    coord_flip()+
    scale_fill_manual(values=c("#FDE725FF",
                               "#481567FF",
                               "#2D708EFF",
                               "#3CBB75FF",
                               "#1f9e89"),
                      na.value = "transparent")+
    xlab(" ")+
    ylab("number of species")+
    theme_minimal()
)


# rest of mod 7 figures ---------------------------------------------------

##*** not to run *** only for extracting AIC and BIC interaction model


# mod7_plot <- lmerTest::lmer(log(cover_percent+1)~
#                               time_since_ab:scale(dissim_SLA)*origin4+
#                               #scale(dissim_leaf_area)*origin4+
#                               time_since_ab:scale(dissim_LDMC)*origin4+
#                               time_since_ab:scale(dissim_height)*origin4+
#                               # log(dissim_leaf_dry_mass+1)*origin4+
#                               # log(dissim_leaf_fresh_mass+1)*origin4+
#                               time_since_ab:log(dissim_seed_mass+1)*origin4+
#                               time_since_ab:scale(phylo_w_dist)*origin4+
#                               time_since_ab:log(fun_w_dist)*origin4+
#                               scale(time_since_ab)*origin4+
#                               time_since_ab:log(mrt1)*origin4+
#                               time_since_ab:BurnTreatment*origin4+#
#                               time_since_ab:nitro_perc*origin4+
#                               #carbon_perc*origin4+
#                               time_since_ab:log(omatter_perc+1)*origin4+#
#                               time_since_ab:light_perc*origin4+
#                               #(1+scale(time_since_ab)|tpl_species)+
#                               (1|tpl_species)+
#                               (1|Year)+
#                               (1|OldField/Transect/Plot),#nested Year?
#                             data=plot_dis_df,
#                             #na.action="na.omit"
# )# filter(origin0!="Unknown")
# 



mod7_plot <- lmerTest::lmer(log(cover_percent+1)~
                              scale(dissim_SLA)*origin4+
                              #scale(dissim_leaf_area)*origin4+
                              scale(dissim_LDMC)*origin4+
                              scale(dissim_height)*origin4+
                              # log(dissim_leaf_dry_mass+1)*origin4+
                              # log(dissim_leaf_fresh_mass+1)*origin4+
                              log(dissim_seed_mass+1)*origin4+
                              scale(phylo_w_dist)*origin4+
                              log(fun_w_dist)*origin4+
                              scale(time_since_ab)*origin4+
                              log(mrt1)*origin4+
                              BurnTreatment*origin4+#
                              nitro_perc*origin4+
                              #carbon_perc*origin4+
                              log(omatter_perc+1)*origin4+#
                              light_perc*origin4+
                              #(1+scale(time_since_ab)|tpl_species)+
                              (1|tpl_species)+
                              (1|Year)+
                              (1|OldField/Transect/Plot),#nested Year?
                            data=plot_dis_df,
                            #na.action="na.omit"
)# filter(origin0!="Unknown")



plot_dis_df%>%
  dplyr::select(dissim_SLA,
                dissim_seed_mass,
                dissim_height,
                dissim_LDMC,
                # dissim_leaf_dry_mass,
                # dissim_leaf_fresh_mass,
                # dissim_leaf_area,
                time_since_ab,
                mrt1,
                #BurnTreatment,
                nitro_perc,
                #carbon_perc,
                omatter_perc,
                light_perc)%>%
  na.omit%>%
  cor()

vif <- check_collinearity(mod7_plot)#run it without interactions
vif
vif_df <- vif %>%as.data.frame()

summary(mod7_plot)
MuMIn::r.squaredGLMM(mod7_plot)
MuMIn::AICc(mod7_plot)
BIC(mod7_plot)


require(emmeans)
em_cat <- emmeans(mod7_plot,  c("origin4"))
contrast(em_cat,  'tukey')

em_int <- emmeans(mod7_plot,  c("origin4", "fun_w_dist"))
contrast(em_int, by="fun_w_dist")
contrast(em_int, by="origin4")

emt7 <- emtrends(mod7_plot, "origin4", var = "fun_w_dist")
emt7          ### estimated slopes of ZcNOF for each level of Old_Lure
pairs(emt7)# i hope it is not the same as contrast contrast(em_cat) 



expl_vars <- c("time_since_ab",
               "mrt1",
               "phylo_w_dist",
               "fun_w_dist",
               "nitro_perc",
               "omatter_perc",
               "light_perc",
               "dissim_SLA",
               "dissim_LDMC",
               "dissim_height",
               "dissim_seed_mass")

slope_df7_plot <- expl_vars %>% 
  purrr::map(function(x){
    emmeans::emtrends(mod7_plot, "origin4", var = x) %>% 
      emmeans::test()%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
      as_tibble() %>% 
      rename(slope = 2) %>%
      mutate(variable = x,
             ci_upper= slope+1.96*SE,
             ci_lower= slope-1.96*SE)
  }) %>% 
  bind_rows()


pairs_df7_plot <- expl_vars %>% 
  purrr::map(function(x){
    emmeans::emtrends(mod7_plot, "origin4", var = x) %>% 
      pairs()%>%
      as_tibble() %>% 
      rename(trend_dif = 2,# add also intro-native column in case they are reversed (native-intro) for some traits
             pvalue_dif=6) %>%
      mutate(variable = x)%>%
      dplyr::select(variable,
                    trend_dif,
                    pvalue_dif)
  }) %>% 
  bind_rows()

em_cat <- emmeans(mod7_plot,  c("origin4", "BurnTreatment"))
burn_treat <- contrast(em_cat, by="origin4", 'tukey')%>%
  as_tibble()%>%
  mutate(variable="BurnTreatment",
         ci_upper= estimate+1.96*SE,
         ci_lower= estimate-1.96*SE)%>%
  rename(slope=estimate)%>%
  dplyr::select(origin4,slope,
                SE,df,z.ratio,
                p.value, variable,
                ci_upper, ci_lower)

pairs_df7_burn <- em_cat%>%
  pairs()%>%
  as_tibble()%>%
  filter(contrast=="Introduced BurnTreatment1 - Rest BurnTreatment1")%>%
  rename(trend_dif = 2,# add also intro-native column in case they are reversed (native-intro) for some traits
         pvalue_dif=6) %>%
  mutate(variable="BurnTreatment")%>%
  dplyr::select(variable,
                trend_dif,
                pvalue_dif)

slope_df7_plot <- rbind(slope_df7_plot,
                        burn_treat)

slope_df7_pairs <- rbind(pairs_df7_plot,
                         pairs_df7_burn)

slope_df7_plot <- slope_df7_plot%>%
  mutate(ef_size=effectsize::z_to_d(z.ratio, 
                                    nrow(mod7_plot@frame))$d,#use nrow(mod@frame)instead nrow(data_table)in case na are authomatically removed in the model
         ef_ci_lower=effectsize::z_to_d(z.ratio, 
                                        nrow(mod7_plot@frame))$CI_low,
         ef_ci_upper=effectsize::z_to_d(z.ratio, 
                                        nrow(mod7_plot@frame))$CI_high,
         r2m=MuMIn::r.squaredGLMM(mod7_plot)[1],
         r2c=MuMIn::r.squaredGLMM(mod7_plot)[2])

slope_df7_plot <- slope_df7_plot%>%
  left_join(slope_df7_pairs, by="variable")


##************************************************

##** not to run** only to check interaction model AIC and BIC

# 
# mod7_trans <- lmerTest::lmer(log(cover_percent+1)~ 
#                                time_since_ab:scale(dissim_SLA)*origin4+
#                                #scale(dissim_leaf_area)*origin4+
#                                time_since_ab:scale(dissim_LDMC)*origin4+
#                                time_since_ab:scale(dissim_height)*origin4+
#                                # log(dissim_leaf_dry_mass+1)*origin4+
#                                # log(dissim_leaf_fresh_mass+1)*origin4+
#                                time_since_ab:log(dissim_seed_mass+1)*origin4+
#                                scale(phylo_w_dist)*origin4+
#                                time_since_ab:log(fun_w_dist)*origin4+
#                                time_since_ab:scale(time_since_ab)*origin4+
#                                time_since_ab:log(mrt1)*origin4+
#                                time_since_ab:BurnTreatment*origin4+#
#                                time_since_ab:nitro_perc*origin4+
#                                #carbon_perc*origin4+
#                                time_since_ab:log(omatter_perc+1)*origin4+#
#                                time_since_ab:light_perc*origin4+
#                                #(1+scale(time_since_ab)|tpl_species)+
#                                (1|tpl_species)+
#                                (1|Year)+
#                                (1|OldField/Transect),#nested Year?
#                              data=trans_dis_df,
#                              #na.action="na.omit"
# )# filter(origin0!="Unknown")
# 
# 

mod7_trans <- lmerTest::lmer(log(cover_percent+1)~ 
                               scale(dissim_SLA)*origin4+
                               #scale(dissim_leaf_area)*origin4+
                               scale(dissim_LDMC)*origin4+
                               scale(dissim_height)*origin4+
                               # log(dissim_leaf_dry_mass+1)*origin4+
                               # log(dissim_leaf_fresh_mass+1)*origin4+
                               log(dissim_seed_mass+1)*origin4+
                               scale(phylo_w_dist)*origin4+
                               log(fun_w_dist)*origin4+
                               scale(time_since_ab)*origin4+
                               log(mrt1)*origin4+
                               BurnTreatment*origin4+#
                               nitro_perc*origin4+
                               #carbon_perc*origin4+
                               log(omatter_perc+1)*origin4+#
                               light_perc*origin4+
                               #(1+scale(time_since_ab)|tpl_species)+
                               (1|tpl_species)+
                               (1|Year)+
                               (1|OldField/Transect),#nested Year?
                             data=trans_dis_df,
                             #na.action="na.omit"
)# filter(origin0!="Unknown")



trans_dis_df%>%
  dplyr::select(dissim_SLA,
                dissim_seed_mass,
                dissim_height,
                dissim_LDMC,
                dissim_leaf_dry_mass,
                dissim_leaf_fresh_mass,
                dissim_leaf_area,
                time_since_ab,
                mrt1,
                #BurnTreatment,
                nitro_perc,
                omatter_perc,
                light_perc)%>%
  na.omit%>%
  cor()

vif <- check_collinearity(mod7_trans)#run it without interactions
vif
vif_df <- vif %>%as.data.frame()

summary(mod7_trans)
MuMIn::r.squaredGLMM(mod7_trans)
MuMIn::AICc(mod7_trans)
BIC(mod7_trans)

require(emmeans)
em_cat <- emmeans(mod7_trans,  c("origin4"))
contrast(em_cat,  'tukey')

em_int <- emmeans(mod7_trans,  c("origin4", "fun_w_dist"))
contrast(em_int, by="fun_w_dist")
contrast(em_int, by="origin4")

emt7 <- emtrends(mod7_trans, "origin4", var = "phylo_w_dist")
emt7          ### estimated slopes of ZcNOF for each level of Old_Lure
pairs(emt7)# i hope it is not the same as contrast contrast(em_cat) 

expl_vars <- c("time_since_ab",
               "mrt1",
               "phylo_w_dist",
               "fun_w_dist",
               "nitro_perc",
               "omatter_perc",
               "light_perc",
               "dissim_SLA",
               "dissim_LDMC",
               "dissim_height",
               "dissim_seed_mass")


slope_df7_trans <- expl_vars %>% 
  purrr::map(function(x){
    emmeans::emtrends(mod7_trans, "origin4", var = x) %>% 
      emmeans::test()%>%#in case of not using test we miss the pvalue but gain the Confindence intervals 95%
      as_tibble() %>% 
      rename(slope = 2) %>%
      mutate(variable = x,
             ci_upper= slope+1.96*SE,
             ci_lower= slope-1.96*SE)
  }) %>% 
  bind_rows()


pairs_df7_trans <- expl_vars %>% 
  purrr::map(function(x){
    emmeans::emtrends(mod7_trans, "origin4", var = x) %>% 
      pairs()%>%
      as_tibble() %>% 
      rename(trend_dif = 2,# add also intro-native column in case they are reversed (native-intro) for some traits
             pvalue_dif=6) %>%
      mutate(variable = x)%>%
      dplyr::select(variable,
                    trend_dif,
                    pvalue_dif)
  }) %>% 
  bind_rows()

em_cat <- emmeans(mod7_trans,  c("origin4", "BurnTreatment"))
burn_treat <- contrast(em_cat, by="origin4", 'tukey')%>%
  as_tibble()%>%
  mutate(variable="BurnTreatment",
         ci_upper= estimate+1.96*SE,
         ci_lower= estimate-1.96*SE)%>%
  rename(slope=estimate)%>%
  dplyr::select(origin4,slope,
                SE,df,z.ratio,
                p.value, variable,
                ci_upper, ci_lower)

pairs_df7_burn <- em_cat%>%
  pairs()%>%
  as_tibble()%>%
  filter(contrast=="Introduced BurnTreatment1 - Rest BurnTreatment1")%>%
  rename(trend_dif = 2,# add also intro-native column in case they are reversed (native-intro) for some traits
         pvalue_dif=6) %>%
  mutate(variable="BurnTreatment")%>%
  dplyr::select(variable,
                trend_dif,
                pvalue_dif)

slope_df7_trans <- rbind(slope_df7_trans,
                         burn_treat)

slope_df7_pairs <- rbind(pairs_df7_trans,
                         pairs_df7_burn)

slope_df7_trans <- slope_df7_trans%>%
  mutate(ef_size=effectsize::z_to_d(z.ratio, 
                                    nrow(mod7_trans@frame))$d,#use nrow(mod@frame)instead nrow(data_table)in case na are authomatically removed in the model
         ef_ci_lower=effectsize::z_to_d(z.ratio, 
                                        nrow(mod7_trans@frame))$CI_low,
         ef_ci_upper=effectsize::z_to_d(z.ratio, 
                                        nrow(mod7_trans@frame))$CI_high,
         r2m=MuMIn::r.squaredGLMM(mod7_trans)[1],
         r2c=MuMIn::r.squaredGLMM(mod7_trans)[2])

slope_df7_trans <- slope_df7_trans%>%
  left_join(slope_df7_pairs, by="variable")




#**** regression plot figures *****##


all_vars <- c("cover_percent",
              "time_since_ab",
              "mrt1",
              "phylo_w_dist",
              "fun_w_dist",
              "nitro_perc",
              "omatter_perc",
              "light_perc",
              "dissim_SLA",
              "dissim_LDMC",
              "dissim_height",
              "dissim_seed_mass",
              "BurnTreatment",
              "origin4",
              "tpl_species",
              "Year",
              "OldField",
              "Transect",
              "Plot")



var_df <- tibble(
  con_var=c("fun_w_dist",
            "phylo_w_dist",
            "dissim_SLA",
            #"dissim_leaf_area",
            "dissim_LDMC",
            "dissim_height",
            "dissim_seed_mass",
            #"dissim_leaf_dry_mass",
            #"dissim_leaf_fresh_mass",
            "time_since_ab",
            "mrt1",
            "BurnTreatment",
            "nitro_perc",
            "omatter_perc",
            "light_perc"
  )
)

var_df <- var_df %>%
  mutate(cat_var="origin4",
         transform=case_when(
           con_var%in%c("fun_w_dist",
                        "mrt1")~"log",
           con_var%in%c("dissim_seed_mass",
                        "omatter_perc")~"log+1",
           con_var%in%c("phylo_w_dist",
                        "dissim_SLA",
                        #"dissim_leaf_area",
                        "dissim_LDMC",
                        "dissim_height",
                        "time_since_ab")~"scale",
         ))

con_list <- unique(var_df$con_var) 
model_list <- c(mod7_plot, mod7_trans)
eff_df_list <- c()
model_sum_list <- c()

for(i in 1:length(model_list)){
  
  effect_var_list <- list()
  
  for(j in 1:length(con_list)){
    
    var_df_sel <- var_df%>%
      filter(con_var==con_list[j])
    
    if(model_list[[i]]@call$data=="plot_dis_df"){
      dataset_mod=plot_dis_df
    }else if(model_list[[i]]@call$data=="trans_dis_df"){
      dataset_mod=trans_dis_df
    }
    
    eff_var <- eff_size_transf(variable_con=unique(var_df_sel$con_var), # variable continuous x variable
                               variable_cat=unique(var_df_sel$cat_var),
                               model=model_list[[i]],#model cat variable, same as real data grouping var
                               transformation=unique(var_df_sel$transform),
                               dataset=dataset_mod)
    
    eff_var <- eff_var%>%
      mutate(dataset=case_when(
        grepl("plot", summary(model_list[[i]])$call[3])~"Neighbourhood",
        grepl("trans", summary(model_list[[i]])$call[3])~"Site"
      ),
      var_name=unique(var_df_sel$con_var))
    
    eff_var$variable_trans <- as.numeric(eff_var$variable_trans)
    eff_var$variable <- as.numeric(eff_var$variable)
    
    effect_var_list[[j]] <- eff_var
    
  }
  
  eff_var_df <- bind_rows(effect_var_list)
  eff_df_list[[i]] <- eff_var_df
  
  
  (param_df <- parameters::parameters(model_list[[i]]))
  param_df <- data.frame(param_df)%>%filter(!is.na(t))
  
  (effect_size_df <- effectsize::t_to_r(param_df$t, 
                                        param_df$df_error))
  
  effect_size_df <- data.frame(cbind(param_df$Parameter,
                                     param_df$Coefficient,
                                     data.frame(effect_size_df)))
  names(effect_size_df)[1] <- "parameter"
  names(effect_size_df)[2] <- "coefficient"
  effect_size_df$r2m <- MuMIn::r.squaredGLMM(model_list[[i]])[1]
  effect_size_df$r2c <- MuMIn::r.squaredGLMM(model_list[[i]])[2]
  effect_size_df$AICc <- MuMIn::AICc(model_list[[i]])
  effect_size_df$pvalue <- round(summary(model_list[[i]])[[10]]%>%
                                   data.frame()%>%
                                   dplyr::select(5), digits=5)[[1]]
  
  effect_size_df <- effect_size_df%>%
    mutate(dataset=case_when(
      grepl("plot", summary(model_list[[i]])$call[3])~"Neighbourhood",
      grepl("trans", summary(model_list[[i]])$call[3])~"Site"
    ))
  
  model_sum_list[[i]] <- effect_size_df
  
  
}

eff_df_all <- bind_rows(eff_df_list)
row.names(eff_df_all) <- NULL

model_sum_df <- bind_rows(model_sum_list)

# prepare random df dataset

plot_random_df <- plot_dis_df%>%
  dplyr::select(c(all_of(all_vars), 
                  species,
                  plot_code))%>%
  drop_na()

plot_random_df$prediction <- predict(mod7_plot)

plot_random_df <- plot_random_df%>%
  group_by(origin4)%>%
  sample_frac(0.025)%>%# select 2.5% rows per each category
  ungroup()%>%
  rename(group=origin4)%>%
  mutate(BurnTreatment=as.numeric(BurnTreatment))%>%
  dplyr::select(species,
                plot_code,
                cover_percent,
                prediction,
                group,
                con_list)%>%
  rename(code=plot_code)%>%
  pivot_longer(!c(species,
                  code,
                  cover_percent,
                  prediction,
                  group),
               names_to="var_name",
               values_to="variable")%>%
  mutate(dataset="Neighbourhood")%>%
  as.data.frame()



trans_random_df <- trans_dis_df%>%
  dplyr::select(c(all_of(all_vars[!all_vars%in%'Plot']), 
                  species,
                  transect_code))%>%
  drop_na()

trans_random_df$prediction <- predict(mod7_trans)


trans_random_df <- trans_random_df%>%
  group_by(origin4)%>%
  sample_frac(0.25)%>%#select 25% rows of each category
  ungroup()%>%
  rename(group=origin4)%>%
  mutate(BurnTreatment=as.numeric(BurnTreatment))%>%# silent in case of no transforming the variables
  dplyr::select(species,
                cover_percent,
                prediction,
                transect_code,
                group,
                con_list)%>%
  rename(code=transect_code)%>%
  pivot_longer(!c(species,
                  code,
                  cover_percent,
                  prediction,
                  group),
               names_to="var_name",
               values_to="variable")%>%
  mutate(dataset="Site")%>%
  as.data.frame()


random_df <- rbind(plot_random_df,
                   trans_random_df)

head(random_df)
head(model_sum_df)
head(eff_df_all)


eff_df_all <- eff_df_all%>%
  mutate(dataset_x_var0= case_when(
    var_name=="fun_w_dist"~paste0("fun_dissim_", dataset),
    var_name=="phylo_w_dist"~paste0("phylo_dist_", dataset),
    var_name=="mrt1"~paste0("mrt_", dataset),
    var_name!=c("fun_w_dist",
                "phylo_w_dist",
                "mrt1")~paste0(var_name, "_", dataset)
  ),
  dataset_x_var=gsub("_", " ", dataset_x_var0))%>%
  dplyr::select(- dataset_x_var0)


random_df <- random_df%>%
  mutate(dataset_x_var0= case_when(
    var_name=="fun_w_dist"~paste0("fun_dissim_", dataset),
    var_name=="phylo_w_dist"~paste0("phylo_dist_", dataset),
    var_name=="mrt1"~paste0("mrt_", dataset),
    var_name!=c("fun_w_dist",
                "phylo_w_dist",
                "mrt1")~paste0(var_name, "_", dataset)
  ),
  dataset_x_var=gsub("_", " ", dataset_x_var0))%>%
  dplyr::select(- dataset_x_var0)


random_df <- random_df%>%
  rename(var_name1=var_name)%>%
  mutate(var_name=case_when(
    var_name1=="fun_w_dist"~"Functional multiv dissim",
    var_name1=="phylo_w_dist"~"Phylogenetic distance",
    var_name1=="mrt1"~"Minimum residence time",
    var_name1=="time_since_ab"~"Succession time",
    var_name1=="nitro_perc"~"N percentage",
    var_name1=="omatter_perc"~"Org matter percentage",
    var_name1=="light_perc"~"Light penetration",
    var_name1=="BurnTreatment"~"Burning",
    var_name1=="dissim_SLA"~"SLA dissim",
    var_name1=="dissim_LDMC"~"LDMC dissim",
    var_name1=="dissim_height"~"Height dissim",
    var_name1=="dissim_seed_mass"~"Seed mass dissim"))

eff_df_all <- eff_df_all%>%
  rename(var_name1=var_name)%>%
  mutate(var_name=case_when(
    var_name1=="fun_w_dist"~"Functional multiv dissim",
    var_name1=="phylo_w_dist"~"Phylogenetic distance",
    var_name1=="mrt1"~"Minimum residence time",
    var_name1=="time_since_ab"~"Succession time",
    var_name1=="nitro_perc"~"N percentage",
    var_name1=="omatter_perc"~"Org matter percentage",
    var_name1=="light_perc"~"Light penetration",
    var_name1=="BurnTreatment"~"Burning",
    var_name1=="dissim_SLA"~"SLA dissim",
    var_name1=="dissim_LDMC"~"LDMC dissim",
    var_name1=="dissim_height"~"Height dissim",
    var_name1=="dissim_seed_mass"~"Seed mass dissim"))



#*** facet wrap plots ***************

random_df <- random_df%>%
  mutate(log_cover=log(cover_percent+1))%>%
  filter(var_name!="Burning")%>%
  filter(dataset=="Neighbourhood"& log_cover<3|
           dataset=="Site"&log_cover<5)

eff_df_all <- eff_df_all%>%
  filter(var_name!="Burning")


eff_df_all_new <- eff_df_all%>%
  filter(var_name1!="omatter_perc"|
           var_name1=="omatter_perc"&
           variable<4)%>%
  filter(var_name1!="light_perc"|
           var_name1=="light_perc"&
           variable<=1)

var_list_plot <- c("Functional multiv dissim",
                   "Phylogenetic distance",
                   "SLA dissim",
                   "LDMC dissim", 
                   "Minimum residence time",
                   "Succession time")


eff_df_all_new <- eff_df_all_new%>%
  filter(!var_name%in%var_list_plot )

random_df <- random_df%>%
  filter(!var_name%in%var_list_plot )

random_df <- random_df%>%
  filter(var_name1!="omatter_perc"|
           var_name1=="omatter_perc"&
           variable<4)%>%
  filter(var_name1!="light_perc"|
           var_name1=="light_perc"&
           variable<=1)

lab_code <- rep(unique(random_df%>%
                         arrange(dataset_x_var)%>%
                         dplyr::select(var_name)%>%
                         pull(var_name)), each=2)

names(lab_code) <- sort(unique(random_df%>%
                                 dplyr::select(dataset_x_var)%>%
                                 pull(dataset_x_var)))


dat_text <- random_df%>%
  group_by(dataset_x_var)%>%
  summarise(minimum=min(variable),
            maximum=max(variable),
            amplitude=diff(range(variable)))%>%
  drop_na()%>%
  mutate(x_pos=minimum+0.05*amplitude,
         y_pos=case_when(
           grepl("Neigh", dataset_x_var)~2.5,
           grepl("Site", dataset_x_var)~4.15
         ))

panel_letters <- c("a)", "b)", "c)", "d)",
                   "e)", "f)", "g)", "h)",
                   "i)", "j)")
dat_text <- dat_text%>%
  cbind(as_tibble(panel_letters))


slope_df7_plot$dataset <- "Neighbourhood"
slope_df7_trans$dataset <- "Site"

slope_df7 <- rbind(slope_df7_plot,
                   slope_df7_trans)



eff_df_all_new1 <- eff_df_all_new%>%
  left_join(slope_df7%>%
              dplyr::select(origin4, p.value,
                            variable, dataset),
            by=c("group"="origin4",
                 "var_name1"="variable",
                 "dataset"="dataset"))

nrow(eff_df_all_new)
nrow(eff_df_all_new1)

eff_df_all_new <- eff_df_all_new1%>%
  mutate(significance=case_when(
    p.value<0.05~"signif",
    p.value>0.05~"no signif"
  ))



(regres7 <- ggplot(random_df)+
    geom_point(aes(x=variable,
                   y=log(cover_percent+1),#prediction, log(cover_percent+1)
                   color=group,
                   alpha=group),
               stroke=0.01)+
    scale_alpha_discrete(range = c(0.5, 0.08))+
    geom_line(data=eff_df_all_new,
              aes(x=variable,
                  y= fit,
                  color=group,
                  linetype=significance),
              size=0.9,alpha=1)+
    scale_linetype_manual(values=c("dashed",
                                   "solid"))+
    geom_ribbon(data=eff_df_all_new,
                aes(x=variable,
                    ymin=lower,
                    ymax=upper,
                    fill=group),
                alpha=0.3)+
    scale_color_manual(values=c("#CFA154",#C18633",
                                "#0E726A"), #5AB2A8"
                                labels=c("Introduced", "Native"))+#"#EEC141","#494884"
    scale_fill_manual(values=c("#CFA154",#C18633",
                               "#0E726A"), #5AB2A8"
                               labels=c("Introduced", "Native"))+ #"#EEC141","#494884"
    facet_wrap(~dataset_x_var, scales="free",
               ncol=2, strip.position = "bottom",
               labeller = labeller(dataset_x_var=lab_code))+
    geom_text(data=dat_text,
              mapping=aes(x=x_pos,
                          y=y_pos,
                          label=value),
              family="serif")+
    ylab ( expression("log(Cover percentage)" )) +
    xlab ( expression(" " )) +
    guides(color=F, size=F,alpha=F,linetype=F,
           fill = guide_legend(override.aes = list(alpha=1)))+
    labs(fill = "Origin")+
    theme_bw()+
    theme(legend.position="bottom",
          legend.title= element_text(size=8),
          legend.text = element_text(size=7),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          axis.text.x = element_text(size=9, color = "black"),
          axis.text.y = element_text(size=9, color = "black"),
          # axis.title.x = element_text(color="black", size=10 
          # ),
          axis.title.y = element_text(
            color="black", size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_blank(), 
          strip.text.x = element_text(size=10, color = "black"),
          strip.placement = "outside",
          text=element_text(family="serif"))
)



library(grid)
library(gridExtra)
d <- unique(random_df$dataset)
grid.table(t(d))
tt3 <- ttheme_default(
  core=list(bg_params = list(fill = "grey70", col="black")), 
  base_family = "serif")
grid.arrange(tableGrob(t(d), theme=tt3))

title_n <- tableGrob("Neighbourhood", theme=tt3)
title_s <- tableGrob("Site", theme=tt3)

title_n$widths <- unit(56, "mm")
title_s$widths <- unit(56, "mm")

hlay <- rbind(c(NA,1,NA,2,NA),
              c(3,3,3,3,3))

multiplot7 <- grid.arrange(title_n, title_s,
                           regres7, 
                           heights = unit(c(0.7,23),'cm'),
                           widths = unit(c(1,5.6,0.5, 5.6, 0.2), "cm"), 
                           #widths = unit(c(1.55,4.8,0.95,4.8,0.3), "cm"), http://127.0.0.1:45401/graphics/plot_zoom_png?width=1134&height=900
                           layout_matrix=hlay)


ggsave("../../results/figures/def0/mod7_regression_rawpoints_supplvars_notrans.pdf",
       plot=multiplot7,
       width = 13.1, height = 25, units="cm")#nopoints, rawpoints, prediction

ggsave("../../results/figures/def0/mod7_regression_rawpoints_supplvars_notrans.png",
       plot=multiplot7, 
       width = 13.1, height = 25, dpi=600,
       units="cm")
















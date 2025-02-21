


##****************************************##
##* Maria Angeles Pérez-Navarro
##* King's College University
##* AlienImpacts project
##* Statistic Analyses and definite figures
##* July 2022
##* **************************************##

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tidyverse)
library(ade4)
library(ape)
library(phytools)
library(MetBrewer)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(DHARMa)
library(ggpointdensity)
library(HH)
library(performance)
library(emmeans)
library(interactions)
#remotes::install_github("mastoffel/partR2") 
library(glmm.hp)
source('plot_required_functions.R')


#*******************************###
#***** 1 load data tables******####
#*******************************###

## plot table 

plot_ctr_df <- read_delim( "../../output/tables/analyses/analysis_df_plot_ctr3d.csv",
                         delim=";", col_names=T)


plot_univ_df <- read_delim("../../output/tables/dissim/dissimilarity_centroid_univariate_plot.csv",
                      delim=";", col_names = T)
 
names(plot_ctr_df)
names(plot_univ_df)

plot_ctr_df <- plot_ctr_df%>%
  mutate(richness=n_intro+n_native)

plot_ctr_df%>%
  filter(genus=="Quercus")

plot_univ_df%>%
  filter(genus=="Quercus")

plot_dis_df <- plot_ctr_df%>%
  left_join(plot_univ_df%>%
              dplyr::select(sp_code,
                            c(setdiff(names(plot_univ_df),
                                      names(plot_ctr_df)))),
            by="sp_code")

nrow(plot_ctr_df)
nrow(plot_univ_df)
nrow(plot_dis_df)

str(plot_dis_df)
plot_dis_df <- plot_dis_df%>%
  mutate(origin0=case_when(
    origin=="Native"~ "Native",
    origin=="Introduced" ~"Introduced",
    !origin%in%c("Introduced", "Native")~ "Unknown"
  ),
  omatter_perc1=case_when(
    omatter_perc<100~omatter_perc,
    omatter_perc>100~NA_real_
  ))%>%
  dplyr::select(-omatter_perc)%>%
  rename(omatter_perc=omatter_perc1)


unique(plot_dis_df$origin)
unique(plot_dis_df$origin0)

plot_dis_df%>%
  dplyr::select(species, origin0)%>%
  distinct()%>%
  group_by(origin0)%>%
  summarise(n=n())

plot_dis_df <- plot_dis_df%>%
  mutate(origin4=case_when(
    tpl_species=="Poa pratensis" ~ "Introduced",
    tpl_species!="Poa pratensis" ~ origin2
  ),
  origin5=case_when(
    tpl_species=="Poa pratensis" ~ "Native",
    tpl_species!="Poa pratensis" ~ origin3
  ),
  origin6=case_when(
    tpl_species=="Poa pratensis" ~ "Introduced",
    tpl_species!="Poa pratensis" ~ origin0 
  ),
  BurnTreatment=as.factor(BurnTreatment))


plot_dis_df_rep <- plot_dis_df%>%
  filter(repres=="Y")

nrow(plot_dis_df)
nrow(plot_dis_df_rep)

mycols <- colorRampPalette(brewer.pal(11, "Spectral"))(26)
mycols1 <- c(met.brewer("Lakota",26),met.brewer("Cross",0))
scales::show_col(c(mycols))
scales::show_col(c(mycols1))

ggplot(data=plot_dis_df, aes(x=longitude,
                             y=latitude,
                             color=as.factor(OldField),
                             shape=Transect))+
  geom_point()+
  scale_color_manual(values = mycols1)



## transect table


trans_ctr_df <- read_delim( "../../output/tables/analyses/analysis_df_dissim_ctr3d_transect.csv",
                            delim=";", col_names=T)


trans_univ_df <- read_delim("../../output/tables/dissim/dissimilarity_centroid_univariate_transect.csv",
                      delim=";", col_names = T)


names(trans_ctr_df)
names(trans_univ_df)

trans_ctr_df <- trans_ctr_df%>%
  mutate(richness=n_intro+n_native)

trans_dis_df <- trans_ctr_df%>%
  left_join(trans_univ_df%>%
              dplyr::select(sp_code,
                            c(setdiff(names(trans_univ_df),
                                      names(trans_ctr_df)))),
            by="sp_code")

nrow(trans_ctr_df)
nrow(trans_univ_df)
nrow(trans_dis_df)

str(trans_dis_df)
trans_dis_df <- trans_dis_df%>%
  mutate(origin0=case_when(
    origin=="Native"~ "Native",
    origin=="Introduced" ~"Introduced",
    !origin%in%c("Introduced", "Native")~ "Unknown"
  ),
  omatter_perc1=case_when(
    omatter_perc<100~omatter_perc,
    omatter_perc>100~NA_real_
  ))%>%
  dplyr::select(-omatter_perc)%>%
  rename(omatter_perc=omatter_perc1)


unique(trans_dis_df$origin)
unique(trans_dis_df$origin0)

trans_dis_df%>%
  dplyr::select(species, origin0)%>%
  distinct()%>%
  group_by(origin0)%>%
  summarise(n=n())

trans_dis_df <- trans_dis_df%>%
  mutate(origin4=case_when(
    tpl_species=="Poa pratensis" ~ "Introduced",
    tpl_species!="Poa pratensis" ~ origin2
  ),
  origin5=case_when(
    tpl_species=="Poa pratensis" ~ "Native",
    tpl_species!="Poa pratensis" ~ origin3
  ),
  origin6=case_when(
    tpl_species=="Poa pratensis" ~ "Introduced",
    tpl_species!="Poa pratensis" ~ origin0 
  ),
  BurnTreatment=as.factor(BurnTreatment))


trans_dis_df_rep <- trans_dis_df%>%
  filter(repres=="Y")

nrow(trans_dis_df)
nrow(trans_dis_df_rep)

mycols <- colorRampPalette(brewer.pal(11, "Spectral"))(26)
mycols1 <- c(met.brewer("Lakota",26),met.brewer("Cross",0))
scales::show_col(c(mycols))
scales::show_col(c(mycols1))

ggplot(data=trans_dis_df, aes(x=longitude,
                        y=latitude,
                        color=as.factor(OldField),
                        shape=Transect))+
  geom_point()+
  scale_color_manual(values = mycols1)

#* add extra information about
#* estim_colon_year1 was estimated with colon delay mode,
#* estim_colon_year2 was estimated with colon delay average,
#* mrt1 was estimated with colon delay mode,
#* mrt2 was estimated with colon delay average,
#* melting_time1 estimated with plot maximum mrt1,
#* melting_time2 estimated with plot maximum mrt2,
#* melting_time1_inf estimated with plot maximum mrt1 and substituting NA -absence of introduced species- as -Inf,
#* melting_time2_inf estimated with plot maximum mrt2 and substituting NA -absence of introduced species- as -Inf,
#* Origin Catford Origin Value -4 categories, I, N, N/I, Unknown, 
#* origin Catford Origin value changing some unknown genus according to NSR-USDA
#* origin2 Introduced vs rest of categories,
#* origin3 Native vs rest of categories,
#* origin0 Introduced/or Native and Unknown
#* Origin4 equal to origin2 but poa pratensis as introduced
#* origin5 equal to origin3 but poa pratensis as native
#* origin6 equal to origin0 but poa pratensis as introduced. It keeps 3 cat: intr, native, unknown

## see cover percent of introduced and natives per year

aa <- plot_dis_df%>%
  mutate(year=as.factor(Year))%>%
  group_by(origin4, year)%>%
  summarise(n=sum(cover_percent))

aa_wider <- aa%>%
  pivot_wider(names_from = origin4, values_from = n)

aa_wider <- aa_wider%>%
  mutate(intro_percent=(Introduced/(Rest+Introduced))*100)

mean(aa_wider$intro_percent)
sd(aa_wider$intro_percent)

bb <- plot_dis_df%>%
  dplyr::select(species, origin0)%>%
  distinct()%>%
  group_by(origin0)%>%
  summarise(n=n())

# check number of rows in case of discarding species with unknown origin

nrow(plot_dis_df)
plot_dis_df%>%
  filter(origin6!="Unknown")%>%
  nrow()
#106248-105436 we only lose 812/106248 rows 0.76%

nrow(trans_dis_df)
trans_dis_df%>%
  filter(origin6!="Unknown")%>%
  nrow()
#13588-13486 we only lose 102/13588 rows 0.75%

##****************************************************###
##***** 2 INTRODUCED vs NATIVE SPECIES ANALYSES******#### 
##****************************************************###

plot_dis_df_nona <- plot_dis_df%>%
  dplyr::select(cover_percent, time_since_ab,
                Duration, Lifeform,
                Pathway, mrt1,
                phylo_w_dist, fun_w_dist,
                BurnTreatment, carbon_perc,
                nitro_perc, omatter_perc,
                light_perc, intr_ratio,
                time_since_burn, melting_time1_inf,
                tpl_species, Year, Plot,
                OldField, Transect,
                origin0, origin2,
                origin3, origin4,
                origin5)%>%
  drop_na()


trans_dis_df_nona <- trans_dis_df%>%
  dplyr::select(cover_percent, time_since_ab,
                Duration, Lifeform,
                Pathway, mrt1,
                phylo_w_dist, fun_w_dist,
                BurnTreatment, carbon_perc,
                nitro_perc, omatter_perc,
                light_perc, intr_ratio,
                time_since_burn, melting_time1_inf,
                tpl_species, Year, 
                OldField, Transect,
                origin0, origin2,
                origin3, origin4,
                origin5)%>%
  drop_na()





##***Final selected model at plot scale **####

## Current code only include the finally selected model in the paper
## Visit Appendix S1 to check different analysed models as modifications of the 
## current one

###*** try not to run ***


mod7_plot <- lmerTest::lmer(log(cover_percent+1)~
                              origin4+
                              scale(dissim_SLA)+
                              #scale(dissim_leaf_area)+
                              scale(dissim_LDMC)+
                              scale(dissim_height)+
                              # log(dissim_leaf_dry_mass+1)+
                              # log(dissim_leaf_fresh_mass+1)+
                              log(dissim_seed_mass+1)+
                              scale(phylo_w_dist)+
                              log(fun_w_dist)+
                              scale(time_since_ab)+
                              log(mrt1)+
                              BurnTreatment+#
                              nitro_perc+
                              #carbon_perc+
                              log(omatter_perc+1)+#
                              light_perc+
                              
                              scale(dissim_SLA):origin4+
                              #scale(dissim_leaf_area):origin4+
                              scale(dissim_LDMC):origin4+
                              scale(dissim_height):origin4+
                              # log(dissim_leaf_dry_mass+1):origin4+
                              # log(dissim_leaf_fresh_mass+1):origin4+
                              log(dissim_seed_mass+1):origin4+
                              scale(phylo_w_dist):origin4+
                              log(fun_w_dist):origin4+
                              scale(time_since_ab):origin4+
                              log(mrt1):origin4+
                              BurnTreatment:origin4+#
                              nitro_perc:origin4+
                              #carbon_perc:origin4+
                              log(omatter_perc+1):origin4+#
                              light_perc:origin4+
                              #(1+scale(time_since_ab)|tpl_species)+
                              (1|tpl_species)+
                              (1|Year)+
                              (1|OldField/Transect/Plot),#nested Year?
                            data=plot_dis_df,
                            #na.action="na.omit"
)# filter(origin0!="Unknown")

plot_dis_df%>%
  dplyr::select(fun_w_dist,
                phylo_w_dist,
                dissim_SLA,
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
  distinct()%>%
  na.omit%>%
  cor()

vif <- check_collinearity(mod7_plot)#run it without interactions
vif
vif_df <- vif %>%as.data.frame()

summary(mod7_plot)

car::Anova(mod7_plot, type="II")
MuMIn::r.squaredGLMM(mod7_plot)
MuMIn::AICc(mod7_plot)
BIC(mod7_plot)
broom.mixed::glance(mod7_plot)
simulationOutput <- simulateResiduals(fittedModel = mod7_plot, n = 250)
plot(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
testZeroInflation(simulationOutput)
hist(resid(mod7_plot))


require(emmeans)
em_cat <- emmeans(mod7_plot,  c("origin4"))
contrast(em_cat,  'tukey')

em_int <- emmeans(mod7_plot,  c("origin4", "fun_w_dist"))
contrast(em_int, by="fun_w_dist")
contrast(em_int, by="origin4")

emt7 <- emtrends(mod7_plot, "origin4", var = "fun_w_dist")
emt7          
test(emt7)
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
                    z.ratio,
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
                z.ratio,
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



##***Final selected model at plot scale **####

## Current code only include the finally selected model in the paper
## Visit Appendix S1 to check different analysed models as modifications of the 
## current one

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
                              BurnTreatment*origin4+
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
  dplyr::select(fun_w_dist,
                phylo_w_dist,
                dissim_SLA,
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
car::Anova(mod7_trans, type="II")
MuMIn::r.squaredGLMM(mod7_trans)
MuMIn::AICc(mod7_trans)
BIC(mod7_trans)
broom.mixed::glance(mod7_trans)


hist(log(plot_dis_df$cover_percent))
plot(mod7_trans)
qqnorm(residuals(mod7_trans))
scatter.smooth(residuals(mod7_trans) ~ fitted(mod7_trans))

simulationOutput <- simulateResiduals(fittedModel = mod7_trans, n = 250)
plot(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
testZeroInflation(simulationOutput)
hist(resid(mod7_trans))


require(emmeans)
em_cat <- emmeans(mod7_trans,  c("origin4"))
contrast(em_cat,  'tukey')

em_int <- emmeans(mod7_trans,  c("origin4", "fun_w_dist"))
contrast(em_int, by="fun_w_dist")
contrast(em_int, by="origin4")


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
                    z.ratio,
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
                z.ratio,
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


# variables slope plot ----------------------------------------------------


slope_df7_plot <- slope_df7_plot%>%
  mutate(scale="Neighbourhood")

slope_df7_trans <- slope_df7_trans%>%
  mutate(scale="Site")

slope_df7 <- rbind(slope_df7_plot,
                   slope_df7_trans)

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

slope_df7 <- slope_df7%>%
  mutate(transform=case_when(
    variable%in%c("fun_w_dist",
                  "seed_mass",
                  "mrt1",
                  "omatter_perc")~"log",
    variable%in%c("phylo_w_dist",
                  "dissim_SLA",
                  "dissim_LDMC",
                  "dissim_height",
                  "time_since_ab")~"scale",
    variable%in%c("BurnTreatment",
                  "nitro_perc",
                  "light_perc")~"none"),## include rest of the variables
    significance=case_when(
      p.value<0.05~"signif",
      p.value>0.05~"no signif"),
    var_name=case_when(
      variable=="fun_w_dist"~"Y_Functional multiv dissim",
      variable=="phylo_w_dist"~"Z_Phylogenetic distance",
      variable=="mrt1"~"We_Colonization time",
      variable=="time_since_ab"~"Wa_Succession time",
      variable=="nitro_perc"~"Wc_N percentage",
      variable=="omatter_perc"~"Wb_Org matter percentage",
      variable=="light_perc"~"Wd_Light penetration",
      variable=="BurnTreatment"~"W_Burning",
      variable=="dissim_SLA"~"Xa_SLA dissim",
      variable=="dissim_LDMC"~"Xd_LDMC dissim",
      variable=="dissim_height"~"Xd_Height dissim",
      variable=="dissim_seed_mass"~"Xb_Seed mass dissim"),
    origin4 = as.character(origin4),
    Origin = case_when(
      origin4=="Introduced"~"Introduced",
      origin4=="Rest"~"Native"
    ),
    originxsig = case_when(
      significance == "no signif"~"no signif",
      significance != "no signif"~ origin4),
    scalexorign= paste0(scale, "_",
                        origin4)
  )


slope_df7 <- slope_df7%>%
  dplyr::select(scale, variable, 
                setdiff(names(slope_df7),
                        c("variable",
                          "scale")))

mycols <- colorRampPalette(brewer.pal(11, "BrBG"))(25)
scales::show_col(c(mycols))


# data.table::fwrite(slope_df,
#                    "../../results/mod_results.csv",
#                    sep=";", dec=".", col.names=T, row.names=F, quote=F)



dat_text <- tibble(
  label= c(paste0("R2=",
                  signif(unique(slope_df7_plot%>%
                                  dplyr::select(r2m)),3)%>%pull()),
           paste0("R2=",
                  signif(unique(slope_df7_trans%>%
                                  dplyr::select(r2m)),3)%>%pull())),
  scale=c("Neighbourhood",
          "Site"),
  x=c(0.015, -0.22),
  y=c(2.45, 2.45)
)


(trend_plot7 <- slope_df7%>%
    arrange(var_name)%>%
    ggplot(aes(x= slope, y = var_name 
               #color= Origin,
               #group=interaction(origin4,var_name)
    ))+
    geom_vline(xintercept= 0, 
               linetype=2, 
               color="grey30", 
               size=0.15)+
    geom_point(aes(#shape=scale, 
      fill = originxsig,
      color=Origin), #fill = originxsig
      size=1.3, #alpha=significance
      stroke=0.5,
      position = position_dodge(0.5, preserve="total"),
      shape=21
    )+
    #scale_shape_manual(values = c(22, 21,23)) +
    # scale_alpha_manual(values= c(0,1))+
    geom_errorbarh(aes(y = var_name, 
                       xmin = ci_lower,
                       xmax = ci_upper,
                       color=Origin),
                   height=0,
                   size=0.25,
                   position = position_dodge(0.5, preserve="total"),
                   # alpha=0.7
    )+
    scale_fill_manual(values= c("#C18633",
                                "white",
                                "#0E726A"#5AB2A8"
    ))+
    scale_colour_manual(values= c("#C18633",
                                  "#0E726A"#5AB2A8"
    ))+
    # geom_text(data=dat_text,
    #           mapping=aes(x=x, y=y, label=label),
    #           hjust=0, size=2.6, family="serif")+
    guides(fill=F, alpha=F)+
    facet_grid(~scale, scales="free")+#var_name~scale
    xlab("Coefficient")+
    ylab(" ")+
    theme_bw()+
    theme(axis.title.y=element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(color="black", size=8),
          axis.title.x = element_text(color="black", size=8),
          axis.ticks= element_line(size=0.4),
          legend.title= element_text(color="black", size=8),
          legend.text=element_text(color="black",size=7),
          legend.key.size = unit(0.4, "mm"),
          legend.key.height = unit(3.5, "mm"),
          legend.key.width = unit(1.2, "mm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          text=element_text(family="serif"),
          legend.box="vertical", 
          panel.border=element_rect(colour="black", size=0.25)
    ))



# ggsave("../../results/figures/def0/mod_slope.pdf",
#        plot=trend_plot,
#        width = 8, height = 8, units="cm")
# 
# ggsave("../../results/figures/def0/mod_slope.png",
#        plot=trend_plot,
#        width = 8, height = 8, dpi=600,
#        units="cm")
# 



# 2.7 effect size plots ---------------------------------------------------



(effsize_plot7 <- slope_df7%>%
   mutate(sig = symnum(pvalue_dif,corr = FALSE,na = FALSE,
                       cutpoints = c(0,0.05, 1),
                       symbols = c("*", " "))
   ) %>% 
   arrange(var_name)%>%
   ggplot(aes(x= ef_size, y = var_name,
              color= Origin
              #group=interaction(origin4,var_name)
   ))+
   geom_vline(xintercept= 0, 
              linetype=2, 
              color="grey30", 
              size=0.15)+
   geom_point(aes(#shape=scale, 
     fill = originxsig), 
     size=1.37, 
     stroke=0.6,
     position = position_dodge(0.6, preserve="total"),
     shape=21
   )+
   #scale_shape_manual(values = c(22, 21,23)) +
   # scale_alpha_manual(values= c(0,1))+
   geom_errorbarh(aes(y = var_name, 
                      xmin = ef_ci_lower,
                      xmax = ef_ci_upper,
                      color=Origin),
                  height=0,
                  size=0.25,
                  position = position_dodge(0.6, preserve="total"),
                  # alpha=0.7
   )+
   scale_fill_manual(values= c("#CFA154",#C18633",
                               "white",
                               "#0E726A"#5AB2A8""#5AB2A8"
   ))+
   scale_colour_manual(values= c("#CFA154",#C18633",
                                 "#0E726A"#5AB2A8"
   ))+
   geom_text(aes(x= -0.2, y = var_name,label=sig),
             vjust=0.7, size=3, family="serif",
             color="black")+
   # geom_text(data=dat_text,
   #           mapping=aes(x=x, y=y, label=label),
   #           hjust=0, size=2.6, family="serif")+
   guides(fill=F, alpha=F)+
   facet_grid(~scale, scales="free")+#var_name~scale
   xlab("Effect size")+
   ylab(" ")+
   theme_bw()+
   theme(axis.title.y=element_blank(),
         legend.position = "bottom",
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.text=element_text(color="black", size=9.5),
         axis.title.x = element_text(color="black", size=10),
         axis.ticks= element_line(size=0.4),
         legend.title= element_text(color="black", size=10),
         legend.text=element_text(color="black",size=9),
         legend.key.size = unit(0.4, "mm"),
         legend.key.height = unit(3.5, "mm"),
         legend.key.width = unit(1.2, "mm"),
         legend.margin=margin(0,0,0,0),
         legend.box.margin=margin(0,0,0,0),
         strip.text=element_text(color="black", size=11),
         text=element_text(family="serif"),
         legend.box="vertical", 
         panel.border=element_rect(colour="black", size=0.25)
         
   ))



# ggsave("../../results/figures/def0/mod_effect_size.pdf",
#        plot=effsize_plot,
#        width = 17.5, height = 12, units="cm")
# 
# ggsave("../../results/figures/def0/mod_effect_size.png",
#        plot=effsize_plot,
#        width = 17.5, height = 12, dpi=600,
#        units="cm")




# 2.7 regression plots ----------------------------------------------------

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
  effect_size_df$BIC <- BIC(model_list[[i]])
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

# load("save_mod_list.Rdata")
# save_mod_list[[13]] <- model_sum_list[[1]]%>%
#   tibble()%>%
#   mutate(model="mod_7")
# save_mod_list[[14]] <- model_sum_list[[2]]%>%
#   tibble()%>%
#   mutate(model="mod_7")
# save_mod_list
# save(save_mod_list,file= "save_mod_list.Rdata")
# 


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
         var_name1=="mrt1"~"Colonization time",
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
    var_name1=="mrt1"~"Colonization time",
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

random_df <- random_df%>%
  mutate(order_code=case_when(
    dataset_x_var=="fun dissim Neighbourhood"~"a",
    dataset_x_var=="fun dissim Site"~"b",
    dataset_x_var=="phylo dist Neighbourhood"~"c",
    dataset_x_var=="phylo dist Site"~"d",
    dataset_x_var=="dissim SLA Neighbourhood"~"e",
    dataset_x_var=="dissim SLA Site"~"f",
    dataset_x_var=="dissim LDMC Neighbourhood"~"g",
    dataset_x_var=="dissim LDMC Site"~ "h",
    dataset_x_var=="time since ab Neighbourhood"~"i",
    dataset_x_var=="time since ab Site"~"j",
    dataset_x_var=="mrt Neighbourhood"~"k",
    dataset_x_var=="mrt Site"~"l",
    dataset_x_var=="nitro perc Neighbourhood"~"m",
    dataset_x_var=="nitro perc Site"~"n"
  ))

eff_df_all <- eff_df_all%>%
  filter(var_name!="Burning")


eff_df_all_new <- eff_df_all%>%
  # filter(var_name1!="omatter_perc"|
  #          var_name1=="omatter_perc"&
  #          variable>0.5&
  #          variable<1.5)%>%
  filter(var_name1!="fun_w_dist"|
           var_name1=="fun_w_dist"&
           variable<7.5)

eff_df_all_new <- eff_df_all_new%>%
  mutate(order_code=case_when(
    dataset_x_var=="fun dissim Neighbourhood"~"a",
    dataset_x_var=="fun dissim Site"~"b",
    dataset_x_var=="phylo dist Neighbourhood"~"c",
    dataset_x_var=="phylo dist Site"~"d",
    dataset_x_var=="dissim SLA Neighbourhood"~"e",
    dataset_x_var=="dissim SLA Site"~"f",
    dataset_x_var=="dissim LDMC Neighbourhood"~"g",
    dataset_x_var=="dissim LDMC Site"~"h",
    dataset_x_var=="time since ab Neighbourhood"~"i",
    dataset_x_var=="time since ab Site"~"j",
    dataset_x_var=="mrt Neighbourhood"~"k",
    dataset_x_var=="mrt Site"~"l",
    dataset_x_var=="nitro perc Neighbourhood"~"m",
    dataset_x_var=="nitro perc Site"~"n"
  ))
  

var_list_plot <- c("Functional multiv dissim",
                   "Phylogenetic distance",
                   "SLA dissim",
                   "LDMC dissim", 
                   "Colonization time",
                   # "N percentage",
                   "Succession time")


eff_df_all_new <- eff_df_all_new%>%
  filter(var_name%in%var_list_plot )

random_df <- random_df%>%
  filter(var_name%in%var_list_plot )

lab_code <- rep(unique(random_df%>%
                         dplyr::select(var_name)%>%
                         pull(var_name)), each=2)

names(lab_code) <- sort(unique(random_df%>%
                                 dplyr::select(order_code)%>%
                                 pull(order_code)))

dat_text <- random_df%>%
  group_by(order_code)%>%
  summarise(minimum=min(variable),
            maximum=max(variable),
            amplitude=diff(range(variable)))%>%
  drop_na()%>%
  mutate(x_pos=minimum+0.05*amplitude,
         y_pos=case_when(
           order_code%in%c("a","c","e","g","i","k")~2.5,
           order_code%in%c("b","d","f","h","j","l")~4.15
         ),
         label=paste0(order_code, ")"))


(regres2 <- ggplot(random_df%>%
                     filter(order_code%in%c("a","b","c","d"))%>%
                     filter(var_name1!="fun_w_dist"|
                              var_name1=="fun_w_dist"&
                              variable<7.5))+
    geom_point(aes(x=variable,
                   y=log(cover_percent+1),#prediction, log(cover_percent+1)
                   color=group,
                   alpha=group),
               stroke=0.01)+
    scale_alpha_discrete(range = c(0.5, 0.08))+
    geom_line(data=eff_df_all_new%>%
                filter(order_code%in%c("a","b","c","d")),
              aes(x=variable,#eff_df_all$variable_trans
                  y= fit,
                  color=group),
              size=0.9,alpha=1)+
    geom_ribbon(data=eff_df_all_new%>%
                  filter(order_code%in%c("a","b","c","d")),
                aes(x=variable,#eff_df_all$variable_trans
                    ymin=lower,
                    ymax=upper,
                    fill=group),
                alpha=0.3)+
    scale_color_manual(values= c("#CFA154",#C18633",
                                 "#0E726A"#5AB2A8"
    ),labels=c("Introduced", "Native"))+#"#EEC141","#494884"
    scale_fill_manual(values= c("#CFA154",#C18633",
                                "#0E726A"#5AB2A8"
    ),labels=c("Introduced", "Native"))+ #"#EEC141","#494884"
    facet_wrap(~order_code, scales="free",
               ncol=2, strip.position = "bottom",
               labeller = labeller(order_code=lab_code))+
    geom_text(data=dat_text%>%
                filter(order_code%in%c("a","b","c","d")),
              mapping=aes(x=x_pos,
                          y=y_pos,
                          label=label),
              family="serif")+
    # facet_wrap(var_name~dataset, scales="free",
    #            ncol=2, strip.position = "bottom",
    #            labeller = labeller(order_code=lab_code))+
    # facet_wrap_custom(var_name~dataset, scales = "free",
    #                   strip.position="bottom",
    #                   ncol = 2, scale_overrides = list(
    #                     scale_override(3, scale_x_continuous(limits = c(-2.2,2.2))),
    #                     scale_override(4, scale_x_continuous(limits = c(-2.2,2.2)))
    #                     )
    #                   ))+
    ylab ( expression("log(Cover percentage)" )) +
    xlab ( expression(" " )) +
    # annotate("text", x=0.9-0.3, y=4+1.8, 
    #          label=paste0("R2=", signif(model_r2_pin[[1]],3)), 
    #          hjust=0, size=2.6, family="serif")+
    # annotate("text", x=0.9-0.3, y=3.5+1.8, 
    #          label=paste0("Pv*=", signif(model_statistic_pin%>%
    #                                        filter(variable=="origin4:light_perc")%>%
    #                                        dplyr::select(pv)%>%pull(),3)), 
    #          hjust=0, size=2.6, family="serif")+
    #xlim(0,1.1)+
    #xlim(0,1)+
    guides(color=F, size=F,alpha=F,
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

title_n$widths <- unit(56.5, "mm")
title_s$widths <- unit(56.5, "mm")

# title_n$widths <- unit(52, "mm")#no points
# title_s$widths <- unit(52, "mm")

hlay <- rbind(c(NA,1,NA,2,NA),
              c(3,3,3,3,3))

multiplot2 <- grid.arrange(title_n, title_s,
                           regres2, 
                           heights = unit(c(0.8,11),'cm'),
                        widths = unit(c(1,5.65,0.5, 5.65, 0.2), "cm"), 
                        #widths = unit(c(1.45,5.2,0.9, 5.2, 0.2), "cm"),# nopoints
                        layout_matrix=hlay)


ggsave("../../results/figures/def0/mod_regression_rawpoints_selvars_notrans.pdf",
       plot=multiplot2,
       width = 13.1, height = 12, units="cm")#width = 13

ggsave("../../results/figures/def0/mod_regression_rawpoints_selvars_notrans.png",
       plot=multiplot2, 
       width = 13.1, height = 12, dpi=600,#width = 13
       units="cm")#nopoints, rawpoints, prediction




dat_text <- random_df%>%
  group_by(order_code)%>%
  summarise(minimum=min(variable),
            maximum=max(variable),
            amplitude=diff(range(variable)))%>%
  drop_na()%>%
  mutate(x_pos=minimum+0.05*amplitude,
         y_pos=case_when(
           order_code%in%c("a","c","e","g","i","k")~2.5,
           order_code%in%c("b","d","f","h","j","l")~4.15
         ),
         label=case_when(
           order_code%in%c("a","b","c","d") ~"x",
           order_code=="e"~"a)",
           order_code=="f"~"b)",
           order_code=="g"~"c)",
           order_code=="h"~"d)",
           order_code=="i"~"e)",
           order_code=="j"~"f)",
           order_code=="k"~"g)",
           order_code=="l"~"h)"
         ))



(regres5 <- ggplot(random_df%>%
                     filter(order_code%in%c("e","f","g","h",
                                            "i","j","k","l"))%>%
                     filter(var_name1!="fun_w_dist"|
                              var_name1=="fun_w_dist"&
                              variable<7.5))+
    geom_point(aes(x=variable,
                   y=log(cover_percent+1),#prediction, log(cover_percent+1)
                   color=group,
                   alpha=group),
               stroke=0.01)+
    scale_alpha_discrete(range = c(0.5, 0.08))+
    geom_line(data=eff_df_all_new%>%
              filter(order_code%in%c("e","f","g","h",
                                     "i","j","k","l")),
              aes(x=variable,#eff_df_all$variable_trans
                  y= fit,
                  color=group),
              size=0.9,alpha=1)+
    geom_ribbon(data=eff_df_all_new%>%
                  filter(order_code%in%c("e","f","g","h",
                                         "i","j","k","l")),
                aes(x=variable,#eff_df_all$variable_trans
                    ymin=lower,
                    ymax=upper,
                    fill=group),
                alpha=0.3)+
    scale_color_manual(values= c("#CFA154",#C18633",
                                          "#0E726A"#5AB2A8"
    ),labels=c("Introduced", "Native"))+#"#EEC141","#494884"
    scale_fill_manual(values= c("#CFA154",#C18633",
                                         "#0E726A"#5AB2A8"
    ),labels=c("Introduced", "Native"))+ #"#EEC141","#494884"
    facet_wrap(~order_code, scales="free",
               ncol=2, strip.position = "bottom",
               labeller = labeller(order_code=lab_code))+
    geom_text(data=dat_text%>%
                filter(order_code%in%c("e","f","g","h",
                                       "i","j","k","l")),
              mapping=aes(x=x_pos,
                          y=y_pos,
                          label=label),
              family="serif")+
    ylab ( expression("log(Cover percentage)" )) +
    xlab ( expression(" " )) +
    guides(color=F, size=F,alpha=F,
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

multiplot5 <- grid.arrange(title_n, title_s,
                           regres5, 
                           heights = unit(c(0.7,20),'cm'),
                           widths = unit(c(1,5.65,0.5, 5.65, 0.2), "cm"), 
                           #widths = unit(c(1.45,5.2,0.9, 5.2, 0.2), "cm"),# nopoints
                           layout_matrix=hlay)


ggsave("../../results/figures/def0/mod_regression_rawpoints_selvars5_notrans.pdf",
       plot=multiplot5,
       width = 13.1, height = 21, units="cm")#width = 13

ggsave("../../results/figures/def0/mod_regression_rawpoints_selvars5_notrans.png",
       plot=multiplot5, 
       width = 13.1, height = 21, dpi=600,#width = 13
       units="cm")#nopoints, rawpoints, prediction




### **** end ****#####


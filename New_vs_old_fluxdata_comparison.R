#Comparison new vs old fluxes


#Description----
#this script is used to compare the new and old (wrong, interpolated) fluxes and assess the relative change in flux estimates (total.flux, LM.flux and HM.flux) as well as the relative change in best.flux used for the paper. 


#Packages------
renv::restore()
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
#clean global environment
rm(list=ls())
#1. CO2 comparison---- 

#Load new and old data:
new_co2<- read.csv("C:/Users/Miguel/Dropbox/GHG_bestflux_calculation/Results/co2_bestflux.csv")
old_co2<- read.csv("C:/Users/Miguel/Dropbox/Cabrera-Brufau_et_al_2026_source_data_and_code/0_sourcedata/co2_bestflux.csv")

#load auxfile
auxfile_co2<- read.csv("C:/Users/Miguel/Dropbox/GHG_bestflux_calculation/Auxfiles/co2_auxfile.csv")

#same columns? 
unique(names(new_co2)==names(old_co2))

#same observations?
unique(new_co2$UniqueID==old_co2$UniqueID)



#Pivot longer to combine
co2fluxes_new<- new_co2 %>% 
  select(UniqueID, best.model,MDF.lim,
         total.flux,total.flux.se,
         LM.flux, LM.flux.se,
         HM.flux, HM.flux.se) %>%
  #Add best.flux column
  mutate(best.flux=case_when(best.model=="total.flux"~total.flux,
                             best.model=="LM"~LM.flux,
                             best.model=="HM"~HM.flux,
                             best.model=="None appropriate"~NA_real_),
         best.flux.se=case_when(best.model=="total.flux"~total.flux.se,
                             best.model=="LM"~LM.flux.se,
                             best.model=="HM"~HM.flux.se,
                             best.model=="None appropriate"~NA_real_)) %>% 
  select(-best.model) %>% 
  left_join(auxfile_co2 %>% select(UniqueID,gas_analiser),by="UniqueID") %>% #add gas analiser for reference
  #Pivot longer to have flux and flux.se separated based on the different estimate: 
  pivot_longer(
    cols = -c(UniqueID, gas_analiser, MDF.lim),
    names_to = c("estimate", ".value"),
    names_pattern = "(.*)\\.(flux(?:\\.se)?)"
  )

co2fluxes_old<-old_co2 %>% 
  select(UniqueID, total.flux,LM.flux,HM.flux,best.model) %>%
  #Add best.flux column
  mutate(best.flux=case_when(best.model=="total.flux"~total.flux,
                             best.model=="LM"~LM.flux,
                             best.model=="HM"~HM.flux,
                             best.model=="None appropriate"~NA_real_)) %>% 
  select(-best.model) %>% 
  left_join(auxfile_co2 %>% select(UniqueID,gas_analiser),by="UniqueID") %>% #add gas analiser for reference
  rename(LM=LM.flux, HM=HM.flux, best=best.flux, total=total.flux) %>% 
  pivot_longer(cols = -c(UniqueID,gas_analiser), names_to = "estimate",values_to = "oldflux")


#Combine fluxes and compute differneces
co2newold<- co2fluxes_new %>% 
  left_join(co2fluxes_old, by=c("UniqueID","estimate","gas_analiser")) %>% 
  mutate(diff_abs=abs(flux)-abs(oldflux),
         diff_percent=diff_abs/abs(oldflux)*100,
         sig_diff=abs(diff_abs)>flux.se)



co2_summarydiff<- co2newold %>% 
  group_by(estimate,gas_analiser) %>% 
  summarise(avg_percent_diff=mean(diff_percent,na.rm=T),
            SD=sd(diff_percent,na.rm=T),
            q5=quantile(diff_percent,0.05,na.rm=T),
            q10=quantile(diff_percent,0.10,na.rm=T),
            q25=quantile(diff_percent,0.25,na.rm=T),
            q50=quantile(diff_percent,0.50,na.rm=T),
            q75=quantile(diff_percent,0.75,na.rm=T),
            q90=quantile(diff_percent,0.90,na.rm=T),
            q95=quantile(diff_percent,0.95,na.rm=T),
            n_fluxes=sum(!is.na(flux)),
            n_sig_diff=sum(sig_diff,na.rm=T)
  ) %>% 
  arrange(gas_analiser,estimate)

#Mean and SD of central 99% distribution
co2_summarydiff_99percent<- co2newold %>% 
  group_by(estimate,gas_analiser) %>% 
  filter(
    between(
      diff_percent,
      quantile(diff_percent, 0.005, na.rm = TRUE),
      quantile(diff_percent, 0.995, na.rm = TRUE)
    )
  ) %>%
  summarise(mean_99 = mean(diff_percent, na.rm = TRUE),
            sd_99= sd(diff_percent,na.rm=T),
            .groups = "drop")%>% 
  arrange(gas_analiser,estimate)


#Boxplots
co2newold %>% 
  select(estimate, diff_percent, gas_analiser) %>% 
  ggplot(aes(x=estimate, y=diff_percent, fill=estimate))+
  geom_boxplot(outliers = F)+
  scale_y_continuous(name="Relative difference new - old flux (%old flux)")+
  facet_wrap(facets=vars(gas_analiser))


#New vs old
co2newold %>% 
  filter(estimate=="best") %>% 
  filter(!is.na(diff_percent)) %>% 
  ggplot(aes(x=abs(flux), y=abs(oldflux), col=sig_diff))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(vars(gas_analiser))


#Percent difference vs absolute value
co2newold %>% 
  filter(estimate=="best") %>% 
  ggplot(aes(x=abs(flux), y=diff_percent, col=estimate))+
  geom_point()+
  # scale_y_log10()+
  scale_x_log10(breaks = c(0.01,0.1,1, 10, 100, 1000),
                labels = c(0.01,0.1,1, 10, 100, 1000))+
  facet_wrap(vars(gas_analiser),scales="free")


#Inspect incubations with substantial difference for best.flux (estimate used in paper). 
#We define substantial as those estimates whose differneces are greater than the estimate Standard error.

str(co2newold)
co2bestflux_5percentdiff <- co2newold %>%
  filter(
    estimate == "best",#estimate chosen as best for paper
    abs(diff_percent) > 5, #discrepancy of 5% from old 
    sig_diff==T #Only "significant" differences, those larger than estimate SE
  )



co2_inspected<- c("s1-ri-r2-4-b-d-10:01", #OK minor change, both estimates below minimum detectable flux (virtually zero-flux)
              "s2-ri-r2-9-v-d-11:31", #OK change is due to interpolation impacts
              "s1-ri-r2-11-v-d-11:50", #OK change is due to interpolation impacts
              "s1-ri-r2-8-v-d-11:05" #OK Minor difference, due to interpolation
              )

co2_bestmodel_different<- c("s3-ri-p2-14-o-d-13:53", #OK
                        "s2-ca-p1-3-v-t-09:06", #OK
                        "s3-da-r1-10-o-d-10:27",#OK LM current, noisy incubation (g.fact just above 4 now)
                        "s4-da-r2-6-v-t-07:46",  #OK LM current, noisy incubation
                        "s4-cu-r1-9-v-t-09:03", #OK LM current, noisy incubation
                        "s3-cu-p1-10-o-d-08:37", #OK LM current, noisy incubation
                        "s4-va-a1-11-b-d-12:08", #OK current HM flux is reasonable (small curvature, not noisy)
                        "s4-du-a1-3-o-d-08:43", #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        "s1-ca-p2-8-b-d-12:27", #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        "s3-ri-p2-9-b-d-12:09", #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        "s4-va-a1-10-b-d-11:54", #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        "s3-ri-r2-7-v-t-10:45",  #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        "s1-ri-r1-8-v-t-11:24", #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        "s1-ca-p2-14-v-d-14:44", #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        "s3-du-r1-13-b-d-10:27"  #OK LM current, both LM and HM could be, HM does not improve MAE by more than 5% now
                        )

co2bestflux_5percentdiff %>%
  filter(!UniqueID%in%c(co2_inspected,co2_bestmodel_different)) %>% 
  arrange(desc(abs(diff_percent))) %>% 
  head()









#2. CH4 comparison---- 

#Load new and old data:
new_ch4<- read.csv("C:/Users/Miguel/Dropbox/GHG_bestflux_calculation/Results/ch4_bestflux.csv")
old_ch4<- read.csv("C:/Users/Miguel/Dropbox/Cabrera-Brufau_et_al_2026_source_data_and_code/0_sourcedata/ch4_bestflux.csv")

#load auxfile
auxfile_ch4<- read.csv("C:/Users/Miguel/Dropbox/GHG_bestflux_calculation/Auxfiles/ch4_auxfile.csv")

#same columns? 
unique(names(new_ch4)==names(old_ch4))

#same observations?
unique(new_ch4$UniqueID==old_ch4$UniqueID)



#Pivot longer to combine
ch4fluxes_new<- new_ch4 %>% 
  select(UniqueID, best.model,MDF.lim,
         total.flux,total.flux.se,
         LM.flux, LM.flux.se,
         HM.flux, HM.flux.se) %>%
  #Add best.flux column
  mutate(best.flux=case_when(best.model=="total.flux"~total.flux,
                             best.model=="LM"~LM.flux,
                             best.model=="HM"~HM.flux,
                             best.model=="None appropriate"~NA_real_),
         best.flux.se=case_when(best.model=="total.flux"~total.flux.se,
                                best.model=="LM"~LM.flux.se,
                                best.model=="HM"~HM.flux.se,
                                best.model=="None appropriate"~NA_real_)) %>% 
  select(-best.model) %>% 
  left_join(auxfile_ch4 %>% select(UniqueID,gas_analiser),by="UniqueID") %>% #add gas analiser for reference
  #Pivot longer to have flux and flux.se separated based on the different estimate: 
  pivot_longer(
    cols = -c(UniqueID, gas_analiser, MDF.lim),
    names_to = c("estimate", ".value"),
    names_pattern = "(.*)\\.(flux(?:\\.se)?)"
  )

ch4fluxes_old<-old_ch4 %>% 
  select(UniqueID, total.flux,LM.flux,HM.flux,best.model) %>%
  #Add best.flux column
  mutate(best.flux=case_when(best.model=="total.flux"~total.flux,
                             best.model=="LM"~LM.flux,
                             best.model=="HM"~HM.flux,
                             best.model=="None appropriate"~NA_real_)) %>% 
  select(-best.model) %>% 
  left_join(auxfile_ch4 %>% select(UniqueID,gas_analiser),by="UniqueID") %>% #add gas analiser for reference
  rename(LM=LM.flux, HM=HM.flux, best=best.flux, total=total.flux) %>% 
  pivot_longer(cols = -c(UniqueID,gas_analiser), names_to = "estimate",values_to = "oldflux")


#Combine fluxes and compute differneces
ch4newold<- ch4fluxes_new %>% 
  left_join(ch4fluxes_old, by=c("UniqueID","estimate","gas_analiser")) %>% 
  mutate(diff_abs=abs(flux)-abs(oldflux),
         diff_percent=diff_abs/abs(oldflux)*100,
         sig_diff=abs(diff_abs)>flux.se)



ch4_summarydiff<- ch4newold %>% 
  group_by(estimate,gas_analiser) %>% 
  summarise(avg_percent_diff=mean(diff_percent,na.rm=T),
            SD=sd(diff_percent,na.rm=T),
            q5=quantile(diff_percent,0.05,na.rm=T),
            q10=quantile(diff_percent,0.10,na.rm=T),
            q25=quantile(diff_percent,0.25,na.rm=T),
            q50=quantile(diff_percent,0.50,na.rm=T),
            q75=quantile(diff_percent,0.75,na.rm=T),
            q90=quantile(diff_percent,0.90,na.rm=T),
            q95=quantile(diff_percent,0.95,na.rm=T),
            n_fluxes=sum(!is.na(flux)),
            n_sig_diff=sum(sig_diff,na.rm=T)
  ) %>% 
  arrange(gas_analiser,estimate)

#Mean and SD of central 99% distribution
ch4_summarydiff_99percent<- ch4newold %>% 
  group_by(estimate,gas_analiser) %>% 
  filter(
    between(
      diff_percent,
      quantile(diff_percent, 0.005, na.rm = TRUE),
      quantile(diff_percent, 0.995, na.rm = TRUE)
    )
  ) %>%
  summarise(mean_99 = mean(diff_percent, na.rm = TRUE),
            sd_99= sd(diff_percent,na.rm=T),
            .groups = "drop")%>% 
  arrange(gas_analiser,estimate)


#Boxplots
ch4newold %>% 
  select(estimate, diff_percent, gas_analiser) %>% 
  ggplot(aes(x=estimate, y=diff_percent, fill=estimate))+
  geom_boxplot(outliers = F)+
  scale_y_continuous(name="Relative difference new - old flux (%old flux)")+
  facet_wrap(facets=vars(gas_analiser))


#New vs old
ch4newold %>% 
  filter(estimate=="best") %>% 
  filter(!is.na(diff_percent)) %>% 
  ggplot(aes(x=abs(flux), y=abs(oldflux), col=sig_diff))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  scale_y_log10(breaks = c(0.01,0.1,1, 10, 100, 1000),
                labels = c(0.01,0.1,1, 10, 100, 1000))+
  scale_x_log10(breaks = c(0.01,0.1,1, 10, 100, 1000),
                labels = c(0.01,0.1,1, 10, 100, 1000))+
  facet_wrap(vars(gas_analiser))


#Percent difference vs absolute value
ch4newold %>% 
  filter(estimate=="best") %>% 
  ggplot(aes(x=abs(flux), y=diff_percent, col=estimate))+
  geom_point()+
  # scale_y_log10()+
  scale_x_log10(breaks = c(0.01,0.1,1, 10, 100, 1000),
                labels = c(0.01,0.1,1, 10, 100, 1000))+
  facet_wrap(vars(gas_analiser),scales="free")


#Inspect incubations with substantial difference for best.flux (estimate used in paper). 
#We define substantial as those estimates whose differneces are greater than the estimate Standard error.

str(ch4newold)
ch4bestflux_5percentdiff <- ch4newold %>%
  filter(
    estimate == "best",#estimate chosen as best for paper
    abs(diff_percent) > 5, #discrepancy of 5% from old 
    sig_diff==T #Only "significant" differences, those larger than estimate SE
  )



ch4_inspected<- c("s4-du-a1-6-b-t-09:32", #OK current LM appropriate, old flux was result of interpolation
              "s3-da-r1-13-o-d-10:54", #OK current LM appropriate, old flux was result of interpolation
              "s2-cu-p1-6-v-d-08:45" #OK current HM appropriate, difference due to interpolation
              )

ch4_bestmodel_different<- c("s2-ri-p2-3-v-t-09:09",  #OK current LM, noisy incubation 
                        "s1-cu-r1-15-v-d-09:44", #OK, current HM, both could be but HM fits best now
                        "s4-ri-p1-4-v-d-07:43", #OK, current LM, noisy incubation, HM does not improve MAE >5%
                        "s4-du-a1-1-o-d-08:13", #OK, current LM, more appropriate
                        "s3-va-p2-11-v-d-11:36", #OK, current LM, more appropriate
                        "s3-va-p2-9-b-t-10:37", #OK, current HM, both could be but HM fits best now
                        "s2-va-p1-5-o-d-09:54", #OK, current LM more appropriate
                        "s3-ri-r2-17-b-t-13:27", #OK current LM more appropriate
                        "s4-ri-a1-15-o-d-09:39",  #OK current LM more appropriate
                        "s1-ri-r2-13-b-t-12:15", #OK, current HM, both could be but HM fits best now
                        "s1-va-p2-1-o-d-08:47", #OK, current HM, both could be but HM fits best now
                        "s3-ri-r2-8-v-t-11:03", #OK, current LM more appropriate
                        "s1-va-p2-1-o-d-08:47", #OK, current HM, both could be but HM fits best now
                        "s3-ri-r2-8-v-t-11:03", #OK, current LM more appropriate
                        "s3-ri-r2-6-b-d-10:30", #OK, current LM more appropriate
                        "s3-ri-r2-16-b-d-13:20",#OK, current LM more appropriate
                        "s2-va-a1-9-o-d-12:42" #OK, current HM more appropriate 
                        )


ch4bestflux_5percentdiff %>%
  filter(!UniqueID%in%c(ch4_inspected,ch4_bestmodel_different)) %>% 
  arrange(desc(abs(diff_percent))) %>% 
  head()


#3.Table and Figures comparison-----
#Here_----------

#Two points to make: 
  #1. Correction is needed: show that there are differences old vs new and that they are not random (los gatos biass).
  #2. Corrected data is generally similar: show that most fluxes are within 5% of old estimate

#Although the differences caused by the correction are small in the overwhelming mayority of cases, there were systematic biasses for some analisers, which could have propagated into meaningful impact on the paper results. 

#The differences in sourcedata is not the main point, what they will be interested in is potential differences in paper results, Keep rawdata differences brief. 






#X percent of CO2 fluxes retained best-flux estimate within 5% of old. 


#Re-do summary of percent change, leaving only relevant fluxes and differences: 
#abs difference is larger than new estimate SE, we focus only on "significant differences". 
#new flux is above MDF, this reduces the number of inflated percent differences (no point comparing fluxes that are "below detection",i.e. non-different than cero based on our method sensitivity)

#Summary percent difference of best-flux estimate for incubations:
 #1. that have fluxes larger than the minimum detectable flux (newflux > MDF.lim)
 #2. whose absolute difference is larger than flux uncertainty, those marked as sig_diff==T (due to diff_absolute> newflux.se)

##ch4-summary----
ch4_difference_relevance<- ch4newold %>% 
  filter(estimate=="best") %>% 
  group_by(gas_analiser) %>%
  summarise(total_incubations=n(),
            incub_with_estimate=sum(!is.na(flux)),
            incub_below_MDF=sum(flux<=MDF.lim, na.rm = T),
            incub_nosigdiff=sum(diff_abs<=flux.se,na.rm=T),
            incub_sigdiff=sum(flux>MDF.lim & diff_abs>flux.se,na.rm = T)
            )

ch4_difference_relevance

#Total % of incubations with flux estimate that retain a value that is within 5% of old estimate: 
ch4_within5percent <- ch4newold %>% 
  filter(estimate=="best",
         !is.na(flux)) %>% 
  summarise(n_incub_with_estimate=n(),
            n_newflux_within5percentofold=sum(abs(diff_percent)<=5),
            percent_newflux_within5percentofold=n_newflux_within5percentofold/n_incub_with_estimate*100)

#Reasons for larger than 5% difference
ch4_largerthan5percent<-ch4newold %>% 
  filter(estimate=="best",
         !is.na(flux),
         abs(diff_percent)>5) %>% 
  summarise(n_total=n(), 
            n_MDF_or_nonsigdiff=sum(abs(flux)<=MDF.lim|!sig_diff),
            n_true_different=n_total-n_MDF_or_nonsigdiff)

#Written summary: considering all incubations, not only those included in the paper (but also those considered not-relevant for comparisons and excluded in the preliminary filtering of data)

paste0("Of the ", 
       ch4_within5percent$n_incub_with_estimate, 
       " incubations with a valid CH4 flux estimate across the full dataset, ",
       ch4_within5percent$n_newflux_within5percentofold,
       " (",
       round(ch4_within5percent$percent_newflux_within5percentofold,2),
       "%) had a new CH4 flux estimate that lied within 5% of their old estimate value. Of the ",
       ch4_within5percent$n_incub_with_estimate-ch4_within5percent$n_newflux_within5percentofold,
       " incubations with larger than 5% difference, ",
       ch4_largerthan5percent$n_MDF_or_nonsigdiff, 
       " either had flux values below the Minimum Detectable Flux (MDF) or had a flux Standard Error greater than the absolute flux difference between the new and old estimates. Thus, only ",
       ch4_largerthan5percent$n_true_different,
       " cases (",
       round(ch4_largerthan5percent$n_true_different/ch4_within5percent$n_incub_with_estimate,2),
       "%) had new flux estimates that were meaningfully different from their old ones.")


##CH4 new vs old biplot-----
ch4newold %>% 
  filter(estimate=="best",
         !is.na(flux)) %>% 
  mutate(incubtype=case_when(abs(diff_percent)<=5~"New within 5% of old",
                             (abs(flux)>MDF.lim)&flux.se>abs(diff_abs)~"Difference within SE",
                             abs(flux)<=MDF.lim~"Below MDF",
                             TRUE~"Meaningfull >5% difference"
                             )) %>% 
  arrange(desc(incubtype)) %>% 
  ggplot(aes(x=abs(oldflux),y=abs(flux), col=incubtype))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1)+
  scale_color_manual(name="Incubation Case",
                     values = c("New within 5% of old"="grey",
                                "Difference within SE"="lightblue",
                                "Below MDF"="lightgreen",
                                "Meaningfull >5% difference"="red"))+
  scale_x_log10(name="Old estimate (abs. value)",limits=c(2e-5,2e4),
                breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000),
                labels=label_log(base = 10))+
  scale_y_log10(name="New estimate (abs. value)", limits=c(2e-5,2e4),
                breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000),
                labels=label_log(base = 10))+
  theme_classic()+
  ggtitle("New vs old absolute CH4 flux best estimates")

#Boxplots of percent differnece, containing only cases with relevant differences: flux>MDF and diff_abs>flux.se
ch4newold %>% 
  filter(estimate=="best",
         !is.na(flux),
         flux>MDF.lim,
         diff_abs>flux.se) %>% 
ggplot(aes(x=gas_analiser, y=diff_percent))+
  geom_boxplot(outliers = F)


#Do not separate gas analisers, only best.flux, to include in response to reviewers 
#create a single table with average percent differences mean sd, quantiles
#boxplot percent difference (no outliers) only for best.flux 
#old vs new fluxes double log scales 1:1 line, only for best flux



##co2-summary----
co2_difference_relevance<- co2newold %>% 
  filter(estimate=="best") %>% 
  group_by(gas_analiser) %>%
  summarise(total_incubations=n(),
            incub_with_estimate=sum(!is.na(flux)),
            incub_below_MDF=sum(flux<=MDF.lim, na.rm = T),
            incub_nosigdiff=sum(diff_abs<=flux.se,na.rm=T),
            incub_sigdiff=sum(flux>MDF.lim & diff_abs>flux.se,na.rm = T)
  )

co2_difference_relevance

#Total % of incubations with flux estimate that retain a value that is within 5% of old estimate: 
co2_within5percent <- co2newold %>% 
  filter(estimate=="best",
         !is.na(flux)) %>% 
  summarise(n_incub_with_estimate=n(),
            n_newflux_within5percentofold=sum(abs(diff_percent)<=5),
            percent_newflux_within5percentofold=n_newflux_within5percentofold/n_incub_with_estimate*100)

#Reasons for larger than 5% difference
co2_largerthan5percent<-co2newold %>% 
  filter(estimate=="best",
         !is.na(flux),
         abs(diff_percent)>5) %>% 
  summarise(n_total=n(), 
            n_MDF_or_nonsigdiff=sum(abs(flux)<=MDF.lim|!sig_diff),
            n_true_different=n_total-n_MDF_or_nonsigdiff)

#Written summary: considering all incubations, not only those included in the paper (but also those considered not-relevant for comparisons and excluded in the preliminary filtering of data)

paste0("Of the ", 
       co2_within5percent$n_incub_with_estimate, 
       " incubations with a valid co2 flux estimate across the full dataset, ",
       co2_within5percent$n_newflux_within5percentofold,
       " (",
       round(co2_within5percent$percent_newflux_within5percentofold,2),
       "%) had a new co2 flux estimate that lied within 5% of their old estimate value. Of the ",
       co2_within5percent$n_incub_with_estimate-co2_within5percent$n_newflux_within5percentofold,
       " incubations with larger than 5% difference, ",
       co2_largerthan5percent$n_MDF_or_nonsigdiff, 
       " either had flux values below the Minimum Detectable Flux (MDF) or had a flux Standard Error greater than the absolute flux difference between the new and old estimates. Thus, only ",
       co2_largerthan5percent$n_true_different,
       " cases (",
       round(co2_largerthan5percent$n_true_different/co2_within5percent$n_incub_with_estimate,2),
       "%) had new flux estimates that were meaningfully different from their old ones.")


#PLot new vs old absolute values, coloring those with:
#flux<= MDF
#flux.se > diff_abs 
##CO2 new vs old biplot-----
co2newold %>% 
  filter(estimate=="best",
         !is.na(flux)) %>% 
  mutate(incubtype=case_when(abs(diff_percent)<=5~"New within 5% of old",
                             (abs(flux)>MDF.lim)&flux.se>abs(diff_abs)~"Difference within SE",
                             abs(flux)<=MDF.lim~"Below MDF",
                             TRUE~"Meaningfull >5% difference"
  )) %>% 
  arrange(desc(incubtype)) %>% 
  ggplot(aes(x=abs(oldflux),y=abs(flux), col=incubtype))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1)+
  scale_color_manual(name="Incubation Case",
                     values = c("New within 5% of old"="grey",
                                "Difference within SE"="lightblue",
                                "Below MDF"="lightgreen",
                                "Meaningfull >5% difference"="red"))+
  scale_x_log10(name="Old estimate (abs. value)",limits=c(5e-5,1e3),
                breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000),
                labels=label_log(base = 10))+
  scale_y_log10(name="New estimate (abs. value)", limits=c(5e-5,1e3),
                breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000),
                labels=label_log(base = 10))+
  theme_classic()+
  ggtitle("New vs old absolute co2 flux best estimates")


#Violin plots----
#Boxplots of percent difference (limited to + or - 5% difference) to show that Li-COR and Picarro differences are centered around zero while Los Gatos has a directional bias, justificating the correction of data in manuscript. 

bind_rows(
  co2newold %>% 
    filter(estimate == "best", !is.na(flux), abs(diff_percent) <= 5) %>% 
    mutate(gas = "CO2"),
  
  ch4newold %>% 
    filter(estimate == "best", !is.na(flux), abs(diff_percent) <= 5) %>% 
    mutate(gas = "CH4")
) %>%
  ggplot(aes(gas_analiser, diff_percent, fill=gas_analiser)) +
  geom_violin(scale = "width") +
  labs(fill="Gas Analizer",
       x="Gas Analizer",
       y="Relative flux difference (%)")+
  facet_wrap(~ gas) +
  theme_bw()+
  ggtitle("Distribution of estimate Relative differences",subtitle = "Only fluxes within 5% of old flux (93% of CH4 data, 98% of CO2 data)")


# ---- Packages and directories ----
library(tidyverse)
library(scales)
library(ggExtra)
library(ggpmisc)
library(ggpubr)

#Get root directory of repository:
path_root <- dirname(rstudioapi::getSourceEditorContext()$path)

#Path to source data: 
path_0_sourcedata <- paste0(path_root,"/0_sourcedata/")

#Path to formatted flux data: 
path_1_paperdata <- paste0(path_root,"/1_paperdata/")

#Path to model outputs:
path_2_modeloutputs <- paste0(path_root,"/2_modeloutputs/")

#Path to save figures and tables: 
path_mainfigures<- paste0(path_root,"/3_Figures_and_tables/")
path_supplementary<- paste0(path_root,"/3_Figures_and_tables/Supplementary/")





#1. Old vs new fluxes-------

#Comparison of old (wrong) and new instantaneous fluxes.
old<- read.csv(file=paste0(path_1_paperdata,"/old_wrong_ChamberData4paper.csv"))
new<- read.csv(file=paste0(path_1_paperdata,"/ChamberData4paper.csv"))



comp<- new %>% 
  rename(new=dailyflux) %>% 
  left_join(old, by=c("plotcode", "season", "casepilot", "status","subsite", "sampling", "strata", "ghgspecies", "unitflux")) %>% 
  rename(old=dailyflux) %>% 
  mutate(absdif=new-old,
         percentdif=100*(new-old)/abs(old))




#Total % of incubations with dailyflux estimate that retain a value that is within 5% of old estimate: 
gwp_within5percent <- comp %>% 
  filter(ghgspecies=="gwp100",
         !is.na(new)) %>% 
  summarise(n_incub_with_estimate=n(),
            n_newflux_within5percentofold=sum(abs(percentdif)<=5),
            percent_newflux_within5percentofold=n_newflux_within5percentofold/n_incub_with_estimate*100)
gwp_within5percent

co2_within5percent <- comp %>% 
  filter(ghgspecies=="co2",
         !is.na(new)) %>% 
  summarise(n_incub_with_estimate=n(),
            n_newflux_within5percentofold=sum(abs(percentdif)<=5),
            percent_newflux_within5percentofold=n_newflux_within5percentofold/n_incub_with_estimate*100)
co2_within5percent


ch4_within5percent <- comp %>% 
  filter(ghgspecies=="ch4",
         !is.na(new)) %>% 
  summarise(n_incub_with_estimate=n(),
            n_newflux_within5percentofold=sum(abs(percentdif)<=5),
            percent_newflux_within5percentofold=n_newflux_within5percentofold/n_incub_with_estimate*100)
ch4_within5percent




#Written summary: considering all incubations, not only those included in the paper (but also those considered not-relevant for comparisons and excluded in the preliminary filtering of data)

paste0("Of the ", 
       co2_within5percent$n_incub_with_estimate, 
       " chamber deployments with a valid daily CO2 flux estimate included in the paper, ",
       co2_within5percent$n_newflux_within5percentofold,
       " (",
       round(co2_within5percent$percent_newflux_within5percentofold,2),
       "%) had a new daily CO2 flux estimate that lied within 5% of their old estimate value.")

paste0("Of the ", 
       ch4_within5percent$n_incub_with_estimate, 
       " chamber deployments with a valid daily CH4 flux estimate included in the paper, ",
       ch4_within5percent$n_newflux_within5percentofold,
       " (",
       round(ch4_within5percent$percent_newflux_within5percentofold,2),
       "%) had a new daily CH4 flux estimate that lied within 5% of their old estimate value.")

paste0("Of the ", 
       gwp_within5percent$n_incub_with_estimate, 
       " incubations with a valid daily combined (CO2+CH4) CO2-eq flux estimate included in the paper, ",
       gwp_within5percent$n_newflux_within5percentofold,
       " (",
       round(gwp_within5percent$percent_newflux_within5percentofold,2),
       "%) had a new daily CO2-eq flux estimate that lied within 5% of their old estimate value.")



##CO2 new vs old biplot-----
comp %>% 
  filter(ghgspecies=="co2",
         !is.na(new)) %>% 
  mutate(incubtype=case_when(abs(percentdif)<=5~"New within 5% of old",
                             TRUE~"New not within 5% of old"
  )) %>% 
  arrange(desc(incubtype)) %>% 
  ggplot(aes(x=abs(old),y=abs(new), col=incubtype))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1)+
  labs(col="Incubation type")+
  scale_x_log10(name="Old estimate (abs. value)",limits=c(5e-6,1e1),
                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100),
                labels=label_log(base = 10))+
  scale_y_log10(name="New estimate (abs. value)", limits=c(5e-6,1e1),
                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100),
                labels=label_log(base = 10))+
  theme_classic()+
  facet_wrap(facets=vars(casepilot), scale="free")+
  ggtitle("New vs old absolute CO2 dailyflux values")


##CH4 new vs old biplot-----
comp %>% 
  filter(ghgspecies=="ch4",
         !is.na(new)) %>% 
  mutate(incubtype=case_when(abs(percentdif)<=5~"New within 5% of old",
                             TRUE~"New not within 5% of old"
  )) %>% 
  arrange(desc(incubtype)) %>% 
  ggplot(aes(x=abs(old),y=abs(new), col=incubtype))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1)+
  labs(col="Incubation type")+
  scale_x_log10(name="Old estimate (abs. value)",limits=c(1e-8,2e3),
                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000),
                labels=label_log(base = 10))+
  scale_y_log10(name="New estimate (abs. value)", limits=c(1e-8,2e3),
                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000),
                labels=label_log(base = 10))+
  theme_classic()+
  facet_wrap(facets=vars(casepilot), scale="free")+
  ggtitle("New vs old absolute CH4 dailyflux values")



##Summary per site-----

co2_summary<- new %>% 
  filter(ghgspecies=="co2") %>% 
  group_by(casepilot, status, subsite,unitflux) %>% 
  summarise(avg=mean(dailyflux,na.rm=T),
            SD=sd(dailyflux, na.rm=T),
            N=sum(!is.na(dailyflux)),
            SE=SD/sqrt(N),
            Q10=quantile(dailyflux, 0.10, na.rm=T),
            Q25=quantile(dailyflux, 0.25, na.rm=T),
            median_Q50=median(dailyflux, na.rm=T),
            Q75=quantile(dailyflux, 0.75, na.rm=T),
            Q90=quantile(dailyflux, 0.10, na.rm=T),
            min=min(dailyflux, na.rm=T),,
            max=max(dailyflux, na.rm=T))


ch4_summary<- new %>% 
  filter(ghgspecies=="ch4") %>% 
  group_by(casepilot, status, subsite,unitflux) %>% 
  summarise(avg=mean(dailyflux,na.rm=T),
            SD=sd(dailyflux, na.rm=T),
            N=sum(!is.na(dailyflux)),
            SE=SD/sqrt(N),
            Q10=quantile(dailyflux, 0.10, na.rm=T),
            Q25=quantile(dailyflux, 0.25, na.rm=T),
            median_Q50=median(dailyflux, na.rm=T),
            Q75=quantile(dailyflux, 0.75, na.rm=T),
            Q90=quantile(dailyflux, 0.10, na.rm=T),
            min=min(dailyflux, na.rm=T),,
            max=max(dailyflux, na.rm=T))





#2. CH4 biogeochemistry------

#This section explores the relationship between CH4 fluxes, biochemical water parameters and inundation patterns. 


#Fluxes in paper
paperdat<- read.csv(paste0(path_1_paperdata,"ChamberData4paper.csv"))


#Chamber deployments 
plotcodes<- read.csv(paste0(path_0_sourcedata, "Plotcode_harmonized_field_obs.csv"))


#Add water-depth (cm) and inundated (T/F) to ch4 dailyfluxes
ch4_daily<- paperdat %>% 
  filter(ghgspecies=="ch4") %>% 
  left_join(plotcodes %>% dplyr::select(plotcode, water_depth), by="plotcode") %>% 
  mutate(inundated=water_depth>0)


#Water dataset: 
#We need sampling-level averages for:
#Chla: microg/L
#TN and TP units: micromol/L
#Conductivity: mS/cm


#Use high-tide samplings for Aveiro, they are the relevant ones for estimate of sulphate-suply via salinity. We might need to manually average and select data for RI due to grouping of high-tide samplings (collected in the channels in between subsites). 


w<- readxl::read_xlsx(path = paste0(path_0_sourcedata, "Supp_Ancillary_Water_Data_RIhightide.xlsx"))

wat<- w %>% 
  dplyr::select(Season,Site,Subsite,Tide,`Label Sample ID`,`Water depth (m)`,`Tw (ºC)`, `Cond (µS/cm)`,Salinity,pH,`DO (%)`,`Ox (mg/l)`,`Total bacteria (cel ml-1) BW`,`DOC (µM)`,`Diss. Total-P (micromol P/L)`,`Diss. Ntotal (micromol N/L)`,`Part. Total-P (micromol P/L)`,`Part. Total -N (micromol N/L)`,`Chl-a with correction for phaeopigments, µg/l`,`Patm (millibar)`) %>% 
  rename(season=Season, casepilot=Site, subsite=Subsite, sampleID=`Label Sample ID`,
         depth=`Water depth (m)`, temp=`Tw (ºC)`,
         conductivity=`Cond (µS/cm)`, salinity=Salinity, ph=pH, DO_percent=`DO (%)`, DO_mgl=`Ox (mg/l)`,
         bact=`Total bacteria (cel ml-1) BW`, doc=`DOC (µM)`, 
         tdp=`Diss. Total-P (micromol P/L)`,tdn=`Diss. Ntotal (micromol N/L)`,
         tpp=`Part. Total-P (micromol P/L)`, tpn=`Part. Total -N (micromol N/L)`,
         chla=`Chl-a with correction for phaeopigments, µg/l`,
         patm_mbar=`Patm (millibar)`) %>% 
  #Remove low-tide samples of Aveiro (the relevant measurement for sulphate suply is the high-tide, low tide never contacts the sites sensu-stricto): 
  filter(!(Tide=="LT"&casepilot=="RI")) %>%
  #subsitute all below detection limit with cero
  mutate(tdp=as.numeric(gsub("<LD", "0", tdp))) %>% 
  #add status: 
  mutate(status=factor(case_when(grepl("A[0-9]",subsite)~"Altered",
                                 grepl("P[0-9]",subsite)~"Preserved",
                                 grepl("R[0-9]",subsite)~"Restored"), levels=c("Altered","Preserved","Restored"), ordered=T)) %>% 
  #Calculate water TN and TP
  mutate(tn=tdn+tpn,
         tp=tdp+tpp)%>%
  dplyr::select(casepilot, season, status, subsite, chla, conductivity, tn,tp) 


#Calculate summary statistics for the relevant variables (seasonal): 
wat_summary_season <- wat%>% 
  mutate(conductivity=conductivity/1000) %>% #Conductivity in mS/cm
  pivot_longer(cols = -c(casepilot,status,season,subsite), names_to = "variable",values_to = "value") %>% 
  group_by(casepilot, variable, status,season,subsite) %>% 
  summarise(avg=mean(value, na.rm=T),
            SD=sd(value, na.rm=T), 
            nobs=sum(!is.na(value)),
            se=SD/sqrt(nobs)) %>% 
  dplyr::select(casepilot, status, season,subsite, variable, avg, se, nobs)

wide_wat_summary_season<- wat_summary_season %>% 
  dplyr::select(-nobs) %>% 
  pivot_wider(names_from = "variable",values_from = c("avg","se")) %>% 
  mutate(subsite_only= subsite,
         subsite= paste(casepilot, subsite_only,sep="-"))



#Calculate summary statistics for the relevant variables (annual): 
wat_summary_year <- wat%>% 
  dplyr::select(-season) %>% 
  mutate(conductivity=conductivity/1000) %>% #Conductivity in mS/cm
  pivot_longer(cols = -c(casepilot,status,subsite), names_to = "variable",values_to = "value") %>% 
  group_by(casepilot, variable, status,subsite) %>% 
  summarise(avg=mean(value, na.rm=T),
            SD=sd(value, na.rm=T), 
            nobs=sum(!is.na(value)),
            se=SD/sqrt(nobs)) %>% 
  dplyr::select(casepilot, status, subsite, variable, avg, se, nobs)

wide_wat_summary_year<- wat_summary_year %>% 
  dplyr::select(-nobs) %>% 
  pivot_wider(names_from = "variable",values_from = c("avg","se")) %>% 
  mutate(subsite_only= subsite,
         subsite= paste(casepilot, subsite_only,sep="-"))



ch4_daily_bch_seasonal<- ch4_daily %>% 
  left_join(wide_wat_summary_season %>% select(-subsite_only), by=c("casepilot", "status", "season", "subsite")) %>% 
  mutate(chamber_typewater=if_else(inundated,"Inundated","Dry"))



#For each subsite visit (total_n=144), calculate mean, median, se, ch4 flux and proportion of chambers with water, add small value to dailyflux to avoid negative and zero-flux values from being excluded during log10 transformation.

ch4_sampling_summary<- ch4_daily_bch_seasonal %>% 
  group_by(casepilot, status, subsite, sampling, ghgspecies, unitflux) %>% 
  mutate(dailyflux=dailyflux+0.001) %>% 
  summarise(mean_ch4=mean(dailyflux, na.rm=T),
            median_ch4=median(dailyflux, na.rm=T),
            q25_ch4=quantile(dailyflux, 0.25, na.rm=T),
            q75_ch4=quantile(dailyflux, 0.75, na.rm=T),
            SD_ch4=sd(dailyflux, na.rm=T),
            n_ch4=sum(!is.na(dailyflux)),
            prop_inundation= sum(inundated)/sum(!is.na(inundated)),
            mean_conductivity=mean(avg_conductivity),
            mean_chla=mean(avg_chla),
            SE_chla=mean(se_chla),
            mean_TP=mean(avg_tp),
            se_TP=mean(se_tp))


##Global CH4-inundation----
#General cross-system correlation CH4 vs inundation  
ch4_sampling_summary %>% 
  ggplot(aes(x=prop_inundation, y=(median_ch4)))+
  geom_point(aes(col=casepilot))+
  stat_cor()+
  scale_y_log10(label=label_log(digits = 2))+
  labs(y=expression(Median ~ CH[4] ~ flux ~ (mmol~m^-2~d^-1)),
       x="Inundation proportion",
       col="Case pilot")+
  geom_smooth(method="lm",se=F, col="black")+
  theme_bw()


##Global CH4-conductivity----
#General cross-system correlation CH4 vs conductivity
ch4_sampling_summary %>% 
  ggplot(aes(x=(mean_conductivity), y=(median_ch4)))+
  geom_point(aes(col=casepilot))+
  stat_cor()+
  scale_y_log10(label=label_log(digits = 2))+
  scale_x_log10()+
  labs(y=expression(Median ~ CH[4] ~ flux ~ (mmol~m^-2~d^-1)),
       x=expression(Mean ~ water ~ conductivity ~ (mS~cm^-1)),
       col="Case pilot")+
  geom_smooth(method="lm",se=F, col="black")+
  theme_bw()


##Global CH4-Chl-a (not included)----
#General cross-sysem correlation CH4 vs Chla
ch4_sampling_summary %>% 
  ggplot(aes(x=(mean_chla), y=(median_ch4)))+
  geom_point(aes(col=casepilot))+
  stat_cor()+
  scale_y_log10(label=label_log(digits = 2))+
  scale_x_log10()+
  labs(y=expression(Median ~ CH[4] ~ flux ~ (mmol~m^-2~d^-1)),
       x=expression(Mean ~ water ~ Chl-a ~ (mg~L^-1)),
       col="Case pilot")+
  geom_smooth(method="lm",se=F, col="black")+
  theme_bw()


##CA & VA CH4-inundation (not included)-----
#Effect of Inundation for sites CA, VA with variable inundation patterns (DA has status-determined inundation, CU sites are completely inundated by definition )

n_labels <- ch4_daily_bch_seasonal %>% 
  filter(casepilot %in% c("CA", "VA")) %>% 
  group_by(casepilot, status, chamber_typewater) %>% 
  mutate(log_dailyfluxdat=log10(dailyflux + 0.001)) %>% 
  summarise(
    n = sum(!is.na(log_dailyfluxdat)),
    y = max(log10(dailyflux + 0.001), na.rm = TRUE) + 0.2,
    .groups = "drop"
  )

ch4_daily_bch_seasonal %>% 
  filter(casepilot %in% c("CA", "VA")) %>% 
  ggplot(aes(
    x = status,
    y = log10(dailyflux + 0.001),
    
    fill=chamber_typewater
  )) +
  geom_boxplot() +
  geom_text(
    data = n_labels,
    aes(
      x = status,
      y = y,
      col = chamber_typewater,
      label = paste0("n=", n)
    ),
    position = position_dodge(width = 0.75),
    size = 3,
    show.legend = FALSE
  ) +
  facet_wrap(vars(casepilot)) +
  stat_compare_means(label = "p.signif") +
  theme_bw()+
  ggtitle("Inundation impact on CH4 fluxes for CA and VA",
          subtitle = "Based on Kruskal-Wallis tests")


ch4_sampling_summary %>% 
  filter(casepilot%in%c("VA","CA")) %>% 
  ggplot(aes(x=prop_inundation, y=(median_ch4),col=casepilot))+
  geom_point()+
  stat_cor()+
  scale_y_log10()+
  labs(y="Median CH4 flux",
       x="Inundation proportion",
       col="Case pilot")+
  geom_smooth(method="lm",se=F)+
  theme_bw()



##CU CH4-TP (not included)-----

#CU only median CH4 vs mean TP
ch4_sampling_summary %>% 
  filter(casepilot=="CU") %>% 
  ggplot(aes(x=(mean_TP), y=(median_ch4)))+
  geom_point(aes(col=status))+
  stat_cor()+
  geom_smooth(method="lm",se=F, col="black")+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()







#3. RI fluxes- literature comparison-----

#Calculate average transparent and dark fluxes for sediment and zostera to compare with those presented in  Bahlmann et al. 2015 Biogeosciences
co2_all<- read.csv(paste0(path_0_sourcedata,"co2_bestflux.csv")) #umol/m2/s
ch4_all<- read.csv(paste0(path_0_sourcedata,"ch4_bestflux.csv")) #nmol/m2/s


#Select best fluxes and format: UniqueID, ghg_best, ghg_se
#CO2
RIco2_best<- co2_all %>% 
  filter(grepl("ri",UniqueID)) %>% 
  mutate(co2_best=case_when(best.model=="LM"~LM.flux,
                            best.model=="HM"~HM.flux,
                            best.model=="total.flux"~total.flux,
                            best.model=="None  appropriate"~NA_real_),
         co2_se=case_when(best.model=="LM"~LM.flux.se,
                          best.model=="HM"~HM.flux.se,
                          best.model=="total.flux"~total.flux.se,
                          best.model=="None  appropriate"~NA_real_),
         co2_model=best.model) %>% 
  dplyr::select(UniqueID, co2_best, co2_se,co2_model)

#CH4
RIch4_best<- ch4_all %>% 
  filter(grepl("ri",UniqueID)) %>% 
  mutate(ch4_best=case_when(best.model=="LM"~LM.flux,
                            best.model=="HM"~HM.flux,
                            best.model=="total.flux"~total.flux,
                            best.model=="None  appropriate"~NA_real_),
         ch4_se=case_when(best.model=="LM"~LM.flux.se,
                          best.model=="HM"~HM.flux.se,
                          best.model=="total.flux"~total.flux.se,
                          best.model=="None  appropriate"~NA_real_),
         ch4_model=best.model) %>% 
  dplyr::select(UniqueID, ch4_best, ch4_se,ch4_model)


#JOIN all ghg best-fluxes per UniqueID
RIghg_best<- RIco2_best %>% 
  merge.data.frame(RIch4_best, by="UniqueID", all=T)

rm(co2_all, RIco2_best, ch4_all, RIch4_best)


#Import harmonized field observations, one row per chamber deployment (including 1 or 2 incubation "UniqueID"s, dark or transparent)

#Harmonized field observations per plotcode (per chamber deployment)
fieldobs_plotcode<- read.csv(paste0(path_0_sourcedata,"plotcode_harmonized_field_obs.csv"))
RIfieldobs_plotcode<- fieldobs_plotcode %>% 
  filter(casepilot=="RI")

#Check correspondence between RIfieldobs_plotcode and RIghg_best
plotharm_IDs<- unique(c(RIfieldobs_plotcode$dark_UniqueID, RIfieldobs_plotcode$transparent_UniqueID))
plotharm_IDs<- plotharm_IDs[!is.na(plotharm_IDs)] #exclude NA (originating from absence of transparent or dark incubation)

ghgbest_IDs<- RIghg_best %>% select(UniqueID) %>% distinct() %>% pull(UniqueID)

#field-obs without ghg_best row: 
plotharm_IDs[which(!(plotharm_IDs%in%ghgbest_IDs))]

#ghg-best-rows without field-obs: 
ghgbest_IDs[which(!(ghgbest_IDs%in%plotharm_IDs))]
#These 27 incubations correspond always to "extra" incubations performed in vegetated areas after removing the vegetation, non representative, implicitly excluded due to not being present in fieldobs_plotcode. 



#Filter valid incubations

#Remove fieldobs from non-valid/non-appropriate incubations: Incubations that
#are not representative of subsite definitions, incubations that are not
#representative of natural conditions

# 1. Remove rising-tide and receding-tide incubations. (non appropriate fluxes for comparisons)
tidal_IDs <- RIfieldobs_plotcode %>%
  filter(grepl("rising-tide|receding-tide", field_observations)) %>%
  select(transparent_UniqueID, dark_UniqueID) %>%
  pivot_longer(cols = everything(), values_to = "UniqueID") %>%
  pull(UniqueID) %>%
  na.omit()

#2.5. Ria d'Aveiro (RI): rising-tide and "after vegetation removal" already
#removed. Inconsistent incubations: tidal-pool incubations in S3-RI-A1 (only
#considered in S3, less than 10% of site), bare incubations in RI-R1 and in
#RI-R2 (wrong extra-incubations, restored sites are 100% vegetated by
#definition)
RI_wrongstrata_IDs<- RIfieldobs_plotcode %>% 
  filter((sampling=="S3-RI-A1"&grepl("tidal-pool", field_observations))|(subsite=="RI-R1"&strata=="bare")|(subsite=="RI-R2"&strata=="bare")) %>%
  select(transparent_UniqueID, dark_UniqueID) %>%
  pivot_longer(cols = everything(), values_to = "UniqueID") %>%
  pull(UniqueID) %>%
  na.omit()


#Combine all IDs that need to be excluded: 
IDs_toexclude<- c(tidal_IDs,RI_wrongstrata_IDs)


#Filter-out the field observations with incubation UniqueIDs that need exclusion: 
valid_fieldobs_plotcode<- RIfieldobs_plotcode %>% 
  filter(!dark_UniqueID%in%IDs_toexclude) %>% 
  filter(!transparent_UniqueID%in%IDs_toexclude)

#Check number and proportion of chamber deployments (plotcodes) excluded from further analysis: 
valid_plotcodes<- valid_fieldobs_plotcode %>% dplyr::select(plotcode) %>% distinct() %>% pull(plotcode)

plotcodes_with_bestflux<- RIghg_best %>% 
  filter(!is.na(co2_best)|!is.na(ch4_best)) %>% 
  dplyr::select(UniqueID) %>% 
  mutate(UniqueID=toupper(UniqueID)) %>% 
  separate(UniqueID, into=c("c1","c2","c3","c4","d1","d2","d3"), sep = "-",remove = T) %>% 
  mutate(plotcode=paste(c1,c2,c3,c4,sep = "-")) %>% 
  dplyr::select(plotcode) %>% distinct() %>% pull(plotcode)

#plotcodes (chamber deployments) with valid fluxes
length(plotcodes_with_bestflux)

n_all_plotcodes_with_bestflux<- length(plotcodes_with_bestflux)
n_valid_plotcodes_with_bestflux<- sum(valid_plotcodes%in%plotcodes_with_bestflux)

#Number of plotcodes with valid fluxes (with at least 1 best-flux estimate) that
#were discarded due to being performed in non-comparable conditions (out of
#site, during high tide, after veg cut only,...)
n_all_plotcodes_with_bestflux-n_valid_plotcodes_with_bestflux

#Percentage
round((n_all_plotcodes_with_bestflux-n_valid_plotcodes_with_bestflux)/n_all_plotcodes_with_bestflux*100,1)

#PLotcodes with valid flux excluded due to tidal behaviour
RIfieldobs_plotcode %>% 
  filter(plotcode%in%plotcodes_with_bestflux) %>% 
  filter((dark_UniqueID%in%tidal_IDs)|(transparent_UniqueID%in%tidal_IDs)) %>% #Tidal behaviour
  dplyr::select(plotcode) %>% distinct() %>% dim() 
#Rest of excluded plotcodes are due to other reasons (performed in non representative conditions, outside site boundaries,...)


rm(tidal_IDs,RI_wrongstrata_IDs)


#Subset RI observations and join to ghgfluxes
RI_obs_long<- valid_fieldobs_plotcode %>% 
  pivot_longer(cols = c(transparent_UniqueID, dark_UniqueID), names_to = "lightcondition", values_to = "UniqueID")

RI_long<- RI_obs_long %>% 
  left_join(RIghg_best, by="UniqueID")


#Summary for RI (transparent or dark, bare or vegetated) all during lowtide
RI_summary<- RI_long %>% 
  filter(strata!="open water") %>% 
  mutate(co2_best=co2_best*60*60*0.001, #from umol m-2 s-1 to mmol m-2 h-1
         ch4_best=ch4_best*60*60*0.001) %>%  #from nmol m-2 s-1 to umol m-2 h-1
  group_by(strata, lightcondition) %>% 
  summarise(co2_mean=mean(co2_best, na.rm=T),
            co2_sd=sd(co2_best, na.rm=T),
            co2_q10=quantile(co2_best,0.1, na.rm=T),
            co2_q25=quantile(co2_best,0.25, na.rm=T),
            co2_q50=quantile(co2_best,0.5, na.rm=T),
            co2_q75=quantile(co2_best,0.75, na.rm=T),
            co2_q90=quantile(co2_best,0.9, na.rm=T),
            ch4_mean=mean(ch4_best, na.rm=T),
            ch4_sd=sd(ch4_best, na.rm=T),
            ch4_q10=quantile(ch4_best,0.1, na.rm=T),
            ch4_q25=quantile(ch4_best,0.25, na.rm=T),
            ch4_q50=quantile(ch4_best,0.5, na.rm=T),
            ch4_q75=quantile(ch4_best,0.75, na.rm=T),
            ch4_q90=quantile(ch4_best,0.9, na.rm=T))

RI_summary
#Our CO2 data agrees reasonably well with that of Bahlmann 2015: 
#Emerged Z.noltii in the light 
#Emerged Z.noltii in the dark

#Emerged sediment in the light
#Emerged sediment in the dark



#Ebullition contribution------

#Reviewers ask to include discussion on the contribution of ebullitive fluxes to total fluxes, especially for DA and CU. 

#Two approaches: 
#1. purely frequentist: proportion of incubations that show ebullitive behaviour. Plotcodes where at least one of the incubations (in case of multiple ones) show ebullition vs total plotcodes (per casepilot and per casepilot-status).

#2. quantitative contribution of flux from plotcodes with ebullitive incubations to total CH4 flux. Not exactly contribution of ebullition to total flux (as incubations with ebullition can also present diffusion), but potentially more informative. 


#In any case, we require: 1. presence/absence of ebullition per plotcode and daily ch4 flux per plotcode. 
#ebullitive inspection: we can obtain it directly form the best.model.flags column of ch4_bestflux.csv

ch4best<- read.csv(paste0(path_0_sourcedata,"/ch4_bestflux.csv"))

#Obtain ebullition presence at each plotcode based on ch4_bestflux.csv 
incu_ebu<- ch4best %>% 
  select(UniqueID, best.model.flags) %>% 
  separate(UniqueID, into = c("s","cp","sub", "plotnum", "strat", "light", "hourminute"),sep = "-", remove = F) %>% 
  mutate(plotcode=toupper(paste(s,cp,sub,plotnum, sep = "-")),
         ebu=grepl("Ebullitive dynamics",best.model.flags)) %>% 
  select(plotcode, ebu) %>% 
  group_by(plotcode) %>% 
  summarise(n_ebullition_incub=sum(ebu),
            ebullitive=n_ebullition_incub>0) %>% 
  select(plotcode, ebullitive)


#Fluxes in paper
paperdat<- read.csv(paste0(path_1_paperdata,"ChamberData4paper.csv"))



#obtain only daily ch4 fluxes in paper and add logical ebullitive column
ch4_dat<- paperdat %>% 
  filter(ghgspecies=="ch4") %>% 
  filter(!is.na(dailyflux)) %>% 
  left_join(incu_ebu, by="plotcode") %>% 
  mutate(dailyflux=if_else(dailyflux<0,0,dailyflux)) %>% 
  mutate(casepilot=factor(casepilot, levels=c("DU","RI","CA","VA","DA","CU"), ordered = T))

#Calculate for each casepilot and for each casepilot and status, the % of ebullitive chamber-measurements
freq_ebullition<- ch4_dat %>% 
  group_by(casepilot, status) %>% 
  summarise(status_n_chambers=sum(!is.na(ebullitive)),
            status_n_ebuchambers=sum(ebullitive),
            status_percent_ebuchambers=status_n_ebuchambers/status_n_chambers*100) %>% 
  ungroup() %>% 
  group_by(casepilot) %>% 
  mutate(casepilot_n_chambers=sum(status_n_chambers),
         casepilot_n_ebuchambers=sum(status_n_ebuchambers),
         casepilot_percent_ebuchambers=casepilot_n_ebuchambers/casepilot_n_chambers*100)


#PRINT summary: % of chamber deployments with evidence of ebullitive dynamics.
freq_ebullition %>% 
  arrange(casepilot, status) %>% 
  select(casepilot, status, status_percent_ebuchambers,casepilot_percent_ebuchambers)



#Now calculate the prevalence of ebullition by measuring how much of the cummulative flux measured at each status and casepilot was obtained from deployments with ebullition (approximation of total ebullitive flux impact)
flux_prop_ebullition<- ch4_dat %>% 
  group_by(casepilot,status, ebullitive) %>% 
  summarise(flux_sum=sum(dailyflux)) %>% 
  ungroup() %>% 
  group_by(casepilot,status) %>% 
  mutate(status_totalflux=sum(flux_sum)) %>% 
  filter(ebullitive==T) %>% 
  rename(status_ebuflux=flux_sum) %>% 
  mutate(status_percentebuflux=status_ebuflux/status_totalflux*100) %>% 
  arrange(casepilot,status)


#Also interesting, summary statistics for average flux depending on presence/absence of ebullition
flux_avg_ebullition<- ch4_dat %>% 
  group_by(casepilot,status, ebullitive, unitflux) %>% 
  summarise(flux_mean=mean(dailyflux),
            flux_sd=sd(dailyflux),
            flux_cv_percent=(flux_sd/flux_mean)*100,
            flux_n=sum(!is.na(dailyflux)),
            flux_q10=quantile(dailyflux,0.10),
            flux_q25=quantile(dailyflux,0.25),
            flux_q50=quantile(dailyflux,0.5),
            flux_q75=quantile(dailyflux,0.75),
            flux_q90=quantile(dailyflux,0.9))

flux_avg_ebullition %>% filter(casepilot%in%c("CA","VA","DA","CU"))

library(e1071)
#Normalized flux distribution according to ebullitive/non-ebullitive chambers:
ch4_dat %>% 
  filter(casepilot%in%c("CA","VA","DA","CU"))%>% 
  group_by(casepilot,status, ebullitive, unitflux) %>% 
  #Standardize per status and casepilot
  mutate(normflux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  ggplot(aes(x=ebullitive, y=normflux))+
  geom_violin(scale = "width")+
  geom_boxplot()+
  geom_text(
    data = . %>% 
      group_by(casepilot, ebullitive) %>% 
      summarise(
        n = n(),
        y = max(normflux, na.rm = TRUE) + 0.05,
        cv= sd(normflux)/mean(normflux)*100,
        skew= skewness(normflux),
        .groups = "drop"
      ),
    aes(x = ebullitive+1.1, y = y, label = paste0("n=", n, ", \ncv=", round(cv,0),"%, \nskew=", round(skew,2))),vjust = 1,hjust=0,size = 3,
    inherit.aes = FALSE
  ) +
  facet_wrap(facets=vars(casepilot))+
  labs(y="Stardardized CH4 flux", x="Chamber with ebullition")+
  ggtitle("Distribution of CH4 flux per type of chamber",
          subtitle="Flux is standardised by case pilot, status and ebullition presence")


#Frequency of ebullition: 
freq_ebullition %>% 
  # filter(casepilot%in%c("CA","VA","DA","CU"))%>% 
  ggplot(aes(x=status, y=status_percent_ebuchambers, fill=status))+
  geom_bar(stat="identity")+
  scale_y_continuous(name="% of chambers with ebullition", limits = c(0,100))+
  facet_wrap(facets=vars(casepilot))+
  ggtitle("Frequency of chamber with ebullition",
          subtitle = "(% total chambers)")


#Relationship between proportion of chambers with ebullition (cp,status,season) and CV of flux
ch4_dat %>% 
  group_by(casepilot, status, season,subsite) %>% 
  mutate(standard_flux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  summarise(percent_cv_flux=sd(standard_flux)/mean(standard_flux)*100,
            flux_skew=skewness(standard_flux),
            flux_sd=sd(standard_flux),
            percent_ebullitive=sum(ebullitive)/sum(!is.na(ebullitive)*100),
            any_ebullitive=sum(ebullitive)>0) %>% 
  ggplot(aes(x=any_ebullitive, y=percent_cv_flux))+
  geom_boxplot()+
  geom_text(
    data = . %>% 
      group_by(any_ebullitive) %>% 
      summarise(n=n(),max_cv=max(percent_cv_flux)),
    aes(x = any_ebullitive, y = 1.05*max_cv, label = paste0("n=", n)),vjust = 1,hjust=0,size = 3,
    inherit.aes = FALSE
  ) +
  stat_compare_means()+
  labs(y="CV of standardised dailyflux (%)", x="Presence of ebullitive chamber in sampling")+
  ggtitle("Flux variability per sampling",
          subtitle = "Each point is a distinct sampling event (n=144)")

#Flux variability does not scale with proportion of ebullition
ch4_dat %>%
  group_by(casepilot, status, season,subsite) %>% 
  mutate(standard_flux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  summarise(percent_cv_flux=sd(standard_flux)/mean(standard_flux)*100,
            percent_ebullitive=sum(ebullitive)/sum(!is.na(ebullitive)*100),
            any_ebullitive=sum(ebullitive)>0) %>% 
  ggplot(aes(x=percent_ebullitive, y=percent_cv_flux, col=any_ebullitive))+
  geom_point()+
  labs(y="CV of standardised dailyflux (%)", x="Proportion of chambers with ebullition")+
  geom_smooth(method="lm")+
  stat_cor()+
  ggtitle("Flux variability vs ebullitive chamber prevalence",
          subtitle = "Each point is a distinct sampling event (n=144)")


#Prevalence of ebullition vs flux magnitude?
ch4_dat %>% 
  group_by(casepilot, status, season,subsite) %>% 
  mutate(standard_flux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  summarise(mean_flux=mean(dailyflux),
            median_flux=median(dailyflux)+0.001,
            sd_flux=sd(dailyflux),
            percent_ebullitive=sum(ebullitive)/sum(!is.na(ebullitive)*100),
            any_ebullitive=sum(ebullitive)>0) %>% 
  ggplot(aes(x=log10(median_flux), y=percent_ebullitive))+
  geom_point(aes(col=any_ebullitive))+
  # labs(y="CV of standardised dailyflux (%)", x="log10(Mean CH4 flux)")+
  # geom_smooth(method="lm")+
  stat_cor()+
  ggtitle("Ebullition prevalence vs median flux magnitude",
          subtitle = "Each point is a distinct sampling event (n=144)")



ebu_vs_mag<-ch4_dat %>% 
  group_by(casepilot, status, season,subsite) %>% 
  mutate(standard_flux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  summarise(mean_flux=mean(dailyflux),
            median_flux=median(dailyflux),
            sd_flux=sd(dailyflux),
            percent_ebullitive=sum(ebullitive)/sum(!is.na(ebullitive)*100),
            any_ebullitive=sum(ebullitive)>0)




#Skewness of flux distribution vs chambers with ebullitive
ch4_dat %>% 
  group_by(casepilot, status, season,subsite) %>% 
  mutate(standard_flux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  summarise(percent_cv_flux=sd(standard_flux)/mean(standard_flux)*100,
            flux_skew=skewness(standard_flux),
            flux_sd=sd(standard_flux),
            percent_ebullitive=sum(ebullitive)/sum(!is.na(ebullitive)*100),
            any_ebullitive=sum(ebullitive)>0) %>% 
  ggplot(aes(x=percent_ebullitive, y=flux_skew, col=any_ebullitive))+
  geom_point()+
  labs(y="Skewness of chamber flux distribution", x="Proportion of chambers with ebullition")+
  geom_smooth(method="lm")+
  stat_cor()+
  ggtitle("Flux variability vs ebullitive chamber prevalence",
          subtitle = "Each point is a distinct sampling event (n=144)")


ch4_dat %>% 
  group_by(casepilot, status, season,subsite) %>% 
  mutate(standard_flux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  summarise(percent_cv_flux=sd(standard_flux)/mean(standard_flux)*100,
            flux_skew=skewness(dailyflux),
            flux_sd=sd(standard_flux),
            percent_ebullitive=sum(ebullitive)/sum(!is.na(ebullitive)*100),
            any_ebullitive=sum(ebullitive)>0) %>% 
  ggplot(aes(x=any_ebullitive, y=flux_skew))+
  geom_boxplot()+
  geom_text(
    data = . %>% 
      group_by(any_ebullitive) %>% 
      summarise(n=n(),y=max(flux_skew)),
    aes(x = any_ebullitive, y = 1.05*y, label = paste0("n=", n)),vjust = 1,hjust=0,size = 3,
    inherit.aes = FALSE
  ) +
  stat_compare_means()+
  labs(y="Skewness of chamber flux distribution", x="Presence of ebullitive chamber in sampling")+
  ggtitle("Skewness of chamber flux distribution",
          subtitle = "Each point is a distinct sampling event (n=144)")




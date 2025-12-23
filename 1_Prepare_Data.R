#Prepare chamber data for modeling

#Author: Miguel Cabrera-Brufau
#Date: September 2025
#Project: Restore4cs



#Description----
#This scrip is used to prepare the datasets that will be used for statistical
#analysis and figures of the main GHG paper that will describe the effects of
#wetland restoration on the fluxes of GHGs using daily GHG exchange rates of
#appropriate incubations derived from in-situ measurements.

#MAIN STEPS of DATA PREPARATION: 
#1. Import all fluxes and combine in single table (UniqueID, ghg_best, ghg_se, ghg_model) 
#2. Import harmonized field observations. 
#3. Discard non-appropriate fluxes (rising/receding tide, after veg cut, appropriate strata per subsite only)
#5. Integrate day-night fluxes,calculate global warming potential GWPco2andch4 (ch4+co2)
#6. Final formatting
#7. Save datasets for further analysis. 

#Naming: 

#UniqueID represents a unique incubation (under transparent or dark conditions)
#associated with a CO2 and a CH4 concentration timeseries

#plotcode represents a chamber deployment (placement of a static chamber in a
#specific location to perform incubations). It can be associated with a single
#UniqueID, or to two UniqueID(s), depending on the strata.

#Inputs to script: 
#Per-UniqueID: co2 and ch4 best fluxes (and se flux)
#Per-UniqueID and per-plotcode harmonized field observations. 

#Outputs: 
#ChamberData4paper.csv: containing valid daily fluxes of co2, ch4 and combined GWP (co2+ch4).


rm(list = ls()) # clear workspace

# ---- Packages ----

#Installs (if needed) and loads required packages:
required_pkgs <- c("tidyverse",
                   "readxl",
                   "lubridate",
                   "zoo",
                   "purrr",
                   "data.table",
                   "tools",
                   "hms",
                   "suncalc")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


# ---- Directories ----
#Get root directory of repository:
path_root <- dirname(rstudioapi::getSourceEditorContext()$path)

#Path to source data:
path_0_sourcedata <- paste0(path_root,"/0_sourcedata/")

#Path to save datasets for main paper: 
path_1_paperdata <- paste0(path_root,"/1_paperdata/")

#Create directory if it does not exist
if (!dir.exists(path_1_paperdata)) {
  dir.create(path_1_paperdata, recursive = TRUE)
}


#__________________------
#Chamber DATA PREP------
{
#1. Import inst. fluxes-----

#All datasets must be identified at the UniqueID level
co2_all<- read.csv(paste0(path_0_sourcedata,"co2_bestflux.csv")) #umol/m2/s
ch4_all<- read.csv(paste0(path_0_sourcedata,"ch4_bestflux.csv")) #nmol/m2/s


#Check proportion of model best-flux choice for co2 and ch4
#CO2:
co2_all %>% 
  dplyr::select(best.model) %>% 
  summarise(n_incubations=sum(!is.na(best.model)),
            n_nonvalid=sum(best.model=="None appropriate"),
            n_valid=n_incubations-n_nonvalid, 
            prop_nonvalid= round(n_nonvalid/n_incubations*100, 1),
            n_LM=sum(best.model=="LM"),
            propvalid_LM=round(n_LM/n_valid*100, 1),
            n_HM=sum(best.model=="HM"),
            propvalid_HM=round(n_HM/n_valid*100, 1)
            )


#CH4:
ch4_all %>% 
  dplyr::select(best.model) %>% 
  summarise(n_incubations=sum(!is.na(best.model)),
            n_nonvalid=sum(best.model=="None appropriate"),
            prop_nonvalid= round(n_nonvalid/n_incubations*100, 1),
            n_valid=n_incubations-n_nonvalid, 
            n_total.flux=sum(best.model=="total.flux"),
            prop_total.flux=round(n_total.flux/n_valid*100, 1),
            n_LM=sum(best.model=="LM"),
            propvalid_LM=round(n_LM/n_valid*100, 1),
            n_HM=sum(best.model=="HM"),
            propvalid_HM=round(n_HM/n_valid*100, 1)
  )
  


#Select best fluxes and format: UniqueID, ghg_best, ghg_se
#CO2
co2_best<- co2_all %>% 
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
ch4_best<- ch4_all %>% 
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
ghg_best<- co2_best %>% 
  merge.data.frame(ch4_best, by="UniqueID", all=T)

rm(co2_all, co2_best, ch4_all, ch4_best)


#2. Import harm. field Obs.------

#Harmonized field observations per UniqueID
fieldobs_uniqueID<-read.csv(paste0(path_0_sourcedata,"UniqueID_harmonized_field_obs.csv"))

#Harmonized field observations per plotcode 
fieldobs_plotcode<- read.csv(paste0(path_0_sourcedata,"plotcode_harmonized_field_obs.csv"))


#3. Filter valid incubations-----

#Remove fieldobs from non-valid/non-appropriate incubations: Incubations that
#are not representative of subsite definitions, incubations that are not
#representative of natural conditions

# 1. Remove rising-tide and receding-tide incubations. (non appropriate fluxes for comparisons)
tidal_IDs<- fieldobs_uniqueID %>% filter(grepl("rising-tide|receding-tide", field_observations)) %>% pull(UniqueID)

# 2. Remove "bare after vegetation removal" (non-appropriate for comparisons)  
vegcut_IDs<- fieldobs_uniqueID %>% filter(grepl("after vegetation removal", field_observations)) %>% pull(UniqueID)

# 3. Remove incubations in non-appropriate strata (according to subsite definition)
#3.1. Camargue (CA): NOTHING TO REMOVE, all incubations are in appropriate strata
CA_wrongstrata_IDs<- c()

#3.2. Curonian lagoon (CU): inconsistent bare representation (bare sediments
#(beach) wrongly considered for restored sites, never for preserved or altered,
#not appropriate for comparisons). Only vegetated and open water areas are part
#of the sites sensu-stricto.
CU_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="CU") %>% filter(strata=="bare") %>% pull(UniqueID)

#3.3. Danube delta (DA): inconsistent strata in some samplings: exclude
#vegetated of S3-DA-A2 (less than 10% of site, not representative), exclude bare
#of S1-DA-R2 (sampled during first campaign but less than 10% of site)
DA_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="DA") %>% 
  filter((sampling=="S3-DA-A2"&strata=="vegetated")|(sampling=="S1-DA-R2"&strata=="bare")) %>% 
  pull(UniqueID)

#3.4. Dutch delta (DU): NOTHING TO REMOVE, all incubations are in appropriate strata
DU_wrongstrata_IDs<- c()

#3.5. Ria d'Aveiro (RI): rising-tide and "after vegetation removal" already
#removed. Inconsistent incubations: tidal-pool incubations in S3-RI-A1 (only
#considered in S3, less than 10% of site), bare incubations in RI-R1 and in
#RI-R2 (wrong extra-incubations, restored sites are 100% vegetated by
#definition)
RI_wrongstrata_IDs<- fieldobs_uniqueID %>% filter(casepilot=="RI") %>% 
  filter((sampling=="S3-RI-A1"&grepl("tidal-pool", field_observations))|(subsite=="RI-R1"&strata=="bare")|(subsite=="RI-R2"&strata=="bare")) %>% pull(UniqueID)

#3.6. Valencian wetland (VA): NOTHING TO REMOVE, all incubations are in appropriate strata
VA_wrongstrata_IDs<-c()

#ONLY valid incubation details
valid_fieldobs_uniqueID<- fieldobs_uniqueID %>% 
  filter(!UniqueID%in%c(tidal_IDs,vegcut_IDs,CA_wrongstrata_IDs,CU_wrongstrata_IDs,DA_wrongstrata_IDs,DU_wrongstrata_IDs,RI_wrongstrata_IDs,VA_wrongstrata_IDs))


#Check number and proportion of chamber deployments (plotcodes) excluded from further analysis: 
valid_plotcodes<- valid_fieldobs_uniqueID %>% dplyr::select(plotcode) %>% distinct() %>% pull(plotcode)

plotcodes_with_bestflux<- ghg_best %>% 
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
fieldobs_uniqueID %>% 
  filter(plotcode%in%plotcodes_with_bestflux) %>% 
  filter(UniqueID%in%tidal_IDs) %>% #Tidal behaviour
  dplyr::select(plotcode) %>% distinct() %>% dim() 
#Rest of excluded plotcodes are due to other reasons (performed in non representative conditions, outside site boundaries,...)


rm(tidal_IDs,vegcut_IDs,CA_wrongstrata_IDs,CU_wrongstrata_IDs,DA_wrongstrata_IDs,DU_wrongstrata_IDs,RI_wrongstrata_IDs,VA_wrongstrata_IDs)

valid_plotcodes<- valid_fieldobs_uniqueID %>% 
  dplyr::select(plotcode, season, casepilot,subsite,sampling, strata) %>% distinct()



#4. Integrate Daily fluxes ------
#From ghg_best (uniqueID) to ghg_daybest (plotcode)


##4.1. Get daylight hours----

#daylight hours are defined as time from sunrise to sunset (using package
#suncalc, median lat/long of each subsite and date of sampling)

# get median coordinates per subsite (including all seasons)
samplings <- valid_fieldobs_uniqueID %>% 
  dplyr::select(subsite, latitude, longitude) %>% 
  group_by(subsite) %>% 
  summarise(subsite_latitude=median(latitude,na.rm=T),
            subsite_longitude=median(longitude,na.rm=T)) %>% 
  ungroup()

#Use dplyr and suncalc to calculate daylight hours for each subsite-date combination (each sampling).
sampling_daylight <- samplings %>%
  merge.data.frame(y=valid_fieldobs_uniqueID %>% dplyr::select(subsite, date,sampling), by="subsite", all=T) %>% 
  distinct() %>% 
  mutate(date=as.Date(date)) %>% 
  rowwise() %>%
  mutate(
    sunlight_times = list(getSunlightTimes(date = date, lat = subsite_latitude, lon = subsite_longitude, keep = c("sunrise","sunset"))),
    daylight_duration = as.numeric(difftime(sunlight_times$sunset, sunlight_times$sunrise, units = "hours"))
  ) %>%
  ungroup() %>%
  mutate(hours_day=daylight_duration,
         hours_night= 24-hours_day) %>% 
  dplyr::select(sampling, hours_day, hours_night)




## 4.2. Scale to daily flux-----

#Daily integration will follow the previously agreed decisions (common for all gas species): 
#Vegetated fluxes will be systematically scaled with daylight hours (T+D)
#Open water fluxes will be directly used (only D)
#Bare fluxes will be scaled according to data availability: 
#Only 1 instantaneous flux --> used directly for daily flux
#T+D available (due to phytobentos)--> scaling with daylight hours (T+D) = daily flux

#get field details of valid plotcodes
valid_plotcode_obs<- fieldobs_plotcode %>%
  filter(plotcode%in%valid_plotcodes$plotcode) 

#Get ghg fluxes for transparent & dark by plotcode (only for valid UniqueIDs)
valid_ghg_transpdark<- valid_fieldobs_uniqueID %>%
  merge.data.frame(ghg_best, by="UniqueID", all.x = T) %>% 
  dplyr::select(plotcode, co2_best, ch4_best, transparent_dark) %>% 
  rename(co2=co2_best, ch4=ch4_best) %>% 
  pivot_longer(cols=c(co2,ch4), names_to = "ghg", values_to = "flux") %>% 
  pivot_wider(names_from = transparent_dark, values_from = flux)


#Integrate to dailyflux
daily_ghg_plotcode<- valid_plotcode_obs %>% 
  #Add daylight hours
  merge.data.frame(sampling_daylight, by=c("sampling"), all = T) %>% 
  dplyr::select(plotcode, season, casepilot, subsite, sampling, plot_num, date, plot_start_time, latitude, longitude, gas_analyzer, strata, water_depth, field_observations, ABG_veg_description, ABG_veg_biomass_gpersquaremeter, hours_day, hours_night) %>% 
  #Add valid_ghg_transpdark
  merge.data.frame(valid_ghg_transpdark, by="plotcode", all.x = T) %>% 
  #Integrate dailyflux: change units to per day (from per s to per hour, then scale with daylight_duration in hours)
  #For all vegetated plots, scale flux with daylight
  mutate(dailyflux=case_when(strata=="vegetated"~(dark*60*60*hours_night)+(transparent*60*60*hours_day),
                             #For bare with both fluxes, scale with daylight
                             (strata=="bare"&!is.na(dark*transparent))~(dark*60*60*hours_night)+(transparent*60*60*hours_day),
                             #For bare with NA in transparent, use dark for all day.
                             (strata=="bare"&is.na(transparent))~dark*60*60*24,
                             #For bare with NA in dark, use transparent for all day
                             (strata=="bare"&is.na(dark))~transparent*60*60*24,
                             #For open water, use always  dark for all day
                             strata=="open water"~dark*60*60*24)) %>% 
  #Discard transparent & dark ghg fluxes
  dplyr::select(-c(dark, transparent)) %>% 
  pivot_wider(names_from = ghg, values_from = dailyflux) %>% 
  #Transform molar units for daily flux: 
  #CO2: umol per m2 per day --> mol per m2 per day
  #CH4: nmol per m2 per day --> mmol per m2 per day
  mutate(co2_mol=co2*1e-6,
         ch4_mmol=ch4*1e-6,
         #Global warming potential GWP100: 
         gwp_co2=co2_mol*44.0095, #molar flux* molar mass-->g/m2/day
         gwp_ch4=(ch4_mmol*1e-3)*16.04246*27,#molar flux * molar mass* GWP100 factor 27
         gwp_co2andch4= gwp_co2+gwp_ch4) %>%  #Sum of gwp of co2 and ch4 only
  dplyr::select(-c(co2, ch4)) %>% 
  #Add modeling variables: 
  mutate(status=substr(x = sampling, start = 7, stop=7),
         status=case_when(status=="A"~"Altered",
                          status=="P"~"Preserved",
                          status=="R"~"Restored"))


#check how many non-vegetated chambers had 2 fluxes (transparent&dark)
valid_plotcode_obs %>% 
  filter(strata!="vegetated") %>% 
  mutate(n_incubations=(!is.na(dark_UniqueID))+!is.na(transparent_UniqueID)) %>% 
  dplyr::select(plotcode,n_incubations) %>% 
  group_by(n_incubations) %>% 
  summarise(n_chambers=sum(!is.na(plotcode)), .groups = "drop") %>% 
  mutate(total_nonvegetated=sum(n_chambers),
         prop=round(n_chambers/total_nonvegetated*100,1))


#Delete objects already used
rm(samplings, sampling_daylight, fieldobs_plotcode, fieldobs_uniqueID, valid_plotcode_obs, valid_plotcodes)





#5. Final format -----

#Produce final dataset to be used in models. Formatting to simplify filtering for modeling. 

#We will use molar fluxes for individual ghgspecies and CO2eq (gCo2eq /m2/d) for
#GWP of combined ghgspecies. Add column with "units" to be able to
#differentiate.

#We format the dataset to be easily used across our modeling approaches: 
#categories are taken as factors, we pivot longer the daily fluxes (response variable)

data4models <- daily_ghg_plotcode %>% 
  dplyr::select(plotcode, season, casepilot, status, subsite, sampling, strata, 
                co2_mol, #mol per m2 per day
                ch4_mmol, #mmol per m2 per day
                gwp_co2andch4) %>% #gCO2eq per m2 per day
  pivot_longer(
    cols = c(co2_mol, ch4_mmol, gwp_co2andch4),
    names_to = "ghgspecies",
    values_to = "dailyflux"
  ) %>% 
  #Add unitflux variable
  mutate(unitflux=case_when(ghgspecies=="co2_mol"~"mol per m2 per day",
                            ghgspecies=="ch4_mmol"~"mmol per m2 per day",
                            ghgspecies=="gwp_co2andch4"~"gCO2eq per m2 per day")) %>% 
  #Remove units from ghgspecies:
  mutate(ghgspecies=case_when(ghgspecies=="co2_mol"~"co2",
                              ghgspecies=="ch4_mmol"~"ch4",
                              ghgspecies=="gwp_co2andch4"~"gwp_co2andch4")) %>% 
  mutate(across(c(season, casepilot, status, subsite, sampling, strata), as.factor))

#7.Save datasets--------

#Save formated and valid in-situ data to be used in paper. 
write.csv(x = data4models, file = paste0(path_1_paperdata,"ChamberData4paper.csv"),row.names = F)
}

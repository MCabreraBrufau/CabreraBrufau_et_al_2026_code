#GHGpaper_figures

#Date: September 2025

#Description----
#This scrip is used to produce the figures and data-tables included in the manuscript

#Inputs:  In situ data + outputs of modelling
#ChamberData4paper.csv
#Supp_Ancillary_Water_Data.xlsx
#EmmeansCLD_RI_chambermodels.csv
#EmmeasnCLD_rest_chambermodels.csv
#Posthoctests_RI_chambermodels.csv
#Posthoctests_rest_chammbermodels.csv

#Outputs: Figures and tables of manuscript.
#NOTE: Figure_1, Table_1, and supplementary Figure_S1 and Figure_S2 were produced manually outside R. 


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#Installs (if needed) and loads required packages:
required_pkgs <- c("tidyverse",
                   "readxl",
                   "lubridate",
                   "zoo",
                   "purrr",
                   "ggpubr",
                   "data.table",
                   "tools",
                   "hms",
                   "suncalc",
                   "cowplot",
                   "ggforce",
                   "flextable",
                   "officer")

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

#Path to formatted flux data: 
path_1_paperdata <- paste0(path_root,"/1_paperdata/")

#Path to model outputs:
path_2_modeloutputs <- paste0(path_root,"/2_modeloutputs/")

#Path to save figures and tables: 
path_mainfigures<- paste0(path_root,"/3_Figures_and_tables/")
path_supplementary<- paste0(path_root,"/3_Figures_and_tables/Supplementary/")

#Create folders if they do not exist:
if (!dir.exists(path_mainfigures)) {
  dir.create(path_mainfigures, recursive = TRUE)
}

if (!dir.exists(path_supplementary)) {
  dir.create(path_supplementary, recursive = TRUE)
}


{

#____________________#------

#0. Import and format------

#Import Chamberdata4paper.csv
data4paper<-read.csv(file = paste0(path_1_paperdata,"ChamberData4paper.csv"))

#Format main data: 
data4paper<- data4paper %>% 
  #Remove NAs
  filter(!is.na(dailyflux)) %>% 
  #Keep only ghg of interest
  filter(ghgspecies%in%c("co2","ch4","gwp_co2andch4")) %>% 
  #Rename to GWPco2andch4
  mutate(ghgspecies=if_else(ghgspecies=="gwp_co2andch4","GWPco2andch4",ghgspecies)) %>% 
  mutate(vegpresence=if_else(strata=="vegetated","Vegetated","Non-vegetated")) %>% 
  #Factor grouping variables:
  mutate(season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         subsite=factor(subsite, ordered = F),
         ghgspecies=factor(ghgspecies, ordered=F),
         strata=factor(strata, levels = c("open water","bare","vegetated"), ordered = F),
         sampling=factor(sampling, ordered = F),
         vegpresence=factor(vegpresence, levels = c("Non-vegetated","Vegetated"), ordered = F))


##Outputs of Simple models------
#Model results for RI casepilot (only status*season as fixed effects)


#Import results of models (R2 and significance)
simplemodel_summary<- read.csv(paste0(path_2_modeloutputs, "Summary_RI_chambermodels.csv"))

simplemodel_significances_alleffects<- simplemodel_summary %>% 
  dplyr::select(dataset, status_pval, season_pval, status.season_pval) %>% 
  separate(dataset, into = c("casepilot","ghgspecies")) %>% 
  pivot_longer(c(status_pval,season_pval,status.season_pval),names_sep = "_",values_to = "pvalue_effect", names_to = c("effect","drop")) %>% 
  dplyr::select(-drop)


#Import Emmeans +SE + cld groupletters 
simplemodel_emmeans_all<-read.csv(paste0(path_2_modeloutputs,"EmmeansCLD_RI_chambermodels.csv"))

#Format Emmeans: filter for RI only & remove CLD group-letters from comparisons
#when all are the same.
RI_simplemodel_emmeans_all<- simplemodel_emmeans_all %>% 
  filter(casepilot=="RI") %>% 
  group_by(ghgspecies, casepilot, comparison) %>% 
  #If there are no differences (same letter for all levels of a given
  #comparison), remove the cld_group letter
  mutate(cld_group = if (n_distinct(cld_group) == 1) "" else cld_group) %>%
  ungroup() %>% 
  mutate(ghgspecies=factor(ghgspecies, ordered=F),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T)) %>% 
  arrange(casepilot, ghgspecies, comparison, season,status)


#separate emmeans per type of comparison
#status
RI_simplemodel_emmeans_status<- RI_simplemodel_emmeans_all %>% 
  filter(comparison=="status")

#season
RI_simplemodel_emmeans_season<- RI_simplemodel_emmeans_all %>% 
  filter(comparison=="season")

#status_within_season, re-filter to remove letters when all status are equal
#under a particular season
RI_simplemodel_emmeans_statuswithinseason<- RI_simplemodel_emmeans_all %>% 
  filter(comparison=="status_within_season") %>% 
  group_by(ghgspecies, casepilot, comparison,season) %>% 
  #If there are no differences (same letter for all levels of a given comparison), remove the cld_group letter 
  mutate(cld_group = if (n_distinct(cld_group) == 1) "" else cld_group)





##Outputs of Complex models------
#Complex models refer to models that consider status*season*vegpresence as fixed
#effects (Appropriate for all casepilot except RI).


#Import results from models (R2 and significances)
complexmodel_summary<- read.csv(paste0(path_2_modeloutputs, "Summary_rest_chambermodels.csv"))

head(complexmodel_summary)

#Get effect significances in long format
complexmodel_significances_alleffects<- complexmodel_summary %>% 
  dplyr::select(dataset, status_pval, season_pval, vegpresence_pval, status.season_pval, status.vegpresence_pval, season.vegpresence_pval, status.season.vegpresence_pval) %>% 
  separate(dataset, into = c("casepilot","ghgspecies")) %>% 
  pivot_longer(c(status_pval, season_pval, vegpresence_pval, status.season_pval, status.vegpresence_pval, season.vegpresence_pval, status.season.vegpresence_pval),names_sep = "_",values_to = "pvalue_effect", names_to = c("effect","drop")) %>% 
  dplyr::select(-drop)

head(complexmodel_significances_alleffects)


#Import Emmeans + cld groupletters (all ghgs Combined into single table)
#Calculated with custom weights for vegpresence. 
complexmodel_emmeans_all<-read.csv(paste0(path_2_modeloutputs,"EmmeansCLD_rest_chambermodels.csv"))

#Remove CLD when all are the same within a given comparison (no differences)
complexmodel_emmeans_all<- complexmodel_emmeans_all %>% 
  group_by(ghgspecies,casepilot,comparison) %>% 
  #If there are no differences (same letter for all status of a given ghgspecies*casepilot combination), remove the cld_group 
  mutate(cld_group = if (n_distinct(cld_group) == 1) "" else cld_group) %>%
  ungroup() %>% 
  mutate(ghgspecies=factor(ghgspecies, ordered=F),
         casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered = T),
         season=factor(season, levels = c("S1","S2","S3","S4"), ordered = T),
         vegpresence=factor(vegpresence, levels = c("Non-vegetated","Vegetated"), ordered = F)) %>% 
  arrange(casepilot, ghgspecies, comparison, season,status,vegpresence)


#separate emmeans per type of comparison:
#status
complexmodel_emmeans_status<- complexmodel_emmeans_all %>% 
  filter(comparison=="status")

#season
complexmodel_emmeans_season<- complexmodel_emmeans_all %>% 
  filter(comparison=="season") 

#status_within_season, need to re-filter to eliminate letters when all are the
#same at the specific season
complexmodel_emmeans_statuswithinseason<- complexmodel_emmeans_all %>% 
  filter(comparison=="status_within_season") %>% 
  group_by(ghgspecies,casepilot,comparison,season) %>% 
  #If there are no differences (same letter for all levels within season),
  #remove the cld_group
  mutate(cld_group = if (n_distinct(cld_group) == 1) "" else cld_group)

#status_within_vegpresence, need to re-filter to eliminate letters when all are
#the same at the specific vegpresence group
complexmodel_emmeans_statuswithinvegpresence<- complexmodel_emmeans_all %>% 
  filter(comparison=="status_within_vegpresence") %>% 
  group_by(ghgspecies,casepilot,comparison,vegpresence) %>% 
  #If there are no differences (same letter for all levels within season),
  #remove the cld_group
  mutate(cld_group = if (n_distinct(cld_group) == 1) "" else cld_group)


##Combine all bestmodel results----
#Combine RI results with rest of casepilot models
#Get emmean (back-transformed) + cld_group for every relevant comparison: 

#status
bestmodel_emmeans_status<- RI_simplemodel_emmeans_status %>% 
  full_join(complexmodel_emmeans_status) %>% 
  dplyr::select(ghgspecies, casepilot, comparison, status, emmean_bt,SE_bt,lower.CL_bt, upper.CL_bt, cld_group)

str(bestmodel_emmeans_status)

#season
bestmodel_emmeans_season<- RI_simplemodel_emmeans_season %>% 
  full_join(complexmodel_emmeans_season) %>% 
  dplyr::select(ghgspecies, casepilot, comparison, season, emmean_bt,SE_bt,lower.CL_bt, upper.CL_bt, cld_group)
str(bestmodel_emmeans_season)

#status_within_season
bestmodel_emmeans_statuswithinseason<- RI_simplemodel_emmeans_statuswithinseason %>% 
  full_join(complexmodel_emmeans_statuswithinseason) %>% 
  dplyr::select(ghgspecies, casepilot, comparison, status,season, emmean_bt,SE_bt, lower.CL_bt, upper.CL_bt, cld_group)

str(bestmodel_emmeans_statuswithinseason)

#status_within_vegpresence, NO need to join RI (no comparison possible)
bestmodel_emmeans_statuswithinvegpresence<- complexmodel_emmeans_statuswithinvegpresence %>% 
  dplyr::select(ghgspecies, casepilot, comparison, status,vegpresence, emmean_bt,SE_bt, lower.CL_bt, upper.CL_bt, cld_group)

str(bestmodel_emmeans_statuswithinvegpresence)

#____________________#------
#1. In-text descriptive stats----

##General sampling stats-----
#TOTAL number of non-na daily fluxes that compose the final dataset presented in the article: 
#n for CO2, CH4 and GWP
data4paper %>% 
 dplyr::select(plotcode,dailyflux, ghgspecies) %>% 
  pivot_wider(names_from = ghgspecies, values_from = dailyflux) %>% 
  summarise(n_co2=sum(!is.na(co2)),
            n_ch4=sum(!is.na(ch4)),
            n_gwp=sum(!is.na(GWPco2andch4)))
  
#Chamber deployments
plotcodes<- read.csv(paste0(path_0_sourcedata, "Plotcode_harmonized_field_obs.csv"))

#Average number of deployments per sampling across sites
plotcodes %>% 
group_by(sampling) %>% 
  summarise(n_chambers=n_distinct(plotcode)) %>% 
  ungroup() %>% 
  summarise(avg_nchambers=mean(n_chambers),
            sd_nchambers=sd(n_chambers),
            median_nchambers=median(n_chambers)) 

##Zostera coverage in RI------
#Chamber-derived proportion of vegetated area in RI per status (excluding
#rising-tide, only considering low-tide chambers in paper)
plotcodes %>%
  filter(plotcode%in%data4paper$plotcode) %>% 
  filter(casepilot=="RI") %>% 
  mutate(status=case_when(grepl("-A[0-9]",subsite)~"Altered",
                          grepl("-P[0-9]",subsite)~"Preserved",
                          grepl("-R[0-9]",subsite)~"Restored")) %>% 
  dplyr::select(strata, status, plotcode,ABG_veg_description) %>% 
  distinct() %>% 
  group_by(status) %>% 
  summarise(n_chambers=sum(!is.na(plotcode)),
            n_zostera=sum(grepl("Zostera",ABG_veg_description)),
            n_vegetated=sum(strata=="vegetated"),
            n_bare=sum(strata=="bare"),
            n_water=sum(strata=="open water")) %>% 
  mutate(percent_bare=round(n_bare/n_chambers*100,2),
         percent_vegetated=round(n_vegetated/n_chambers*100,2),
         percent_zostera=round(n_zostera/n_chambers*100,2))


##Inundation pattents-----
#Frequency of inundated (water_depth>0) chambers between preserved, altered and restored sites.

#In CA, to support different flooding patterns affecting CH4 emissions.
#Overall proportion: 
plotcodes %>% 
  filter(casepilot=="CA") %>% 
  mutate(inundated=if_else(water_depth>0,T,F),
         status=case_when(grepl("P[0-9]",subsite)~"Preserved",
                          grepl("A[0-9]",subsite)~"Altered",
                          grepl("R[0-9]",subsite)~"Restored")) %>% 
  group_by(status,season) %>% 
  mutate(season_chambers=sum(!is.na(inundated)),
         season_floodprop=sum(inundated)/sum(!is.na(inundated))) %>% 
  ungroup() %>% 
  group_by(status) %>% 
  mutate(yearly_chambers=sum(!is.na(inundated)), 
         yearly_floodprop=sum(inundated)/sum(!is.na(inundated))) %>% 
  group_by(status,season) %>% 
  summarise(
    season_chambers=unique(season_chambers),
    season_floodprop=unique(season_floodprop),
    yearly_chambers=unique(yearly_chambers),
    yearly_floodprop=unique(yearly_floodprop))

#In VA: 
plotcodes %>% 
  filter(casepilot=="VA") %>% 
  mutate(inundated=if_else(water_depth>0,T,F),
         status=case_when(grepl("P[0-9]",subsite)~"Preserved",
                          grepl("A[0-9]",subsite)~"Altered",
                          grepl("R[0-9]",subsite)~"Restored")) %>% 
  group_by(status,season) %>% 
  mutate(season_chambers=sum(!is.na(inundated)),
         season_floodprop=sum(inundated)/sum(!is.na(inundated))) %>% 
  ungroup() %>% 
  group_by(status) %>% 
  mutate(yearly_chambers=sum(!is.na(inundated)), 
         yearly_floodprop=sum(inundated)/sum(!is.na(inundated))) %>% 
  group_by(status,season) %>% 
  summarise(
    season_chambers=unique(season_chambers),
    season_floodprop=unique(season_floodprop),
    yearly_chambers=unique(yearly_chambers),
    yearly_floodprop=unique(yearly_floodprop))




##GHG importance on CO2eq------

#In data4paper, units are: Co2 (molperm2perday), ch4(millimolper m2 per day),
#GWPco2andch4(gCO2eq per m2 per day)

#Co2 to CO2 eq: mutliply by 44 to get gCO2eq
#CH4 to CO2eq: divide by 1000 & multiply by 16 to get gCH4, then multiply by 28 to get CO2eq. 

ghgprop<- data4paper %>% 
 dplyr::select(season, casepilot, subsite, sampling, status,plotcode, ghgspecies, dailyflux) %>% 
  pivot_wider(names_from = ghgspecies, values_from = dailyflux) %>% 
  mutate(co2=abs(co2*44),
         ch4=abs((ch4/1000)*16*28),
         Co2eq_absolutesum=abs(co2)+abs(ch4),
         co2_prop=co2/Co2eq_absolutesum,
         ch4_prop=ch4/Co2eq_absolutesum) %>% 
  group_by(casepilot) %>% 
  mutate(avg_cp_ch4prop=mean(ch4_prop,na.rm=T),
         sd_cp_ch4prop=sd(ch4_prop, na.rm=T),
         n_cp_ch4prop=sum(!is.na(ch4_prop))) %>% 
  group_by(casepilot, status) %>% 
  mutate(avg_cpstatus_ch4prop=mean(ch4_prop,na.rm=T),
         sd_cpstatus_ch4prop=sd(ch4_prop, na.rm=T),
         n_cpstatus_ch4prop=sum(!is.na(ch4_prop)))

#CH4 proportion of CO2eq in Preserved sites
ghgprop %>% 
 dplyr::select(casepilot, status, avg_cpstatus_ch4prop,sd_cpstatus_ch4prop,n_cpstatus_ch4prop,
         avg_cp_ch4prop, sd_cp_ch4prop, n_cp_ch4prop) %>% 
  distinct() %>% 
  filter(status=="Preserved") %>% 
  transmute(casepilot, status, avg_cpstatus_ch4prop, sd_cpstatus_ch4prop, se_cpstatus_ch4prop=sd_cpstatus_ch4prop/sqrt(n_cpstatus_ch4prop))


#CH4 proportion across casepilots (including all Conservation status)
ghgprop %>% 
  ungroup() %>% 
 dplyr::select(casepilot, avg_cp_ch4prop, sd_cp_ch4prop, n_cp_ch4prop) %>% 
  distinct()


#____________________--------

#2. Main Figures------

##Fig. 1 (Map)-----

#Map produced outside R-enviroment


##Fig. 2 (Preserved Fluxes)-----
#A 3-panel vertical plot. 1 panel per gas. Boxplot (without outliers). Linear scale. 
#For CH4, inset with DU,RI,CA,VA zoomed in. 

base_co2_ident<-
  data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("co2")) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux*1000)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(Daily~CO[2]~flux~(mmol~m^-2~d^-1)),
       x=NULL,
       fill=paste0("Status"))



# Main plot
base_ch4_ident_main <- data4paper %>%
  filter(status == "Preserved", ghgspecies %in% c("ch4")) %>%
  ggplot(aes(x = casepilot, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = FALSE, fill = "#00BA38", size = 0.7) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  ) +
  guides(color = "none", fill = "none") +
  labs(
    y = expression(Daily~CH[4] ~ flux ~ (mmol~m^-2~d^-1)),
    x = NULL,
    fill = "Status"
  )

# Inset plot
inset_plot <- data4paper %>%
  filter(status == "Preserved", ghgspecies %in% c("ch4"),
         casepilot %in% c("DU", "RI", "CA", "VA")) %>%
  ggplot(aes(x = casepilot, y = dailyflux)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.2, outliers = FALSE, fill = "#00BA38", size = 0.7) +
  theme_bw(base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 6, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

#Create the grob of your inset plot
inset_grob <- ggplotGrob(inset_plot)

#define position of inset and geom_rect (in main-plot scale)
xmin<- 1
xmax<- 4
ymin<- 5
ymax<- 20

#Extract and store the y-range (without outliers) of the main plot (they are
#overriden when using annotation_custom
y_range <- range(ggplot_build(base_ch4_ident_main)$layout$panel_params[[1]]$y.range)


#Combine plot with inset (using positions and limits defined above) using annotation_custom
base_ch4_ident_comb <- base_ch4_ident_main +
  annotation_custom(
    grob = inset_grob,
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax
  ) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "black", fill = NA) +
  coord_cartesian(ylim = c(y_range[1], y_range[2]))


  
base_gwp_ident<-data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("GWPco2andch4")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"))) %>% 
  ggplot(aes(x=casepilot, y=(dailyflux)))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(width=0.2,outliers = F, fill="#00BA38",size = 0.7)+
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))+
  guides(color = "none", fill="none")+
  labs(y= expression(Daily~CO[2]*eq~flux~(g~CO[2]*eq~m^-2~d^-1)),
       x=paste0("Case pilot"),
       fill=paste0("Status"))


#COMBINE all 3 panels (with CH4 inset)
baseline_3sp_ident_withinset <- plot_grid(
  base_co2_ident,
  base_ch4_ident_comb,
  base_gwp_ident,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("(a)", "(b)", "(c)"),label_x =0.165,label_y = 0.975,label_size = 13
  # label_size = 12,
  # label_x = 0.02,  # move label closer to left
  # label_y = 1      # top-aligned
)

#Save:
ggsave(plot = baseline_3sp_ident_withinset, 
       filename = "Figure_2_preserved_fluxes.png",
       path = path_mainfigures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 200, 
       width = 90)

##In-text Avg preserved fluxes----

data4paper %>% 
  filter(status=="Preserved", ghgspecies%in%c("co2","ch4","GWPco2andch4")) %>% 
  group_by(casepilot, ghgspecies) %>% 
  #Transform to consistent unit (co2 and ch4 in mmol, gwp in gCO2eq)
  mutate(dailyflux=if_else(ghgspecies=="co2",dailyflux*1000, dailyflux),
         unitflux=if_else(ghgspecies=="co2","mmol per m2 per day",  unitflux)) %>% 
  summarise(avg_flux=mean(dailyflux, na.rm=T),
            sd_flux=sd(dailyflux, na.rm = T),
            median_flux=median(dailyflux, na.rm = T),
            q1= quantile(dailyflux, 0.25, na.rm=T),
            q3= quantile(dailyflux, 0.75, na.rm=T),
            iqr = IQR(dailyflux, na.rm = T),
            unit_flux=unique(unitflux)) %>% 
  arrange(ghgspecies, median_flux)



##Function for Status plots------

#Preserved-Altered-Restored fluxes (boxplots, no without data outside 1.5*IQR
#"outliers") and with group differences (EMMeans as diamonds, letters from
#significant post-hoc differences), for each casepilot (vertical pannels).

#We create a single plot template function (for which we will change the
#data-origin and yaxis label for each ghgspecies). It takes 4 arguments: GHG
#("co2", "ch4","GWPco2andch4"),  y_label (the expression call to use for the y
#axis), y_multiplier: (multiplier of dailyflux to use, 1000 for CO2 in mmol, 1
#for CH4mmol, 1 for GWP in gCO2eq), drawboxplot_outliers (T/F)

plot_ghg_faceted <- function(GHG,
                             y_label, 
                             y_multiplier = 1,
                             drawboxplot_outliers=F) {
  
  #filter data for GHG
  data<- data4paper %>% filter(ghgspecies==GHG)
  emmeans_data<- bestmodel_emmeans_status %>% filter(ghgspecies==GHG) 
  panel_label<- bestmodel_emmeans_status %>% filter(ghgspecies==GHG) %>% 
    distinct(casepilot) %>% 
    arrange(casepilot) %>% 
    mutate(label=c("(a)","(b)","(c)","(d)","(e)","(f)"), x=-Inf, y=Inf)
  
  #Produce plot: 
  ggplot(data, aes(x = status, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(width = 0.2, outliers = drawboxplot_outliers, aes(fill = status), size = 0.7) +
    #Add panel letters
    geom_text(data = panel_label,  aes(x, y, label = label),
      inherit.aes = FALSE,  hjust = -0.2, vjust = 1.3,
      fontface = "bold",
      size = 4
    )+
    # Add emmeans points
    geom_point(data = emmeans_data, aes(x = status, y = emmean_bt * y_multiplier),
               shape = 23, size = 3, fill = "black") +
    
    # Add group letters
    geom_text(data = emmeans_data, aes(x = status, label = cld_group, y = Inf),
              vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
              size = 5, fontface = "bold") +
    theme_bw() +
    scale_fill_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    scale_color_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
    theme(
      axis.text = element_text(face = "bold")
    ) +
    guides(color = "none", fill = "none") +
    labs(
      y = y_label,
      x = "Conservation status",
      fill = "Status"
    ) +
    facet_grid(rows = vars(casepilot), scales = "free")
}


##Fig. 3 (CO2 status)-----

#Produce plot using status template and subset data for CO2: 
plot_ghg_faceted(GHG="co2",
                 y_label = expression(Daily~CO[2]~flux~(mmol~m^-2~d^-1)),
                 y_multiplier = 1000)


#Save plot: 
ggsave(filename = "Figure_3_CO2perStatus.png",
       path = path_mainfigures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)


##Fig. 4 (CH4 status)-----

#Produce plot using status template and subset data for CH4: 
plot_ghg_faceted(GHG="ch4",
                 y_label = expression(Daily~CH[4]~flux~(mmol~m^-2~d^-1)),
                 y_multiplier = 1)

#Save plot: 
ggsave(filename = "Figure_4_CH4perStatus.png",
       path = path_mainfigures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)


##Fig. 5 (CO2eq status)-----

#Produce plot using status template and subset data for CO2eq: 
plot_ghg_faceted(GHG="GWPco2andch4",
                 y_label = expression(Daily~CO[2]*eq~flux~(g~CO[2]*eq~m^-2~d^-1)),
                 y_multiplier = 1)

#Save plot: 
ggsave(filename = "Figure_5_CO2eqperStatus.png",
       path = path_mainfigures,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 90)



#____________________--------
#3. Supplementary Figures------

#Fig. S1 (wetland pictures)------
#Produced manually outside R

#Fig. S2 (chamber types) ------
#Produced manually outside R


#Create function for status per season flux figures: 

##Plot both data and emmeans across status*season for each casepilot. 
##Color for status, x-position for season. 
##Points for emmeans and group-letters for status_within_season comparisons
plot_status_per_season <- function(GHG, y_label, 
                                   y_multiplier = 1, 
                                   dodge_width = 0.4) {
  #Filter datasets for GHG:
  data<- data4paper %>% filter(ghgspecies==GHG)%>% 
    mutate(season=case_when(season=="S1"~"Autumn",
                            season=="S2"~"Winter",
                            season=="S3"~"Spring",
                            season=="S4"~"Summer")) %>% 
    mutate(season=factor(season, levels=c("Autumn", "Winter","Spring","Summer"), ordered = T))
  
  emmeans_data<- bestmodel_emmeans_statuswithinseason %>% filter(ghgspecies==GHG)%>% 
    mutate(season=case_when(season=="S1"~"Autumn",
                            season=="S2"~"Winter",
                            season=="S3"~"Spring",
                            season=="S4"~"Summer")) %>% 
    mutate(season=factor(season, levels=c("Autumn", "Winter","Spring","Summer"), ordered = T))
  
  panel_label<- bestmodel_emmeans_statuswithinseason %>% filter(ghgspecies==GHG) %>% 
    ungroup() %>%     dplyr::select(casepilot) %>% 
    distinct(casepilot) %>% 
    arrange(casepilot) %>% 
    mutate(label=c("(a)","(b)","(c)","(d)","(e)","(f)"), x=-Inf, y=Inf)
  
  
  #Produce plots:
  ggplot(data, aes(x = season, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # Boxplots with dodge
    geom_boxplot(aes(fill = status), width = 0.2, outliers = F,
                 position = position_dodge(width = dodge_width), size = 0.7) +
    # Emmeans points with same dodge
    geom_point(data = emmeans_data,
               aes(x = season, y = emmean_bt * y_multiplier, group = status),
               shape = 23, size = 2, fill = "black",
               position = position_dodge(width = dodge_width)) +
    #Add panel letters
    geom_text(data = panel_label,  aes(x, y, label = label),
              inherit.aes = FALSE,  hjust = -0.2, vjust = 1.3,
              fontface = "bold",
              size = 4
    )+
    # Add group letters with same dodge
    geom_text(data = emmeans_data, 
              aes(x = season, label = cld_group, y = Inf, group=status),
              vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
              size = 4.5, fontface = "bold",
              position = position_dodge(width = dodge_width)) +
    #expand y-scale to avoid overlap between group_letters and wiskers
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
    theme_bw() +
    scale_fill_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    theme(
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5)
    ) +
    guides(color = "none") +
    labs(
      y = y_label,
      x = "Season",
      fill = "Status"
    ) +
    facet_grid(rows = vars(casepilot), scales = "free") 
  
}



#Fig. S3: (CO2 status&season)------ 
#Plot CO2: 
plot_status_per_season(GHG = "co2",
                       y_label = expression(Daily~CO[2]~flux~(mmol~m^-2~d^-1)),
                       y_multiplier = 1000, #to change from molar to milli-molar units
                       dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Figure_S3_CO2perStatusandSeason.png",
       path = path_supplementary,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)


#Fig. S4: (CH4 status&season)------
#Plot CH4: 
plot_status_per_season(GHG = "ch4",
                       y_label = expression(Daily~CH[4]~flux~(mmol~m^-2~d^-1)),
                       y_multiplier = 1,
                       dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Figure_S4_CH4perStatusandSeason.png",
       path = path_supplementary,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)

#Fig. S5: (CO2eq status&season)------
#Plot CO2eq: 
plot_status_per_season(GHG = "GWPco2andch4",
                       y_label = expression(Daily~CO[2]*eq~flux~(g*CO[2~eq]~m^-2~d^-1)),
                       y_multiplier = 1,
                       dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Figure_S5_CO2eqperStatusandSeason.png",
       path = path_supplementary,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)


#Create function to plot status per vegpresence:
  #Plot both data and emmeans across status*vegpresence for each casepilot that can support it (not RI). 
  #Color for status, x-position for vegpresence 
  #Points for emmeans and group-letters for status_within_vegpresence comparisons
plot_status_per_vegpresence <- function(GHG, y_label, 
                                        y_multiplier = 1, 
                                        dodge_width = 0.4) {
  #Filter datasets for GHG:
  data<- data4paper %>% filter(ghgspecies==GHG) %>% filter(casepilot!="RI")
  emmeans_data<- complexmodel_emmeans_statuswithinvegpresence %>% filter(ghgspecies==GHG) %>% filter(casepilot!="RI")
  panel_label<- complexmodel_emmeans_statuswithinvegpresence %>% filter(ghgspecies==GHG) %>% 
    ungroup() %>%     dplyr::select(casepilot) %>% 
    distinct(casepilot) %>% 
    arrange(casepilot) %>% 
    mutate(label=c("(a)","(b)","(c)","(d)","(e)"), x=-Inf, y=Inf)
  
  
  #Produce plots:
  ggplot(data, aes(x = vegpresence, y = dailyflux * y_multiplier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # Boxplots with dodge
    geom_boxplot(aes(fill = status), width = 0.2, outliers = F,
                 position = position_dodge(width = dodge_width), size = 0.7) +
    #Add panel letters
    geom_text(data = panel_label,  aes(x, y, label = label),
              inherit.aes = FALSE,  hjust = -0.2, vjust = 1.3,
              fontface = "bold",
              size = 4)+
    # Emmeans points with same dodge
    geom_point(data = emmeans_data,
               aes(x = vegpresence, y = emmean_bt * y_multiplier, group = status),
               shape = 23, size = 2, fill = "black",
               position = position_dodge(width = dodge_width)) +
    # Add group letters with same dodge
    geom_text(data = emmeans_data, 
              aes(x = vegpresence, label = cld_group, y = Inf, group=status),
              vjust = 1.2, hjust = 0.5, color = "black", inherit.aes = FALSE,
              size = 4.5, fontface = "bold",
              position = position_dodge(width = dodge_width)) +
    #expand y-scale to avoid overlap between group_letters and wiskers
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2), 0)) +
    theme_bw() +
    scale_fill_manual(values = c(
      "Preserved" = "#009E73",
      "Altered"   = "#D55E00",
      "Restored"  = "#56B4E9"
    )) +
    theme(
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5)
    ) +
    guides(color = "none") +
    labs(
      y = y_label,
      x = "Vegetation presence",
      fill = "Status"
    ) +
    facet_grid(rows = vars(casepilot), scales = "free") 
  
}



#Fig. S6: (CO2 status&vegpresence)------ 
#CO2: 
plot_status_per_vegpresence(GHG = "co2",
                            y_label = expression(Daily~CO[2]~flux~(mmol~m^-2~d^-1)),
                            y_multiplier = 1000, #to change from molar to milli-molar units
                            dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Figure_S6_CO2perStatusandVegpresence.png",
       path = path_supplementary,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)


#Fig. S7: (CH4 status&vegpresence)------ 
#CH4: 
plot_status_per_vegpresence(GHG = "ch4",
                            y_label = expression(Daily~CH[4]~flux~(mmol~m^-2~d^-1)),
                            y_multiplier = 1,
                            dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Figure_S7_CH4perStatusandVegpresence.png",
       path = path_supplementary,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)


#Fig. S8: (CO2eq status&vegpresence)------ 
#CO2eq: 
plot_status_per_vegpresence(GHG = "GWPco2andch4",
                            y_label = expression(Daily~CO[2]*eq~flux~(g*CO[2~eq]~m^-2~d^-1)),
                            y_multiplier = 1,
                            dodge_width = 0.6)

#Save plot: 
ggsave(filename = "Figure_S8_CO2eqperStatusandVegpresence.png",
       path = path_supplementary,
       device = "png",
       dpi = 400,
       units = "mm",
       height = 229,
       width = 180)










#____________________--------
#4. Tables--------

#Main tables: 
  #1. Site descriptions (manually created)

#Supplementary tables output word files: 
  #1. Model summary (R-based)
  #2. Emmeans status (R-based)
  #3. Contrasts post-hoc status (R-based, reordered contrasts for clarity)
  #4. Water sample summary (R-based)


#Table 1: Site descriptions------
#Manually filed. 


#Table S1: Models summary------
#Supplementary table with model summaries 18 models (6CP x 3 GHG). 

#Format model summary into table with effect significance, R2c and R2m
#Arrange by co2,ch4,GWP, then casepilot (DU,RI,CA,VA,DA,CU) , then effect (status,season,interaction).

#RI model: 
formated_summary_simplemodels<- simplemodel_summary %>% 
  dplyr::select(dataset, transformation,formula, family,nobs, status_pval, season_pval, status.season_pval, homoced_R2m, homoced_R2c) %>% 
  pivot_longer(cols=c(status_pval, season_pval, status.season_pval), names_to = c("effect","drop"),names_sep = "_", values_to = "p_value") %>% 
  separate(dataset, into = c("casepilot", "ghgspecies")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         ghgspecies=factor(ghgspecies, levels = c("co2","ch4","GWPco2andch4"), ordered = T),
         effect=factor(effect, levels = c("status","season","status.season"), ordered = T)) %>% 
  rename(R2c=homoced_R2c, R2m=homoced_R2m) %>% 
  dplyr::select(ghgspecies, casepilot, transformation,formula,family, nobs, R2m, R2c, effect, p_value) %>% 
  filter(casepilot=="RI") %>% 
  mutate(R2m=round(R2m,digits = 3), R2c=round(R2c, digits = 3)) %>% 
  arrange(ghgspecies, casepilot, effect) %>% 
  mutate(rounded_p_value=if_else(condition = p_value<0.001,
                                 true= "< 0.001",
                                 false=as.character(round(p_value, digits = 3))))


#FOR complex models
formated_summary_complexmodels<- complexmodel_summary %>% 
  dplyr::select(dataset, transformation,formula, family,nobs, 
                status_pval, season_pval,vegpresence_pval, status.season_pval, status.vegpresence_pval, season.vegpresence_pval,status.season.vegpresence_pval, 
                homoced_R2m, homoced_R2c) %>% 
  pivot_longer(cols=c(status_pval, season_pval,vegpresence_pval, status.season_pval, status.vegpresence_pval, status.season.vegpresence_pval), names_to = c("effect","drop"),names_sep = "_", values_to = "p_value") %>% 
  separate(dataset, into = c("casepilot", "ghgspecies")) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         ghgspecies=factor(ghgspecies, levels = c("co2","ch4","GWPco2andch4"), ordered = T),
         effect=factor(effect, levels = c("status","season","vegpresence","status.season","status.vegpresence","season.vegpresence", "status.season.vegpresence"), ordered = T)) %>% 
  rename(R2c=homoced_R2c, R2m=homoced_R2m) %>% 
  dplyr::select(ghgspecies, casepilot, transformation,formula,family, nobs, R2m, R2c, effect, p_value) %>% 
  mutate(R2m=round(R2m,digits = 3), R2c=round(R2c, digits = 3)) %>% 
  arrange(ghgspecies, casepilot, effect)%>% 
  mutate(rounded_p_value=if_else(condition = p_value<0.001,
                                 true= "< 0.001",
                                 false=as.character(round(p_value, digits = 3))))

#JOIN summaries of all models and save
combined_model_summary<- rbind(formated_summary_simplemodels,formated_summary_complexmodels) %>% 
  arrange(ghgspecies, casepilot, effect)
rm(formated_summary_complexmodels, formated_summary_simplemodels)

# write.csv(combined_model_summary, 
#           file=paste0(main_figures, "Formated_summary_models.csv"), row.names = F)


#Format summary table into word-table:
wordtable_model_summary<- combined_model_summary %>% 
 dplyr::select(ghgspecies, casepilot, transformation, formula, family, nobs, R2m, R2c, effect, rounded_p_value) %>% 
  mutate(effect = factor(gsub("\\."," : ",effect), levels = c("status","season","vegpresence","status : season","status : vegpresence","season : vegpresence","status : season : vegpresence"), ordered = T)) %>% 
  arrange(ghgspecies, casepilot, effect) %>% 
  mutate(transformation=case_when(transformation=="best_pseudo_log"~"pseudo-log",
                                  transformation=="log_x"~"log", 
                                  transformation=="arcsinh_x"~"arcsinh",
                                  transformation=="yeojohnson"~"Yeo-Johnson"),
         formula=gsub("dailyflux_trans","Flux",formula),
         ghgspecies=case_when(ghgspecies=="co2"~"CO2", 
                              ghgspecies=="ch4"~"CH4",
                              ghgspecies=="GWPco2andch4"~ "CO2eq"),
         formula= paste0("Call: ",formula, 
                         ", \nDistribution: ",family,
                         ", \nTransformation: ",transformation), 
         dataset=paste(casepilot, ghgspecies, sep = " - ")) %>% 
  dplyr::select(dataset, formula, nobs, R2m, R2c, effect, rounded_p_value)


ft<- flextable(wordtable_model_summary) %>% 
  merge_v( j= c("dataset","formula","nobs","R2m","R2c"), combine = T) %>% 
  set_header_labels(
    dataset = "Dataset",
    formula = "Best-Supported Model",
    nobs = "N",
    R2m = "R2m",
    R2c = "R2c",
    effect = "Effect",
    rounded_p_value = "p-Value"
  ) %>% 
  bold(part = "header") %>%
  fontsize(size=10, part = "body") %>% 
  theme_vanilla() %>% 
  autofit()


doc <- read_docx() %>%
  body_add_par("Table S1. Model structures and effects significance", style = "heading 1") %>%
  body_add_flextable(ft) %>% 
  body_end_section_landscape()

print(doc, 
      target=paste0(path_supplementary,"Table_S1_ModelSummaries.docx"))



#Table S2: Emmeans status-----
#Table with full emmeans for each status, casepilot and ghg

formated_emmeans<- bestmodel_emmeans_status %>%
  #change units co2  to mmol 
  mutate(emmean=if_else(ghgspecies=="co2", emmean_bt*1000,emmean_bt),
         SE=if_else(ghgspecies=="co2", SE_bt*1000,SE_bt),
         lower.CL=if_else(ghgspecies=="co2", lower.CL_bt*1000,lower.CL_bt),
         upper.CL=if_else(ghgspecies=="co2", upper.CL_bt*1000,upper.CL_bt)) %>% 
  mutate(emmean_SE=sprintf("%.3g ± %.3g", emmean, SE),# 3 significant figures (normal or scientific notation)
         CI95= sprintf("%.3g to %.3g", lower.CL, upper.CL)) %>% # 3 significant figures (normal or scientific notation)
  dplyr::select(c(ghgspecies, casepilot, status, emmean_SE, CI95)) %>% 
  mutate(ghgspecies=if_else(ghgspecies=="GWPco2andch4", "gwp",ghgspecies)) %>% 
  pivot_wider(names_from = ghgspecies, values_from = c("emmean_SE", "CI95")) %>% 
  dplyr::select(casepilot, status, emmean_SE_co2, CI95_co2, emmean_SE_ch4, CI95_ch4,emmean_SE_gwp, CI95_gwp)



ft_emmeans<- formated_emmeans %>% 
  flextable() %>% 
  add_header_row(
    values = c("Case pilot","Status","CO2 flux (mmol m-2 d-1)","CH4 flux (mmol m-2 d-1)","GWP flux (g CO2eq. m-2 d-1)"),  
    colwidths = c(1,1,2,2,2)                          # number of columns each value spans
  ) %>%
  set_header_labels(
    casepilot = "Case pilot",
    status = "Status",
    emmean_SE_co2 = "Mean ± SE",
    CI95_co2 = "95% CI",
    emmean_SE_ch4 = "Mean ± SE",
    CI95_ch4 = "95% CI",
    emmean_SE_gwp = "Mean ± SE",
    CI95_gwp = "95% CI") %>%
  merge_v(j="casepilot") %>% 
  merge_v(j="casepilot", part="header") %>% 
  merge_v(j="status", part="header") %>% 
  bold(part = "header") %>%
  fontsize(size=10, part = "header") %>% 
  fontsize(size=9, part = "body") %>%
  theme_vanilla() %>%
  autofit()


#Create word doc with emmeans table.

doc <- read_docx() %>%
  body_add_par("Table S2. Model-derived Estimated marginal means of GHG fluxes across status") %>%
  body_add_flextable(ft_emmeans) %>% 
  body_end_section_landscape()

print(doc,
      target=paste0(path_supplementary,"Table_S2_StatusEMMs.docx"))



#Table S3: Contrasts status-----
#Supplementary tables with full contrasts for status: 
#case pilot    Contrast    Difference +- SE    95%CI   p.value
simplemodel_contrasts<- read.csv(file = paste0(path_2_modeloutputs,"Posthoctests_RI_chambermodels.csv"))

formatedcontrasts_RI<- simplemodel_contrasts %>% 
  #Select RI status comparisons: 
  filter(casepilot=="RI", comparison=="status") %>% 
  #Change the order of contrast to show flux-change after restoration
  #(Restored-Altered), avoided emissions through conservation
  #(Preserved-Altered), and "restoration debt" (Restored-Preserved).
  #Sign-inversion and swap of CL order
  mutate(difference=-estimate_bt, SE=SE_bt, lower.CL=-upper.CL_bt, upper.CL=-lower.CL_bt) %>% 
  mutate(Fghg=case_when(ghgspecies=="co2"~"CO2",
                        ghgspecies=="ch4"~"CH4",
                        ghgspecies=="GWPco2andch4"~"CO2eq"),
         contrast_reordered=case_when(contrast == "Altered - Restored"~"Restored - Altered",
                                      contrast == "Altered - Preserved"~"Preserved - Altered",
                                      contrast == "Preserved - Restored"~"Restored - Preserved")) %>% 
  mutate(dataset=paste0(casepilot, " - ",Fghg)) %>% 
  dplyr::select(casepilot, ghgspecies,dataset, contrast_reordered, difference, SE, lower.CL, upper.CL,p.value)


complexmodel_contrasts<- read.csv(file=paste0(path_2_modeloutputs,"Posthoctests_rest_chambermodels.csv"))

formatedcontrasts_rest<- complexmodel_contrasts %>% 
  #Select status comparison for all casepilots (except RI, not in this dataset)
  filter(comparison=="status") %>% 
  #Change the order of contrast to show flux-change after restoration
  #(Restored-Altered), avoided emissions through conservation
  #(Preserved-Altered), and restoration debt (Restored-Preserved).
  #Sign-inversion and swap of CL order
  mutate(difference=-estimate_bt, SE=SE_bt, lower.CL=-upper.CL_bt, upper.CL=-lower.CL_bt) %>% 
  mutate(Fghg=case_when(ghgspecies=="co2"~"CO2",
                        ghgspecies=="ch4"~"CH4",
                        ghgspecies=="GWPco2andch4"~"CO2eq"),
         contrast_reordered=case_when(contrast == "Altered - Restored"~"Restored - Altered",
                                      contrast == "Altered - Preserved"~"Preserved - Altered",
                                      contrast == "Preserved - Restored"~"Restored - Preserved")) %>% 
  mutate(dataset=paste0(casepilot, " - ",Fghg)) %>% 
  dplyr::select(casepilot, ghgspecies,dataset, contrast_reordered, difference, SE, lower.CL, upper.CL,p.value)



formatedcontrasts_all<-  rbind(formatedcontrasts_rest,formatedcontrasts_RI) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         ghgspecies=factor(ghgspecies, levels=c("co2","ch4","GWPco2andch4"), ordered = T)) %>% 
  #Change units of CO2 contrast from mol units to mmol m-2 d-1, ch4 is already
  #in mmol m-2d-1, GWP is in gCO2eq m-2d-1
  mutate(difference= if_else(ghgspecies=="co2", difference*1000, difference), 
         SE=if_else(ghgspecies=="co2", SE*1000, SE),
         lower.CL=if_else(ghgspecies=="co2", lower.CL*1000,lower.CL),
         upper.CL=if_else(ghgspecies=="co2", upper.CL*1000, upper.CL), 
         flux_units=case_when(ghgspecies=="co2"~"mmol m-2 d-1",
                              ghgspecies=="ch4"~"mmol m-2 d-1",
                              ghgspecies=="GWPco2andch4"~"g CO2 eq.m-2 d-1"
         )) %>% 
  arrange(casepilot,ghgspecies) %>%  
  mutate(difference_SE=sprintf("%.3g ± %.3g", difference, SE),# 3 significant figures (normal or scientific notation)
         CI95= sprintf("%.3g to %.3g", lower.CL, upper.CL),# 3 significant figures (normal or scientific notation)
         pvalue=if_else(p.value<0.001,"< 0.001",as.character(round(p.value,digits = 3)))) %>% 
  dplyr::select(casepilot, ghgspecies, dataset, contrast_reordered, difference_SE,CI95, pvalue) %>% 
  rename(contrast=contrast_reordered)


#We need a different column with units
#3 tables (1perGHG), with columns dataset, contrast, Difference (units), Pvalue

#within the Difference column we will have 3 actual columns (estimate +- SE, 95% CI, Pvalue)

ft_contrasts_co2<- formatedcontrasts_all %>%
  filter(ghgspecies=="co2") %>%
  dplyr::select(casepilot, contrast, difference_SE, CI95, pvalue) %>% 
  flextable() %>% 
  add_header_row(
    values = c("Case pilot","Contrast","Difference (mmol CO2 m-2 d-1)"),  
    colwidths = c(1,1,3)                          # number of columns each value spans
  ) %>%
  set_header_labels(
    casepilot = "Case pilot",
    contrast = "Contrast",
    difference_SE = "Estimate ± SE",
    CI95 = "95% CI",
    pvalue = "P-value"
  ) %>% 
  merge_v(j="casepilot") %>% 
  merge_v(j="casepilot", part="header") %>% 
  merge_v(j="contrast", part="header") %>% 
  bold(part = "header") %>%
  fontsize(size=10, part = "header") %>% 
  fontsize(size=10, part = "body") %>%
  theme_vanilla() %>%
  autofit()



ft_contrasts_ch4<- formatedcontrasts_all %>%
  filter(ghgspecies=="ch4") %>%
  dplyr::select(casepilot, contrast, difference_SE, CI95, pvalue) %>% 
  flextable() %>% 
  add_header_row(
    values = c("Case pilot","Contrast","Difference (mmol CH4 m-2 d-1)"),  
    colwidths = c(1,1,3)                          # number of columns each value spans
  ) %>%
  set_header_labels(
    casepilot = "Case pilot",
    contrast = "Contrast",
    difference_SE = "Estimate ± SE",
    CI95 = "95% CI",
    pvalue = "P-value"
  ) %>% 
  merge_v(j="casepilot") %>% 
  merge_v(j="casepilot", part="header") %>% 
  merge_v(j="contrast", part="header") %>% 
  bold(part = "header") %>%
  fontsize(size=10, part = "header") %>% 
  fontsize(size=10, part = "body") %>%
  theme_vanilla() %>%
  autofit()


ft_contrasts_gwp<- formatedcontrasts_all %>%
  filter(ghgspecies=="GWPco2andch4") %>%
  dplyr::select(casepilot, contrast, difference_SE, CI95, pvalue) %>% 
  flextable() %>% 
  add_header_row(
    values = c("Case pilot","Contrast","Difference (g CO2eq. m-2 d-1)"),  
    colwidths = c(1,1,3)                          # number of columns each value spans
  ) %>%
  set_header_labels(
    casepilot = "Case pilot",
    contrast = "Contrast",
    difference_SE = "Estimate ± SE",
    CI95 = "95% CI",
    pvalue = "P-value"
  ) %>% 
  merge_v(j="casepilot") %>% 
  merge_v(j="casepilot", part="header") %>% 
  merge_v(j="contrast", part="header") %>% 
  bold(part = "header") %>%
  fontsize(size=10, part = "header") %>% 
  fontsize(size=10, part = "body") %>%
  theme_vanilla() %>%
  autofit()


doc <- read_docx() %>%
  body_add_par("Table S3. Post-hoc contrasts between Conservation status classes for CO2, CH4 and CO2eq.", style = "heading 1") %>%
  body_add_par("Table S3a. CO2 Post-hoc contrasts.", style = "Normal") %>%
  body_add_flextable(ft_contrasts_co2) %>% 
  body_add_par("Table S3b. CH4 Post-hoc contrasts.", style = "Normal") %>%
  body_add_flextable(ft_contrasts_ch4) %>% 
  body_add_par("Table S3c. CO2eq Post-hoc contrasts.", style = "Normal") %>%
  body_add_flextable(ft_contrasts_gwp) 


print(doc,
      target=paste0(path_supplementary,"Table_S3_posthocContrasts_co2_ch4_co2eq.docx"))



#Table S4: Ancillary Water summary-------

#Use ancillary water data to extract few relevant variables useful for interpretation: 
#Chla: microg/L
#TN and TP units: micromol/L
#Conductivity: mS/cm

#Round to 2 decimals in all cases. 
#Format to avg-+SE (n) 

w<- readxl::read_xlsx(path = paste0(path_0_sourcedata, "Supp_Ancillary_Water_Data.xlsx"))

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
  #Remove high tide samples: 
  filter(is.na(Tide)|Tide=="LT") %>% 
  #subsitute all below detection limit with cero
  mutate(tdp=as.numeric(gsub("<LD", "0", tdp))) %>% 
  #add status: 
  mutate(status=factor(case_when(grepl("A[0-9]",subsite)~"Altered",
                                 grepl("P[0-9]",subsite)~"Preserved",
                                 grepl("R[0-9]",subsite)~"Restored"), levels=c("Altered","Preserved","Restored"), ordered=T)) %>% 
  #Calculate water TN and TP
  mutate(tn=tdn+tpn,
         tp=tdp+tpp)%>%
  dplyr::select(casepilot, status, chla, conductivity, tn,tp) 

#Calculate summary statistics for the relevant variables: 
wat_summary <- wat%>% 
  mutate(conductivity=conductivity/1000) %>% #Conductivity in mS/cm
  pivot_longer(cols = -c(casepilot,status), names_to = "variable",values_to = "value") %>% 
  group_by(casepilot, variable, status) %>% 
  summarise(avg=mean(value, na.rm=T),
            SD=sd(value, na.rm=T), 
            nobs=sum(!is.na(value)),
            se=SD/sqrt(nobs)) %>% 
  dplyr::select(casepilot, status, variable, avg, se, nobs)

head(wat_summary)

#Format numbers
formated_watersummary<- wat_summary %>%
  mutate(avg=round(avg, 2),
         se=round(se,2)) %>% 
  mutate(avg_se=sprintf("%.2f ± %.2f", avg, se)) %>%  
  dplyr::select(casepilot, status, nobs,variable, avg_se) %>% 
  group_by(casepilot, status) %>% 
  mutate(nobs=max(nobs)) %>% 
  pivot_wider(names_from = variable, values_from = avg_se) %>% 
  mutate(casepilot=factor(casepilot, levels = c("DU","RI","CA","VA","DA","CU"), ordered = T),
         status=factor(status, levels = c("Preserved","Altered","Restored"), ordered=T)) %>% 
  arrange(casepilot, status)

#Format table 
ft_water<- formated_watersummary %>% 
  flextable() %>% 
  add_header_row(
    values = c("Case pilot","Status","N","Chl-a (µg L-1)","EC (mS cm-1)","Total N (µM)","Total P (µM)"),  
    colwidths = c(1,1,1,1,1,1,1)                          # number of columns each value spans
  ) %>%
  set_header_labels(
    casepilot = "Case pilot",
    status = "Status",
    nobs= "N",
    chla = "Mean ± SE",
    conductivity = "Mean ± SE",
    tn = "Mean ± SE",
    tp = "Mean ± SE") %>%
  merge_v(j="casepilot") %>% 
  merge_v(j="casepilot", part="header") %>% 
  merge_v(j="status", part="header") %>% 
  merge_v(j="nobs", part="header") %>% 
  bold(part = "header") %>%
  fontsize(size=10, part = "header") %>% 
  fontsize(size=9, part = "body") %>%
  theme_vanilla() %>%
  autofit()


#Create word doc with table.
doc <- read_docx() %>%
  body_add_par("Table S4. Summary of water parameters") %>%
  body_add_flextable(ft_water) 

print(doc,
      target=paste0(path_supplementary,"Table_S4_ancillaryWaterSummary.docx"))



#____________________--------
}

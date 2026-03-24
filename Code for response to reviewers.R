

library(tidyverse)
library(scales)





old<- read.csv(file="C:/Users/Miguel/Documents/GitHub/CabreraBrufau_et_al_2026_code/1_paperdata/ChamberData4paper.csv")

new<- read.csv(file="C:/Users/Miguel/Dropbox/UPDATED_Cabrera-Brufau_et_al_2026_source_data_and_code/1_paperdata/ChamberData4paper.csv")



comp<- new %>% 
  rename(new=dailyflux) %>% 
  left_join(old, by=c("plotcode", "season", "casepilot", "status","subsite", "sampling", "strata", "ghgspecies", "unitflux")) %>% 
  rename(old=dailyflux) %>% 
  mutate(absdif=new-old,
         percentdif=100*(new-old)/abs(old))




#Total % of incubations with dailyflux estimate that retain a value that is within 5% of old estimate: 
gwp_within5percent <- comp %>% 
  filter(ghgspecies=="gwp_co2andch4",
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
       " incubations with a valid daily CO2-eq flux estimate included in the paper, ",
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



##GWP new vs old biplot-----
comp %>% 
  filter(ghgspecies=="gwp_co2andch4",
         !is.na(new)) %>% 
  mutate(incubtype=case_when(abs(percentdif)<=5~"New within 5% of old",
                             TRUE~"New not within 5% of old"
  )) %>% 
  arrange(desc(incubtype)) %>% 
  ggplot(aes(x=abs(old),y=abs(new), col=incubtype))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1)+
  labs(col="Incubation type")+
  scale_x_log10(name="Old estimate (abs. value)",limits=c(1e-3,2e3),
                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000),
                labels=label_log(base = 10))+
  scale_y_log10(name="New estimate (abs. value)", limits=c(1e-3,2e3),
                breaks=c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000),
                labels=label_log(base = 10))+
  theme_classic()+
  facet_wrap(facets=vars(casepilot), scale="free")+
  ggtitle("New vs old absolute CO2eq dailyflux values")





#Boxplots of percent differnece, containing only cases with relevant differences: flux>MDF and diff_abs>flux.se
comp %>% 
  filter(ghgspecies=="co2",
         !is.na(new)) %>% 
  ggplot(aes(x=casepilot, y=absdif, col=status))+
  geom_boxplot(outliers = T)






comp %>% 
  filter(casepilot=="DA") %>% 
  ggplot(aes(x=status, y=absdif))+
  geom_boxplot()


comp %>% 
  filter(casepilot=="DA") %>% 
  pivot_longer(cols=c(new,old), names_to = "newold",values_to = "dailyflux") %>% 
  ggplot(aes(x=newold,y=dailyflux, col=status))+
  geom_boxplot()+
  geom_violin()+
  facet_wrap(facets=vars(ghgspecies),scales="free")




#Ebullition contribution------
library(tidyverse)
library(ggExtra)
library(ggpmisc)
library(ggpubr)
#Reviewers ask to include discussion on the contribution of ebullitive fluxes to total fluxes, especially for DA and CU. 


#Two approaches: 
#1. purely frequentist: proportion of incubations that show ebullitive behaviour. Plotcodes where at least one of the incubations (in case of multiple ones) show ebullition vs total plotcodes (per casepilot and per casepilot-status).


#2. quantitative contribution of flux from plotcodes with ebullitive incubations to total CH4 flux. Not exactly contribution of ebullition to total flux (as incubations with ebullition can also present diffusion), but potentially more informative. 


#In any case, we require: 1. presence/absence of ebullition per plotcode and daily ch4 flux per plotcode. 
#ebullitive inspection: we can obtain it directly form the best.model.flags column of ch4_bestflux.csv

ch4best<- read.csv("C:/Users/Miguel/Dropbox/UPDATED_Cabrera-Brufau_et_al_2026_source_data_and_code/0_sourcedata/ch4_bestflux.csv")

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
paperdat<- read.csv("C:/Users/Miguel/Dropbox/UPDATED_Cabrera-Brufau_et_al_2026_source_data_and_code/1_paperdata/ChamberData4paper.csv")



#obtain only daily ch4 fluxes in paper and add logical ebullitive column
ch4_dat<- paperdat %>% 
  filter(ghgspecies=="ch4") %>% 
  filter(!is.na(dailyflux)) %>% 
  left_join(incu_ebu, by="plotcode") %>% 
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
          subtitle="Flux is standardised by casepilot and status")


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

ch4_dat %>% 
  group_by(casepilot, status, season,subsite) %>% 
  mutate(standard_flux=(dailyflux-min(dailyflux))/(max(dailyflux) - min(dailyflux))) %>% 
  summarise(percent_cv_flux=sd(standard_flux)/mean(standard_flux)*100,
            flux_skew=skewness(standard_flux),
            flux_sd=sd(standard_flux),
            percent_ebullitive=sum(ebullitive)/sum(!is.na(ebullitive)*100),
            any_ebullitive=sum(ebullitive)>0) %>% 
  ggplot(aes(x=percent_ebullitive, y=percent_cv_flux, col=any_ebullitive))+
  geom_point()+
  labs(y="CV of standardised dailyflux (%)", x="Proportion of chambers with ebullition")+
  geom_smooth(method="lm")+
  stat_cor()+
  ggtitle("Flux variability vs ebullitive chamber prevalence",
          subtitle = "Each point is a distinct sampling event (n=144)")




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




#PLot flux distribution per casepilot, status and ebullition 
ggplot(ch4_dat, aes(x=ebullitive, y=dailyflux+0.20, col=status))+
  geom_boxplot()+
  facet_wrap(facets=vars(casepilot),scales="free")+
  scale_y_log10()

ch4_dat %>% 
  group_by(casepilot,status, ebullitive) %>% 
  summarise(flux_sum=sum(dailyflux)) %>% 
  ungroup() %>% 
  group_by(casepilot,status) %>% 
  mutate(prop_flux_sum=flux_sum/sum(flux_sum)) %>% 
ggplot(aes(x=status, y=flux_sum, fill=ebullitive))+
  geom_bar(stat = "identity")+
  facet_wrap(facets=vars(casepilot),scales="free")+
  ggtitle("Sum of fluxes per chamber type")


ch4_dat %>% 
  group_by(casepilot,status, ebullitive) %>% 
  summarise(flux_sum=sum(dailyflux)) %>% 
  ungroup() %>% 
  group_by(casepilot,status) %>% 
  mutate(prop_flux_sum=flux_sum/sum(flux_sum)) %>% 
  ggplot(aes(x=status, y=prop_flux_sum*100, fill=ebullitive))+
  geom_bar(stat = "identity")+
  labs(y="% of total flux)")+
  facet_wrap(facets=vars(casepilot),scales="free")+
  ggtitle("Contribution of chambers with ebullition to total flux")



# Load packages
library(ggplot2)
library(dplyr)
library(emmeans)
library(multcompView)
library(multcomp)

# 1. Transform response to match log scale
ch4_dat <- ch4_dat %>%
  mutate(log_flux = log10(dailyflux + 1))

# 2. Fit model (including interaction for faceting variable)
model <- lm(log_flux ~ status * casepilot, data = ch4_dat)

# 3. Estimated marginal means
emm <- emmeans(model, ~ status | casepilot)

# 4. Compact letter display (CLD)
cld_df <- cld(emm, Letters = letters)

# Clean up letters (remove spaces)
cld_df <- cld_df %>%
  mutate(.group = gsub(" ", "", .group)) %>% 
  mutate(casepilot=factor(casepilot, levels=c("DU","RI","CA","VA","DA","CU"), ordered = T))

# 5. Get y-positions for letters (per facet)
y_pos <- ch4_dat %>%
  group_by(casepilot, status) %>%
  summarise(y = max(dailyflux, na.rm = TRUE), .groups = "drop")

# Merge positions with CLD
cld_df <- left_join(cld_df, y_pos, by = c("casepilot", "status"))

# Add offset so letters sit above plots
cld_df <- cld_df %>%
  mutate(y = y * 1.5)

# 6. Plot
ggplot(ch4_dat, aes(x = status, y = dailyflux + 1, col = status)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.2) +
  facet_wrap(vars(casepilot), scales = "free") +
  scale_y_log10() +
  geom_text(
    data = cld_df,
    aes(x = status, y = y, label = .group),
    inherit.aes = FALSE,
    color = "black"
  ) +
  labs(y = "Daily flux (log scale)", x = "Status")





#Alternative Weighting EMMEANs----
#Example of "custom weighting scheme" using emmeans to layer the weighting: emmean is proportional to vegpresence within each season, but each season has the same weight in final emmean. 


library(lme4)
library(emmeans)
library(dplyr)

set.seed(123)

# Base grid
dat <- expand.grid(
  status = c("A","B"),
  season = c("spring","summer","fall","winter"),
  site = 1:20
)

# Assign unequal vegpresence probabilities per status × season
veg_probs <- expand.grid(
  status = c("A","B"),
  season = c("spring","summer","fall","winter")
) %>%
  mutate(p_T = runif(n(), 0.3, 0.8))

dat <- dat %>%
  left_join(veg_probs, by = c("status","season")) %>%
  mutate(
    vegpresence = ifelse(runif(n()) < p_T, "T", "F")
  ) %>%
  dplyr::select(-p_T)

# Simulate response
dat$flux <- rnorm(nrow(dat),
                  mean = ifelse(dat$status=="A", 5, 7) +
                    ifelse(dat$vegpresence=="T", 2, -1) +
                    as.numeric(factor(dat$season)),
                  sd = 1)

# Fit mixed model
model <- lmer(flux ~ status * season * vegpresence + (1|site), data = dat)

emm_season <- emmeans(
  model,
  ~ status | season,
  weights = "proportional"
)
emm_season
annual_emm <- emmeans(
  emm_season,
  ~ status,
  weights = "equal"
)

annual_emm

equal_emm<- emmeans(
  model,
  ~ status,
  weights = "equal"
)

prop_emm<-  emmeans(
  model,
  ~ status,
  weights = "proportional"
)

annual_emm
equal_emm
prop_emm

dat %>%
  group_by(status, season, vegpresence) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(status, season) %>%
  mutate(prop = n/sum(n))

pairs(annual_emm)
pairs(equal_emm)
pairs(prop_emm)



emmeans(model, ~ status | season * vegpresence, weights = "show.levels")


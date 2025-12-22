#GHGpaper_modelchamberdata.R

#Author: Miguel Cabrera-Brufau
#Date: September 2025
#Project: Restore4cs


#Description----
#This scrip is used to model the effect of restoration in each casepilot. Using net daily GHG exchange rates of appropriate incubations, data preparation is in GHGpaper_prepChamberData.R script. 


#DECISSIONS: 
  #Pooling of non-vegetated strata: due to seasonal variability and site-specific differences, using the 3 original strata classes is impossible (would lead to to rank-deficient models for some combinations). Instead we group them into vegpresence: vegetated or not-vegetated
  
#Modelling approach for DU,CA,VA,DA,CU: GLMMtmb with general formula dailyflux~status*season*vegpresence + (1|subsite)
  #Modelling approach for RI: reduced model without vegpresence (to avoid rank-deficient due to restored sites being 100% vegetated)

  #Contrast: Set contrast options to "contr. sum", so that we are not using any one level of a factor as the reference level, but rather the tests will asses whether there is an overall average effect of one factor across all levels of the other factors. For example, is there an effect of status across all levels of season?


#STEPS: 
#1. Import and format data: remove Nas, as.factor levels, 
#2. Transform data: Using BestNormalize package, find and apply for each casepilot-GHGspecies, the transformation that maximizes normality. 
#3. Optimize and select modelling structure (family-distribution) for each casepilot*GHGspecies combo (co2, ch4, GWPco2andch4).
#4. Calculate model outputs: 
  
  #Model summaries
  #Model residual-tests
  #Estimated marginal means (EMMs) for desired levels 
  #Post-hoc contrasts between EMMs 


rm(list = ls()) # clear workspace
cat("/014") # clear console


# ---- Packages ----
#For data-handling & plotting
library(tidyverse)
library(lubridate)
library(zoo)
library(ggpubr)
library(rlang)
library(ggtext)#for plot captions
require(purrr)
require(data.table)
require(tools)
library(hms)
library(suncalc)


#For modelling: 

library(bestNormalize) #For data-transformations
library(numDeriv) #For back-transformation of emmeans and contrasts SEs
library(scales)#For pseudolog transformation
library(glmmTMB) #For glmm modelling
library(DHARMa) #For residual evaluation
library(performance)#Checking singluar models
library(rstatix) #Outlier exploration & Anova of models
library(emmeans) #To estimate EMMs
library(multcomp)# for emmeans post-hoc comparisons
library(multcompView) #for clc (letters groups)



# ---- Directories ----

#Get root directory of repository:
path_root <- dirname(rstudioapi::getSourceEditorContext()$path)

#Path with formated data for modelling: 
path_1_paperdata <- paste0(path_root,"/1_paperdata/")

#Path to save model outputs:
path_2_modeloutputs <- paste0(path_root,"/2_modeloutputs/")

#Subfolder to save extra plots:
extraplots_path <- paste0(path_2_modeloutputs,"/Exploratory_figures/")

#Create folders if they do not exists: 
if (!dir.exists(path_2_modeloutputs)) {
  dir.create(path_2_modeloutputs, recursive = TRUE)
}

if (!dir.exists(extraplots_path)) {
  dir.create(extraplots_path, recursive = TRUE)
}



#0.Contrast options-----
#NOTES on contrasts: in R by default, contrast is "contr.treatment" which uses the first level of each factor as the "reference" level and all subsequent as "treatment" levels. This implies that, with default contrast type, when looking at the main effects of my model, what is shown for each main effect is whether there is a significant effect of 1 factor (eg. status) at the reference level of all other factors (i.e. season). What we want is to asses whether there is an overall average effect of status across all levels of season. For this purpose we need to set contrasts to "contr. sum", and always call Anova (model, type="III").

#Set contrasts to contr.sum (for the whole R session)
options(contrasts = c("contr.sum", "contr.poly"))


#0. Import and format--------

#Import data and format: 
data4models<- read.csv(paste0(path_1_paperdata,"ChamberData4paper.csv"))

#Format and filter:
data4models<- data4models %>% 
  filter(ghgspecies%in%c("co2","ch4","gwp_co2andch4")) %>% 
  #Rename to GWPco2andch4
  mutate(ghgspecies=if_else(ghgspecies=="gwp_co2andch4","GWPco2andch4",ghgspecies)) %>% 
  #Remove NAs
  filter(!is.na(dailyflux)) %>% 
  mutate(season=factor(season, ordered = F),
         status=factor(status, ordered = F),
         strata=factor(strata, ordered = F)) %>% 
  #Pool non-vegetated strata:
  mutate(vegpresence=if_else(strata=="vegetated","Vegetated","Non-vegetated")) %>% 
  mutate(vegpresence=factor(vegpresence, ordered=F))
  

#Notes -----
#Evaluate model structure and assumptions (residuals), not significance of effects nor even explained variability (we expect very little variance explained by status in many cases). 

#1st step is to decide the approach based on our data structure and distribution. 

#Fixed decissions: 
#Each case-pilot is modelled independently (different data-distributions and transformation needs)
#Need to include subsite as random effect (accounts for repeated samplings and site-specific intercepts)

#Approach: use Generalized Linear Mixed Models (glmmTMB)
#They can account for non-normal data (even after transformation) by using different distribution families


#STEPS in model optimization: 
#1. Decide transformation: that which maximizes normality of data
#2. Decide modelling family and structure (Gaussian vs T-family): Gaussian preferred unless residuals fail


#0. Custom Functions -------
#Here functions to access relevant results from models contained inside a named list. Used for ease of formatting and to avoid repetition. Will be used to summarise information.  

#Formula to get marginal and conditional R2s from model list:
get_R2s <- function(model_list) {
  library(glmmTMB)
  library(performance)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    
    r2_vals<- NULL
    # Try to extract R2 values
    r2_vals <- tryCatch({
      MuMIn::r.squaredGLMM(model) #Using MuMin to be able to get marginal R2 even with singularities
    }, error = function(e) {
      warning(paste("R2 extraction failed for", cp))
      return(c(NA_real_, NA_real_))
      
    })
    
    tibble(
      dataset = cpghg,
      transformation=transformation,
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      homoced_R2m = r2_vals[1,1],
      homoced_R2c = r2_vals[1,2]
    )
  })
}



#Function to get various pseudoR2 estimates and extract structure of each model
#calculate various pseudo-R2 metrics (log-Likelyhood based) to get the "improvement in model fit relative to a null model (intercept + random effects), capturing the combined explanatory contribution of fixed effects".
#We have 4 estimates: 
#Efron pseudoR2: squared correlation between observed and predicted values. 
#The other 3 pseudoR2 come from Log-likelihood comparisons between actual and null model.

get_pseudoR2s<- function(model_list){
  
  map_df(names(model_list), function(cpghg){
    model<- model_list[[cpghg]]
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    
    
    obs<- model$frame[,1]
    pred<- predict(model, type="response")
    #Efron_pseudoR2 is equal to the squared correlation between predicted and observed values
    Efron_pseudoR2<- 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)
    
    #rcompanion calcualtes 3 different pseudoR2 based on comparision of actual vs null model (mantaining the structure, but removing the predictors from the model formula)
    res<-rcompanion::nagelkerke(model)  
    
    tibble(
      dataset = cpghg,
      transformation = transformation,
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      pseudoR2_Efron = Efron_pseudoR2,
      pseudoR2_McFadden = res$Pseudo.R.squared.for.model.vs.null[1],
      pseudoR2_Cox_and_snell=res$Pseudo.R.squared.for.model.vs.null[2],
      pseudoR2_Nagelkerke = res$Pseudo.R.squared.for.model.vs.null[3]
    )
  })
}


#Function to extract model structure and all fixed effects significance (regardless of their naming, can differ in order or presence across the models in model_list)
get_all_anova_results <- function(model_list) {
  library(car)
  library(glmmTMB)
  library(dplyr)
  library(purrr)
  
  map(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    # transformation <- table_trans %>% filter(dataset == cpghg) %>% pull(trans_name)
    
    # Run Type III ANOVA
    anova_res <- tryCatch(car::Anova(model, type = "III"), error = function(e) NULL)
    
    # Extract p-values and rename terms
    pval_tibble <- if (!is.null(anova_res)) {
      terms <- rownames(anova_res)
      pvals <- anova_res[, "Pr(>Chisq)"]
      names(pvals) <- paste0(gsub(":", ".", terms), "_pval")
      as_tibble(as.list(pvals))
    } else {
      tibble()
    }
    
    # Combine with metadata
    tibble(
      dataset = cpghg,
      # transformation = transformation,
      formula = deparse(formula(model)),
      family = model$modelInfo$family$family,
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      nobs = dim(model$frame)[1]
    ) %>%
      bind_cols(pval_tibble)
  }) %>%
    bind_rows()
}


#Summary of Model assumptions: 
#We want to structure the results in a table: model structure,  distribution family, dispformula, 
#add pvalue for: normality of residuals, heteroscedasticity (levene test).

summarize_dharma_diagnostics <- function(model_list) {
  library(DHARMa)
  library(dplyr)
  library(purrr)
  
  map_df(names(model_list), function(cpghg) {
    model <- model_list[[cpghg]]
    data  <- model$frame
    transformation<- table_trans %>% filter(dataset==cpghg) %>% pull(trans_name)
    # Simulate residuals
    sim_res <- tryCatch({
      simulateResiduals(fittedModel = model, plot = FALSE)
    }, error = function(e) {
      warning(paste("DHARMa simulation failed for", cp))
      return(NULL)
    })
    
    if (is.null(sim_res)) {
      return(tibble(
        dataset = cpghg,
        transformation=transformation, 
        uniformity_pval = NA_real_,
        dispersion_pval = NA_real_,
        hetero_status_pval = NA_real_,
        hetero_season_pval = NA_real_
      ))
    }
    
    # Run tests
    tibble(
      dataset = cpghg,
      transformation=transformation, 
      formula = deparse(formula(model)),
      dispformula = deparse(model$modelInfo$allForm$dispformula),
      family = model$modelInfo$family$family,
      convergence=check_convergence(model),
      singular=check_singularity(model,),
      uniformity_pval = testUniformity(sim_res)$p.value,
      dispersion_pval = testDispersion(sim_res)$p.value,
      hetero_status_pval = testCategorical(sim_res,catPred =data$status,  plot = F)$homogeneity$`Pr(>F)`[1],
      hetero_season_pval = testCategorical(sim_res,catPred =data$season,  plot = F)$homogeneity$`Pr(>F)`[1]
    )
  })
}


#Function to Produce comparative histograms and qqplots of untransformed and best-Normalize transformed data. 
plot_bestNormalize <- function(bn_obj, n_bins = 30, title = NULL) {
  library(ggplot2)
  library(gridExtra)
  
  # Extract original and transformed data
  x_orig <- bn_obj$x
  x_trans <- predict(bn_obj)
  
  # Get transformation name from class attribute
  trans_name <- class(bn_obj$chosen_transform)[1]
  
  # Make combined title (if title is provided)
  plot_title <- paste0("Transformation: ", trans_name, if (!is.null(title)) paste0(" - ", title) else "")
  
  # Perform KS tests
  ks_orig <- ks.test(x_orig, "pnorm", mean(x_orig), sd(x_orig))
  ks_trans <- ks.test(x_trans, "pnorm", mean(x_trans), sd(x_trans))
  
  # Create data frame for plotting
  df <- data.frame(
    value = c(x_orig, x_trans),
    type = factor(rep(c("Original", "Transformed"), each = length(x_orig)), levels = c("Original", "Transformed"))
  )
  
  # Histogram
  p_hist <- ggplot(df, aes(x = value, fill = type)) +
    geom_histogram(bins = n_bins, alpha = 0.6, position = "identity") +
    facet_wrap(~type, scales = "free") +
    labs(title = paste("Histogram:", plot_title), x = "Value", y = "Count") +
    theme_minimal() +
    scale_fill_manual(values = c("Original" = "steelblue", "Transformed" = "darkgreen")) +
    theme(legend.position = "none")
  
  # QQ plot data
  qq_data <- data.frame(
    sample = c(x_orig, x_trans),
    type = factor(rep(c("Original", "Transformed"), each = length(x_orig)), levels = c("Original", "Transformed"))
  )
  
  # KS p-values for annotation
  ks_labels <- data.frame(
    type = c("Original", "Transformed"),
    label = c(
      paste0("KS p-value: ", signif(ks_orig$p.value, 4)),
      paste0("KS p-value: ", signif(ks_trans$p.value, 4))
    ),
    x = c(Inf, Inf),
    y = c(-Inf, -Inf)
  )
  
  # QQ plot with annotation
  p_qq <- ggplot(qq_data, aes(sample = sample)) +
    stat_qq(color = "black") +
    stat_qq_line(color = "red", linetype = "dashed") +
    facet_wrap(~type, scales = "free") +
    geom_text(data = ks_labels, aes(x = x, y = y, label = label), hjust = 1.1, vjust = -0.5, inherit.aes = FALSE) +
    labs(title = paste("QQ Plot:", plot_title), x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  
  # Combine plots
  gridExtra::grid.arrange(p_hist, p_qq, nrow = 2)
}



#Function that plots observed vs predicted values for a model (or pair of models fitted under the same data), coloring points via color_var
plot_obs_vs_pred_models <- function(model1, data, model2 = NULL, color_var = "vegpresence") {
  # Get name of the data object
  data_name <- deparse(substitute(data))
  
  # Extract and clean fixed effect formula (RHS only)
  model1_name <- gsub("^dailyflux_trans\\s*~\\s*", "", deparse(formula(model1)))
  
  df1 <- data.frame(
    obs = model1$frame$dailyflux_trans,
    pred = predict(model1, type = "response"),
    color = data[[color_var]],
    model = model1_name
  )
  
  if (!is.null(model2)) {
    model2_name <- gsub("^dailyflux_trans\\s*~\\s*", "", deparse(formula(model2)))
    
    df2 <- data.frame(
      obs = model2$frame$dailyflux_trans,
      pred = predict(model2, type = "response"),
      color = data[[color_var]],
      model = model2_name
    )
    
    df_combined <- bind_rows(df1, df2) %>%
      mutate(model = factor(model, levels = c(model1_name, model2_name)))
  } else {
    df_combined <- df1 %>%
      mutate(model = factor(model, levels = model1_name))
  }
  
  # Compute R² and label positions
  r2_labels <- df_combined %>%
    group_by(model) %>%
    summarise(
      r2 = summary(lm(obs ~ pred))$r.squared,
      label = paste0("R² = ", round(r2, 2)),
      .groups = "drop"
    )
  
  # Get range of observed values
  obs_range <- range(df_combined$obs, na.rm = TRUE)
  
  # Plot
  ggplot(df_combined, aes(x = pred, y = obs)) +
    geom_point(aes(col = color), size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_smooth(method = "lm", col = "black", se = FALSE) +
    geom_smooth(method = "lm", aes(col = color), se = FALSE) +
    facet_wrap(~model) +
    geom_text(
      data = r2_labels,
      aes(x = min(obs_range), y = max(obs_range), label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 4
    ) +
    scale_x_continuous(limits = obs_range) +
    scale_y_continuous(limits = obs_range) +
    labs(
      title = paste(data_name, "Observed vs Predicted"),
      subtitle = paste("Best-fit lines and R² from linear regression in model scale\nColored by", color_var),
      x = "Predicted",
      y = "Observed",
      color = color_var
    ) +
    
    theme_bw() +
    theme(
      plot.caption = element_textbox_simple(),
      legend.position = "bottom",
      legend.box.margin = margin(t = -5, b = -5),  # reduce top/bottom margin
      legend.spacing = unit(0.2, "lines")          # reduce spacing between legend items
    )
}

# #Example usage:
# plot_obs_vs_pred_models(model1=ca_simplemodel_co2,
#                         data=ca_co2, 
#                         model2=ca_best_co2,
#                         color_var = "status")


##Pseudolog-------
#Adapted from scales package

#Notes on pseudolog (~signed log) transformation: https://stratosida.github.io/regression-regrets/Pseudo_log_explainer.html#finding-a-parameter-that-best-achieves-normality
#Pseudolog tranformation needs tunning of its parameter to each data-set (similar to yeo-johnson). To incorporate this transformation into the BestNormalize, we need to create several functions. 

#Pseudolog transformation allows and maintains NAs, allows and mantains negative and positive values 
#Range of potential sigmas are data-driven: based on 10th-90th percentile ranges. 

# 1. Constructor: fits pseudo-log transformation with optimal sigma
#best_pseudo_log with data-driven optimization of sigma using 10th and 90th percentile range
best_pseudo_log <- function(x, base = 10, sigma_range = NULL, standardize = TRUE) {
  tryCatch({
    # Clean input
    x_clean <- x[!is.na(x) & is.finite(x)]
    if (length(x_clean) < 3) stop("Not enough valid data points.")
    
    # Compute data-driven sigma_range if not provided (using 10th-90th percentile range)
    if (is.null(sigma_range)) {
      q <- quantile(abs(x_clean), probs = c(0.1, 0.9), na.rm = TRUE)
      lower <- max(q[1] / 10, .Machine$double.eps)  # Avoid zero or near-zero
      upper <- max(q[2], lower * 10)
      sigma_range <- 10^seq(log10(lower), log10(upper), length.out = 100)#100 potential sigmas
    }
    
    # Evaluate normality across candidate sigma values
    stats <- sapply(sigma_range, function(sigma) {
      x.t <- tryCatch(
        pseudo_log_trans(base = base, sigma = sigma)$transform(x_clean),
        error = function(e) rep(NA, length(x_clean))
      )
      if (length(unique(x.t[!is.na(x.t)])) < 3) return(NA)
      tryCatch(stats::shapiro.test(x.t[!is.na(x.t)])$statistic, error = function(e) NA)
    })
    
    # Select best sigma
    if (all(is.na(stats))) stop("No valid transformation found.")
    best_sigma <- sigma_range[which.max(stats)]
    
    # Apply transformation
    trans <- pseudo_log_trans(base = base, sigma = best_sigma)
    x.t <- trans$transform(x)
    
    # Optionally standardize
    mu <- mean(x.t, na.rm = TRUE)
    sd <- sd(x.t, na.rm = TRUE)
    if (standardize) x.t <- (x.t - mu) / sd
    
    # Pearson normality statistic
    ptest <- nortest::pearson.test(x.t[!is.na(x.t)])
    norm_stat <- unname(ptest$statistic / ptest$df)
    
    # Return result
    val <- list(
      x = x,
      x.t = x.t,
      base = base,
      sigma = best_sigma,
      trans = trans,
      mean = mu,
      sd = sd,
      standardize = standardize,
      n = length(x_clean),
      norm_stat = norm_stat
    )
    class(val) <- c("best_pseudo_log", class(val))
    return(val)
  }, error = function(e) {
    structure(list(norm_stat = NA), class = "best_pseudo_log")
  })
}

# 2. Predict method: applies transformation or inverse using fitted parameters
#Predict method with NA handling:
predict.best_pseudo_log <- function(object, newdata = NULL, inverse = FALSE, ...) {
  if (is.null(newdata) && !inverse) newdata <- object$x
  if (is.null(newdata) && inverse) newdata <- object$x.t
  
  na_idx <- is.na(newdata)
  
  if (inverse) {
    if (object$standardize) {
      newdata[!na_idx] <- newdata[!na_idx] * object$sd + object$mean
    }
    newdata[!na_idx] <- object$trans$inverse(newdata[!na_idx])
  } else {
    newdata[!na_idx] <- object$trans$transform(newdata[!na_idx])
    if (object$standardize) {
      newdata[!na_idx] <- (newdata[!na_idx] - object$mean) / object$sd
    }
  }
  
  return(unname(newdata))
}



# 3. Print method (optional): displays transformation details
print.best_pseudo_log <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      "pseudo_log(x) Transformation with", x$n, "nonmissing obs.:\n",
      "Relevant statistics:\n",
      "- base =", x$base, "\n",
      "- sigma =", x$sigma, "\n",
      "- mean (before standardization) =", x$mean, "\n",
      "- sd (before standardization) =", x$sd, "\n")
}

# 4. Register custom transformation structure to use within bestNormalize
pseudolog_transform <- list(
  best_pseudo_log = best_pseudo_log,
  predict.best_pseudo_log = predict.best_pseudo_log,
  print.best_pseudo_log = print.best_pseudo_log
)


#__________------
#1.Transformations------

#Data needs transformation for compliance with model assumptions. 
#We will use a completely objective decision for transformation: using BestNormalize function (with pseudolog tranformation as a potential option), we will chose the transformation that results in the most-normal data distribution. 

#USE bestNormalize (allowing for pseudo-log trans) to obtain for each dataset (ghg*casepilot combo) the most-Normal transformation: 
rm(table_trans, table_trans_i)
for (cp in unique(data4models$casepilot)) {
  for (ghg in  unique(data4models$ghgspecies)) {
    
    x<- data4models %>% filter(casepilot==cp&ghgspecies==ghg) %>%
      pull(dailyflux)
    bn_result<- bestNormalize(x = x, new_transforms = pseudolog_transform, standardize = TRUE, warn = TRUE, k = 5,allow_orderNorm = F)
    #calculate normality test (shapiro)
    x_trans <- bn_result$chosen_transform$x.t
    pval <- shapiro.test(x_trans)$p.value
    
    table_trans_i <- tibble(
      dataset = paste(cp, ghg,sep = "_"),
      casepilot = cp,
      ghgspecies = ghg,
      nobs= length(x),
      trans_name = attr(x = bn_result$chosen_transform,which = "class")[1],
      trans_Pearson_p=bn_result$chosen_transform$norm_stat, #Pearson's P (the lower, the more-normal data, only used for ranking our transformation options)
      shapiro_pval=pval, #P value of Normality
      bn_result = list(bn_result)  # wrap BestNormalize object in list to store as list-column
    )
    
    
    #Create table_trans for first round of loop, then apend to table
    if(cp==unique(data4models$casepilot)[1]&
       ghg==unique(data4models$ghgspecies)[1]){
      table_trans<- table_trans_i
    }else{ table_trans<- bind_rows(table_trans, table_trans_i)}
  }
}

#View all bestNormalizing transformations:
table_trans

rm(cp, bn_result, table_trans_i, ghg, x, x_trans, pval)


#Automatically apply the best transformation to each dataset within data4models using the bestNormalize objects stored in table_trans 

# Create a named list of BestNormalize objects for easy lookup
bn_list <- table_trans %>%
  dplyr::select(dataset, bn_result) %>%
  deframe()

# Add a key column to data4models for matching
data4models <- data4models %>%
  mutate(dataset = paste(casepilot, ghgspecies, sep = "_"))

# Apply the corresponding BestNormalize transformation
data4models <- data4models %>%
  rowwise() %>%
  mutate(
    trans_name = table_trans$trans_name[match(dataset, paste(table_trans$casepilot, table_trans$ghgspecies, sep = "_"))],
    dailyflux_trans = if (!is.null(bn_list[[dataset]])) predict(bn_list[[dataset]], dailyflux) else NA_real_
  ) %>%
  ungroup()


#2. Model fitting-------

#GENERAL MODEL potential OPTIONS: 

#For RI, we cannot include vegpresence in model (would lead to deficient models, restored is 100%vegetated)USE ONLY BASIC FORMULA: dailyflux_trans~season*status + (1|subsite)

#For REST, allow vegpresence as fixed effect, decide best model based on residuals and fit. 

#Always use transformed data (most parsimonious option).
#IF residuals are not normal, try using t_family instead of gaussian. 


#CO2______-------
##2.1. CA_co2 (ok) -------

#Subset data and check transformation: 
ca_co2<- data4models %>% filter(casepilot=="CA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = gaussian(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Good residuals, keep


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                   data = ca_co2,
                   family = gaussian(),
                   dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: does not converge. 


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                      data = ca_co2,
                      family = t_family(),
                      dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#GOOD residuals, keep


#Compare valid models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
ca_simplemodel_co2<- m1_gaus_nostrata
ca_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_co2, 
                        model2 = ca_complexmodel_co2,
                        data=ca_co2,
                        color_var = "vegpresence")

ggsave(filename = "CA_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_co2 (ok) -------

#Subset data and check transformation: 
cu_co2<- data4models %>% filter(casepilot=="CU"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = cu_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#GOOD residuals, keep


#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_co2<- m3_t_nostrata
cu_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_co2, 
                        model2 = cu_complexmodel_co2,
                        data=cu_co2,
                        color_var = "vegpresence")

ggsave(filename = "CU_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.3. DA_co2 (ok) -------

#Subset data and check transformation: 
da_co2<- data4models %>% filter(casepilot=="DA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#Best residuals

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
da_simplemodel_co2<- m3_t_nostrata
da_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_co2, 
                        model2 = da_complexmodel_co2,
                        data=da_co2,
                        color_var = "vegpresence")

ggsave(filename = "DA_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_co2 (ok) -------

#Subset data and check transformation: 
du_co2<- data4models %>% filter(casepilot=="DU"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: GOOD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: Model does not converge


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = du_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#DECISSION: good residuals

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
du_simplemodel_co2<- m1_gaus_nostrata
du_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_co2, 
                        model2 = du_complexmodel_co2,
                        data=du_co2,
                        color_var = "vegpresence")

ggsave(filename = "DU_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_co2 (ok) -------

#Cannot include vegpresence (inherently related to status, leads to deficient model)

#Subset data and check transformation: 
ri_co2<- data4models %>% filter(casepilot=="RI"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: GOOD residuals


#SKIP t-family option: gaussian is already good
if (F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)

#DECISSION: 
}


#Compare models (SKIP!)
# anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

#Save best models (simple and complex best)
ri_simplemodel_co2<- m1_gaus_nostrata
# ri_complexmodel_co2<- NA  NO complex model possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_co2, 
                        # model2 = ri_complexmodel_co2,
                        data=ri_co2,
                        color_var = "vegpresence")

ggsave(filename = "RI_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_co2 (ok) -------

#Subset data and check transformation: 
va_co2<- data4models %>% filter(casepilot=="VA"&ghgspecies=="co2")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="co2") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_co2,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: BAD residuals (bad model)


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_co2,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: Model does not converge


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = va_co2,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#DECISSION: good residuals

#Compare models 

anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_co2<- m1_gaus_nostrata
va_complexmodel_co2<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_co2, 
                        model2 = va_complexmodel_co2,
                        data=va_co2,
                        color_var = "vegpresence")

ggsave(filename = "VA_co2_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#CH4______-------

##2.1. CA_ch4 (ok) -------

#Subset data and check transformation: 
ca_ch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ca_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ca_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: GOOD residuals, keep


#SKIP t-family options, gaussians already good
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: 


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = ca_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
check_overdispersion(m4_t_vegpresence)
#DECISSION: 
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does not significantly improve fit. We will use it regardless for consistency and to check for vegetation-transport effects

#Save best models (simple and complex best)
ca_simplemodel_ch4<- m1_gaus_nostrata
ca_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_ch4, 
                        model2 = ca_complexmodel_ch4,
                        data=ca_ch4,
                        color_var = "vegpresence")

ggsave(filename = "CA_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_ch4 (ok) -------

#Subset data and check transformation: 
cu_ch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Good-enough residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: Good residuals


#SKIP, GAUSSIAN ALREADY GOOD
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = cu_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-9)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#DECISSION: 
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_ch4<- m1_gaus_nostrata
cu_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_ch4, 
                        model2 = cu_complexmodel_ch4,
                        data=cu_ch4,
                        color_var = "vegpresence")

ggsave(filename = "CU_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.3. DA_ch4 (ok) -------

#Subset data and check transformation: 
da_ch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#Better residuals, although they present issues

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
da_simplemodel_ch4<- m3_t_nostrata
da_complexmodel_ch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_ch4, 
                        model2 = da_complexmodel_ch4,
                        data=da_ch4,
                        color_var = "vegpresence")

ggsave(filename = "DA_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_ch4 (ok) -------

#Subset data and check transformation: 
du_ch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: GOOD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence) #OK
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #Singular
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: GOOD residuals


#SKIP, gaussian already good: 
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = du_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#DECISSION:

}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does not significantly improve fit (not worth the extra complexity), we will still use it for consistency and to evaluate vegetation effects.

#Save best models (simple and complex best)
du_simplemodel_ch4<- m1_gaus_nostrata
du_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_ch4, 
                        model2 = du_complexmodel_ch4,
                        data=du_ch4,
                        color_var = "vegpresence")

ggsave(filename = "DU_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_ch4 (ok) -------

#Subset data and check transformation: 
ri_ch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: BAD residuals


#SKIP! cannot include vegpresence as fixed with status
if(F){
  #Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
  m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                                data = ri_ch4,
                                family = gaussian(),
                                dispformula = ~1)
  #Evaluate model:
  check_convergence(m2_gaus_vegpresence)
  check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
  r2(m2_gaus_vegpresence, tolerance = 1e-10)
  Anova(m2_gaus_vegpresence)
  summary(m2_gaus_vegpresence)
  res<- simulateResiduals(m2_gaus_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m2_gaus_vegpresence)
  #DECISSION: 
}



  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = ri_ch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  #DECISSION: GOOD residuals
  
  #SKIP, cannot include vegpresence as fixed effect
  if (F){
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ri_ch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  summary(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
  summary(em)
  pairs(em)
  #DECISSION: 
}


#Compare models 
# anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

  #Save best models (simple and complex best)
ri_simplemodel_ch4<- m3_t_nostrata
# ri_complexmodel_ch4<- NA #NOt possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_ch4, 
                        # model2 = ri_complexmodel_ch4,
                        data=ri_ch4,
                        color_var = "vegpresence")

ggsave(filename = "RI_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_ch4 (ok) -------

#Subset data and check transformation: 
va_ch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="ch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="ch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_ch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: GOOD residuals 


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_ch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: GOOD residuals


#SKIP t-family options, Gaussian already good
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_ch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
em<- emmeans(m3_t_nostrata, ~ status, weights = "equal")
summary(em)
pairs(em)
#DECISSION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = va_ch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence) #OK
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-10)
Anova(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
em<- emmeans(m4_t_vegpresence, ~ status, weights = "proportional")
summary(em)
pairs(em)
#DECISSION:
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_ch4<- m1_gaus_nostrata
va_complexmodel_ch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_ch4, 
                        model2 = va_complexmodel_ch4,
                        data=va_ch4,
                        color_var = "vegpresence")

ggsave(filename = "VA_ch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




#GWP______-------

##2.1. CA_GWPco2andch4 (ok) -------

#Subset data and check transformation: 
ca_GWPco2andch4<- data4models %>% filter(casepilot=="CA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ca_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ca_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = ca_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  #DECISSION: Good residuals
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ca_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-9)
  Anova(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  #DECISSION: Good  residuals
  

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence does significantly improve results

#Save best models (simple and complex best)
ca_simplemodel_GWPco2andch4<- m1_gaus_nostrata
ca_complexmodel_GWPco2andch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_GWPco2andch4, 
                        model2 = ca_complexmodel_GWPco2andch4,
                        data=ca_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "CA_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_GWPco2andch4 (ok) -------

#Subset data and check transformation: 
cu_GWPco2andch4<- data4models %>% filter(casepilot=="CU"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: BAD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = cu_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  summary(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  #DECISSION: Good enough residuals
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = cu_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-9)
  Anova(m4_t_vegpresence)
  summary(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  #DECISSION: GOOD residuals

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_GWPco2andch4<- m3_t_nostrata
cu_complexmodel_GWPco2andch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_GWPco2andch4, 
                        model2 = cu_complexmodel_GWPco2andch4,
                        data=cu_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "CU_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.3. DA_GWPco2andch4 (ok) -------

#Subset data and check transformation: 
da_GWPco2andch4<- data4models %>% filter(casepilot=="DA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_GWPco2andch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-11) 
Anova(m3_t_nostrata)
summary(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION: Better residuals, although they still have issues


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_GWPco2andch4,
                           family = t_family(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m4_t_vegpresence)
check_singularity(m4_t_vegpresence,tolerance = 1e-8)
r2(m4_t_vegpresence, tolerance = 1e-11)
Anova(m4_t_vegpresence)
summary(m4_t_vegpresence)
res<- simulateResiduals(m4_t_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m4_t_vegpresence)
#Better residuals, although they present issues

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
da_simplemodel_GWPco2andch4<- m3_t_nostrata
da_complexmodel_GWPco2andch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_GWPco2andch4, 
                        model2 = da_complexmodel_GWPco2andch4,
                        data=da_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "DA_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_GWPco2andch4 (ok) -------

#Subset data and check transformation: 
du_GWPco2andch4<- data4models %>% filter(casepilot=="DU"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: GOOD-enough residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = du_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  #DECISSION: DOES NOT CONVERGE
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = du_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  #DECISSION: good-enough

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does significantly improve fit (worth the extra complexity)

#Save best models (simple and complex best)
du_simplemodel_GWPco2andch4<- m1_gaus_nostrata
du_complexmodel_GWPco2andch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_GWPco2andch4, 
                        model2 = du_complexmodel_GWPco2andch4,
                        data=du_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "DU_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_GWPco2andch4 (ok) -------

#Subset data and check transformation: 
ri_GWPco2andch4<- data4models %>% filter(casepilot=="RI"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
summary(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: GOOD residuals


#SKIP! cannot include vegpresence as fixed with status, and gaus is already good
if(F){
  #Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
  m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                                data = ri_GWPco2andch4,
                                family = gaussian(),
                                dispformula = ~1)
  #Evaluate model:
  check_convergence(m2_gaus_vegpresence) #OK
  check_singularity(m2_gaus_vegpresence,tolerance = 1e-8) #
  r2(m2_gaus_vegpresence, tolerance = 1e-10)
  Anova(m2_gaus_vegpresence)
  summary(m2_gaus_vegpresence)
  res<- simulateResiduals(m2_gaus_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m2_gaus_vegpresence)
  #DECISSION: 

#SKIP t_options: 

#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_GWPco2andch4,
                        family = t_family(),
                        dispformula = ~1)
#Evaluate basic model:
check_convergence(m3_t_nostrata)
check_singularity(m3_t_nostrata,tolerance = 1e-8) 
r2(m3_t_nostrata, tolerance=1e-10) 
Anova(m3_t_nostrata)
res<- simulateResiduals(m3_t_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m3_t_nostrata)
#DECISSION:

#SKIP, cannot include vegpresence as fixed effect
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ri_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  summary(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  #DECISSION: 
}

#Compare models: gaussian simple is already the best 
# anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

#Save best models (simple and complex best)
ri_simplemodel_GWPco2andch4<- m1_gaus_nostrata
# ri_complexmodel_GWPco2andch4<- NA # Complex not possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_GWPco2andch4, 
                        # model2 = ri_complexmodel_GWPco2andch4,
                        data=ri_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "RI_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_GWPco2andch4 (ok) -------

#Subset data and check transformation: 
va_GWPco2andch4<- data4models %>% filter(casepilot=="VA"&ghgspecies=="GWPco2andch4")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="GWPco2andch4") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_GWPco2andch4,
                           family = gaussian(),
                           dispformula = ~1)
#Evaluate basic model:
check_convergence(m1_gaus_nostrata)
check_singularity(m1_gaus_nostrata,tolerance = 1e-8)
r2(m1_gaus_nostrata, tolerance = 1e-10) 
Anova(m1_gaus_nostrata)
res<- simulateResiduals(m1_gaus_nostrata)
plotQQunif(res)
plotResiduals(res)
check_residuals(m1_gaus_nostrata)
#DECISSION: BAD residuals 


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_GWPco2andch4,
                              family = gaussian(),
                              dispformula = ~1)
#Evaluate model:
check_convergence(m2_gaus_vegpresence)
check_singularity(m2_gaus_vegpresence,tolerance = 1e-8)
r2(m2_gaus_vegpresence, tolerance = 1e-10)
Anova(m2_gaus_vegpresence)
summary(m2_gaus_vegpresence)
res<- simulateResiduals(m2_gaus_vegpresence)
plotQQunif(res)
plotResiduals(res)
check_residuals(m2_gaus_vegpresence)
#DECISSION: BAD residuals 


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = va_GWPco2andch4,
                          family = t_family(),
                          dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m3_t_nostrata)
  check_singularity(m3_t_nostrata,tolerance = 1e-8) 
  r2(m3_t_nostrata, tolerance=1e-10) 
  Anova(m3_t_nostrata)
  res<- simulateResiduals(m3_t_nostrata)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m3_t_nostrata)
  #DECISSION: DOES NOT CONVERGE
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = va_GWPco2andch4,
                             family = t_family(),
                             dispformula = ~1)
  #Evaluate basic model:
  check_convergence(m4_t_vegpresence)
  check_singularity(m4_t_vegpresence,tolerance = 1e-8)
  r2(m4_t_vegpresence, tolerance = 1e-10)
  Anova(m4_t_vegpresence)
  res<- simulateResiduals(m4_t_vegpresence)
  plotQQunif(res)
  plotResiduals(res)
  check_residuals(m4_t_vegpresence)
  #DECISSION:GOOD residuals

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_GWPco2andch4<- m1_gaus_nostrata
va_complexmodel_GWPco2andch4<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_GWPco2andch4, 
                        model2 = va_complexmodel_GWPco2andch4,
                        data=va_GWPco2andch4,
                        color_var = "vegpresence")

ggsave(filename = "VA_GWPco2andch4_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


#3. Create bestmodel-list -------

#Save lists of models: 
  #1. simple_model_list (only status*season models), only RI casepilot, 
  #2. complex_model_list (status*season*vegpresence), rest of casepilots 
#All models  named as: casepilot_ghgspecies


#List of simple models for each casepilot*ghgspecies: only appropriate for RI 
simplemodel_list_allghg<- list(
  "RI_co2" = ri_simplemodel_co2,
  "RI_ch4" = ri_simplemodel_ch4,
  "RI_GWPco2andch4" = ri_simplemodel_GWPco2andch4
)

#List of complex models for each casepilot*ghgspecies: all but RI 
complexmodel_list_allghg<- list(
  "CA_co2" = ca_complexmodel_co2,
  "CU_co2" = cu_complexmodel_co2, 
  "DA_co2" = da_complexmodel_co2,
  "DU_co2" = du_complexmodel_co2,
  "VA_co2" = va_complexmodel_co2,
  "CA_ch4" = ca_complexmodel_ch4,
  "CU_ch4" = cu_complexmodel_ch4, 
  "DA_ch4" = da_complexmodel_ch4,
  "DU_ch4" = du_complexmodel_ch4,
  "VA_ch4" = va_complexmodel_ch4,
  "CA_GWPco2andch4" = ca_complexmodel_GWPco2andch4,
  "CU_GWPco2andch4" = cu_complexmodel_GWPco2andch4, 
  "DA_GWPco2andch4" = da_complexmodel_GWPco2andch4,
  "DU_GWPco2andch4" = du_complexmodel_GWPco2andch4,
  "VA_GWPco2andch4" = va_complexmodel_GWPco2andch4
)




## 3.1. Overall fit -----
#Structure, pseudoR2 of best model:
simplemodel_fit<- get_pseudoR2s(simplemodel_list_allghg)
complexmodel_fit<- get_pseudoR2s(complexmodel_list_allghg)

#R2c and R2m of homocedastic model:
simplemodel_homoc_r2<- get_R2s(simplemodel_list_allghg)
complexmodel_homoc_r2<- get_R2s(complexmodel_list_allghg)

#DHARMa Residual diagnostics of best model:
simplemodel_resid_diag_summary <- summarize_dharma_diagnostics(simplemodel_list_allghg)
complexmodel_resid_diag_summary <- summarize_dharma_diagnostics(complexmodel_list_allghg)


##3.2. Significance of effects-----
#Effect of status will first be assessed via Anova (best_model,type=3): is there an effect averaging across all seasons? yes/no, is the effect dependent on season?, is the effect depenent on vegpresence? 

#Obtain significance for main effects of best_models
simplemodel_results_anova<-get_all_anova_results(simplemodel_list_allghg)
complexmodel_results_anova<-get_all_anova_results(complexmodel_list_allghg)


all_simplemodel_outputs<- simplemodel_results_anova %>% 
  full_join(simplemodel_fit) %>% 
  full_join(simplemodel_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))

all_complexmodel_outputs<- complexmodel_results_anova %>% 
  full_join(complexmodel_fit) %>% 
  full_join(complexmodel_homoc_r2 %>% dplyr::select(dataset, homoced_R2m, homoced_R2c))


#Save tables with the models summary outputs: 

#Simple models (RI)
write.csv(x = all_simplemodel_outputs, file = paste0(path_2_modeloutputs,"Summary_RI_chambermodels.csv"),row.names = F)

write.csv(x = simplemodel_resid_diag_summary, file = paste0(path_2_modeloutputs,"Residualtests_RI_chambermodels.csv"),row.names = F)


#Complex models
write.csv(x = all_complexmodel_outputs, file = paste0(path_2_modeloutputs,"Summary_rest_chambermodels.csv"),row.names = F)

write.csv(x = complexmodel_resid_diag_summary, file = paste0(path_2_modeloutputs,"Residualtests_rest_chambermodels.csv"),row.names = F)


##3.3. EMEANs & post-hoc----

#When calculating emmeans: I have to use the full variance-covariance structure to get accurate SEs. 
#If the model is gaussian, we use t-test with Degrees of freedom available (via Satterthwaite approximation)
#IF the model is t_family, it will have Inf df, not an issue but forces to use z-test for post-hocs instead (No DF calculation possible)

#CLDs (group letters) should be assigned to emmeans based on the above pairwise tests, do not repeat the test, but assign letters based on already calcualted tests via (multcompLetters function). Additionally, re-code letters so that they give emmean ranking info (letter "a" denotes the smallest-emmean significant group, letter"b" the signficant group with the next lowest emmean

#Calculate relevant emmeans and post-hocs tests: separately for simple model list (status*season, only for RI) and complex model list (status*season*vegpresence, for rest of casepilots)



#DEDICATED FUNCTIONS:

##function for emmeans------
#=========================================================-
# FUNCTION: custom_weighted_emmeans_backtrans
#---------------------------------------------------------_
# Calculates weighted emmeans (with custom weights) from a glmmTMB model.
# Returns both model-scale and back-transformed emmeans with proper SEs and CIs (using full covariance matrix).
# Back-transformation is done via the supplied BestNormalize object or custom inverse.
#=========================================================-

custom_weighted_emmeans_backtrans <- function(
    model_object,       # glmmTMB model
    trans_object,       # BestNormalize object (or similar, must support predict(inverse=TRUE))
    grouping_vars,      # character vector of fixed effects to group by
    custom_weights = FALSE, # logical: use custom_weight_df or not
    custom_weight_df = NULL, # optional: data frame with all level combinations + prop_weight
    conf_level = 0.95   # confidence level for CIs
) {
  library(emmeans)
  library(dplyr)
  library(numDeriv)
  
  #-------------------------------#
  # Helper: critical value
  #-------------------------------#
  get_crit_value <- function(df, conf_level) {
    alpha <- 1 - conf_level
    if (is.infinite(df)) qnorm(1 - alpha / 2) else qt(1 - alpha / 2, df)
  }
  
  #-------------------------------#
  # Extract factor variables
  #-------------------------------#
  all_terms <- attr(terms(model_object), "term.labels")
  factor_vars <- all_terms[
    sapply(all_terms, function(v) is.factor(model.frame(model_object)[[v]]))
  ]
  
  # Fallback: if none detected, warn user
  if (length(factor_vars) == 0)
    stop("No factor variables detected in the model.")
  
  #-------------------------------#
  # Compute full emmeans object
  #-------------------------------#
  full_emms_object <- emmeans(model_object, specs = factor_vars)
  vc <- vcov(full_emms_object)
  emm_df <- as.data.frame(full_emms_object)
  
  #-------------------------------#
  # Apply weighting scheme
  #-------------------------------#
  if (custom_weights) {
    if (is.null(custom_weight_df))
      stop("If custom_weights = TRUE, you must provide custom_weight_df.")
    
    # Merge with custom weights
    emm_df <- emm_df %>%
      left_join(custom_weight_df, by = factor_vars)
    
    if (any(is.na(emm_df$prop_weight)))
      stop("Mismatch between factor levels and custom_weight_df.")
    
  } else {
    # Use equal weights within each grouping
    emm_df <- emm_df %>%
      group_by(across(all_of(grouping_vars))) %>%
      mutate(prop_weight = 1 / n()) %>%
      ungroup()
  }
  
  # Create unique labels for vcov indexing
  emm_df <- emm_df %>%
    mutate(row_label = do.call(paste, c(across(all_of(factor_vars)), sep = " ")))
  
  #-------------------------------#
  # Group by user-specified variables
  #-------------------------------#
  grouped_df <- emm_df %>%
    group_by(across(all_of(grouping_vars))) %>%
    group_split()
  
  results_list <- list()
  
  #-------------------------------#
  # Loop over groups
  #-------------------------------#
  for (group in grouped_df) {
    group_label <- group %>%
      dplyr::select(all_of(grouping_vars)) %>%
      slice(1)
    
    w <- group$prop_weight / sum(group$prop_weight)
    mu <- group$emmean
    vc_sub <- vc[group$row_label, group$row_label, drop = FALSE]
    
    # Model-scale weighted mean and SE
    emmean_og <- sum(w * mu)
    SE_og <- sqrt(as.numeric(t(w) %*% vc_sub %*% w))
    
    df_val <- unique(group$df)
    if (length(df_val) != 1) stop("Non-unique df within group.")
    crit <- get_crit_value(df_val, conf_level)
    
    lower.CL_og <- emmean_og - crit * SE_og
    upper.CL_og <- emmean_og + crit * SE_og
    
    #-------------------------------#
    # Back-transform using delta method
    #-------------------------------#
    inv_fun <- function(x) predict(trans_object, x, inverse = TRUE)
    deriv_fun <- function(x) numDeriv::grad(inv_fun, x)
    
    emmean_bt <- inv_fun(emmean_og)
    deriv_val <- deriv_fun(emmean_og)
    SE_bt <- abs(deriv_val) * SE_og
    
    lower.CL_bt <- emmean_bt - crit * SE_bt
    upper.CL_bt <- emmean_bt + crit * SE_bt
    
    #-------------------------------#
    # Store results
    #-------------------------------#
    res <- group_label %>%
      mutate(
        emmean_og = emmean_og,
        SE_og = SE_og,
        df = df_val,
        lower.CL_og = lower.CL_og,
        upper.CL_og = upper.CL_og,
        emmean_bt = emmean_bt,
        SE_bt = SE_bt,
        lower.CL_bt = lower.CL_bt,
        upper.CL_bt = upper.CL_bt
      )
    
    collapsed_vars <- setdiff(factor_vars, grouping_vars)
    for (v in collapsed_vars) res[[v]] <- NA
    
    results_list <- append(results_list, list(res))
  }
  
  #-------------------------------#
  # Combine all results
  #-------------------------------#
  final_df <- bind_rows(results_list) %>%
    dplyr::select(all_of(grouping_vars),
                  emmean_og, SE_og, lower.CL_og, upper.CL_og,
                  emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt,
                  df, everything())
  
  return(final_df)
}


##function for contrats--------
#=========================================================-
# FUNCTION: custom_pairwise_contrasts_fullvcov
#---------------------------------------------------------_
# Calculates weighted pairwise contrasts using full covariance
# structure among emmeans, not assuming independence.
# Returns both model-scale and back-transformed contrasts.
#=========================================================-

custom_pairwise_contrasts_fullvcov <- function(
    model_object,              # glmmTMB model
    trans_object,              # BestNormalize or similar transformation object
    compare_var,               # factor for which contrasts are computed
    group_vars = character(0), # optional grouping variables
    use_custom_weights = FALSE,# whether to use custom weights
    custom_weight_df = NULL,   # optional custom weights df (must include prop_weight)
    conf_level = 0.95,         # CI level
    return_emmeans = FALSE     # return weighted emmeans if TRUE
) {
  library(emmeans)
  library(dplyr)
  library(purrr)
  library(numDeriv)
  
  #---------------------------------------------#
  # Helper function: critical value for CI
  #---------------------------------------------#
  get_crit_value <- function(df, conf_level) {
    alpha <- 1 - conf_level
    if (is.infinite(df)) qnorm(1 - alpha / 2) else qt(1 - alpha / 2, df)
  }
  
  #---------------------------------------------#
  # STEP 1. Extract factor variables from model
  #---------------------------------------------#
  all_terms <- attr(terms(model_object), "term.labels")
  factor_vars <- all_terms[
    sapply(all_terms, function(v) is.factor(model.frame(model_object)[[v]]))
  ]
  if (length(factor_vars) == 0)
    stop("No factor variables detected in model.")
  
  #---------------------------------------------#
  # STEP 2. Compute emmeans + variance-covariance matrix
  #---------------------------------------------#
  full_emms <- emmeans(model_object, specs = factor_vars)
  vc <- vcov(full_emms)
  emm_df <- as.data.frame(full_emms)
  
  #---------------------------------------------#
  # STEP 3. Apply weights (custom or equal)
  #---------------------------------------------#
  if (use_custom_weights) {
    if (is.null(custom_weight_df))
      stop("If use_custom_weights = TRUE, custom_weight_df must be provided.")
    
    emm_df <- emm_df %>%
      left_join(custom_weight_df, by = factor_vars)
    
    if (any(is.na(emm_df$prop_weight)))
      stop("Mismatch between factor levels and custom_weight_df.")
  } else {
    # Use equal weights across all combinations of the factor_vars
    emm_df <- emm_df %>%
      group_by(across(all_of(factor_vars))) %>%
      mutate(prop_weight = 1 / n()) %>%
      ungroup()
  }
  
  # Create row labels for vcov indexing
  emm_df <- emm_df %>%
    mutate(row_label = do.call(paste, c(across(all_of(factor_vars)), sep = " ")))
  
  #---------------------------------------------#
  # STEP 4. Compute weighted emmeans (model scale)
  #---------------------------------------------#
  grouping_vars_all <- unique(c(group_vars, compare_var))
  grouped_df <- emm_df %>%
    group_by(across(all_of(grouping_vars_all))) %>%
    group_split()
  
  weighted_emmeans <- list()
  
  for (group in grouped_df) {
    group_label <- group %>% dplyr::select(all_of(grouping_vars_all)) %>% slice(1)
    
    w <- group$prop_weight / sum(group$prop_weight)
    mu <- group$emmean
    vc_sub <- vc[group$row_label, group$row_label, drop = FALSE]
    
    emmean_og <- sum(w * mu)
    SE_og <- sqrt(as.numeric(t(w) %*% vc_sub %*% w))
    
    df_val <- unique(group$df)
    if (length(df_val) != 1) stop("Non-unique df within group.")
    
    weighted_emmeans <- append(weighted_emmeans, list(
      group_label %>%
        mutate(weights = list(w),
               row_label = list(group$row_label),
               emmean_og = emmean_og,
               SE_og = SE_og,
               df = df_val)
    ))
  }
  
  emmean_df <- bind_rows(weighted_emmeans)
  
  #---------------------------------------------#
  # STEP 5. Compute pairwise contrasts
  #---------------------------------------------#
  contrast_groups <- if (length(group_vars) > 0) {
    emmean_df %>% group_by(across(all_of(group_vars))) %>% group_split()
  } else {
    list(emmean_df)
  }
  
  results <- purrr::map_dfr(contrast_groups, function(group_df) {
    levels <- unique(group_df[[compare_var]])
    combs <- combn(levels, 2, simplify = FALSE)
    m <- length(combs)  # number of pairwise comparisons for Sidak correction
    
    df_val <- unique(group_df$df)
    if (length(df_val) != 1) stop("Non-unique df within group.")
    crit <- get_crit_value(df_val, conf_level)
    
    purrr::map_dfr(combs, function(pair) {
      g1 <- group_df %>% filter(!!sym(compare_var) == pair[1])
      g2 <- group_df %>% filter(!!sym(compare_var) == pair[2])
      
      #---------------------------------------------#
      # Model-scale contrast (using full vcov)
      #---------------------------------------------#
      w1 <- unlist(g1$weights)
      w2 <- unlist(g2$weights)
      L <- numeric(nrow(vc))
      idx1 <- match(unlist(g1$row_label), rownames(vc))
      idx2 <- match(unlist(g2$row_label), rownames(vc))
      L[idx1] <- w1
      L[idx2] <- -w2
      
      estimate_og <- g1$emmean_og - g2$emmean_og
      SE_og <- sqrt(as.numeric(t(L) %*% vc %*% L))
      stat <- estimate_og / SE_og
      
      # test type + p-value
      if (is.infinite(df_val)) {
        test_used <- "Z.test"
        p_value_raw <- 2 * pnorm(-abs(stat))
      } else {
        test_used <- "T.test"
        p_value_raw <- 2 * pt(-abs(stat), df = df_val)
      }
      
      # Sidak correction for multiple comparisons
      p_value_sidak <- 1 - (1 - p_value_raw)^m
      
      lower.CL_og <- estimate_og - crit * SE_og
      upper.CL_og <- estimate_og + crit * SE_og
      
      #---------------------------------------------#
      # Back-transformation (delta method)
      #---------------------------------------------#
      inv_fun <- function(x) predict(trans_object, x, inverse = TRUE)
      deriv_fun <- function(x) numDeriv::grad(inv_fun, x)
      
      mu1_bt <- inv_fun(g1$emmean_og)
      mu2_bt <- inv_fun(g2$emmean_og)
      d1 <- deriv_fun(g1$emmean_og)
      d2 <- deriv_fun(g2$emmean_og)
      
      estimate_bt <- mu1_bt - mu2_bt
      SE_bt <- sqrt((d1^2 * g1$SE_og^2) + (d2^2 * g2$SE_og^2))
      lower.CL_bt <- estimate_bt - crit * SE_bt
      upper.CL_bt <- estimate_bt + crit * SE_bt
      
      #---------------------------------------------#
      # Return results
      #---------------------------------------------#
      res <- tibble(
        contrast = paste0(pair[1], " - ", pair[2]),
        !!compare_var := NA,
        estimate_og = estimate_og,
        SE_og = SE_og,
        lower.CL_og = lower.CL_og,
        upper.CL_og = upper.CL_og,
        estimate_bt = estimate_bt,
        SE_bt = SE_bt,
        lower.CL_bt = lower.CL_bt,
        upper.CL_bt = upper.CL_bt,
        df = df_val,
        stat.ratio = stat,
        test_used = test_used,
        # p.value_raw = p_value_raw, #omit raw p-value
        p.value = p_value_sidak
      )
      
      if (length(group_vars) > 0)
        res <- bind_cols(group_df[1, group_vars, drop = FALSE], res)
      
      return(res)
    })
  })
  
  #---------------------------------------------#
  # STEP 6. Organize output
  #---------------------------------------------#
  out <- results %>%
    relocate(all_of(group_vars), .before = everything()) %>%
    relocate(contrast, .before = everything()) %>% 
    dplyr::select(-all_of(compare_var))
  
  if (return_emmeans) {
    return(list(emmeans = emmean_df, contrasts = out))
  } else {
    return(out)
  }
}



#Simple models (ok)-------

#Only applicable to RI

#Emmeans are extracted, pairwise tests done (based on model distribution family), CLDs are assigned to significantly different groups (and re-ordered based on emmean ranking). WE save both pairwise post-hoc tests in original (model scale) and back-transformed scales; AND emmeanCLDs in original and back-transformed scales. 

#Equal Weights used for all EMMs: this assumes that all seasons and all status have the same weight for each other (appropriate). 


#Initialize simple model results list
simple_comparison_list<- list()

# Loop to extract all emmeans and perform pairwise tests for appropriate comparisons for every casepilot*ghg simple model. 

for (dataset in names(simplemodel_list_allghg)) {
  #get model
  cp_model <- simplemodel_list_allghg[[dataset]]
  #get transformation object
  cp_trans_obj<- bn_list[[dataset]]
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  # Obtain emmean object for each comparison: status, season, status_within_season. Using custom function for consistency. 
  status_emmeans <- custom_weighted_emmeans_backtrans(
    model_object = cp_model, trans_object = cp_trans_obj,
    grouping_vars = c("status"),
    custom_weights = F)%>% 
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA)
  
  season_emmeans <- custom_weighted_emmeans_backtrans(
    model_object = cp_model, trans_object = cp_trans_obj,
    grouping_vars = c("season"),
    custom_weights = F)%>% 
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA)
  
  statuswithinseason_emmeans <- custom_weighted_emmeans_backtrans(
    model_object = cp_model, trans_object = cp_trans_obj,
    grouping_vars = c("status","season"),
    custom_weights = F) %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  
  #Calculate contrasts for each comparison and add identifying columns:   
  #Status comparison:
  status_contrasts <- custom_pairwise_contrasts_fullvcov(
    model_object = cp_model, trans_object = cp_trans_obj,
    compare_var = "status",
    use_custom_weights = FALSE)%>% 
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season= NA)
  
  season_contrasts <- custom_pairwise_contrasts_fullvcov(
    model_object = cp_model, trans_object = cp_trans_obj,
    compare_var = "season",
    use_custom_weights = FALSE)%>% 
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status= NA)
  
  statuswithinseason_contrasts <- custom_pairwise_contrasts_fullvcov(
    model_object = cp_model, trans_object = cp_trans_obj,
    compare_var = "status",group_vars = "season",
    use_custom_weights = FALSE)%>% 
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name)
  
  
  #Join all emmeans and all contrasts:
  all_emmeans<- status_emmeans %>%
    full_join(season_emmeans) %>%
    full_join(statuswithinseason_emmeans)
  
  all_contrasts<- status_contrasts %>%
    full_join(season_contrasts) %>%
    full_join(statuswithinseason_contrasts)
  
  #Store in named list: 
  simple_comparison_list[[dataset]] <- list(
    emmeans_og = all_emmeans,
    posthoc_comparisons = all_contrasts)
  
}


#Remove within-loop objects
rm(all_contrasts,all_emmeans, cp_model,cp_trans_obj, casepilot_name, ghgspecies,
   status_emmeans, season_emmeans, statuswithinseason_emmeans,
   status_contrasts,season_contrasts,statuswithinseason_contrasts)


#Get pairwise posthoc tests (in model scale): 
#T-test or Z-test depending on model distribution family used (automatically assigned by contrasts (method="pairwise)). p.value is was adjusted for multiple comparisons sidak 
simplemodel_posthoc_tests <- purrr::map_dfr(simple_comparison_list, "posthoc_comparisons") %>%
  #Identify the model_distribution based on df estimation:
  mutate(model_distribution=if_else(df==Inf, "t_family", "gaussian")) %>% 
  dplyr::select(casepilot, ghgspecies, model_distribution,comparison, test_used, contrast, season,
                estimate_og, SE_og, lower.CL_og, upper.CL_og,
                estimate_bt, SE_bt, lower.CL_bt, upper.CL_bt,
                df, stat.ratio, p.value)

#Save pairwise posthoc tests as csv.
#Save post-hoc comparisons in model scale and back-transformed (to go to supplementary table)
write.csv(x = simplemodel_posthoc_tests, file = paste0(path_2_modeloutputs,"Posthoctests_RI_chambermodels.csv"),row.names = F)


#Obtain CLDs (letter-groups) based on post-hoc tests for every comparison:
CLD_letters <- simplemodel_posthoc_tests %>%
  #leave only (1) variables that identify unique comparisons, (2)contrast column, (3) pvalue column
  #season is kept to account for status_within_season comparison. 
  dplyr::select(casepilot, ghgspecies, comparison, season, contrast, p.value) %>% 
  #Group by comparison identifiers
  group_by(casepilot, ghgspecies, comparison, season) %>%
  #Remove any spaces from contrast column (mutcompLetters expects levels to be only separated by a hyphen "-")
  mutate(contrast=gsub(" ","", contrast)) %>% 
  #Nest data for each group
  nest() %>%
  # Step 3: Apply multcompLetters to each group
  mutate(letters = purrr::map(data, function(group_df) {
    # Extract levels and p-values 
    contrast_matrix <- group_df %>%
      dplyr::select(contrast, p.value) %>%
      deframe()
    
    # Apply multcompLetters
    multcompLetters(contrast_matrix)$Letters %>%
      enframe(name = "level", value = "cld_group")
  })) %>%
  # Step 4: Unnest results
  dplyr::select(-data) %>%
  unnest(letters) %>% 
  #Reformat level to fit the appropriate columns: season or status
  mutate(season=if_else(level%in%c("S1","S2","S3","S4"),level,season),
         status=if_else(level%in%c("Altered","Preserved","Restored"),level, NA),
         #Mantain level to identify nested comparisons 
         seasonlevel=if_else(comparison=="status_within_season", season, NA)) %>% 
  dplyr::select(casepilot, ghgspecies, comparison,status,season,seasonlevel, cld_group)


#Get emmeans (in model scale) from loop list 
simplemodel_emmeans<- purrr::map_dfr(simple_comparison_list, "emmeans_og")

#Add CLD and back transform emmeans
simplemodel_emmeansCLD <- simplemodel_emmeans %>%
  #ADD CLD group-letters (re-coding them so that "a" always identifies the significant group with lowest emmean, b the next lowest emmean, and so on...)
  left_join(CLD_letters, by = c("status", "comparison", "ghgspecies", "casepilot", "season")) %>%
  #Separate each comparison group to do the letter re-coding
  group_split(casepilot, ghgspecies, comparison, seasonlevel) %>%
  map_dfr(function(group_df) {
    # Step 1: Order by emmean
    group_df <- group_df %>% arrange(emmean_og)
    # Step 2: Extract and map letters
    letter_list <- str_split(group_df$cld_group, "")
    unique_letters <- unique(unlist(letter_list))
    new_letters <- letters[seq_along(unique_letters)]
    letter_map <- setNames(new_letters, unique_letters)
    # Step 3: Apply mapping
    group_df <- group_df %>%
      mutate(cld_group = map_chr(letter_list, ~ paste0(letter_map[.x], collapse = "")))
    return(group_df)
  }) %>% 
  #remove grouping season variable
  dplyr::select(-seasonlevel) %>% 
  mutate(weights_used="equal") %>% #add weights info for all (all equal)
  #Format final emmean_CLD: select columns
  dplyr::select(casepilot, ghgspecies, comparison, status, season,
                df, weights_used,
                emmean_og, SE_og,lower.CL_og, upper.CL_og,
                cld_group,
                emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt) %>% 
  #Format final emmean_CLD: arrange values
  arrange(casepilot,ghgspecies, comparison, season,status )

#Save emmeans and groupletters (back-transformed):  
write.csv(x = simplemodel_emmeansCLD, file = paste0(path_2_modeloutputs,"EmmeansCLD_RI_chambermodels.csv"),row.names = F)





#Complex models (ok)------

#We need to take into account the proportion of vegpresence in the field to estimate the status, season and status_within_season emmeans. For status_within_vegpresence, equal weights should be applied (to give equal importance to each season). 

#Overall status: weights should be used to account for different vegpresence in different status and  season, but all seasons should have the same combined weight for status emmeans and comparisons (to give the same importance to each season regardless of how small deviations in number of observations between seasons)

#RATIONALE FOR CUSTOM WEIGHTS:
#WE NEED TO APPLY weights to scale the vegpresence effects according to their proportion in the field. This applies for comparisions: status, season, and status_within_season. For status_within_vegpresence comparisons, equal weights should be applied (to give equal importance to each season). 
#Weights should be calculated as the proportion of veg/noveg in every combination of status*season so they sum 1 (and give equal importance to all seasons).
#The function "custom_weighted_emmeans" already re-normalizes the weigths as needed (for status and for season emmeans, no need to re-normalize for status*season. It calculates the associated weighted SE taking into account the variance-covariance structure of the model.

#Status comparison: each status estimate is calculated taking into account its particular vegpresence composition while giving equal importance to all seasons.

#Season comparison: each seasonal estimate is calculated taking into account the seasonally variable vegpresence composition, while giving equal importance to all status (overall seasonal effect across all status).

#Status_within_season comparison: each seasonal estimate of each status takes into account its particular vegpresence composition. 

#For status_within_vegpresence, we are only averaging across season and we want to have equal impact of all seasons (overall yearly effect), use weights="equal".


##Calculate weights------
#First, calculate vegpresence proportions at every level of casepilot, status and season 

#With weights calculated from chamber distribution, emmeans are proportional to the data used to fit the model, with the exception that different N across seasons does not influence the importance of each season in the status emmean (all season treated equally).
{
casepilot_weights<- data4models %>% 
  dplyr::select(plotcode, casepilot, season, status, vegpresence, subsite) %>% 
  #remove duplicates (same plotcode, different ghgspecies)
  distinct() %>% 
  dplyr::group_by(casepilot, season, status, vegpresence, subsite) %>%
  #Calculate vegpresence deployment counts for each (subsite) sampling and pivot_wider
  summarise(n_vegpresence=sum(!is.na(vegpresence)), .groups = "drop") %>% 
  pivot_wider(names_from = vegpresence, values_from = n_vegpresence,values_fill = 0) %>% 
  #Calculate vegpresence proportion for every (subsite) sampling.
  mutate(sum_all=`Non-vegetated`+Vegetated) %>% 
  mutate(Vegetated=Vegetated/sum_all,
         `Non-vegetated`=`Non-vegetated`/sum_all) %>% 
  dplyr::select(-sum_all) %>% 
  pivot_longer(cols = c(`Non-vegetated`,Vegetated), names_to = "vegpresence", values_to = "proportions") %>% 
  #Calculate the average of each status (giving equal weight to the two subsites)
  dplyr::group_by(casepilot, status, season, vegpresence) %>% 
  summarise(prop_weight=mean(proportions, na.rm = T), .groups = "drop")

#For CURONIAN: Calculate overall veg/noveg of curonian, across status and seasons. No seasonal variability was observed in vegetation extent, additionally due to initial incorrect restored site boundaries, the proportions of vegpresence in these sites are not representative of actual site composition. We use constant average proportion across all status and seasons. 
cu_weights<- data4models %>% 
  filter(casepilot=="CU") %>% 
  group_by(plotcode, casepilot,vegpresence) %>% 
  distinct() %>% #to remove duplicates (same plotcode, different ghgspecies)
  group_by(casepilot, vegpresence) %>% 
  #calculate overall vegpresence counts (across all seasons and all subsites)
  summarise(n_vegpresence=sum(!is.na(vegpresence)), .groups = "drop") %>% 
  pivot_wider(names_from = vegpresence, values_from = n_vegpresence,values_fill = 0) %>% 
  #Calculate vegpresence proportion 
  mutate(sum_all=`Non-vegetated`+Vegetated) %>% 
  mutate(Vegetated=Vegetated/sum_all,
         `Non-vegetated`=`Non-vegetated`/sum_all) %>% 
  dplyr::select(-sum_all)

#Override CU composition: constant actual composition, differences in chamber deployment are due to systematic sampling biass in Restored sites
casepilot_weights<-casepilot_weights %>% 
  mutate(prop_weight=if_else(casepilot=="CU"&vegpresence=="Vegetated", cu_weights %>% pull(`Vegetated`),prop_weight),
         prop_weight=if_else(casepilot=="CU"&vegpresence=="Non-vegetated",cu_weights %>% pull(`Non-vegetated`),prop_weight))
}


##Weighted-emmean-loop------

#Calculate emmeans for each relevant comparison, using weights when appropriate:

#Intialize list for complex models
complex_comparison_list<- list()

# Loop to extract all emmeans and perform pairwise tests for appropriate comparisons for every casepilot complexbestmodel. 
for (dataset in names(complexmodel_list_allghg)) {
  #get model
  cp_model <- complexmodel_list_allghg[[dataset]]
  #get transformation object
  cp_trans_obj<- bn_list[[dataset]]
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  #Prepare casepilot custom weights: 
  weight_df<- casepilot_weights %>% filter(casepilot==casepilot_name) %>% dplyr::select(-casepilot)
  
  #Obtain weighted emmeans for comparisons: status, season, status_within_season
  #Obtain equal-weighted emmeans for comparison: status_within_vegpresence
  status_emmeans <- custom_weighted_emmeans_backtrans(
    model_object = cp_model, trans_object = cp_trans_obj,
    grouping_vars = c("status"),
    custom_weights = T, custom_weight_df = weight_df) %>% 
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA, vegpresence = NA, weights_used = "custom")
  
  season_emmeans <- custom_weighted_emmeans_backtrans(
    model_object = cp_model, trans_object = cp_trans_obj,
    grouping_vars = c("season"),
    custom_weights = T, custom_weight_df = weight_df) %>% 
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           status = NA, vegpresence = NA, weights_used = "custom")
  
  statuswithinseason_emmeans <- custom_weighted_emmeans_backtrans(
    model_object = cp_model, trans_object = cp_trans_obj,
    grouping_vars = c("status","season"),
    custom_weights = T, custom_weight_df = weight_df) %>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name, 
           vegpresence = NA, weights_used = "custom")
  
  statuswithinvegpresence_emmeans <- custom_weighted_emmeans_backtrans(
    model_object = cp_model, trans_object = cp_trans_obj,
    grouping_vars = c("status","vegpresence"),
    custom_weights = F) %>% #Using equal weights (same importance to all seasons)
    mutate(comparison = "status_within_vegpresence",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name, 
           season = NA, weights_used = "equal")
  
  #Calculate contrasts for each comparison and add identifying columns:   
  #Status comparison:
  status_contrasts <- custom_pairwise_contrasts_fullvcov(
    model_object = cp_model, trans_object = cp_trans_obj,
    compare_var = "status",
    use_custom_weights = T, custom_weight_df = weight_df)%>% 
    mutate(comparison = "status",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           season = NA, vegpresence = NA, weights_used = "custom")
  
  season_contrasts <- custom_pairwise_contrasts_fullvcov(
    model_object = cp_model, trans_object = cp_trans_obj,
    compare_var = "season",
    use_custom_weights = T, custom_weight_df = weight_df)%>% 
    mutate(comparison = "season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name,
           vegpresence = NA, weights_used = "custom")
  
  statuswithinseason_contrasts <- custom_pairwise_contrasts_fullvcov(
    model_object = cp_model, trans_object = cp_trans_obj,
    compare_var = "status", group_vars = "season",
    use_custom_weights = T, custom_weight_df = weight_df)%>%
    mutate(comparison = "status_within_season",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name, 
           vegpresence = NA, weights_used = "custom")
  
  statuswithinvegpresence_contrasts <- custom_pairwise_contrasts_fullvcov(
    model_object = cp_model, trans_object = cp_trans_obj,
    compare_var = "status", group_vars = "vegpresence",
    use_custom_weights = F)%>%#Using equal weights (same importance to all seasons)
    mutate(comparison = "status_within_vegpresence",
           ghgspecies = ghgspecies,
           casepilot = casepilot_name, 
           season = NA, weights_used = "equal")
  
  #Join all emmeans and all contrasts:
  all_emmeans<- status_emmeans %>%
    full_join(season_emmeans) %>%
    full_join(statuswithinseason_emmeans) %>% 
    full_join(statuswithinvegpresence_emmeans)
  
  all_contrasts<- status_contrasts %>%
    full_join(season_contrasts) %>%
    full_join(statuswithinseason_contrasts) %>% 
    full_join(statuswithinvegpresence_contrasts)
  
  #Store in named list: 
  complex_comparison_list[[dataset]] <- list(
    emmeans_og = all_emmeans,
    posthoc_comparisons = all_contrasts)
  
}

#Remove within-loop objects
rm(all_contrasts,all_emmeans, cp_model, cp_trans_obj, casepilot_name, ghgspecies,weight_df,
   status_emmeans, season_emmeans, statuswithinseason_emmeans,statuswithinvegpresence_emmeans,
   status_contrasts,season_contrasts,statuswithinseason_contrasts,statuswithinvegpresence_contrasts)


#Get pairwise posthoc tests (in model scale and back-transformed): 
#T-test or Z-test depending on model distribution family used (automatically assigned by contrasts (method="pairwise)). p.value is adjusted for multiple comparisons sidak 
complexmodel_customweight_posthoc_tests <- purrr::map_dfr(complex_comparison_list, "posthoc_comparisons") %>%
  #Identify the test used (and model_distribution) based on df estimation:
  mutate(model_distribution=if_else(df==Inf, "t_family", "gaussian")) %>% 
  dplyr::select(casepilot, ghgspecies, model_distribution, comparison, test_used, 
                contrast,season, vegpresence, 
                estimate_og, SE_og, lower.CL_og, upper.CL_og, 
                estimate_bt, SE_bt, lower.CL_bt, upper.CL_bt, 
                df, stat.ratio, p.value)


#Save pairwise posthoc tests as csv.
#Save post-hoc comparisons in model scale (to go to supplementary table)
write.csv(x = complexmodel_customweight_posthoc_tests, file = paste0(path_2_modeloutputs,"Posthoctests_rest_chambermodels.csv"),row.names = F)



#Obtain CLDs (letter-groups) based on post-hoc tests for every comparison:
CLD_letters <- complexmodel_customweight_posthoc_tests %>%
  #leave only (1) variables that identify unique comparisons, (2)contrast column, (3) p.value column
  #season and vegpresence are kept to account for status_within_season and status_within_vegpresence comparison. 
  dplyr::select(casepilot, ghgspecies, comparison, season, vegpresence, contrast, p.value) %>%
  #Group by comparison identifiers
  group_by(casepilot, ghgspecies, comparison, season, vegpresence) %>%
  #Remove any spaces from contrast column (mutcompLetters expects levels separated by a hyphen "-")
  mutate(contrast=gsub(" ","", contrast)) %>% 
  #Nest data for each group
  nest() %>%
  # Step 3: Apply multcompLetters to each group
  mutate(letters = purrr::map(data, function(group_df) {
    # Extract levels and p-values 
    contrast_matrix <- group_df %>%
      dplyr::select(contrast, p.value) %>%
      deframe()
    # Apply multcompLetters
    multcompLetters(contrast_matrix)$Letters %>%
      enframe(name = "level", value = "cld_group")
  })) %>%
  # Step 4: Unnest results
  dplyr::select(-data) %>%
  unnest(letters) %>% 
  #Reformat level to fit the appropriate columns: season or status
  mutate(season=if_else(level%in%c("S1","S2","S3","S4"),level,season),
         status=if_else(level%in%c("Altered","Preserved","Restored"),level, NA),
         #Mantain level to identify nested comparisons 
         seasonlevel=if_else(comparison=="status_within_season",season,NA),
         vegpresencelevel=if_else(comparison=="status_within_vegpresence",vegpresence,NA)) %>% 
  dplyr::select(casepilot, ghgspecies, comparison,status,season,vegpresence, seasonlevel, vegpresencelevel, cld_group)


#Get emmeans (in model scale) from loop list 
complexmodel_customweight_emmeans<- purrr::map_dfr(complex_comparison_list, "emmeans_og")

#Add CLD and back transform emmeans
complexmodel_customweight_emmeansCLD <- complexmodel_customweight_emmeans %>%
  #ADD CLD group-letters (re-coding them so that "a" always identifies the significant group with lowest emmean, b the next lowest emmean, and so on...)
  left_join(CLD_letters, by = c("status", "comparison", "ghgspecies", "casepilot", "season","vegpresence")) %>%
  #Separate each comparison group to do the letter re-coding
  group_split(casepilot, ghgspecies, comparison, seasonlevel,vegpresencelevel) %>%
  map_dfr(function(group_df) {
    # Step 1: Order by emmean
    group_df <- group_df %>% arrange(emmean_og)
    # Step 2: Extract and map letters
    letter_list <- str_split(group_df$cld_group, "")
    unique_letters <- unique(unlist(letter_list))
    new_letters <- letters[seq_along(unique_letters)]
    letter_map <- setNames(new_letters, unique_letters)
    # Step 3: Apply mapping
    group_df <- group_df %>%
      mutate(cld_group = map_chr(letter_list, ~ paste0(letter_map[.x], collapse = "")))
    return(group_df)
  }) %>% 
  #remove unnecesary grouping variables
  dplyr::select(-c(seasonlevel,vegpresencelevel)) %>% 
  #Format final emmean_CLD: select columns
  dplyr::select(casepilot, ghgspecies, comparison, status, season,vegpresence, weights_used,
                df, emmean_og, SE_og,lower.CL_og, upper.CL_og,
                cld_group,
                emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt) %>% 
  #Format final emmean_CLD: arrange values
  arrange(casepilot,ghgspecies, comparison,season,vegpresence, status)

#Save emmeans and groupletters (back-transformed):  
write.csv(x = complexmodel_customweight_emmeansCLD, file = paste0(path_2_modeloutputs,"EmmeansCLD_rest_chambermodels.csv"),row.names = F)


#_____________----


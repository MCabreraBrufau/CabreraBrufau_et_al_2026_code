#GHGpaper_modelchamberdata.R

#Date: September 2025


#Description----
#ALTERNATIVE APPROACH: 
#using separate and independent models for co2,ch4,gwp100 and gwp20, is not ideal: models on built on co2eq directly lose gas-specfic behaviours, and, due to dataset-specific transformations and different modeling quality, results derived from GWP100 and GWP20 models do not agree well with each other. 
#Instead of modeling co2eq, we will try to produce reliable paired estimates of co2 and ch4 from independent models (via paired bootstraping with resampling). Once we have the paired emmean distributions calculated, we can derive gwp100 emmeans distributions and gwp20 emmeans distributions, which can be tested and reflect the exact same behaviour modeledd by co2 and ch4 models, but with different relative importance for the different timescales.

#Model summaries: 1 for co2 1 for ch4, from full paired dataset 
#Emmeans CO2 and CH4: obtained via bootstraping, back-transformed
#Contrasts CO2 and CH4: obtained via back-transformed emmeans 
#p-values calculated in model scale (contrasts from direct model emmmeans)

#CO2eq(GWP100, GWP20)
#NO model for them (their results emerge from the separate CO2 and CH4 models)
#Emmeans: calculated from paired bootstraped emmeans distributions of ch4 and co2 (in real scale)
#Contrasts: calculated from real scale emmeans distributions
#P-values: calculated in real scale (cannot be in model scale) from contrasts

#steps: 
#1. filter dataset for paired co2 and ch4 observations (required)
#2. transform for normality (as usual)
#3. build models and derive effects and fit

#4. Bootstrap: 
#4.1. Create resampled datasets (some obs missing, others duplicated) that:
    #A. mantain subsite grouping (subiste A1 data cannot be mixed with subsite A2 when resampling)
    #B. result in valid datasets (where all factors from model are represented. Matrix build with levels of season*status*vegpresence must be filled), otherwise reject simulated dataset and try again. 

#4.2. For each simulated dataset, refit models (co2, ch4), calculate EMMEANs (weighted appropriately) for all comparisons of interest. 

#5. Calculate for co2 and for ch4: 
  #5.1. Summary of back-transformed emmeans for reporting in paper, (mean, SE, CI95)
  #5.2. Conrasts of model-scale emmeans, derive significance of comparisons for paper (empirical p-value )
  #5.3. Contrasts of back-transformed emmeans. Summarised for reporting in paper (mean, SE, CI95 ) 

#6. CO2eq (gwp100 and gwp20): 
  #6.1. Calculate emmean distributions of CO2eq from paired distributions of back-transformed co2 and ch4 emmeans. Using gwp100 and gwp20 for ch4. 
  #6.2. Create summary of emmeans CO2eq for paper (mean, SE, CI95)
  #6.3. Calculate contrasts and p-values, summary for paper (mean, SE, CI95, empirical p-value)


#Caveats/Limitations: 
  #we would be discarding some observations (non-paired obs of co2 or ch4)
  #Slightly more complex to explain
  #May be less powerfull for CO2eq comparisons

#BUT: 
  #statisitically sound approach, 
  #we mantain covariation/dependence of co2 and ch4 via paired bootstrap, 
  #we do not obscure gas-specific behaviours, we model each gas independently without forcing a common distribution, common significant effects or common directional response (all things that a single gwp model built on individually-summed co2eq fluxes would be doing). 
  #Restults for GWP100 and GWP20 will be consistent, just scaled by different GWP factors, but rely on the same estimates as CO2 and CH4 models. 




#This scrip is used to model the effect of restoration in each casepilot. Using
#net daily GHG exchange rates of appropriate incubations, data preparation is in
#1_Prepare_Data.R script.


#DECISIONS: 
  #Pooling of non-vegetated strata: due to seasonal variability and
  #site-specific differences, using the 3 original sampling strata classes is impossible
  #(would lead to to rank-deficient models for some combinations). Instead we
  #group them into vegpresence: vegetated or not-vegetated
  
#Modelling approach for DU,CA,VA,DA,CU: GLMMtmb with general formula dailyflux~status*season*vegpresence + (1|subsite)
  #Modelling approach for RI: reduced model without vegpresence (to avoid rank-deficient due to restored sites being 100% vegetated)

  #Contrast: Set contrast options to "contr. sum", so that we are not using any
  #one level of a factor as the reference level, but rather the tests will asses
  #whether there is an overall average effect of one factor across all levels of
  #the other factors. For example, is there an effect of status across all
  #levels of season?


#STEPS: 
#0. Import and format data: remove Nas, as.factor levels, 
#1. Transform data: Using BestNormalize package, find and apply for each casepilot-GHGspecies, the transformation that maximizes normality. 
#2. Optimize and select modelling structure (family-distribution) for each casepilot*GHGspecies combo (co2, ch4, gwp100 and gwp20).
#3. Calculate model outputs: 
  
  #Model summaries
  #Model residual-tests
  #Estimated marginal means (EMMs) for desired levels 
  #Post-hoc contrasts between EMMs 


rm(list = ls()) # clear workspace


# ---- Packages ----
#Installs (if needed) and loads required packages:
required_pkgs <- c("tidyverse",
                   "lubridate",
                   "zoo",
                   "ggpubr",
                   "rlang",
                   "ggtext",
                   "purrr",
                   "data.table",
                   "tools",
                   "hms",
                   "suncalc",
                   "bestNormalize",
                   "numDeriv",
                   "scales",
                   "glmmTMB",
                   "DHARMa",
                   "performance",
                   "rstatix",
                   "emmeans",
                   "progressr",
                   "future.apply",
                   "multcomp",
                   "multcompView"
                   )

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}



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

#Run up to model-fittng (included)


#0.Contrast options-----
#NOTES on contrasts: in R by default, contrast is "contr.treatment" which uses
#the first level of each factor as the "reference" level and all subsequent as
#"treatment" levels. This implies that, with default contrast type, when looking
#at the main effects of my model, what is shown for each main effect is whether
#there is a significant effect of 1 factor (eg. status) at the reference level
#of all other factors (i.e. season). What we want is to asses whether there is
#an overall average effect of status across all levels of season. For this
#purpose we need to set contrasts to "contr. sum", and always call Anova (model,
#type="III").

#Set contrasts to contr.sum (for the whole R session)
options(contrasts = c("contr.sum", "contr.poly"))


#0. Import and format--------

#Import data and format: 
data4models<- read.csv(paste0(path_1_paperdata,"ChamberData4paper.csv"))

#Format and filter:
data4models<- data4models %>% 
  filter(ghgspecies%in%c("co2","ch4")) %>% 
  #Remove NAs
  filter(!is.na(dailyflux)) %>% 
  mutate(season=factor(season, ordered = F),
         status=factor(status, ordered = F),
         strata=factor(strata, ordered = F)) %>% 
  #Pool non-vegetated strata:
  mutate(vegpresence=if_else(strata=="vegetated","Vegetated","Non-vegetated")) %>% 
  mutate(vegpresence=factor(vegpresence, ordered=F))
  

#Notes -----
#Evaluate model structure and assumptions (residuals), not significance of
#effects nor even explained variability (we expect very little variance
#explained by status in many cases).

#1st step is to decide the approach based on our data structure and distribution. 

#Fixed decisions: 
#Each case-pilot is modeled independently (different data-distributions and transformation needs)
#Need to include subsite as random effect (accounts for repeated samplings and site-specific intercepts)

#Approach: use Generalized Linear Mixed Models (glmmTMB)
#They can account for non-normal data (even after transformation) by using different distribution families


#STEPS in model optimization: 
#1. Decide transformation: that which maximizes normality of data
#2. Decide modelling family and structure (Gaussian vs T-family): Gaussian preferred unless residuals fail


#0. Custom Functions -------
#Here functions to access relevant results from models contained inside a named
#list. Used for ease of formatting and to avoid repetition. Will be used to summarise information.

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
#calculate various pseudo-R2 metrics (log-Likelyhood based) to get the
#"improvement in model fit relative to a null model (intercept + random
#effects), capturing the combined explanatory contribution of fixed effects".
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


#Function to extract model structure and all fixed effects significance
#(regardless of their naming, can differ in order or presence across the models
#in model_list)
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


#Function to Produce comparative histograms and QQplots of untransformed and best-Normalize transformed data. 
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

#Notes on pseudolog (~signed log) transformation:
#https://stratosida.github.io/regression-regrets/Pseudo_log_explainer.html#finding-a-parameter-that-best-achieves-normality
#Pseudolog tranformation needs tunning of its parameter to each data-set
#(similar to yeo-johnson). To incorporate this transformation into the
#BestNormalize, we need to create several functions.

#Pseudolog transformation allows and maintains NAs, allows and mantains negative and positive values 
#Range of potential sigmas are data-driven: based on 10th-90th percentile ranges. 

# 1. Constructor: fits pseudo-log transformation with optimal sigma
#best_pseudo_log with data-driven optimization of sigma using 10th and 90th percentile range
best_pseudo_log <- function(x, base = 10, sigma_range = NULL, standardize = TRUE) {
  tryCatch({
    # Clean input
    x_clean <- x[!is.na(x) & is.finite(x)]
    if (length(x_clean) < 3) stop("Not enough valid data points.")
    
    # Compute data-driven sigma_range if not provided (using 5th-95th percentile range)
    if (is.null(sigma_range)) {
      q <- quantile(abs(x_clean), probs = c(0.05, 0.95), na.rm = TRUE)
      lower <- max(q[1] / 10, .Machine$double.eps)  # Avoid zero or near-zero
      upper <- max(q[2], lower * 10)
      sigma_range <- 10^seq(log10(lower), log10(upper), length.out = 500)#500 potential sigmas
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
#We will use a completely objective decision for transformation: using
#BestNormalize function (with pseudolog transformation as a potential option),
#we will chose the transformation that results in the most-normal data
#distribution.

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


#Automatically apply the best transformation to each dataset within data4models
#using the bestNormalize objects stored in table_trans

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

#For RI, we cannot include vegpresence in model (would lead to deficient models,
#restored is 100%vegetated)USE ONLY BASIC FORMULA: dailyflux_trans~season*status
#+ (1|subsite)

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
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: Good residuals, keep


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                   data = ca_co2,
                   family = gaussian(),
                   control = glmmTMBControl(
                     optCtrl = list(iter.max = 5000, eval.max = 5000)
                   ),
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
#DECISION: BAD residuals


#T_family:
if(F){
  #skip, gaussian nonveg already good
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_co2,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: 
}

#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                      data = ca_co2,
                      family = t_family(),
                      control = glmmTMBControl(
                        optCtrl = list(iter.max = 5000, eval.max = 5000)
                      ),
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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_co2,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_co2,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = cu_co2,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_co2,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_co2,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: Better residuals, although they still fail


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_co2,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: GOOD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_co2,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_co2,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: DOES not converge


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = du_co2,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: good residuals



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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: GOOD residuals


#SKIP t-family option: gaussian is already good
if (F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_co2,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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

#DECISION: 
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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: BAD residuals (bad model)


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_co2,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_co2,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: Bad residuals (BAD model)


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = va_co2,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: good residuals

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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ca_ch4,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: GOOD residuals, keep


#SKIP t-family options, gaussians already good
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_ch4,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: 


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = ca_ch4,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: 
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does not significantly improve fit. We will use it
#regardless for consistency and to check for vegetation-transport effects

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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: Good residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_ch4,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: Good residuals


#SKIP, GAUSSIAN ALREADY GOOD
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_ch4,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = cu_ch4,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: 
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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_ch4,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_ch4,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION: Bad residuals

#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_ch4,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
anova(m1_gaus_nostrata,m3_t_nostrata)
anova(m2_gaus_vegpresence, m4_t_vegpresence)
anova(m3_t_nostrata, m4_t_vegpresence)
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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: GOOD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_ch4,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: GOOD residuals


#SKIP, gaussian already good: 
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_ch4,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION:


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = du_ch4,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION:

}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does not significantly improve fit (not worth the extra
#complexity), we will still use it for consistency and to evaluate vegetation
#effects.

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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: BAD residuals


#SKIP! cannot include vegpresence as fixed with status
if(F){
  #Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
  m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                                data = ri_ch4,
                                family = gaussian(),
                                control = glmmTMBControl(
                                  optCtrl = list(iter.max = 5000, eval.max = 5000)
                                ),
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
  #DECISION: 
}



  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = ri_ch4,
                          family = t_family(),
                          control = glmmTMBControl(
                            optCtrl = list(iter.max = 5000, eval.max = 5000)
                          ),
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
  #DECISION: GOOD residuals
  
  #SKIP, cannot include vegpresence as fixed effect
  if (F){
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ri_ch4,
                             family = t_family(),
                             control = glmmTMBControl(
                               optCtrl = list(iter.max = 5000, eval.max = 5000)
                             ),
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
  #DECISION: 
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
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: GOOD residuals 


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_ch4,
                              family = gaussian(),
                              control = glmmTMBControl(
                                optCtrl = list(iter.max = 5000, eval.max = 5000)
                              ),
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
#DECISION: Significant outliers detected 


#T_family:
if(F){
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_ch4,
                        family = t_family(),
                        control = glmmTMBControl(
                          optCtrl = list(iter.max = 5000, eval.max = 5000)
                        ),
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
#DECISION:
}

#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = va_ch4,
                           family = t_family(),
                           control = glmmTMBControl(
                             optCtrl = list(iter.max = 5000, eval.max = 5000)
                           ),
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
#DECISION: Good residuals

#Compare models 
anova(m2_gaus_vegpresence, m4_t_vegpresence)# T-family significantly better than gaussian
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_ch4<- m1_gaus_nostrata
va_complexmodel_ch4<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_ch4, 
                        model2 = va_complexmodel_ch4,
                        data=va_ch4,
                        color_var = "vegpresence")

ggsave(filename = "VA_ch4_Observed_VS_predicted.png", 
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
  "RI_ch4" = ri_simplemodel_ch4
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
  "VA_ch4" = va_complexmodel_ch4
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
#Effect of status will first be assessed via Anova (best_model,type=3): is there
#an effect averaging across all seasons? yes/no, is the effect dependent on
#season?, is the effect depenent on vegpresence?

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



#cleanup of Global environment: 
#we no longer need individual datasets, individual models, or summary tables.
rm(ca_co2,cu_co2,da_co2,du_co2,ri_co2,va_co2,
   ca_ch4,cu_ch4,da_ch4,du_ch4,ri_ch4,va_ch4,
   ca_simplemodel_co2,ca_simplemodel_ch4,
   cu_simplemodel_co2,cu_simplemodel_ch4,
   da_simplemodel_co2,da_simplemodel_ch4,
   du_simplemodel_co2,du_simplemodel_ch4,
   ri_simplemodel_co2,ri_simplemodel_ch4,
   va_simplemodel_co2,va_simplemodel_ch4,
   ca_complexmodel_co2,ca_complexmodel_ch4,
   cu_complexmodel_co2,cu_complexmodel_ch4,
   da_complexmodel_co2,da_complexmodel_ch4,
   du_complexmodel_co2,du_complexmodel_ch4,
   #No RI complexmodel to remove
   va_complexmodel_co2,va_complexmodel_ch4,
   m1_gaus_nostrata,m2_gaus_vegpresence,m3_t_nostrata,m4_t_vegpresence,res,
   complexmodel_fit, complexmodel_homoc_r2, complexmodel_resid_diag_summary, complexmodel_results_anova,
   simplemodel_fit, simplemodel_homoc_r2, simplemodel_resid_diag_summary, simplemodel_results_anova)


#Bootstrap EMMEANS--------
#We need to produce the paired bootstrap (allowing only valid datasets that contain all levels of factors used for modelling: all levels of status, season and vegpresence)


#1. set parallel plan ----
plan(multisession, workers = parallel::detectCores() - 1)
Nboot<- 5000 #Define number of bootstrap iterations

#2. Define functions-----

#Pair models: 
# Groups models like "CA_co2", "CA_ch4" → list(CA = list(co2, ch4))
pair_models <- function(model_list) {
  
  split_names <- sub("_(co2|ch4)$", "", names(model_list))
  paired <- split(model_list, split_names)
  
  # Force correct order: co2 first, ch4 second
  paired <- lapply(paired, function(x) {
    x[c(grep("co2", names(x)), grep("ch4", names(x)))]
  })
  
  return(paired)
}


#function for simplemodels emmeans (status, season, status_within_season, all with equal weights)
compute_emm_simple <- function(mod) {
  # Build reference grid
  rg <- ref_grid(mod)
  # Compute emmeans
  list(
    status = suppressMessages(emmeans(rg, ~ status, weights = "equal")),
    season = suppressMessages(emmeans(rg, ~ season, weights = "equal")),
    inseason_status = suppressMessages(emmeans(rg, ~ status | season, weights = "equal"))
  )
}


#function for complexmodels emmeans (status, season,status_within_season, status_within_vegpresence), using appropriate weights (from seasonal vegpresence distribution)
compute_emm_complex <- function(mod, weights_df) {
  
  # Build reference grid
  rg <- ref_grid(mod)
  rg_df <- as.data.frame(rg)
  
  # Join custom weights to the reference grid of the model
  rg_df <- rg_df %>%
    left_join(weights_df,
              by = c("status", "season", "vegpresence"))
  
  # Check missing weights (IMPORTANT)
  if (any(is.na(rg_df$prop_weight))) {
    stop("Missing weights after join")
  }
  
  # Attach weights
  rg@grid$w <- rg_df$prop_weight
  
  #2. Calculate emmeans
  
  # STATUS (equal across seasons for overall annual effect, but taking into account seasonally variable vegetation presence)
  emm_status <- suppressMessages(emmeans(rg, ~ status | season, weights = "proportional"))
  emm_status <- suppressMessages(emmeans(emm_status, ~ status, weights = "equal"))
  
  # SEASON (equal across status for overall seasonal effect, but taking into account status-variable vegetation presence)
  emm_season <- suppressMessages(emmeans(rg, ~ season | status, weights = "proportional"))
  emm_season <- suppressMessages(emmeans(emm_season, ~ season, weights = "equal"))
  
  # STATUS WITHIN SEASON (takes into account seasonally variable vegpresence)
  emm_ss <- suppressMessages(emmeans(rg, ~ status | season, weights = "proportional"))
  
  # STATUS WITHIN VEGPRESENCE (gives equal weight to all seasons to look at the overall annual effect)
  emm_sv <- suppressMessages(emmeans(rg, ~ status | vegpresence, weights = "equal"))
  
  return(list(
    status = emm_status,
    season = emm_season,
    inseason_status = emm_ss,
    inveg_status = emm_sv
  ))
}

#Function to extract only relevant parameters from emmeans objects: 
extract_emm_df <- function(emm_list, gas_label) {
  
  do.call(rbind, lapply(names(emm_list), function(term_name) {
    
    df <- as.data.frame(emm_list[[term_name]])
    
    df$term <- term_name
    df$ghgspecies <- gas_label
    
    df
  }))
}


#DEFINE function to bootstrap (parametric bootstraping paried for co2 and ch4)
one_boot <- function(df_case, model_pair, model_type, weights_df = NULL, max_tries = 50) {
  
  # helper: safe model fit capturing warnings
  safe_fit <- function(expr) {
    
    warn <- NULL
    
    result <- withCallingHandlers(
      tryCatch(expr, error = function(e) return(NULL)),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    list(model = result, warnings = warn)
  }
  
  for (i in 1:max_tries) {
    
    # --- 1. SPLIT ORIGINAL DATA (STRUCTURE PRESERVED) ---
    df_co2 <- df_case[df_case$ghgspecies == "co2", ]
    df_ch4 <- df_case[df_case$ghgspecies == "ch4", ]
    
    
    # --- 2. PARAMETRIC SIMULATION ---
    y_co2_sim <- tryCatch(simulate(model_pair[[1]])[[1]], error = function(e) NULL)
    y_ch4_sim <- tryCatch(simulate(model_pair[[2]])[[1]], error = function(e) NULL)
    
    if (is.null(y_co2_sim) || is.null(y_ch4_sim)) next
    
    
    # attach simulated responses (pairing preserved automatically)
    df_co2_sim <- df_co2
    df_ch4_sim <- df_ch4
    
    df_co2_sim$dailyflux_trans <- y_co2_sim
    df_ch4_sim$dailyflux_trans <- y_ch4_sim
    
    
    # --- 3. REFIT MODELS ---
    fit_co2 <- safe_fit(update(model_pair[[1]], data = df_co2_sim))
    fit_ch4 <- safe_fit(update(model_pair[[2]], data = df_ch4_sim))
    
    mod_co2 <- fit_co2$model
    mod_ch4 <- fit_ch4$model
    
    if (is.null(mod_co2) || is.null(mod_ch4)) next
    
    
    # --- 4. FILTER BAD FITS ---
    bad_warn_patterns <- c("rank-deficient", "convergence", "Hessian", "non-positive-definite")
    
    if (any(grepl(paste(bad_warn_patterns, collapse = "|"), fit_co2$warnings))) next
    if (any(grepl(paste(bad_warn_patterns, collapse = "|"), fit_ch4$warnings))) next
    
    
    # --- 5. COMPUTE EMMs ---
    emm_co2 <- tryCatch({
      if (model_type == "simple") {
        compute_emm_simple(mod_co2)
      } else {
        compute_emm_complex(mod_co2, weights_df)
      }
    }, error = function(e) return(NULL))
    
    emm_ch4 <- tryCatch({
      if (model_type == "simple") {
        compute_emm_simple(mod_ch4)
      } else {
        compute_emm_complex(mod_ch4, weights_df)
      }
    }, error = function(e) return(NULL))
    
    if (is.null(emm_co2) || is.null(emm_ch4)) next
    
    
    # --- 6. EXTRACT OUTPUT ---
    df_co2_emm <- extract_emm_df(emm_co2, "co2")
    df_ch4_emm <- extract_emm_df(emm_ch4, "ch4")
    
    df_out <- rbind(df_co2_emm, df_ch4_emm)
    
    df_out <- df_out %>%
      dplyr::select(any_of(c("ghgspecies","term","status","season","vegpresence","emmean")))
    
    
    # --- 7. RETURN SUCCESS ---
    return(df_out)
  }
  
  # --- FAILED ---
  return(NULL)
}


#Bootstrap one casepilot:
run_boot_case <- function(df, case, model_pair, model_type,
                          B = 1000, weights_df = NULL) {
  
  df_case <- df[df$casepilot == case, ]
  
  with_progress({
    
    p <- progressor(steps = B)
    
    results <- future_lapply(seq_len(B), function(b) {
      
      res <- tryCatch({
        
        out <- one_boot(
          df_case = df_case,
          model_pair = model_pair,
          model_type = model_type,
          weights_df = weights_df
        )
        
        # assign correct bootstrap iteration
        if (!is.null(out)) {
          out$iteration <- b
        }
        
        out
        
      }, error = function(e) {
        return(NULL)
      })
      
      p(sprintf("case %s - iter %d", case, b))
      
      return(res)
      
    }, 
    future.seed = TRUE,
    future.packages = c("glmmTMB", "emmeans", "dplyr")
    )
    
  })
  
  # Remove failed iterations
  results <- Filter(Negate(is.null), results)
  
  return(results)
}


#3.1. Simplemodels-----

#ACTUAL CALCULATION FOR SIMPLEMODELS
simple_pairs <- pair_models(simplemodel_list_allghg)

boot_simple <- list()

for (case in names(simple_pairs)) {
  
  cat("\nRunning SIMPLE bootstrap for case:", case, "\n")
  
  boot_simple[[case]] <- run_boot_case(
    df = data4models,
    case = case,
    model_pair = simple_pairs[[case]],
    model_type = "simple",
    B = Nboot #Bootstrap N
  )
}




#3.2. Complexmodels-----
#We need to take into account the proportion of vegpresence in the field to
#estimate the status, season and status_within_season emmeans. For
#status_within_vegpresence, equal weights should be applied (to give equal
#importance to each season).

#Overall status: weights should be used to account for different vegpresence in
#different status and  season, but all seasons should have the same combined
#weight for status emmeans and comparisons (to give the same importance to each
#season regardless of how small deviations in number of observations between
#seasons)

#RATIONALE FOR CUSTOM WEIGHTS:
#WE NEED TO APPLY weights to scale the vegpresence effects according to their
#proportion in the field. This applies for comparisions: status, season, and
#status_within_season. For status_within_vegpresence comparisons, equal weights
#should be applied (to give equal importance to each season). Weights should be
#calculated as the proportion of veg/noveg in every combination of status*season
#so they sum 1 (and give equal importance to all seasons). The function
#"custom_weighted_emmeans" already re-normalizes the weigths as needed (for
#status and for season emmeans, no need to re-normalize for status*season. It
#calculates the associated weighted SE taking into account the
#variance-covariance structure of the model.

#Status comparison: each status estimate is calculated taking into account its
#particular vegpresence composition while giving equal importance to all
#seasons.

#Season comparison: each seasonal estimate is calculated taking into account the
#seasonally variable vegpresence composition, while giving equal importance to
#all status (overall seasonal effect across all status).

#Status_within_season comparison: each seasonal estimate of each status takes
#into account its particular vegpresence composition.

#For status_within_vegpresence, we are only averaging across season and we want
#to have equal impact of all seasons (overall yearly effect), use
#weights="equal".


##Calculate weights------
#First, calculate vegpresence proportions at every level of casepilot, status and season 

#With weights calculated from chamber distribution, emmeans are proportional to
#the data used to fit the model, with the exception that different N across
#seasons does not influence the importance of each season in the status emmean
#(all season treated equally).
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
  
  #For CURONIAN: Calculate overall veg/noveg of curonian, across status and
  #seasons. No seasonal variability was observed in vegetation extent,
  #additionally due to initial incorrect restored site boundaries, the proportions
  #of vegpresence in these sites are not representative of actual site
  #composition. We use constant average proportion across all status and seasons.
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
  
  #Override CU composition: constant actual composition, differences in chamber
  #deployment are due to systematic sampling biass in Restored sites
  casepilot_weights<-casepilot_weights %>% 
    mutate(prop_weight=if_else(casepilot=="CU"&vegpresence=="Vegetated", cu_weights %>% pull(`Vegetated`),prop_weight),
           prop_weight=if_else(casepilot=="CU"&vegpresence=="Non-vegetated",cu_weights %>% pull(`Non-vegetated`),prop_weight))
}

#Build a list of weights to be called inside the loop (based on casepilot model name)
weights_list<- split(casepilot_weights, casepilot_weights$casepilot)



##Boot emmeans-----
complex_pairs <- pair_models(complexmodel_list_allghg)

boot_complex <- list()

for (case in names(complex_pairs)) {
  
  cat("\nRunning COMPLEX bootstrap for case:", case, "\n")
  
  boot_complex[[case]] <- run_boot_case(
    df = data4models,
    case = case,
    model_pair = complex_pairs[[case]],
    model_type = "complex",
    B = Nboot,
    weights_df = weights_list[[case]]
  )
}

# boot_simple and boot_complex are nested lists storing bootstrap results by casepilot. 
# The first level contains casepilot names (e.g., "RI", "CA", "DA"), and each case contains a list of bootstrap iterations. 
# Each iteration is a data frame where rows correspond to estimated marginal means (EMMs) for a given gas (co2/ch4), term (e.g., "status", "season", "inseason_status", and additionally "inveg_status" for complex models), and factor levels (status, season, and vegpresence when applicable). 
# The columns store the identifiers (ghgspecies, term, status, season, iteration) and the model-scale estimate (emmean). 
# Together, these lists represent the distribution of EMMs across bootstrap samples, preserving pairing between CO2 and CH4 within each iteration for downstream contrasts, back-transformation, and CO2eq calculations.



#4. CONTRASTS AND P-values-----

# For each casepilot:
  
# Model-scale contrasts + empirical p-values
# Back-transform EMMs
# Summarise EMMs (original scale)
# Compute contrasts on back-transformed scale + summaries

#to combine boot_simple and boot_complex lists into single data-frames on which to perform the calculations.
combine_boot <- function(boot_list) {
  
  map2_df(
    boot_list,
    names(boot_list),
    ~ bind_rows(.x) %>%
      mutate(casepilot = .y)
  )
}

boot_simple_df  <- combine_boot(boot_simple)
boot_complex_df <- combine_boot(boot_complex)


#1. CONTRASTS (MODEL SCALE + empirical p-values)

#helper: pairwise contrasts (optionally grouped, value can be direct "emmean" or "bt_emmean" for convenience)
compute_pairwise <- function(df, value_col = "emmean",
                             group_vars = NULL,
                             level_var = "status") {
  
  df %>%
    group_by(across(all_of(c("casepilot", "iteration", "ghgspecies", group_vars)))) %>%
    summarise(
      contrast_df = list({
        
        levs <- unique(.data[[level_var]])
        vals <- .data[[value_col]]
        
        combs <- combn(levs, 2, simplify = FALSE)
        
        do.call(rbind, lapply(combs, function(cmb) {
          
          i1 <- which(levs == cmb[1])
          i2 <- which(levs == cmb[2])
          
          data.frame(
            contrast = paste(cmb[1], "-", cmb[2]),
            diff = vals[i1] - vals[i2]
          )
        }))
        
      }),
      .groups = "drop"
    ) %>%
    tidyr::unnest(contrast_df)
}

#helper: empirical p-value
p_empirical <- function(x) {
  2 * min(mean(x >= 0), mean(x <= 0))
}

#summarise contrasts (optional pvalue)
summarise_contrasts <- function(df, compute_p = TRUE) {
  
  # Step 1: summarise per contrast
  df_sum <- df %>%
    group_by(across(any_of(c("casepilot", "ghgspecies", "contrast",
                             "season", "vegpresence")))) %>%
    summarise(
      mean = mean(diff, na.rm = TRUE),
      se   = sd(diff, na.rm = TRUE),
      lwr  = quantile(diff, 0.025, na.rm = TRUE),
      upr  = quantile(diff, 0.975, na.rm = TRUE),
      n_boot = sum(!is.na(diff)),
      p    = if (compute_p) p_empirical(diff) else NA_real_,
      .groups = "drop"
    )
  
  # Step 2: apply Holm correction across ALL contrasts
  # within each casepilot × ghgspecies combination
  df_sum <- df_sum %>%
    group_by(across(any_of(c("casepilot", "ghgspecies")))) %>%
    mutate(
      p_adj = if (compute_p) p.adjust(p, method = "holm") else NA_real_
    ) %>%
    ungroup()
  
  return(df_sum)
}


#Function to summarise emmean distribution (or emmean_bt), 
summarise_emm <- function(df, value_col = "emmean") {

df %>%
  group_by(across(any_of(c("casepilot", "ghgspecies", "term",
                           "status", "season", "vegpresence")))) %>%
  summarise(
    mean = mean(.data[[value_col]], na.rm = TRUE),
    se   = sd(.data[[value_col]], na.rm = TRUE),
    lwr  = quantile(.data[[value_col]], 0.025, na.rm = TRUE),
    upr  = quantile(.data[[value_col]], 0.975, na.rm = TRUE),
    n_boot = sum(!is.na(.data[[value_col]])),
    .groups = "drop"
  )
}

#Function to back transform emmeans distribution objects: boot_complex_df and boot_simple_df
backtransform_emmeans <- function(x, bn_obj) {
  predict(bn_obj, newdata = x, inverse = TRUE)
}

#4.1. Simplemodels -----

##A) model-scale cotrasts&pvalue-----
# STATUS
contr_status_RI_modscale <- boot_simple_df %>%
  filter(term == "status") %>%
  compute_pairwise(value_col = "emmean",group_vars = NULL, level_var = "status") %>% 
  summarise_contrasts(compute_p = T)

# SEASON
contr_season_RI_modscale <- boot_simple_df %>%
  filter(term == "season") %>%
  compute_pairwise(value_col = "emmean",group_vars = NULL, level_var = "season") %>% 
  summarise_contrasts(compute_p = T)

# STATUS WITHIN SEASON
contr_ss_RI_modscale <- boot_simple_df %>%
  filter(term == "inseason_status") %>%
  compute_pairwise(value_col = "emmean",group_vars = "season", level_var = "status") %>% 
  summarise_contrasts(compute_p = T)

##B) Bt_emmeans and Bt_contrasts------

#Back-transform emmeans
bt_boot_simple_df <- boot_simple_df %>%
  rowwise() %>%
  mutate(
    emmean_bt = {
      key <- paste0(unique(casepilot), "_", ghgspecies)
      backtransform_emmeans(emmean, bn_list[[key]])
    }
  ) %>%
  ungroup()

#Save back-transformed emmeans distributions (just in case): 
write.csv(bt_boot_simple_df,  file = paste0(path_2_modeloutputs,"bt_boot_simple_df.csv"),  row.names = FALSE)

#Back-transformed emmeans summary:
bt_sum_emmean_RI<- summarise_emm(bt_boot_simple_df,value_col = "emmean_bt")


#Calculate bt_scale contrasts (for tables in paper), add pvalue from model-scale test
# STATUS
bt_contr_status_RI<- bt_boot_simple_df %>%
  filter(term == "status") %>%
  compute_pairwise(value_col = "emmean_bt",group_vars = NULL, level_var = "status") %>% 
  summarise_contrasts(compute_p = F)  %>% dplyr::select(-c(p,p_adj)) %>% 
  left_join(contr_status_RI_modscale %>% dplyr::select(-c(mean, se, lwr,upr,n_boot))) %>% 
  mutate(comparison="status")

# SEASON
bt_contr_season_RI <- bt_boot_simple_df %>%
  filter(term == "season") %>%
  compute_pairwise(value_col = "emmean_bt",group_vars = NULL, level_var = "season") %>% 
  summarise_contrasts(compute_p = F) %>% dplyr::select(-c(p,p_adj)) %>% 
  left_join(contr_season_RI_modscale %>% dplyr::select(-c(mean, se, lwr,upr,n_boot)))%>% 
  mutate(comparison="season")

# STATUS WITHIN SEASON
bt_contr_ss_RI <- bt_boot_simple_df %>%
  filter(term == "inseason_status") %>%
  compute_pairwise(value_col = "emmean",group_vars = "season", level_var = "status") %>% 
  summarise_contrasts(compute_p = F) %>% dplyr::select(-c(p,p_adj)) %>% 
  left_join(contr_ss_RI_modscale %>% dplyr::select(-c(mean, se, lwr,upr,n_boot)))%>% 
  mutate(comparison="inseason_status")

#Format and combine contrasts: 
bt_contr_RI<- bind_rows(bt_contr_status_RI,bt_contr_season_RI,bt_contr_ss_RI) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, contrast, season, mean, se, lwr, upr, n_boot, p, p_adj)






#4.2. Complexmodels -----

##A) model-scale cotrasts&pvalue-----
# STATUS
contr_status_rest_modscale <- boot_complex_df %>%
  filter(term == "status") %>%
  compute_pairwise(value_col = "emmean",group_vars = NULL, level_var = "status") %>% 
  summarise_contrasts(compute_p = T)

# SEASON
contr_season_rest_modscale <- boot_complex_df %>%
  filter(term == "season") %>%
  compute_pairwise(value_col = "emmean",group_vars = NULL, level_var = "season") %>% 
  summarise_contrasts(compute_p = T)

# STATUS WITHIN SEASON
contr_ss_rest_modscale <- boot_complex_df %>%
  filter(term == "inseason_status") %>%
  compute_pairwise(value_col = "emmean",group_vars = "season", level_var = "status") %>% 
  summarise_contrasts(compute_p = T)

# STATUS WITHIN VEGPRESENCE
contr_sv_rest_modscale <- boot_complex_df %>% 
  filter(term == "inveg_status") %>%
  compute_pairwise(value_col = "emmean",group_vars = "vegpresence", level_var = "status") %>% 
  summarise_contrasts(compute_p = T)


##B) Bt_emmeans and Bt_contrasts------

#Back-transform emmeans
bt_boot_complex_df <- boot_complex_df %>%
  rowwise() %>%
  mutate(
    emmean_bt = {
      key <- paste0(unique(casepilot), "_", ghgspecies)
      backtransform_emmeans(emmean, bn_list[[key]])
    }
  ) %>%
  ungroup()

#Save back-transformed emmeans distributions (just in case): 
write.csv(bt_boot_complex_df,  file = paste0(path_2_modeloutputs,"bt_boot_complex_df.csv"),  row.names = FALSE)


#Back-transformed emmeans summary:
bt_sum_emmean_rest<- summarise_emm(bt_boot_complex_df,value_col = "emmean_bt")


#Calculate bt_scale contrasts (for tables in paper), add pvalue from model-scale test
# STATUS
bt_contr_status_rest<- bt_boot_complex_df %>%
  filter(term == "status") %>%
  compute_pairwise(value_col = "emmean_bt",group_vars = NULL, level_var = "status") %>% 
  summarise_contrasts(compute_p = F)  %>% dplyr::select(-c(p,p_adj)) %>% 
  left_join(contr_status_rest_modscale %>% dplyr::select(-c(mean, se, lwr,upr,n_boot))) %>% 
  mutate(comparison="status")

# SEASON
bt_contr_season_rest <- bt_boot_complex_df %>%
  filter(term == "season") %>%
  compute_pairwise(value_col = "emmean_bt",group_vars = NULL, level_var = "season") %>% 
  summarise_contrasts(compute_p = F) %>% dplyr::select(-c(p,p_adj)) %>% 
  left_join(contr_season_rest_modscale %>% dplyr::select(-c(mean, se, lwr,upr,n_boot)))%>% 
  mutate(comparison="season")

# STATUS WITHIN SEASON
bt_contr_ss_rest <- bt_boot_complex_df %>%
  filter(term == "inseason_status") %>%
  compute_pairwise(value_col = "emmean",group_vars = "season", level_var = "status") %>% 
  summarise_contrasts(compute_p = F) %>% dplyr::select(-c(p,p_adj)) %>% 
  left_join(contr_ss_rest_modscale %>% dplyr::select(-c(mean, se, lwr,upr,n_boot)))%>% 
  mutate(comparison="inseason_status")

# STATUS WITHIN VEGPRESENCE
bt_contr_sv_rest <- bt_boot_complex_df %>%
  filter(term == "inveg_status") %>%
  compute_pairwise(value_col = "emmean",group_vars = "vegpresence", level_var = "status") %>% 
  summarise_contrasts(compute_p = F) %>% dplyr::select(-c(p,p_adj)) %>% 
  left_join(contr_sv_rest_modscale %>% dplyr::select(-c(mean, se, lwr,upr,n_boot)))%>% 
  mutate(comparison="inveg_status")

#Format and combine contrasts: 
bt_contr_rest<- bind_rows(bt_contr_status_rest,bt_contr_season_rest,bt_contr_ss_rest,bt_contr_sv_rest) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, contrast, season,vegpresence, mean, se, lwr, upr,n_boot,p,p_adj)



#4.3. CO2eq calc-----
#We now calulate the CO2eq metric (gwp100, gwp20) from the paired bootstrapped emmeans.
#we need: 
  #1. Emmeans summary
  #2. Contrast summary (with pvalue, all from real-scale)

#GWP factors for ch4:
gwp100factor<- 27
gwp20factor<- 79.7

#Co2 is in mol per m2 per day: use (co2*44.0095)  to transform to g per m2 per day
#ch4 is in mmol per m2 per day: use (ch4*1e-3*16.04246) to transform to g per m2 per day


##A) BUILD CO2eq DATASET (from bt scale) ----
#FOR simplemodels:

# reshape to wide to combine gases
bt_simple_wide <- bt_boot_simple_df %>%
  dplyr::select(casepilot, iteration, term, status, season, ghgspecies, emmean_bt) %>%
  tidyr::pivot_wider(names_from = ghgspecies, values_from = emmean_bt)

# compute CO2eq
gwp_simple_df <- bt_simple_wide %>%
  mutate(
    gwp100 = (co2*44.0095) + ((ch4*1e-3*16.04246) * gwp100factor),
    gwp20  =  (co2*44.0095) + ((ch4*1e-3*16.04246) * gwp20factor)
  ) %>%
  tidyr::pivot_longer(
    cols = c(gwp100, gwp20),
    names_to = "ghgspecies",
    values_to = "emmean_bt"
  )

#Save back-transformed emmeans distributions (just in case): 
write.csv(gwp_simple_df,  file = paste0(path_2_modeloutputs,"gwp_simple_df.csv"),  row.names = FALSE)


#FOR complexmodels:
bt_complex_wide <- bt_boot_complex_df %>%
  dplyr::select(casepilot, iteration, term, status, season, vegpresence, ghgspecies, emmean_bt) %>%
  tidyr::pivot_wider(names_from = ghgspecies, values_from = emmean_bt)

gwp_complex_df <- bt_complex_wide %>%
  mutate(
    gwp100 = (co2*44.0095) + ((ch4*1e-3*16.04246) * gwp100factor),
    gwp20  =  (co2*44.0095) + ((ch4*1e-3*16.04246) * gwp20factor)
  ) %>%
  tidyr::pivot_longer(
    cols = c(gwp100, gwp20),
    names_to = "ghgspecies",
    values_to = "emmean_bt"
  )


#Save back-transformed emmeans distributions (just in case): 
write.csv(gwp_complex_df,  file = paste0(path_2_modeloutputs,"gwp_complex_df.csv"),  row.names = FALSE)

##B) EMMEAN ------
gwp_sum_emmean_RI <- summarise_emm(gwp_simple_df, value_col = "emmean_bt")
gwp_sum_emmean_rest <- summarise_emm(gwp_complex_df, value_col = "emmean_bt")

##C) CONTRASTS-----
#Contrasts and pvalues for co2eq are calculated directly in the real scale (no model-scale for these metrics)

#STATUS simple
gwp_contr_status_simple <- gwp_simple_df %>%
  filter(term == "status") %>%
  compute_pairwise(value_col = "emmean_bt", group_vars = NULL, level_var = "status") %>%
  summarise_contrasts(compute_p = TRUE) %>%
  mutate(comparison = "status")

#SEASON simple 
gwp_contr_season_simple <- gwp_simple_df %>%
  filter(term == "season") %>%
  compute_pairwise(value_col = "emmean_bt", group_vars = NULL, level_var = "season") %>%
  summarise_contrasts(compute_p = TRUE) %>%
  mutate(comparison = "season")

#STATUS WITHIN SEASON simple
gwp_contr_ss_simple <- gwp_simple_df %>%
  filter(term == "inseason_status") %>%
  compute_pairwise(value_col = "emmean_bt", group_vars = "season", level_var = "status") %>%
  summarise_contrasts(compute_p = TRUE) %>%
  mutate(comparison = "inseason_status")

#Combine Simple contrasts: 
gwp_contr_RI <- bind_rows(
  gwp_contr_status_simple,
  gwp_contr_season_simple,
  gwp_contr_ss_simple) %>%
  dplyr::select(casepilot, ghgspecies, comparison, contrast, season, mean, se, lwr, upr,n_boot, p, p_adj)


#STATUS complex
gwp_contr_status_complex <- gwp_complex_df %>%
  filter(term == "status") %>%
  compute_pairwise(value_col = "emmean_bt", group_vars = NULL, level_var = "status") %>%
  summarise_contrasts(compute_p = TRUE) %>%
  mutate(comparison = "status")

#SEASON complex
gwp_contr_season_complex <- gwp_complex_df %>%
  filter(term == "season") %>%
  compute_pairwise(value_col = "emmean_bt", group_vars = NULL, level_var = "season") %>%
  summarise_contrasts(compute_p = TRUE) %>%
  mutate(comparison = "season")

#STATUS WITHIN SEASON complex
gwp_contr_ss_complex <- gwp_complex_df %>%
  filter(term == "inseason_status") %>%
  compute_pairwise(value_col = "emmean_bt", group_vars = "season", level_var = "status") %>%
  summarise_contrasts(compute_p = TRUE) %>%
  mutate(comparison = "inseason_status")

#STATUS WITHIN VEGPRESENCE complex
gwp_contr_sv_complex <- gwp_complex_df %>%
  filter(term == "inveg_status") %>%
  compute_pairwise(value_col = "emmean_bt", group_vars = "vegpresence", level_var = "status") %>%
  summarise_contrasts(compute_p = TRUE) %>%
  mutate(comparison = "inveg_status")

#COMBINE complex:
gwp_contr_rest <- bind_rows(
  gwp_contr_status_complex,
  gwp_contr_season_complex,
  gwp_contr_ss_complex,
  gwp_contr_sv_complex) %>%
  dplyr::select(casepilot, ghgspecies, comparison, contrast, season, vegpresence, mean, se, lwr, upr,n_boot, p, p_adj)


#5. Format and save------

#COMBINE OUTPUTS FROM DIFFERENT GASSES. 
all_contr_RI <- bind_rows(gwp_contr_RI, bt_contr_RI)
all_contr_rest <- bind_rows(gwp_contr_rest, bt_contr_rest)

all_emmean_RI <- bind_rows(bt_sum_emmean_RI, gwp_sum_emmean_RI)
all_emmean_rest <- bind_rows(bt_sum_emmean_rest, gwp_sum_emmean_rest)


##5.1. Simplemodels------

#Format contrasts RI: 
posthoctests_RI_chambermodels<- all_contr_RI %>% 
  rename(estimate_bt=mean, SE_bt=se,lower.CL_bt=lwr, upper.CL_bt=upr, boot_n= n_boot, p.value= p_adj) %>% 
  mutate(comparison=case_when(comparison=="inseason_status"~"status_within_season",
                              comparison=="inveg_status"~"status_within_vegpresence",
                              comparison=="status"~"status",
                              comparison=="season"~"season",
                              TRUE~NA)) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, contrast, season, 
                estimate_bt, SE_bt, lower.CL_bt, upper.CL_bt, boot_n, p.value) %>% 
  mutate(
    # Desired ordering
    ghgspecies = factor(ghgspecies, levels = c("co2", "ch4", "gwp100", "gwp20")),
    comparison = factor(comparison, levels = c("status","season","status_within_season")),
    season = factor(season, levels = c("S1", "S2", "S3", "S4"))
  ) %>%
  arrange(casepilot,ghgspecies,comparison,season)
  


#Obtain CLDs (letter-groups) based on post-hoc tests for every comparison:
#we need to substitute the NAs in status and season with "." to allow valid join afterwards. 
# Obtain CLDs (letter-groups)
simpleCLD_letters <- posthoctests_RI_chambermodels %>%
  # Keep only needed columns
  dplyr::select(casepilot, ghgspecies, comparison, season, contrast, p.value) %>%
  # Clean contrast + replace NA with placeholder
  mutate(
    contrast = gsub(" ", "", contrast),
    #Add seasonlevel for status_within_season comparison
    seasonlevel=if_else(comparison=="status_within_season", season, ".")
  ) %>%
  # Group by comparison context
  group_by(casepilot, ghgspecies, comparison, seasonlevel) %>%
  nest() %>%
  # Apply multcompLetters per group
  mutate(
    letters = purrr::map(data, function(group_df) {
      contrast_matrix <- group_df %>%
        dplyr::select(contrast, p.value) %>%
        deframe()
      multcompLetters(contrast_matrix)$Letters %>%
        enframe(name = "level", value = "cld_group")
    })
  ) %>%
  # Unnest
  dplyr::select(-data) %>%
  unnest(letters) %>%
  # Recover status and season
  mutate(
    status = stringr::str_extract(level, "Altered|Preserved|Restored"),
    season = if_else(
      comparison == "season",
      level,      # season levels are in "level"
      seasonlevel # otherwise comes from grouping
    ),
    # Replace NA status with placeholder (for join compatibility)
    status = ifelse(is.na(status), ".", status)
  ) %>%
  # Final columns (JOIN-READY)
  dplyr::select(casepilot, ghgspecies, comparison, status, season, seasonlevel,
                cld_group)


#Format EMMs
EMMs_RI_formated<- all_emmean_RI %>% 
  rename(emmean_bt=mean, SE_bt=se,lower.CL_bt=lwr,upper.CL_bt=upr,boot_n=n_boot) %>% 
  mutate(comparison=case_when(term=="inseason_status"~"status_within_season",
                              term=="inveg_status"~"status_within_vegpresence",
                              term=="status"~"status",
                              term=="season"~"season",
                              TRUE~NA),
         seasonlevel=if_else(comparison=="status_within_season", season, ".")) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, status, season, seasonlevel,
                boot_n, emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt)


#Add CLD to formated emmeans
EMMsCLD_RI_formated <- EMMs_RI_formated %>%
  #ADD CLD group-letters (re-coding them so that "a" always identifies the
  #significant group with lowest emmean, b the next lowest emmean, and so on...)
  left_join(
    simpleCLD_letters,
    by = c("casepilot","ghgspecies","comparison","status","season","seasonlevel")
  ) %>%
  #Separate each comparison group to do the letter re-coding
  group_split(casepilot, ghgspecies, comparison, seasonlevel) %>%
  map_dfr(function(group_df) {
    # Step 1: Order by emmean
    group_df <- group_df %>% arrange(emmean_bt)
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
  #Alphabetical order of mixed "ba" "ab" cld_groups
  mutate(cld_group=sapply(strsplit(cld_group, ""), function(x) paste(sort(x), collapse = ""))) %>% 
  #Format final emmean_CLD: select columns
  dplyr::select(casepilot, ghgspecies, comparison, status, season, boot_n, emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt, cld_group) %>% 
  #Format final emmean_CLD: arrange values
  mutate(
    # Desired ordering
    ghgspecies = factor(ghgspecies, levels = c("co2", "ch4", "gwp100", "gwp20")),
    comparison = factor(comparison, levels = c("status","season","status_within_season")),
    season = factor(season, levels = c("S1", "S2", "S3", "S4")),
    status = factor(status, levels = c("Altered", "Preserved", "Restored"))
  ) %>%
  arrange(casepilot,ghgspecies,comparison,season,status)


#Save Final outputs: 
EMMsCLD_RI_formated
#save as EmmeansCLD_RI_chambermodels_boot.csv
write.csv(x = EMMsCLD_RI_formated, 
          file = paste0(path_2_modeloutputs,"EmmeansCLD_RI_chambermodels_boot.csv"),
          row.names = F)


posthoctests_RI_chambermodels
#save as Posthoctests_RI_chambermodels_boot.csv
write.csv(x = posthoctests_RI_chambermodels, 
          file = paste0(path_2_modeloutputs,"Posthoctests_RI_chambermodels_boot.csv"),
          row.names = F)






##5.2. Complexmodels----


#Format contrasts REST: 
posthoctests_rest_chambermodels<- all_contr_rest %>% 
  rename(estimate_bt=mean, SE_bt=se,lower.CL_bt=lwr, upper.CL_bt=upr, boot_n= n_boot, p.value= p_adj) %>% 
  mutate(comparison=case_when(comparison=="inseason_status"~"status_within_season",
                              comparison=="inveg_status"~"status_within_vegpresence",
                              comparison=="status"~"status",
                              comparison=="season"~"season",
                              TRUE~NA)) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, contrast, season, vegpresence, 
                estimate_bt, SE_bt, lower.CL_bt, upper.CL_bt, boot_n, p.value) %>% 
  # Desired ordering
  mutate(
    ghgspecies = factor(ghgspecies, levels = c("co2", "ch4", "gwp100", "gwp20")),
    comparison = factor(comparison, levels = c("status","season","status_within_season","status_within_vegpresence")),
    season = factor(season, levels = c("S1", "S2", "S3", "S4")),
    vegpresence=factor(vegpresence, levels = c("Non-vegetated","Vegetated"))
  ) %>%
  arrange(casepilot,ghgspecies,comparison,season,vegpresence)


#Obtain CLDs (letter-groups) based on post-hoc tests for every comparison:
#we need to substitute the NAs in status and season with "." to allow valid join afterwards. 
complexCLD_letters <- posthoctests_rest_chambermodels %>%
  #leave only (1) variables that identify unique comparisons, (2)contrast column, (3) pvalue column
  #seasonlevel and vegpresencelevel is kept to account for status_within_season comparison. 
  dplyr::select(casepilot, ghgspecies, comparison, season,vegpresence, contrast, p.value) %>% 
  #Remove any spaces from contrast column (mutcompLetters expects levels to be only separated by a hyphen "-")
  mutate(contrast=gsub(" ","", contrast),
         #Add seasonlevel for status_within_season comparison
         seasonlevel=if_else(comparison=="status_within_season", season, "."),
         #Add vegpresenelevel for status_within_vegpresence comparison
         vegpresencelevel=if_else(comparison=="status_within_vegpresence", vegpresence, ".")) %>% 
  #Group by comparison identifiers
  group_by(casepilot, ghgspecies, comparison, seasonlevel, vegpresencelevel) %>%
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
  #Reformat level to fit the appropriate columns: season, status, vegpresence
  mutate(
    status = stringr::str_extract(level, "Altered|Preserved|Restored"),
    season = if_else(comparison=="season", level, seasonlevel),
    vegpresence = vegpresencelevel,
    # Replace NA status with placeholder (for join compatibility)
    status = ifelse(is.na(status), ".", status)
    ) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, status, season, vegpresence, seasonlevel, vegpresencelevel, cld_group)


#Format EMMs
EMMs_rest_formated<- all_emmean_rest %>% 
  rename(emmean_bt=mean, SE_bt=se,lower.CL_bt=lwr,upper.CL_bt=upr,boot_n=n_boot) %>% 
  mutate(comparison=case_when(term=="inseason_status"~"status_within_season",
                              term=="inveg_status"~"status_within_vegpresence",
                              term=="status"~"status",
                              term=="season"~"season",
                              TRUE~NA)) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, status, season,vegpresence, boot_n, emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt)



#Add CLD to formated emmeans
EMMsCLD_RI_formated <- EMMs_RI_formated %>%
  #ADD CLD group-letters (re-coding them so that "a" always identifies the
  #significant group with lowest emmean, b the next lowest emmean, and so on...)
  left_join(
    simpleCLD_letters,
    by = c("casepilot","ghgspecies","comparison","status","season","seasonlevel")
  ) %>%
  #Separate each comparison group to do the letter re-coding
  group_split(casepilot, ghgspecies, comparison, seasonlevel) %>%
  map_dfr(function(group_df) {
    # Step 1: Order by emmean
    group_df <- group_df %>% arrange(emmean_bt)
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
  #Alphabetical order of mixed "ba" "ab" cld_groups
  mutate(cld_group=sapply(strsplit(cld_group, ""), function(x) paste(sort(x), collapse = ""))) %>% 
  #Format final emmean_CLD: select columns
  dplyr::select(casepilot, ghgspecies, comparison, status, season, boot_n, emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt, cld_group) %>% 
  #Format final emmean_CLD: arrange values
  mutate(
    # Desired ordering
    ghgspecies = factor(ghgspecies, levels = c("co2", "ch4", "gwp100", "gwp20")),
    comparison = factor(comparison, levels = c("status","season","status_within_season","status_within_vegpresence")),
    season = factor(season, levels = c("S1", "S2", "S3", "S4")),
    status = factor(status, levels = c("Altered", "Preserved", "Restored"))
  ) %>%
  arrange(casepilot,ghgspecies,comparison,season,status)


#Add CLD to formated emmeans
EMMsCLD_rest_formated <- EMMs_rest_formated %>%
  #ADD CLD group-letters (re-coding them so that "a" always identifies the
  #significant group with lowest emmean, b the next lowest emmean, and so on...)
  left_join(complexCLD_letters, by = c("casepilot","ghgspecies","comparison","status","season","vegpresence")) %>%
  #Separate each comparison group to do the letter re-coding
  group_split(casepilot, ghgspecies, comparison, seasonlevel,vegpresencelevel) %>%
  map_dfr(function(group_df) {
    # Step 1: Order by emmean
    group_df <- group_df %>% arrange(emmean_bt)
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
  #Alphabetical order of mixed "ba" "ab" cld_groups
  mutate(cld_group=sapply(strsplit(cld_group, ""), function(x) paste(sort(x), collapse = ""))) %>% 
  #Format final emmean_CLD: select columns
  dplyr::select(casepilot, ghgspecies, comparison, status, season, vegpresence, boot_n, emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt, cld_group) %>% 
  #Format final emmean_CLD: arrange values
  mutate(
    # Desired ordering
    ghgspecies = factor(ghgspecies, levels = c("co2", "ch4", "gwp100", "gwp20")),
    comparison = factor(comparison, levels = c("status","season","status_within_season","status_within_vegpresence")),
    season = factor(season, levels = c("S1", "S2", "S3", "S4")),
    status = factor(status, levels = c("Altered", "Preserved", "Restored")),
    vegpresence = factor(vegpresence, levels=c("Non-vegetated","Vegetated"))
  ) %>%
  arrange(casepilot,ghgspecies,comparison,season,vegpresence,status)




#Save Final outputs: 
EMMsCLD_rest_formated
#save as EmmeansCLD_rest_chambermodels_boot.csv
write.csv(x = EMMsCLD_rest_formated, 
          file = paste0(path_2_modeloutputs,"EmmeansCLD_rest_chambermodels_boot.csv"),
          row.names = F)


posthoctests_rest_chambermodels
#save as Posthoctests_rest_chambermodels_boot.csv
write.csv(x = posthoctests_rest_chambermodels, 
          file = paste0(path_2_modeloutputs,"Posthoctests_rest_chambermodels_boot.csv"),
          row.names = F)


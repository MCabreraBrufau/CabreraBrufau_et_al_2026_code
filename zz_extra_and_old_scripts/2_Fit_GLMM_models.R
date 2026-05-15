#GHGpaper_modelchamberdata.R

#Date: September 2025


#Description----
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
{

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
  filter(ghgspecies%in%c("co2","ch4","gwp100", "gwp20")) %>% 
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
#DECISION: does not converge. 


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
#DECISION: Bad residuals


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
#DECISION: BAD residuals


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
#DECISION: Better residuals, although they still fail


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
#DECISION: Bad residuals


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
#DECISION: BAD residuals



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
#DECISION: Better residuals, although they still fail


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
#DECISION: GOOD residuals


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
#DECISION: BAD residuals


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
#DECISION: Model does not converge


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
#DECISION: Model does not converge


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
#DECISION: Good-enough residuals


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
#DECISION: Good residuals


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
#DECISION:


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
#DECISION: GOOD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_ch4,
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
#DECISION: GOOD residuals


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
#DECISION:


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
#DECISION:


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
#DECISION:
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


#GWP100______-------

##2.1. CA_gwp100 (ok) -------

#Subset data and check transformation: 
ca_gwp100<- data4models %>% filter(casepilot=="CA"&ghgspecies=="gwp100")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="gwp100") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ca_gwp100,
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
#DECISION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ca_gwp100,
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
#DECISION: BAD residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = ca_gwp100,
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
  #DECISION: Does not converge
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ca_gwp100,
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
  #DECISION: Good  residuals
  

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence does significantly improve results

#Save best models (simple and complex best)
ca_simplemodel_gwp100<- m1_gaus_nostrata
ca_complexmodel_gwp100<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_gwp100, 
                        model2 = ca_complexmodel_gwp100,
                        data=ca_gwp100,
                        color_var = "vegpresence")

ggsave(filename = "CA_gwp100_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_gwp100 (ok) -------

#Subset data and check transformation: 
cu_gwp100<- data4models %>% filter(casepilot=="CU"&ghgspecies=="gwp100")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="gwp100") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_gwp100,
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
#DECISION: BAD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_gwp100,
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
#DECISION: BAD residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = cu_gwp100,
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
  #DECISION: Good enough residuals
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = cu_gwp100,
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
  #DECISION: GOOD residuals

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_gwp100<- m3_t_nostrata
cu_complexmodel_gwp100<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_gwp100, 
                        model2 = cu_complexmodel_gwp100,
                        data=cu_gwp100,
                        color_var = "vegpresence")

ggsave(filename = "CU_gwp100_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.3. DA_gwp100 (ok) -------

#Subset data and check transformation: 
da_gwp100<- data4models %>% filter(casepilot=="DA"&ghgspecies=="gwp100")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="gwp100") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_gwp100,
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
#DECISION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_gwp100,
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
#DECISION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_gwp100,
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
#DECISION: Better residuals, although they still have issues


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_gwp100,
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
da_simplemodel_gwp100<- m3_t_nostrata
da_complexmodel_gwp100<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_gwp100, 
                        model2 = da_complexmodel_gwp100,
                        data=da_gwp100,
                        color_var = "vegpresence")

ggsave(filename = "DA_gwp100_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_gwp100 (ok) -------

#Subset data and check transformation: 
du_gwp100<- data4models %>% filter(casepilot=="DU"&ghgspecies=="gwp100")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="gwp100") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_gwp100,
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
#DECISION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_gwp100,
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
#DECISION: GOOD-enough residuals


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = du_gwp100,
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
  #DECISION: DOES NOT CONVERGE
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = du_gwp100,
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
  #DECISION: good-enough

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does significantly improve fit (worth the extra complexity)

#Save best models (simple and complex best)
du_simplemodel_gwp100<- m1_gaus_nostrata
du_complexmodel_gwp100<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_gwp100, 
                        model2 = du_complexmodel_gwp100,
                        data=du_gwp100,
                        color_var = "vegpresence")

ggsave(filename = "DU_gwp100_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_gwp100 (ok) -------

#Subset data and check transformation: 
ri_gwp100<- data4models %>% filter(casepilot=="RI"&ghgspecies=="gwp100")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="gwp100") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_gwp100,
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
#DECISION: GOOD residuals


#SKIP! cannot include vegpresence as fixed with status, and gaus is already good
if(F){
  #Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
  m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                                data = ri_gwp100,
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
  #DECISION: 

#SKIP t_options: 

#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ri_gwp100,
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
#DECISION:

#SKIP, cannot include vegpresence as fixed effect
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ri_gwp100,
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
  #DECISION: 
}

#Compare models: gaussian simple is already the best 
# anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

#Save best models (simple and complex best)
ri_simplemodel_gwp100<- m1_gaus_nostrata
# ri_complexmodel_gwp100<- NA # Complex not possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_gwp100, 
                        # model2 = ri_complexmodel_gwp100,
                        data=ri_gwp100,
                        color_var = "vegpresence")

ggsave(filename = "RI_gwp100_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_gwp100 (ok) -------

#Subset data and check transformation: 
va_gwp100<- data4models %>% filter(casepilot=="VA"&ghgspecies=="gwp100")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="gwp100") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_gwp100,
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
#DECISION: BAD residuals 


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_gwp100,
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
#DECISION: BAD residuals 


  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = va_gwp100,
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
  #DECISION: DOES NOT CONVERGE
  
  
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = va_gwp100,
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
  #DECISION:GOOD residuals

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_gwp100<- m1_gaus_nostrata
va_complexmodel_gwp100<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_gwp100, 
                        model2 = va_complexmodel_gwp100,
                        data=va_gwp100,
                        color_var = "vegpresence")

ggsave(filename = "VA_gwp100_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)



#GWP20______-------

##2.1. CA_gwp20 (ok) -------

#Subset data and check transformation: 
ca_gwp20<- data4models %>% filter(casepilot=="CA"&ghgspecies=="gwp20")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CA"&
           ghgspecies=="gwp20") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ca_gwp20,
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
#DECISION: Good residuals, keep



#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = ca_gwp20,
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
#DECISION: Good residuals, keep

#Skip t-family, gaussian already good:
if(F){
#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = ca_gwp20,
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
#DECISION: 

#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = ca_gwp20,
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
#DECISION: 
}

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does significantly improve results

#Save best models (simple and complex best)
ca_simplemodel_gwp20<- m1_gaus_nostrata
ca_complexmodel_gwp20<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = ca_simplemodel_gwp20, 
                        model2 = ca_complexmodel_gwp20,
                        data=ca_gwp20,
                        color_var = "vegpresence")

ggsave(filename = "CA_gwp20_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)




##2.2. CU_gwp20 (ok) -------

#Subset data and check transformation: 
cu_gwp20<- data4models %>% filter(casepilot=="CU"&ghgspecies=="gwp20")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="CU"&
           ghgspecies=="gwp20") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = cu_gwp20,
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
#DECISION: BAD residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = cu_gwp20,
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
#DECISION: BAD residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = cu_gwp20,
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
#DECISION: Good enough residuals


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = cu_gwp20,
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
#DECISION: GOOD residuals

#Compare models 
anova(m3_t_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
cu_simplemodel_gwp20<- m3_t_nostrata
cu_complexmodel_gwp20<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = cu_simplemodel_gwp20, 
                        model2 = cu_complexmodel_gwp20,
                        data=cu_gwp20,
                        color_var = "vegpresence")

ggsave(filename = "CU_gwp20_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.3. DA_gwp20 (ok) -------

#Subset data and check transformation: 
da_gwp20<- data4models %>% filter(casepilot=="DA"&ghgspecies=="gwp20")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DA"&
           ghgspecies=="gwp20") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = da_gwp20,
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
#DECISION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = da_gwp20,
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
#DECISION: BAD residuals



#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = da_gwp20,
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
#DECISION: Better residuals, although they still have issues


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = da_gwp20,
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
da_simplemodel_gwp20<- m3_t_nostrata
da_complexmodel_gwp20<- m4_t_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = da_simplemodel_gwp20, 
                        model2 = da_complexmodel_gwp20,
                        data=da_gwp20,
                        color_var = "vegpresence")

ggsave(filename = "DA_gwp20_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)


##2.1. DU_gwp20 (ok) -------

#Subset data and check transformation: 
du_gwp20<- data4models %>% filter(casepilot=="DU"&ghgspecies=="gwp20")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="DU"&
           ghgspecies=="gwp20") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = du_gwp20,
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
#DECISION: Bad residuals


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = du_gwp20,
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
#DECISION: GOOD-enough residuals


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = du_gwp20,
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
#DECISION: Bad residuals


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = du_gwp20,
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
#DECISION: good-enough

#Compare models 
anova(m1_gaus_nostrata,m2_gaus_vegpresence)
#Adding vegpresence does significantly improve fit (worth the extra complexity)

#Save best models (simple and complex best)
du_simplemodel_gwp20<- m1_gaus_nostrata
du_complexmodel_gwp20<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = du_simplemodel_gwp20, 
                        model2 = du_complexmodel_gwp20,
                        data=du_gwp20,
                        color_var = "vegpresence")

ggsave(filename = "DU_gwp20_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. RI_gwp20 (ok) -------

#Subset data and check transformation: 
ri_gwp20<- data4models %>% filter(casepilot=="RI"&ghgspecies=="gwp20")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="RI"&
           ghgspecies=="gwp20") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = ri_gwp20,
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
#DECISION: GOOD residuals


#SKIP! cannot include vegpresence as fixed with status, and gaus is already good
if(F){
  #Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
  m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                                data = ri_gwp20,
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
  #DECISION: 
  
  #SKIP t_options: 
  
  #T_family:
  #Most basic: t_family, only status and season 
  m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                          data = ri_gwp20,
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
  #DECISION:
  
  #SKIP, cannot include vegpresence as fixed effect
  #Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
  m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                             data = ri_gwp20,
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
  #DECISION: 
}

#Compare models: gaussian simple is already the best 
# anova(m1_gaus_nostrata,m3_t_nostrata)

#CANNOT HAVE VEGPRESENCE AS EFFECT FOR RIA DE AVEIRO, most basic model is best (and only allowed)

#Save best models (simple and complex best)
ri_simplemodel_gwp20<- m1_gaus_nostrata
# ri_complexmodel_gwp20<- NA # Complex not possible


#PLot observed vs predicted ONLY option: 
plot_obs_vs_pred_models(model1 = ri_simplemodel_gwp20, 
                        # model2 = ri_complexmodel_gwp20,
                        data=ri_gwp20,
                        color_var = "vegpresence")

ggsave(filename = "RI_gwp20_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

##2.1. VA_gwp20 (ok) -------

#Subset data and check transformation: 
va_gwp20<- data4models %>% filter(casepilot=="VA"&ghgspecies=="gwp20")

#See effect of transformation:
table_trans %>%
  filter(casepilot=="VA"&
           ghgspecies=="gwp20") %>%
  pull(bn_result) %>%
  pluck(1) %>% 
  plot_bestNormalize()


#Most basic: gaussian, only status and season 
m1_gaus_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                           data = va_gwp20,
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
#DECISION: BAD residuals, bad model 


#Vegetation presence: gaussian, status, season and veg presence as fixed, subsite as random
m2_gaus_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                              data = va_gwp20,
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
#DECISION: BAD residuals 


#T_family:
#Most basic: t_family, only status and season 
m3_t_nostrata<- glmmTMB(formula = dailyflux_trans~status*season + (1|subsite), 
                        data = va_gwp20,
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
#DECISION: DOES NOT CONVERGE


#Vegetation presence: t_family, status, season and veg presence as fixed, subsite as random
m4_t_vegpresence<- glmmTMB(formula = dailyflux_trans~status*season*vegpresence + (1|subsite), 
                           data = va_gwp20,
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
#DECISION:GOOD residuals

#Compare models 
anova(m1_gaus_nostrata,m4_t_vegpresence)
#Adding vegpresence is significantly better (better fit worth complexity increase)

#Save best models (simple and complex best)
va_simplemodel_gwp20<- m1_gaus_nostrata
va_complexmodel_gwp20<- m2_gaus_vegpresence


#PLot observed vs predicted of two options: 
plot_obs_vs_pred_models(model1 = va_simplemodel_gwp20, 
                        model2 = va_complexmodel_gwp20,
                        data=va_gwp20,
                        color_var = "vegpresence")

ggsave(filename = "VA_gwp20_Observed_VS_predicted.png", 
       width = 160,height = 120,units = "mm",device = "png",dpi = 400,
       path = extraplots_path
)

}

#3. Create bestmodel-list -------

#Save lists of models: 
  #1. simple_model_list (only status*season models), only RI casepilot, 
  #2. complex_model_list (status*season*vegpresence), rest of casepilots 
#All models  named as: casepilot_ghgspecies


#List of simple models for each casepilot*ghgspecies: only appropriate for RI 
simplemodel_list_allghg<- list(
  "RI_co2" = ri_simplemodel_co2,
  "RI_ch4" = ri_simplemodel_ch4,
  "RI_gwp100" = ri_simplemodel_gwp100,
  "RI_gwp20" = ri_simplemodel_gwp20
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
  "CA_gwp100" = ca_complexmodel_gwp100,
  "CU_gwp100" = cu_complexmodel_gwp100, 
  "DA_gwp100" = da_complexmodel_gwp100,
  "DU_gwp100" = du_complexmodel_gwp100,
  "VA_gwp100" = va_complexmodel_gwp100,
  "CA_gwp20" = ca_complexmodel_gwp20,
  "CU_gwp20" = cu_complexmodel_gwp20, 
  "DA_gwp20" = da_complexmodel_gwp20,
  "DU_gwp20" = du_complexmodel_gwp20,
  "VA_gwp20" = va_complexmodel_gwp20
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
rm(ca_co2,ca_ch4,ca_gwp100,ca_gwp20,
   cu_co2,cu_ch4,cu_gwp100,cu_gwp20,
   da_co2,da_ch4,da_gwp100,da_gwp20,
   du_co2,du_ch4,du_gwp100,du_gwp20,
   ri_co2,ri_ch4,ri_gwp100,ri_gwp20,
   va_co2,va_ch4,va_gwp100,va_gwp20,
   ca_simplemodel_co2,ca_simplemodel_ch4,ca_simplemodel_gwp100,ca_simplemodel_gwp20,
   cu_simplemodel_co2,cu_simplemodel_ch4,cu_simplemodel_gwp100,cu_simplemodel_gwp20,
   da_simplemodel_co2,da_simplemodel_ch4,da_simplemodel_gwp100,da_simplemodel_gwp20,
   du_simplemodel_co2,du_simplemodel_ch4,du_simplemodel_gwp100,du_simplemodel_gwp20,
   ri_simplemodel_co2,ri_simplemodel_ch4,ri_simplemodel_gwp100,ri_simplemodel_gwp20,
   va_simplemodel_co2,va_simplemodel_ch4,va_simplemodel_gwp100,va_simplemodel_gwp20,
   ca_complexmodel_co2,ca_complexmodel_ch4,ca_complexmodel_gwp100,ca_complexmodel_gwp20,
   cu_complexmodel_co2,cu_complexmodel_ch4,cu_complexmodel_gwp100,cu_complexmodel_gwp20,
   da_complexmodel_co2,da_complexmodel_ch4,da_complexmodel_gwp100,da_complexmodel_gwp20,
   du_complexmodel_co2,du_complexmodel_ch4,du_complexmodel_gwp100,du_complexmodel_gwp20,
   #No RI complexmodel to remove
   va_complexmodel_co2,va_complexmodel_ch4,va_complexmodel_gwp100,va_complexmodel_gwp20,
   m1_gaus_nostrata,m2_gaus_vegpresence,m3_t_nostrata,m4_t_vegpresence,res,
   complexmodel_fit, complexmodel_homoc_r2, complexmodel_resid_diag_summary, complexmodel_results_anova,
   simplemodel_fit, simplemodel_homoc_r2, simplemodel_resid_diag_summary, simplemodel_results_anova)



## 3.3. EMEANS & post-hoc (bootstraped)------

#Calculate relevant emmeans and post-hocs tests: separately for simple model
#list (status*season, only for RI) and complex model list
#(status*season*vegpresence, for rest of casepilots)

#For simplemodels, use equal weights (same importance to all seasons for status comparison and same importance of statuses to assess seasonal trends)

#For complexmodels, use custom weights to account for seasonally-variable and for status-variable presence of vegetation. See complex models section below. 


#For more reliable results, instead of relying on the models results and uncertainty directly for statistical tests, we obtain EMMs, contrasts and p-values using bootstraping with re-fiting. 


#For each fitted GLMM, we implement the following:
# 1. Bootstrap simulation and refitting:
# Generate multiple (1000) simulated response datasets refit the model to each dataset to propagate uncertainty.
# 2. Weighted EMM estimation (model scale):
# For each bootstrap replicate, compute EMMs using a reference grid with custom weights reflecting the observed proportion of vegetation presence within each status × season combination, ensuring results represent field conditions rather than sampling imbalance. For simple models (no vegetation included), use equal weights. 
# 3. Back-transformation to original scale:
# Back-transform the bootstrap distributions of EMMs to the original scale to retain interpretability of contrasts while correctly propagating uncertainty.
# 4. Summarisation of EMM distributions:
# Summarise bootstrap EMMs using the empirical mean, standard deviation (SE), and percentile-based confidence intervals (2.5–97.5%).
# 5. Contrast calculation and uncertainty:
# Compute contrasts from back-transformed paired bootstrap EMMs (preserving covariance), applying the same weighting scheme, and summarise using mean, SE, and confidence intervals.
# 6. Empirical inference and visualization support:
# Derive two-sided p-values from bootstrap contrast distributions (in model-scale, with same weighting scheme) and construct compact letter displays (CLDs) to facilitate graphical comparison.

#Empirical p-values are derived from the bootstrap contrast distributions as the proportion of estimates crossing zero: 

#For each contrast (eg. altered - restored), get the proportion of values above and below zero, and use the minimum for the p value. This is an empirical p-value: if true contrast is at zero, we would expect the distribution to be centered at zero (with same proportion of values above and below it), if true contrast is not cero, we would expect the distribution to be mostly above or below it (the proportion of values in the minority tells us how "far" from zero the distribution is). p-values are adjusted with "holm" method

#Functions are defined the first time they are needed, then re-used for complex-models

#Set number of bootstrap simulations
nsim <- 1000

#NOTE: each model takes ~3 minutes with 1000 simulations. 
#3min*24models ~ 75 minutes of total computing time!


#SIMPLE MODELS------

##A) Bootstrap emmean-----

#Loop to bootstrap emmeans from simple models (using equal weights): 

# PARALLEL SETUP
workers <- max(1, parallel::detectCores() - 2)
plan(multisession, workers = workers)

# Enable progressr
handlers(global = TRUE)
handlers("txtprogressbar")  # or "progress" for nicer display

# INITIALIZE STORAGE
simple_emmean_list <- vector("list", length(simplemodel_list_allghg))
names(simple_emmean_list) <- names(simplemodel_list_allghg)

# LOOP OVER MODELS
for (m in seq_along(simplemodel_list_allghg)) {
  
  model_name <- names(simplemodel_list_allghg)[m]
  cat("\n==============================\n")
  cat("Bootstrapping model", m, "of", length(simplemodel_list_allghg), ":", model_name, "\n")
  cat("==============================\n")
  
  mod <- simplemodel_list_allghg[[m]]
  
  dat <- model.frame(mod)
  response_name <- all.vars(formula(mod))[1]
  
  # Factor levels (for consistent ordering)
  status_levels <- levels(dat$status)
  season_levels <- levels(dat$season)
  
  inseason_status_labels <- as.vector(
    outer(season_levels, status_levels, paste, sep = "_")
  )
  
  # Pre-simulate responses
  sim_matrix <- simulate(mod, nsim = nsim)
  
  # PARALLEL BOOTSTRAP WITH PROGRESS
  results <- with_progress({
    
    p <- progressor(steps = nsim)
    
    future_lapply(
      1:nsim,
      function(i) {
        
        # Load packages inside worker
        library(glmmTMB)
        library(emmeans)
        
        sim_y <- sim_matrix[[i]]
        
        # Create local copy of data
        dat_sim <- dat
        dat_sim[[response_name]] <- sim_y
        
        # Refit model
        mod_sim <- try(update(mod, data = dat_sim), silent = TRUE)
        if (inherits(mod_sim, "try-error")) {
          p()  # still update progress
          return(NULL)
        }
        
        # Build reference grid
        rg <- ref_grid(mod_sim)
        
        # Compute emmeans
        em_status <- suppressMessages(emmeans(rg, ~ status))
        em_season <- suppressMessages(emmeans(rg, ~ season))
        em_ss <- suppressMessages(emmeans(rg, ~ status | season))
        
        # Extract numeric values
        
        # STATUS
        status_vals <- as.data.frame(em_status)$emmean
        
        # SEASON
        season_vals <- as.data.frame(em_season)$emmean
        
        # STATUS WITHIN SEASON
        ss_df <- as.data.frame(em_ss)
        ss_df$label <- paste(ss_df$season, ss_df$status, sep = "_")
        ss_df <- ss_df[match(inseason_status_labels, ss_df$label), ]
        
        # Update progress
        p()
        
        return(list(
          status = status_vals,
          season = season_vals,
          ss = ss_df$emmean
        ))
      },
      future.seed = TRUE,
      future.packages = c("glmmTMB", "emmeans")
    )
  })
  
  # COLLECT RESULTS INTO MATRICES
  emm_status_boot <- matrix(NA, nsim, length(status_levels),
                            dimnames = list(NULL, status_levels))
  
  emm_season_boot <- matrix(NA, nsim, length(season_levels),
                            dimnames = list(NULL, season_levels))
  
  emm_inseason_status_boot <- matrix(NA, nsim, length(inseason_status_labels),
                                     dimnames = list(NULL, inseason_status_labels))
  
  for (i in seq_along(results)) {
    
    res <- results[[i]]
    if (is.null(res)) next
    
    emm_status_boot[i, ] <- res$status
    emm_season_boot[i, ] <- res$season
    emm_inseason_status_boot[i, ] <- res$ss
  }
  
  # STORE RESULTS
  simple_emmean_list[[model_name]] <- list(
    status = emm_status_boot,
    season = emm_season_boot,
    inseason_status = emm_inseason_status_boot
  )
}



##B) Back-transform emmeans----- 

#Function for back-transforming with bn_object (created at beginning of script for each dataset)
backtransform_emmeans <- function(emm_matrix, bn_obj) {
  
  # Apply inverse transform column-wise
  bt_matrix <- apply(emm_matrix, 2, function(x) {
    predict(bn_obj, newdata = x, inverse = TRUE)
  })
  
  # Ensure matrix structure is preserved (apply can simplify)
  bt_matrix <- as.matrix(bt_matrix)
  
  # Restore column names
  colnames(bt_matrix) <- colnames(emm_matrix)
  
  return(bt_matrix)
}


#Apply back-transformation: 

bt_simple_emmean_list <- vector("list", length(simple_emmean_list))
names(bt_simple_emmean_list) <- names(simple_emmean_list)

for (name in names(simple_emmean_list)) {
  
  # Get corresponding transformation object
  bn_obj <- bn_list[[name]]
  
  # Get emmeans object
  model_res <- simple_emmean_list[[name]]
  
  # Apply back-transformation to each component
  bt_simple_emmean_list[[name]] <- list(
    status = backtransform_emmeans(model_res$status, bn_obj),
    season = backtransform_emmeans(model_res$season, bn_obj),
    inseason_status = backtransform_emmeans(model_res$inseason_status, bn_obj)
  )
}


##C) Emmean summary-----

#Create function to summarise EMMs distribution for a given model and parameter
summarise_boot_matrix <- function(mat, model_name, parameter_name) {
  
  # Convert matrix to long format
  df_long <- as.data.frame(mat) %>%
    mutate(iter = row_number()) %>%
    pivot_longer(
      cols = -iter,
      names_to = "level",
      values_to = "value"
    )
  
  # Summarise per level
  df_summary <- df_long %>%
    group_by(level) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      SE = sd(value, na.rm = TRUE),
      lowerCI = quantile(value, 0.025, na.rm = TRUE),
      upperCI = quantile(value, 0.975, na.rm = TRUE),
      n = sum(!is.na(value)),
      .groups = "drop"
    ) %>%
    mutate(
      model = model_name,
      parameter = parameter_name
    ) %>%
    dplyr::select(model, parameter, level, mean, SE, lowerCI, upperCI, n)
  
  return(df_summary)
}


#Summarise EMMs for simplemodels
simpleemean_summary_list <- list()

for (model_name in names(bt_simple_emmean_list)) {
  
  model_res <- bt_simple_emmean_list[[model_name]]
  
  # Summarise each parameter
  summary_status <- summarise_boot_matrix(
    mat = model_res$status,
    model_name = model_name,
    parameter_name = "status"
  )
  
  summary_season <- summarise_boot_matrix(
    mat = model_res$season,
    model_name = model_name,
    parameter_name = "season"
  )
  
  summary_inseason_status <- summarise_boot_matrix(
    mat = model_res$inseason_status,
    model_name = model_name,
    parameter_name = "inseason_status"
  )
  
  # Combine
  simpleemean_summary_list[[model_name]] <- bind_rows(
    summary_status,
    summary_season,
    summary_inseason_status
  )
}

# Unlist into data frame
simpleemmean_summary_df <- bind_rows(simpleemean_summary_list)



##D) Contrasts------

# HELPER FUNCTION: compute contrasts from bootstrap matrix
# mat: [nsim x levels]
# contrast_pairs: list of pairs like c("A","B")
compute_contrasts <- function(mat, contrast_pairs) {
  
  # Preallocate matrix
  res <- matrix(NA, nrow = nrow(mat), ncol = length(contrast_pairs))
  
  # Fill matrix
  for (i in seq_along(contrast_pairs)) {
    pair <- contrast_pairs[[i]]
    res[, i] <- mat[, pair[1]] - mat[, pair[2]]
  }
  
  # Assign names immediately (no repair needed)
  colnames(res) <- sapply(contrast_pairs, function(x) {
    paste0(x[1], " - ", x[2])
  })
  
  return(res)
}



# MAIN LOOP OVER MODELS

contrast_summary_list <- list()

for (model_name in names(bt_simple_emmean_list)) {
  
  #Get emmeans distributions
  #Back-transformed emmeans (for interpretation effect size)
  model_res <- bt_simple_emmean_list[[model_name]]
  
  #OG-scale emmeans (for inference of significance, p-values)
  model_res_og <- simple_emmean_list[[model_name]] 
  
  #Ensure consistent column ordering: 
  # STATUS
  model_res$status <- model_res$status[, status_levels, drop = FALSE]
  model_res_og$status <- model_res_og$status[, status_levels, drop = FALSE]
  
  # SEASON
  model_res$season <- model_res$season[, season_levels, drop = FALSE]
  model_res_og$season <- model_res_og$season[, season_levels, drop = FALSE]
  
  # STATUS WITHIN SEASON
  model_res$inseason_status <- model_res$inseason_status[, inseason_status_labels, drop = FALSE]
  model_res_og$inseason_status <- model_res_og$inseason_status[, inseason_status_labels, drop = FALSE]
  
  # 1) STATUS CONTRASTS (fixed order)
  status_pairs <- list(
    c("Altered", "Preserved"),
    c("Altered", "Restored"),
    c("Preserved", "Restored")
  )
  
  status_contr <- compute_contrasts(model_res$status, status_pairs)
  status_contr_og <- compute_contrasts(model_res_og$status, status_pairs)
  
  # 2) SEASON CONTRASTS (all pairwise combinations)
  season_levels <- colnames(model_res$season)
  season_pairs <- combn(season_levels, 2, simplify = FALSE)
  
  season_contr <- compute_contrasts(model_res$season, season_pairs)
  season_contr_og <- compute_contrasts(model_res_og$season, season_pairs)
  
  # 3) STATUS WITHIN SEASON CONTRASTS
  ss_mat <- model_res$inseason_status
  ss_mat_og <- model_res_og$inseason_status
  
  # Extract seasons and statuses from column names
  col_split <- strsplit(colnames(ss_mat), "_")
  
  seasons <- unique(sapply(col_split, `[`, 1))
  statuses <- unique(sapply(col_split, `[`, 2))
  
  #Init contrast lists for ss 
  ss_contr_list <- list()
  ss_contr_og_list <- list()
  
  for (s in seasons) {
    
    # Build full column names for this season
    cols <- paste0(s, "_", statuses)
    
    # Subset matrix (keep order of statuses!)
    sub_mat <- ss_mat[, cols, drop = FALSE]
    sub_mat_og <- ss_mat_og[, cols, drop = FALSE]
    
    # Now build contrasts using FULL column names
    contrast_pairs_full <- lapply(status_pairs, function(pair) {
      c(paste0(s, "_", pair[1]),
        paste0(s, "_", pair[2]))
    })
    
    # Compute contrasts
    sub_contr <- compute_contrasts(sub_mat, contrast_pairs_full)
    sub_contr_og <- compute_contrasts(sub_mat_og, contrast_pairs_full)
    
    #Assign contrasts    
    ss_contr_list[[s]] <- sub_contr
    ss_contr_og_list[[s]] <- sub_contr_og
  }
  
  # Combine all seasons
  ss_contr <- do.call(cbind, ss_contr_list)
  ss_contr_og <- do.call(cbind, ss_contr_og_list)
  
  
#Function to summarise contrasts from back-transformed Emmeans (to report in tables as effect sizes)
#With functionality to re-use for complexmodels 
  summarise_contr <- function(mat, param_name) {
    df <- as.data.frame(mat) %>%
      mutate(iter = row_number()) %>%
      pivot_longer(-iter, names_to = "contrast", values_to = "value") %>%
      group_by(contrast) %>%
      summarise(
        # --- central tendency ---
        mean = mean(value, na.rm = TRUE),
        
        # --- bootstrap standard error ---
        SE = sd(value, na.rm = TRUE),
        
        # --- percentile CI ---
        lowerCI = quantile(value, 0.025, na.rm = TRUE),
        upperCI = quantile(value, 0.975, na.rm = TRUE),
        
        # --- number of bootstrap samples ---
        n = sum(!is.na(value)),
        .groups = "drop"
      )
    #Add metadata
    df <- df %>%
      mutate(
        model = model_name,
        parameter = param_name
      ) %>%
      dplyr::select(model, parameter, contrast, mean, SE, lowerCI, upperCI, n)
    return(df)
  }
  

#Function to calculate empirical pvalue from contrasts derived from model-scale emmeans. 
#With functionality to re-use for complexmodels 
  pval_contr_og <- function(mat, param_name) {
    
    df <- as.data.frame(mat) %>%
      mutate(iter = row_number()) %>%
      pivot_longer(-iter, names_to = "contrast", values_to = "value") %>%
      group_by(contrast) %>%
      summarise(
        # --- bootstraped empirical p-value (two-sided) ---
        p = 2 * min(mean(value <= 0, na.rm = TRUE),
                    mean(value >= 0, na.rm = TRUE)),
        .groups = "drop"
      )
    
    # --- adjust p-values---
    #IF statements for different comparisons:
    #inseason_status only compares status within each level of season 
    if (param_name == "inseason_status") {
      # Extract season from contrast name
      df <- df %>%
        mutate(season = sub("^([^_]+)_.*", "\\1", contrast)) %>%
        group_by(season) %>%
        mutate(p_adj = p.adjust(p, method = "holm")) %>%
        ungroup() %>%
        dplyr::select(-season)
      
    } else if(param_name == "inveg_status"){
      #inveg_status only compares status within each level of vegpresence
      # Extract vegpresence from contrast name
      df <- df %>%
        mutate(vegpresence = sub("^([^_]+)_.*", "\\1", contrast)) %>%
        group_by(vegpresence) %>%
        mutate(p_adj = p.adjust(p, method = "holm")) %>%
        ungroup() %>%
        dplyr::select(-vegpresence)
      
    }else {
      # Default: adjust across all contrasts (for status and season comparisons)
      df <- df %>%
        mutate(p_adj = p.adjust(p, method = "holm"))
    }
    
    #Add metadata
    df <- df %>%
      mutate(
        model = model_name,
        parameter = param_name
      ) %>%
      dplyr::select(model, parameter, contrast, p, p_adj)
    
    return(df)
  }
  
  #Perform summary of contrasts (back-transformed) and add p-value (model-scale)
  #status contrasts and pvalues
  sum_status_contr<- summarise_contr(status_contr, "status")
  pval_status_contr<- pval_contr_og(status_contr_og, "status")
  status_contrast_sumpval<- sum_status_contr %>% left_join(pval_status_contr,
                                                           by = c("model","parameter","contrast"))
  
  #season contrasts and pvalues
  sum_season_contr<- summarise_contr(season_contr, "season")
  pval_season_contr<- pval_contr_og(season_contr_og, "season")
  season_contrast_sumpval<- sum_season_contr %>% left_join(pval_season_contr,
                                                           by = c("model","parameter","contrast"))
  
  #inseasonstatus contrasts and pvalues
  sum_ss_contr<- summarise_contr(ss_contr, "inseason_status")
  pval_ss_contr<- pval_contr_og(ss_contr_og, "inseason_status")
  ss_contrast_sumpval<- sum_ss_contr %>% left_join(pval_ss_contr,
                                                   by = c("model","parameter","contrast"))
  
  #Assign contrasts and pvalues
  contrast_summary_list[[model_name]] <- bind_rows(
    status_contrast_sumpval,
    season_contrast_sumpval,
    ss_contrast_sumpval
  )
}

#Unlist into data frame
simplecontrast_summary_df <- bind_rows(contrast_summary_list)



##E) CLD and format-----
#Add CLD to emmeans and format final tables (EMMs+CLD and Contrasts)

#Format contrasts: 
posthoctests_RI_chambermodels<- simplecontrast_summary_df %>% 
  separate(model, into = c("casepilot","ghgspecies")) %>% 
  rename(estimate_bt=mean, rawcontrast=contrast, SE_bt=SE,lower.CL_bt=lowerCI,upper.CL_bt=upperCI,boot_n=n,p.value= p_adj) %>% 
  mutate(comparison=case_when(parameter=="inseason_status"~"status_within_season",
                              parameter=="inveg_status"~"status_within_vegpresence",
                              parameter=="status"~"status",
                              parameter=="season"~"season",
                              TRUE~NA),
         season = stringr::str_extract(rawcontrast, "S[1-4]"),
         contrast = stringr::str_remove(stringr::str_remove(rawcontrast, "S[1-4]_"), "S[1-4]_")
  ) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, contrast,season, estimate_bt,SE_bt,lower.CL_bt,upper.CL_bt,boot_n,p.value)


#Obtain CLDs (letter-groups) based on post-hoc tests for every comparison:
simpleCLD_letters <- posthoctests_RI_chambermodels %>%
  #leave only (1) variables that identify unique comparisons, (2)contrast column, (3) pvalue column
  #seasonlevel is kept to account for status_within_season comparison. 
  dplyr::select(casepilot, ghgspecies, comparison, season, contrast, p.value) %>% 
  #Remove any spaces from contrast column (mutcompLetters expects levels to be only separated by a hyphen "-")
  mutate(contrast=gsub(" ","", contrast),
         #Add seasonlevel for status_within_season comparison
         seasonlevel=if_else(comparison=="status_within_season", season, NA)) %>% 
  #Group by comparison identifiers
  group_by(casepilot, ghgspecies, comparison, seasonlevel) %>%
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
  mutate(
    status = stringr::str_extract(level, "Altered|Preserved|Restored"),
    season = coalesce(seasonlevel,  stringr::str_extract(level, "S[1-4]"))
    )%>%
  dplyr::select(casepilot, ghgspecies, comparison, status, season, seasonlevel, cld_group)


#Format EMMs
EMMs_RI_formated<- simpleemmean_summary_df %>% 
  separate(model, into = c("casepilot","ghgspecies")) %>% 
  rename(emmean_bt=mean, SE_bt=SE,lower.CL_bt=lowerCI,upper.CL_bt=upperCI,boot_n=n) %>% 
  mutate(comparison=case_when(parameter=="inseason_status"~"status_within_season",
                              parameter=="inveg_status"~"status_within_vegpresence",
                              parameter=="status"~"status",
                              parameter=="season"~"season",
                              TRUE~NA),
         season = str_extract(level, "S[1-4]"),
         status = str_extract(level, "Altered|Preserved|Restored")) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, status, season, boot_n, emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt)


#Add CLD to formated emmeans
EMMsCLD_RI_formated <- EMMs_RI_formated %>%
  #ADD CLD group-letters (re-coding them so that "a" always identifies the
  #significant group with lowest emmean, b the next lowest emmean, and so on...)
  left_join(simpleCLD_letters, by = c("casepilot","ghgspecies","comparison","status","season")) %>%
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





#COMPLEX MODELS------

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


##A) Bootstrap emmean-----


#First, Build function to calculate propperly weighted emmeans from each model within loop
compute_weighted_emmeans <- function(mod, weights_df) {
  
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



# INITIALIZE OUTPUT
complex_emmean_list <- vector("list", length(complexmodel_list_allghg))
names(complex_emmean_list) <- names(complexmodel_list_allghg)


# LOOP OVER MODELS
for (m in seq_along(complexmodel_list_allghg)) {
  
  model_name <- names(complexmodel_list_allghg)[m]
  
  cat("\n==============================\n")
  cat("Bootstrapping model", m, "of", length(complexmodel_list_allghg), ":", model_name, "\n")
  cat("==============================\n")
  
  # Model + data
  mod <- complexmodel_list_allghg[[model_name]]
  dat <- model.frame(mod)
  response_name <- all.vars(formula(mod))[1]
  
  # Get weights
  cp <- sub("_.*", "", model_name)
  weights_df <- weights_list[[cp]]
  
  # Factor levels (REFERENCE ORDER)
  status_levels <- levels(dat$status)
  season_levels <- levels(dat$season)
  vegpresence_levels <- levels(dat$vegpresence)
  
  # Labels (STRICT ORDERING)
  inseason_status_labels <- as.vector(
    outer(season_levels, status_levels, paste, sep = "_")
  )
  
  inveg_status_labels <- as.vector(
    outer(vegpresence_levels, status_levels, paste, sep = "_")
  )
  
  # Pre-simulate responses (BIG speed gain)
  sim_matrix <- simulate(mod, nsim = nsim)

  # PARALLEL BOOTSTRAP WITH PROGRESS
  results <- with_progress({
    
    p <- progressor(steps = nsim)
    
    future_lapply(
      1:nsim,
      function(i) {
        
        library(glmmTMB)
        library(emmeans)
        library(dplyr)
        
        sim_y <- sim_matrix[[i]]
        
        # local copy (safe for parallel)
        dat_sim <- dat
        dat_sim[[response_name]] <- sim_y
        
        # refit
        mod_sim <- try(update(mod, data = dat_sim), silent = TRUE)
        if (inherits(mod_sim, "try-error")) {
          p()
          return(NULL)
        }
        
        # weighted emmeans (using function above)
        emm_list <- try(
          compute_weighted_emmeans(mod_sim, weights_df),
          silent = TRUE
        )
        
        if (inherits(emm_list, "try-error")) {
          p()
          return(NULL)
        }
        
        # FORMAT + ENFORCE ORDER
        # STATUS
        emm_status_df <- as.data.frame(emm_list$status) %>%
          arrange(status)
        
        # SEASON
        emm_season_df <- as.data.frame(emm_list$season) %>%
          arrange(season)
        
        # STATUS × SEASON
        emm_ss_df <- as.data.frame(emm_list$inseason_status) %>%
          mutate(label = paste(season, status, sep = "_"))
        
        emm_ss_df <- emm_ss_df[match(inseason_status_labels, emm_ss_df$label), ]
        
        # STATUS × VEGPRESENCE
        emm_sv_df <- as.data.frame(emm_list$inveg_status) %>%
          mutate(label = paste(vegpresence, status, sep = "_"))
        
        emm_sv_df <- emm_sv_df[match(inveg_status_labels, emm_sv_df$label), ]
        
        # progress update
        p()
        
        return(list(
          status = emm_status_df$emmean,
          season = emm_season_df$emmean,
          ss = emm_ss_df$emmean,
          sv = emm_sv_df$emmean
        ))
        
      },
      future.seed = TRUE,
      future.packages = c("glmmTMB", "emmeans", "dplyr")
    )
  })

  # COLLECT RESULTS
  emm_status_boot <- matrix(NA, nsim, length(status_levels),
                            dimnames = list(NULL, status_levels))
  
  emm_season_boot <- matrix(NA, nsim, length(season_levels),
                            dimnames = list(NULL, season_levels))
  
  emm_inseason_status_boot <- matrix(NA, nsim, length(inseason_status_labels),
                                     dimnames = list(NULL, inseason_status_labels))
  
  emm_inveg_status_boot <- matrix(NA, nsim, length(inveg_status_labels),
                                  dimnames = list(NULL, inveg_status_labels))
  
  for (i in seq_along(results)) {
    
    res <- results[[i]]
    if (is.null(res)) next
    
    emm_status_boot[i, ] <- res$status
    emm_season_boot[i, ] <- res$season
    emm_inseason_status_boot[i, ] <- res$ss
    emm_inveg_status_boot[i, ] <- res$sv
  }
  
  # STORE
  complex_emmean_list[[model_name]] <- list(
    status = emm_status_boot,
    season = emm_season_boot,
    inseason_status = emm_inseason_status_boot,
    inveg_status = emm_inveg_status_boot
  )
}


##B) Back-transform emmeans----- 
#Using function defined in simple models section.

#Init list: 
bt_complex_emmean_list <- vector("list", length(complex_emmean_list))
names(bt_complex_emmean_list) <- names(complex_emmean_list)

for (name in names(complex_emmean_list)) {
  
  # Get corresponding transformation object
  bn_obj <- bn_list[[name]]
  
  # Get emmeans object
  model_res <- complex_emmean_list[[name]]
  
  # Apply back-transformation to each component
  bt_complex_emmean_list[[name]] <- list(
    status = backtransform_emmeans(model_res$status, bn_obj),
    season = backtransform_emmeans(model_res$season, bn_obj),
    inseason_status = backtransform_emmeans(model_res$inseason_status, bn_obj),
    inveg_status = backtransform_emmeans(model_res$inveg_status, bn_obj)
  )
}


##C) Emmean summary-----
##Using function defined in simple models section

#Initialize list for storage
complexemmean_summary_list <- list()

#Compute summary
for (model_name in names(bt_complex_emmean_list)) {
  
  model_res <- bt_complex_emmean_list[[model_name]]
  
  # Summarise each parameter
  summary_status <- summarise_boot_matrix(
    mat = model_res$status,
    model_name = model_name,
    parameter_name = "status"
  )
  
  summary_season <- summarise_boot_matrix(
    mat = model_res$season,
    model_name = model_name,
    parameter_name = "season"
  )
  
  summary_inseason_status <- summarise_boot_matrix(
    mat = model_res$inseason_status,
    model_name = model_name,
    parameter_name = "inseason_status"
  )
  
  summary_inveg_status <- summarise_boot_matrix(
    mat = model_res$inveg_status,
    model_name = model_name,
    parameter_name = "inveg_status"
  )
  
  # Combine
  complexemmean_summary_list[[model_name]] <- bind_rows(
    summary_status,
    summary_season,
    summary_inseason_status, 
    summary_inveg_status
  )
}

# Final tidy dataframe
complexemmean_summary_df <- bind_rows(complexemmean_summary_list)



##D) Contrasts------
#Using functions "compute_contrast", "summarise_contr" and "pval_contr_og" defined in simple model section

#Init list
contrast_summary_list <- list()

for (model_name in names(bt_complex_emmean_list)) {
  
  #Get emmeans distributions
  #Back-transformed emmeans (for interpretation effect size)
  model_res <- bt_complex_emmean_list[[model_name]]
  
  #OG-scale emmeans (for inference of significance, p-values)
  model_res_og <- complex_emmean_list[[model_name]] 
  
  #Ensure consistent column ordering: 
  # STATUS
  model_res$status <- model_res$status[, status_levels, drop = FALSE]
  model_res_og$status <- model_res_og$status[, status_levels, drop = FALSE]
  
  # SEASON
  model_res$season <- model_res$season[, season_levels, drop = FALSE]
  model_res_og$season <- model_res_og$season[, season_levels, drop = FALSE]
  
  # STATUS WITHIN SEASON
  model_res$inseason_status <- model_res$inseason_status[, inseason_status_labels, drop = FALSE]
  model_res_og$inseason_status <- model_res_og$inseason_status[, inseason_status_labels, drop = FALSE]
  
  # STATUS WITHIN VEGPRESENCE
  model_res$inveg_status <- model_res$inveg_status[, inveg_status_labels, drop = FALSE]
  model_res_og$inveg_status <- model_res_og$inveg_status[, inveg_status_labels, drop = FALSE]
  
  
  # 1) STATUS CONTRASTS (fixed order)
  status_pairs <- list(
    c("Altered", "Preserved"),
    c("Altered", "Restored"),
    c("Preserved", "Restored")
  )
  
  status_contr <- compute_contrasts(model_res$status, status_pairs)
  status_contr_og <- compute_contrasts(model_res_og$status, status_pairs)
  
  
  # 2) SEASON CONTRASTS (all pairwise combinations)
  
  season_levels <- colnames(model_res$season)
  season_pairs <- combn(season_levels, 2, simplify = FALSE)
  
  season_contr <- compute_contrasts(model_res$season, season_pairs)
  season_contr_og <- compute_contrasts(model_res_og$season, season_pairs)
  
  
  # 3) STATUS WITHIN SEASON CONTRASTS
  ss_mat <- model_res$inseason_status
  ss_mat_og <- model_res_og$inseason_status
  
  # Extract seasons and statuses from column names
  col_split <- strsplit(colnames(ss_mat), "_")
  seasons <- unique(sapply(col_split, `[`, 1))
  statuses <- unique(sapply(col_split, `[`, 2))
  
  #Init contrast lists for ss 
  ss_contr_list <- list()
  ss_contr_og_list <- list()
  
  for (s in seasons) {
    
    # Build full column names for this season
    cols <- paste0(s, "_", statuses)
    
    # Subset matrix (keep order of statuses!)
    sub_mat <- ss_mat[, cols, drop = FALSE]
    sub_mat_og <- ss_mat_og[, cols, drop = FALSE]
    
    # Now build contrasts using FULL column names
    contrast_pairs_full <- lapply(status_pairs, function(pair) {
      c(paste0(s, "_", pair[1]),
        paste0(s, "_", pair[2]))
    })
    
    # Compute contrasts
    sub_contr <- compute_contrasts(sub_mat, contrast_pairs_full)
    sub_contr_og <- compute_contrasts(sub_mat_og, contrast_pairs_full)
    
    #Assign contrasts
    ss_contr_list[[s]] <- sub_contr
    ss_contr_og_list[[s]] <- sub_contr_og
  }
  
  # Combine all seasons
  ss_contr <- do.call(cbind, ss_contr_list)
  ss_contr_og <- do.call(cbind, ss_contr_og_list)
  
  
  # 4) STATUS WITHIN VEGPRESENCE CONTRASTS
  sv_mat <- model_res$inveg_status
  sv_mat_og <- model_res_og$inveg_status
  
  # Extract vegpresences and statuses from column names
  col_split <- strsplit(colnames(sv_mat), "_")
  vegpresences <- unique(sapply(col_split, `[`, 1))
  statuses <- unique(sapply(col_split, `[`, 2))
  
  sv_contr_list <- list()
  sv_contr_og_list <- list()
  
  for (v in vegpresences) {
    
    # Build full column names for this vegpresence
    cols <- paste0(v, "_", statuses)
    
    # Subset matrix (keep order of statuses defined before!)
    sub_mat <- sv_mat[, cols, drop = FALSE]
    sub_mat_og<- sv_mat_og[, cols, drop = FALSE]
    
    # Now build contrasts using FULL column names
    contrast_pairs_full <- lapply(status_pairs, function(pair) {
      c(paste0(v, "_", pair[1]),
        paste0(v, "_", pair[2]))
    })
    
    # Compute contrasts
    sub_contr <- compute_contrasts(sub_mat, contrast_pairs_full)
    sub_contr_og <- compute_contrasts(sub_mat_og, contrast_pairs_full)
    
    #Assign contrasts
    sv_contr_list[[v]] <- sub_contr
    sv_contr_og_list[[v]] <- sub_contr_og
  }
  
  # Combine all vegpresences
  sv_contr <- do.call(cbind, sv_contr_list)
  sv_contr_og <- do.call(cbind, sv_contr_og_list)
  
  # APPLY SUMMARISATION
  #Perform summary of contrasts (back-transformed) and add p-value (model-scale)
  #status contrasts and pvalues
  sum_status_contr<- summarise_contr(status_contr, "status")
  pval_status_contr<- pval_contr_og(status_contr_og, "status")
  status_contrast_sumpval<- sum_status_contr %>% left_join(pval_status_contr,
                                                           by = c("model","parameter","contrast"))
  
  #season contrasts and pvalues
  sum_season_contr<- summarise_contr(season_contr, "season")
  pval_season_contr<- pval_contr_og(season_contr_og, "season")
  season_contrast_sumpval<- sum_season_contr %>% left_join(pval_season_contr,
                                                           by = c("model","parameter","contrast"))
  
  #inseasonstatus contrasts and pvalues
  sum_ss_contr<- summarise_contr(ss_contr, "inseason_status")
  pval_ss_contr<- pval_contr_og(ss_contr_og, "inseason_status")
  ss_contrast_sumpval<- sum_ss_contr %>% left_join(pval_ss_contr,
                                                   by = c("model","parameter","contrast"))
  
  #invegstatus contrasts and pvalues
  sum_sv_contr<- summarise_contr(sv_contr, "inveg_status")
  pval_sv_contr<- pval_contr_og(sv_contr_og, "inveg_status")
  sv_contrast_sumpval<- sum_sv_contr %>% left_join(pval_sv_contr,
                                                   by = c("model","parameter","contrast"))
  
  #Combine results:
  contrast_summary_list[[model_name]] <- bind_rows(
    status_contrast_sumpval,
    season_contrast_sumpval,
    ss_contrast_sumpval,
    sv_contrast_sumpval
  )
  
}

#Unlist into data frame:
complexcontrast_summary_df <- bind_rows(contrast_summary_list)



##E)CLD and format-----
#Add CLD to emmeans and format final tables (EMMs+CLD and Contrasts)

#Format contrasts: 
posthoctests_rest_chambermodels<- complexcontrast_summary_df %>% 
  separate(model, into = c("casepilot","ghgspecies")) %>% 
  rename(estimate_bt=mean, rawcontrast=contrast, SE_bt=SE,lower.CL_bt=lowerCI,upper.CL_bt=upperCI,boot_n=n,p.value= p_adj) %>% 
  mutate(comparison=case_when(parameter=="inseason_status"~"status_within_season",
                              parameter=="inveg_status"~"status_within_vegpresence",
                              parameter=="status"~"status",
                              parameter=="season"~"season",
                              TRUE~NA),
         season = stringr::str_extract(rawcontrast, "S[1-4]"),
         vegpresence = stringr::str_extract(rawcontrast, "Non-vegetated|Vegetated"),
         #remove season and vegpresence grouping from contrast
         contrast = gsub("S[1-4]_|Non-vegetated_|Vegetated_","",rawcontrast)
         ) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, contrast,season, vegpresence,estimate_bt,SE_bt,lower.CL_bt,upper.CL_bt,boot_n,p.value)


#Obtain CLDs (letter-groups) based on post-hoc tests for every comparison:
complexCLD_letters <- posthoctests_rest_chambermodels %>%
  #leave only (1) variables that identify unique comparisons, (2)contrast column, (3) pvalue column
  #seasonlevel and vegpresencelevel is kept to account for status_within_season comparison. 
  dplyr::select(casepilot, ghgspecies, comparison, season,vegpresence, contrast, p.value) %>% 
  #Remove any spaces from contrast column (mutcompLetters expects levels to be only separated by a hyphen "-")
  mutate(contrast=gsub(" ","", contrast),
         #Add seasonlevel for status_within_season comparison
         seasonlevel=if_else(comparison=="status_within_season", season, NA),
         #Add vegpresenelevel for status_within_vegpresence comparison
         vegpresencelevel=if_else(comparison=="status_within_vegpresence", vegpresence, NA)) %>% 
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
    season = coalesce(seasonlevel,  stringr::str_extract(level, "S[1-4]")),
    vegpresence = coalesce(vegpresencelevel, stringr::str_extract(level,"Non-vegetated|Vegetated" ))
  )%>%
  dplyr::select(casepilot, ghgspecies, comparison, status, season, vegpresence, seasonlevel, vegpresencelevel, cld_group)


#Format EMMs
EMMs_rest_formated<- complexemmean_summary_df %>% 
  separate(model, into = c("casepilot","ghgspecies")) %>% 
  rename(emmean_bt=mean, SE_bt=SE,lower.CL_bt=lowerCI,upper.CL_bt=upperCI,boot_n=n) %>% 
  mutate(comparison=case_when(parameter=="inseason_status"~"status_within_season",
                              parameter=="inveg_status"~"status_within_vegpresence",
                              parameter=="status"~"status",
                              parameter=="season"~"season",
                              TRUE~NA),
         season = str_extract(level, "S[1-4]"),
         status = str_extract(level, "Altered|Preserved|Restored"),
         vegpresence = str_extract(level, "Non-vegetated|Vegetated")) %>% 
  dplyr::select(casepilot, ghgspecies, comparison, status, season,vegpresence, boot_n, emmean_bt, SE_bt, lower.CL_bt, upper.CL_bt)


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





#_____________----
#OLD CODE BELOW-------
##3.3. EMEANs & post-hoc----

#When calculating emmeans: I have to use the full variance-covariance structure
#to get accurate SEs. If the model is gaussian, we use t-test with Degrees of
#freedom available (via Satterthwaite approximation) IF the model is t_family,
#it will have Inf df, not an issue but forces to use z-test for post-hocs
#instead (No DF calculation possible)

#CLDs (group letters) should be assigned to emmeans based on the above pairwise
#tests, do not repeat the test, but assign letters based on already calcualted
#tests via (multcompLetters function). Additionally, re-code letters so that
#they give emmean ranking info (letter "a" denotes the smallest-emmean
#significant group, letter"b" the signficant group with the next lowest emmean

#Calculate relevant emmeans and post-hocs tests: separately for simple model
#list (status*season, only for RI) and complex model list
#(status*season*vegpresence, for rest of casepilots)



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

#Emmeans are extracted, pairwise tests done (based on model distribution
#family), CLDs are assigned to significantly different groups (and re-ordered
#based on emmean ranking). WE save both pairwise post-hoc tests in original
#(model scale) and back-transformed scales; AND emmeanCLDs in original and
#back-transformed scales.

#Equal Weights used for all EMMs: this assumes that all seasons and all status
#have the same weight for each other (appropriate).


#Initialize simple model results list
simple_comparison_list<- list()

# Loop to extract all emmeans and perform pairwise tests for appropriate
# comparisons for every casepilot*ghg simple model.

for (dataset in names(simplemodel_list_allghg)) {
  #get model
  cp_model <- simplemodel_list_allghg[[dataset]]
  #get transformation object
  cp_trans_obj<- bn_list[[dataset]]
  
  #extract casepilot name and ghgspecies from model-list names
  casepilot_name<- sub("_.*", "", dataset)
  ghgspecies<- sub(paste0(casepilot_name,"_"),"",dataset)
  
  # Obtain emmean object for each comparison: status, season,
  # status_within_season. Using custom function for consistency.
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
#T-test or Z-test depending on model distribution family used (automatically
#assigned by contrasts (method="pairwise)). p.value is was adjusted for multiple
#comparisons sidak
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
  #ADD CLD group-letters (re-coding them so that "a" always identifies the
  #significant group with lowest emmean, b the next lowest emmean, and so on...)
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
  arrange(casepilot,ghgspecies, comparison, season,status)

#Save emmeans and groupletters (back-transformed):  
write.csv(x = simplemodel_emmeansCLD, file = paste0(path_2_modeloutputs,"EmmeansCLD_RI_chambermodels.csv"),row.names = F)





#Complex models (ok)------

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


##Weighted-emmean-loop------

#Calculate emmeans for each relevant comparison, using weights when appropriate:

#Intialize list for complex models
complex_comparison_list<- list()

# Loop to extract all emmeans and perform pairwise tests for appropriate
# comparisons for every casepilot complexbestmodel.
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
#T-test or Z-test depending on model distribution family used (automatically
#assigned by contrasts (method="pairwise)). p.value is adjusted for multiple
#comparisons sidak
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
  #ADD CLD group-letters (re-coding them so that "a" always identifies the
  #significant group with lowest emmean, b the next lowest emmean, and so on...)
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

}
#_____________----

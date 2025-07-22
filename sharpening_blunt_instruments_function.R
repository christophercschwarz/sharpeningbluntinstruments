################################################################################
#               Sharpening Blunt Instruments: Bootstrap Inference              #
#                                                                              #
# Author:  Christopher Schwarz                                                 #
# Date:    07/19/2025                                                          #
# Purpose: Function for bootstraped inference of non-linear first stage IV     #
# Edits:   07/20/2025: Better handling of control variables                    #
################################################################################

library(mgcv)
library(future)
library(future.apply)
library(progressr)

################################################################################
#                                Progress Handler                              #
################################################################################

handlers(global = TRUE)
handlers(list(handler_progress(
  format   = ":spin :current/:total batches complete [:bar] :percent in :elapsedfull ETA: :eta",
  width    = 60,
  complete = "+",
  clear = FALSE
)))

################################################################################
#       Non-linear first stage two stage residual inclusion; single shot       #
################################################################################

# Controls is a list of control variables, the first being for the second stage
#  and the subsequent entries being for each of the first stage regressions.

nlfstsri_core <- function(data, outcome, endogenous, exogenous, controls = NULL, save_fs = TRUE){
  
  # Ensure controls list is valid
  if (is.null(controls) || length(controls) == 0) {
    controls <- vector("list", length(endogenous))
  } else if (length(controls) != length(endogenous) + 1) {
    stop("Length of controls must match number of endogenous variables plus one.")
  }
  
  # Step 1: Estimate first-stage models and get residuals
  first_stage_resids <- list()
  
  if(save_fs){
    
    fs_mods = list()
    
  }
  
  for (i in seq_along(endogenous)) {
    
    endo_var <- endogenous[i]
    ctrl_vars <- controls[[i + 1]]
  
    formula_str <- paste0(endo_var, " ~ ")
    
    for(inst in exogenous){
      
      formula_str <- paste(formula_str, paste0("+ te(",inst,")"))
      
    }
    
    if (!is.null(ctrl_vars)) {
      formula_str <- paste0(formula_str, " + ", paste(ctrl_vars, collapse = " + "))
    }
    
    model <- gam(as.formula(formula_str), data = data)
    first_stage_resids[[endo_var]] <- resid(model)
    
    if(save_fs){
      
      fs_mods[[i]] <- model
      
    }
    
  }
  
  # Step 2: Construct second-stage data matrix and formula
  resids_mat <- as.data.frame(first_stage_resids)
  colnames(resids_mat) <- paste0(endogenous,"_resid")
  data_with_resids <- cbind(data,resids_mat)
  
  formula_str <- paste0(outcome," ~ ")
  formula_str <- paste0(formula_str, paste(endogenous, collapse = " + "))
  formula_str <- paste0(formula_str, " + ", paste(paste0(endogenous,"_resid"), collapse = " + "))
  
  if (!is.null(unlist(controls))) {
    
    formula_str <- paste0(formula_str, " + ", paste(controls[[1]], collapse = " + "))
    
  }
  
  # Estimate second stage
  second_stage_lm <- gam(as.formula(formula_str), data = data_with_resids)
  
  # Save and return
  if(save_fs){
    
    out <- list(second_stage_lm = second_stage_lm,
                first_stages = fs_mods)
    
  }else{
    
    out <- second_stage_lm
    
  }
  
  return(out)
  
}

################################################################################
#       Non-linear first stage two stage residual inclusion; bootstrapper      #
################################################################################

# Controls is a list of control variables, the first being for the second stage
#  and the subsequent entries being for each of the first stage regressions.

nlfstsri_boot <- function(data, n_boot, outcome, endogenous, exogenous, controls = NULL, save_fs = TRUE, cores = parallel::detectCores()){
  
  initialization <- nlfstsri_core(data, outcome, endogenous, exogenous, controls, save_fs)
  
  plan(multisession, workers = cores)
  
  with_progress({
    p <- progressor(along = 1:n_boot)
    out <- future_lapply(1:n_boot, function(x) {
      p()
      nlfstsri_core(data[sample(1:nrow(data),replace = TRUE),],
                    outcome,
                    endogenous,
                    exogenous,
                    controls,
                    save_fs = FALSE)$coefficients
    }, future.seed = TRUE)
  })
  
  plan(sequential)
  
  out <- list(initialization = initialization,
              second_stage_boots = do.call("rbind",out))

  return(out)

}

################################################################################
#          Diagnostic: First Stage Precision Loss from Non-Linear Fits         #
################################################################################

system_kappa <- function(nlfstsri_output){
  
  lapply(nlfstsri_output$initialization$first_stages, function(first_stage){
    
    all_terms <- predict(first_stage, type = "terms")
    endo_term <- all_terms[,grepl("te(",colnames(all_terms),fixed = TRUE)][,1]
    term_name <- grep("te(",colnames(all_terms),fixed = TRUE,value = TRUE)[1]
    endo_term <- as.data.frame(endo_term)
    colnames(endo_term) <- names(first_stage$model)[1]
    endo_term
    
  }) -> first_stage_fits
  
  Sigma <- cor(do.call("cbind",first_stage_fits))
  
  kappa_val <- log10(kappa(Sigma, exact = TRUE))
  
  return(kappa_val)
  
}


################################################################################
#           Diagnostic: First Stage Correlations from Non-Linear Fits          #
################################################################################

system_corrs <- function(nlfstsri_output){
  
  lapply(nlfstsri_output$initialization$first_stages, function(first_stage){
    
    all_terms <- predict(first_stage, type = "terms")
    endo_term <- all_terms[,grepl("te(",colnames(all_terms),fixed = TRUE)][,1]
    term_name <- grep("te(",colnames(all_terms),fixed = TRUE,value = TRUE)[1]
    endo_term <- as.data.frame(endo_term)
    colnames(endo_term) <- names(first_stage$model)[1]
    endo_term
    
  }) -> first_stage_fits
  
  Sigma <- cor(do.call("cbind",first_stage_fits))
  
  return(Sigma)
  
}

# ################################################################################
# #                               Simulated Example                              #
# ################################################################################
# 
# #####################################
# # Source of Exogenous Variation
# nobs <- 10000
# z <- rnorm(nobs)
# 
# #####################################  
# # 3 x variables, one unobserved and two a function of z
# covm <- matrix(c(1,0.5,0.5,
#                  0.5,1,-0.5,
#                  0.5,-0.5,1),nrow=3)
# 
# X_mat <- data.frame(mvrnorm(nobs,rep(0,3),covm,empirical = T),z)
# 
# z1 <- z 
# z2 <- z + z^2
# 
# X_mat$X1 <- X_mat$X1 + z1
# X_mat$X2 <- X_mat$X2 + z2
# 
# X_mat$y <- with(X_mat,X1 + X2 +X3 + rnorm(nobs,0,3))
# 
# #######################################
# # Run
# 
# nlfstsri_boot(data = as.data.frame(X_mat),
#               n_boot = 1000,
#               outcome = "y",
#               endogenous = c("X1","X2"),
#               exogenous = "z",
#               controls = NULL) -> test_run
# 
# #######################################
# # Diagnose; lower is better. Ideally below 1.
# 
# system_kappa(test_run)
# 
# #######################################
# # Example post-processing
# 
# # Point estimates and non-bootstrapped results
# summary(test_run$initialization$second_stage_lm)
# 
# # Simple standard deviation of bootstrap
# apply(test_run$second_stage_boots, 2, FUN = function(x){sd(x)})
# 
# # 95% quantile intervals
# apply(test_run$second_stage_boots, 2, FUN = function(x){quantile(x, c(0.025, 0.975))})


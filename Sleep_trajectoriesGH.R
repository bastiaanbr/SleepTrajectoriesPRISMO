#------------- SCRIPT FOR THE TRAJECTORY ANALYSES --------------#
#
# Description:  This script calculates and plots the trajectories of SCL-90 insomnia symptoms scale for the post-deployment time points
#
# Authors:      B.Bruinsma
# Date:         January 2024 - edited December 2024 for Github, edit may 2025 for revision manuscript
# Version:      1.3
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                          Settings & Dependencies
#------------------------------------------------------------------------------#
setwd("Folder location")

# Define path of the Rproject to get and save files
save_loc = "Folder location"

# data import
library(haven)

# data manipulation
library(tidyverse) # includes dply, tidyr and ggplot2

# visualizations
library(Cairo)

# trajectory analysis
library(tidySEM) # for SEM/ trajectory analysis
library(OpenMx) # Required for tidySEM

#------------------------------------------------------------------------------#
#                              Data Collection
#------------------------------------------------------------------------------#
# Load the RData file
load("LocationToFile.RData")


#------------------------------------------------------------------------------#
#                       Trajectory analysis SCL Sleep
#------------------------------------------------------------------------------#

# Set this condition to TRUE if you want to re-run the LCGA. It will take a LONG time!
 run_everything <- TRUE
 
 if(run_everything){
   
  set.seed(123456)
  #dat[["id"]] <- NULL
  names(Scl_sleep_boxcox) <- paste0("scl", 1:6)
  Scl_sleep_boxcox <- data.frame(Scl_sleep_boxcox)
  write.csv(Scl_sleep_boxcox, paste("lcga_dat_scl.csv", sep = ""), row.names = FALSE)

  # Linear models
  #----------------#
  set.seed(12345)
  res_step <- mx_growth_mixture(
    model =
      "i =~ 1*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
  step =~ 0*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
  s =~ 0*scl1 + 0*scl2 + 1*scl3 + 2*scl4 + 4*scl5 + 20*scl6
  scl1 ~~ vscl1*scl1
  scl2 ~~ vscl2*scl2
  scl3 ~~ vscl3*scl3
  scl4 ~~ vscl4*scl4
  scl5 ~~ vscl5*scl5
  scl6 ~~ vscl6*scl6
  i ~~ 0*i
  step ~~ 0*step
  s ~~ 0*s

  i ~~ 0*step
  i ~~ 0*s

  step ~~ 0*s", classes = 1:7,
    data = Scl_sleep_boxcox)

  # In case of convergence problems in the first model:
  #res_postDep[[4]] <- mxTryHardWideSearch(res_step[[4]], extraTries = 300)

  saveRDS(res_postDep, paste0(paste("res_post_deployment", sep = ""), Sys.Date(), ".RData"))

  # fit statistics
  res <- c(res_post_deployment)
  class(res) <- c("mixture_list", "list")
  names(res) <- c(paste0("step", 1:length(res_step)))
  tab_fit <- table_fit(res)
  write.csv(tab_fit, paste("tab_fit_res.csv", sep = ""), row.names = FALSE)

 knitr::kable(tab_fit, digits = 2, caption = "Fit of LCGA models")

 # Quadratic models
 #------------------#
 set.seed(12345)
 res_step_quadratic <- mx_growth_mixture(
   model = "
    # Intercept factor (i)
    i =~ 1*scl1 + 1*scl2 + 1*scl3 + 1*scl4 + 1*scl5 + 1*scl6

    # Linear slope factor (s)
    s =~ 0*scl1 + 0*scl2 + 1*scl3 + 2*scl4 + 4*scl5 + 20*scl6

    # Quadratic slope factor (q)
    q =~ 0*scl1 + 0*scl2 + 1*scl3 + 4*scl4 + 16*scl5 + 400*scl6

    # Step function: 0 at scl1, 1 from scl2 onwards
    step =~ 0*scl1 + 1*scl2 + 1*scl3 + 1*scl4 + 1*scl5 + 1*scl6

    # Residual variances
    scl1 ~~ vscl1*scl1
    scl2 ~~ vscl2*scl2
    scl3 ~~ vscl3*scl3
    scl4 ~~ vscl4*scl4
    scl5 ~~ vscl5*scl5
    scl6 ~~ vscl6*scl6

    # Variances and covariances among the growth factors (LCGA: variances fixed to 0)
    i ~~ 0*i
    s ~~ 0*s
    q ~~ 0*q
    step ~~ 0*step

    # Covariances (fixed to 0 for LCGA)
    i ~~ 0*s
    i ~~ 0*q
    i ~~ 0*step
    s ~~ 0*q
    s ~~ 0*step
    q ~~ 0*step
  ",
   classes = 1:7,
   data = Scl_sleep_boxcox
 )
 
 # Save the results
 saveRDS(res_step_quadratic, paste0("res_step_quadratic_", Sys.Date(), ".RData"))

# Fit statistics
   res <- c(res_quadratic)
   class(res) <- c("mixture_list", "list")
   names(res) <- c(paste0("step", 1:length(res_quadratic)))
   tab_fit <- table_fit(res)
   knitr::kable(tab_fit, digits = 2, caption = "Fit of quadratic LCGA models")
   write.csv(tab_fit, paste("filename.csv", sep = ""), row.names = FALSE)



# Figure for 2-class model
#-----------------------------#
# 
p <- tidySEM::plot_growth(x=res_step[[5]], rawdata = FALSE, alpha_range = c(0, .05))
brks <- seq(0, 1, length.out = 5)
labs <- round(invbc(scales::rescale(brks, from = c(0, 1), to = out$rng_bc), lambda))
p <- p + scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = labs) + ylab("SCL (rescaled from Box-Cox)")
ggsave(paste("filename.svg", sep = ""), p, device = "svg", width = 210, height = 120, units = "mm")
print(p)



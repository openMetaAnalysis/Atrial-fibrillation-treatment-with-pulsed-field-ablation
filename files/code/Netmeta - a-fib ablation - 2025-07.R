# This file is available at http://openmetaanalysis.github.io/afib-rx-pulsed-field
# Author:rbadgett@kumc.edu
# Permissions:
#* Code GNU GPLv3 https://choosealicense.com/licenses/gpl-3.0/
#* Images CC BY-NC-SA 4.0 https://creativecommons.org/licenses/by-nc-sa/4.0/
# Optimized for coding with R Studio document outline view
# Last edited 2025-07-01

#== Startup ======
library(tcltk) # For interactions and troubleshooting, part of base package so no install needed.
#* Set working directory -----
if (Sys.getenv("RSTUDIO") != "1"){
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])  
  script_path <- dirname(script_path)
  setwd(script_path)
}else{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
getwd()

# Packages/libraries -----------
library(netmeta)
library(meta)
library(dplyr)
library (grid)

# Functions ------
function_plot_print <- function (plotname, plotwidth, plotheight, imagetype = "png") {
  
  #plotname <- gsub("[:\\s\n?!']", "", plotname)
  plotname <- gsub(":|\\s|\\n|\\?|\\!|\\'", "", plotname)
  
  current.date <- as.character(strftime(Sys.time(), format="%Y-%m-%d", tz="", usetz=FALSE))
  
  rstudioapi::savePlotAsImage(
    paste(plotname, ' -- ', current.date, '.', imagetype, sep=''),
    format = imagetype, width = plotwidth, height = plotheight)
}

# _________________________________----------

# Data creation ===============================

# Two arms using new RCTs only ------
data_afib_ablation_new <- read.table(textConnection('
studlab,  event1, n1, event2, n2, event3, n3, treat1, treat2, treat3
"Advent (Reddy), 2023",   51,   305,    48, 302, NA,    NA,   "Pulsed field", "Cryoablation", NA
"Single Shot (Reichlin), 2025",   39,   105,    53, 106, NA,    NA,   "Pulsed field", "Cryoablation", NA
'), header=TRUE, sep=",",strip.white=TRUE)

meta_afib_ablation_new <- meta::metabin(event1, n1 , event2, n2, sm = "RR", common = FALSE, data = data_afib_ablation_new)
forest(meta_afib_ablation_new)

# Network meta-analysis ----------------

#* afib_ablation_new pairwise ------
meta_afib_ablation_new$TE.random
meta_afib_ablation_new$seTE.random

pw_afib_ablation_new <- pairwise(list(treat1, treat2, treat3),
                                 event = list(event1, event2, event3),
                                 n     = list(n1,     n2,     n3),
                                 studlab = studlab,
                                 data  = data_afib_ablation_new,
                                 sm    = "RR")         # ← use RR, not OR, for consistency

#* afib_ablation_turagam pairwise ------
# Published RR from Turagam
rr    <- 0.60
ci    <- c(0.50, 0.72)
logRR <- log(rr)
se    <- (log(ci[2]) - log(ci[1])) / (2 * 1.96)

data_afib_ablation_new <- read.table(textConnection('
studlab,  event1, n1, event2, n2, event3, n3, treat1, treat2, treat3
"Turagam, 2021",   NA,   NA,    NA, NA, NA,    NA,   "Cryoablation", "Medication", NA
'), header=TRUE, sep=",",strip.white=TRUE)

pw_afib_ablation_turagam <- data.frame(
  studlab   = "Turagam 2021",
  treat1    = "Cryoablation",
  treat2    = "Medication",
  TE     = log(rr),
  seTE   = (log(ci[2]) - log(ci[1])) / (2*1.96),
  event1 = NA,
  n1 = NA,
  event2 = NA,
  n2 = NA
)

pw_all <- bind_rows(pw_afib_ablation_new, pw_afib_ablation_turagam)   # columns unmatched → NA

net_all <- netmeta(
  pw_all,
  common          = FALSE,        # random-effects
  reference.group = "Medication",
  sm              = "RR"
)

meta::forest(net_all)

# Three arms with Turagam -------
#* Cryotherapy meta-analysis from Turagam -------
rr    <- 0.60
ci    <- c(0.50, 0.72)
logRR <- log(rr)
se    <- (log(ci[2]) - log(ci[1])) / (2 * 1.96)

# Correct format: Two rows, one per arm
turagam_cryo <- data.frame(
  studlab = rep("Turagam 2021", 2),
  treat   = c("Cryoablation", "Medication"),
  TE      = c(logRR, 0),      # Reference treatment (Medication) has logRR=0
  seTE    = c(se, NA)         # Reference treatment has seTE=NA
)

# Now create the pairwise object
pw_turagam <- pairwise(
  studlab = studlab,
  treat   = treat,
  TE      = TE,
  seTE    = seTE,
  sm      = "RR",
  data    = turagam_cryo
)

pw_all <- bind_rows(pw_afib_new, pw_turagam)   # columns unmatched → NA
# class(pw_all) <- c("pairwise", class(pw_all))  # restore the label

net_all <- netmeta(
  pw_all,
  common          = FALSE,        # random-effects
  reference.group = "Medication",
  sm              = "RR"
)

meta::forest(net_all)
grid.text('Network meta-analysis', 0.5, 0.95, gp=gpar(cex=1.4))
grid.text('Notes:', 0.08, 0.25, hjust=0, gp=gpar(cex=1, font=2))
grid.text("Medication data from Turagam, 2021 PMID 33909022", 
          0.08, 0.20, hjust=0, gp=gpar(cex=0.9))
grid.text("Results include mixed methods of monitoring for dysrhythmia.", 
          0.08, 0.15, hjust=0, gp=gpar(cex=0.9)) 
rights <- paste0("Copyleft: Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).\nVersion: ", Sys.Date(), ". rbadgett@kumc.edu")
grid.text(rights, 1, 0.06, hjust=1, gp=gpar(cex=0.75))  #grid.text(Footer, 0.10, 0.075, hjust = 0, gp=gpar(cex=1))

function_plot_print("forest plot of network analysis", 600,300, "png")

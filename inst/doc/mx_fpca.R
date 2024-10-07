## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  figs.out = "../figures"
)
invisible(suppressPackageStartupMessages(library(tidyverse)))

## ----setup--------------------------------------------------------------------
library(mxfda)
library(tidyverse)
library(ggpubr)

## -----------------------------------------------------------------------------
data("ovarian_FDA")
ovarian_FDA

## -----------------------------------------------------------------------------
plot(ovarian_FDA, y = "fundiff", what = "uni g", sampleID = "patient_id") +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  theme_minimal() 
  

## -----------------------------------------------------------------------------
ovarian_FDA <- run_fpca(ovarian_FDA, 
                        metric = "uni g", 
                        r = "r", 
                        value = "fundiff",
                        pve = .95)


## -----------------------------------------------------------------------------
summary(ovarian_FDA)

## ----eval = FALSE-------------------------------------------------------------
#  ovarian_FDA@functional_pca

## ----fpc_plots, fig.width = 10------------------------------------------------

p1 = plot(ovarian_FDA, what = 'uni g fpca', pc_choice = 1)
p2 = plot(ovarian_FDA, what = 'uni g fpca', pc_choice = 2)

ggarrange(p1, p2, nrow = 1, ncol = 2)

## ----refund.shiny, eval = FALSE-----------------------------------------------
#  G_fpca = extract_fpca_object(ovarian_FDA,
#                               what = "uni g fpca")
#  
#  library(refund.shiny)
#  plot_shiny(G_fpca)
#  

## -----------------------------------------------------------------------------
data(lung_df)

clinical = lung_df %>%
  select(image_id, patient_id, patientImage_id, gender, age, survival_days, survival_status, stage) %>%
  distinct()

spatial = lung_df %>%
  select(-image_id, -gender, -age, -survival_days, -survival_status, -stage)

mxFDAobject = make_mxfda(metadata = clinical,
                         spatial = spatial,
                         subject_key = "patient_id",
                         sample_key = "patientImage_id"
                         )

mxFDAobject = extract_summary_functions(mxFDAobject, 
                                        extract_func = univariate,
                                        summary_func = Kest,
                                        r_vec = seq(0, 100, by = 1),
                                        edge_correction = "iso",
                                        markvar = "immune",
                                        mark1 = "immune")

## -----------------------------------------------------------------------------
plot(mxFDAobject, y = "fundiff", what = "uni k", sampleID = "patientImage_id") +
  geom_hline(yintercept = 0, color = "red", linetype = 2) 

## -----------------------------------------------------------------------------
mxFDAobject <- run_mfpca(mxFDAobject, 
                         metric = "uni k", 
                         r = "r", 
                         value = "fundiff",
                         pve = .99)

mxFDAobject

## -----------------------------------------------------------------------------
p = plot(mxFDAobject, what = 'uni k mfpca', level1 = 1, level2 = 1)

ggarrange(plotlist = p, nrow = 1, ncol = 2)



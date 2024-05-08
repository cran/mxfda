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

## -----------------------------------------------------------------------------
data(lung_df)

## ----object-------------------------------------------------------------------
clinical = lung_df %>%
  select(image_id, patient_id, patientImage_id, gender, age, survival_days, survival_status, stage) %>%
  distinct()

spatial = lung_df %>%
  select(-image_id, -gender, -age, -survival_days, -survival_status, -stage)

mxFDAobject = make_mxfda(metadata = clinical,
                         spatial = spatial,
                         subject_key = "patient_id",
                         sample_key = "patientImage_id")

## -----------------------------------------------------------------------------
class(mxFDAobject)

## ----univariate_k-------------------------------------------------------------
mxFDAobject = extract_summary_functions(mxFDAobject,
                                        extract_func = univariate,
                                        summary_func = Kest,
                                        r_vec = seq(0, 100, by = 1),
                                        edge_correction = "iso",
                                        markvar = "immune",
                                        mark1 = "immune")

## -----------------------------------------------------------------------------
mxFDAobject@univariate_summaries$Kest

## -----------------------------------------------------------------------------
plot(mxFDAobject, y = "fundiff", what = "uni k") +
  geom_hline(yintercept = 0, color = "red", linetype = 2)

## -----------------------------------------------------------------------------
lung_df = lung_df %>%
  mutate(phenotype = case_when(phenotype_cd8 == "CD8+" ~ "T-cell",
                               phenotype_cd14 == "CD14+" ~ "macrophage",
                               TRUE ~ "other"),
         phenotype = factor(phenotype))


## -----------------------------------------------------------------------------
spatial = lung_df %>%
  select(-image_id, -gender, -age, -survival_days, -survival_status, -stage)

mxFDAobject = make_mxfda(metadata = clinical,
                         spatial = spatial,
                         subject_key = "patient_id",
                         sample_key = "patientImage_id")

## -----------------------------------------------------------------------------
mxFDAobject = extract_summary_functions(mxFDAobject,
                summary_func = Gcross,
                extract_func = bivariate,
                r_vec = seq(0, 100, by = 1),
                edge_correction = "rs",
                markvar = "phenotype",
                mark1 = "T-cell",
                mark2 = "macrophage")

## -----------------------------------------------------------------------------
plot(mxFDAobject, y = "fundiff", what = "bi g") +
  geom_hline(yintercept = 0, color = "red", linetype = 2)

## ----summary------------------------------------------------------------------
mxFDAobject

## -----------------------------------------------------------------------------
#Step 1
spatialTIME_spatial_df = spatial %>% 
  select(-phenotype) %>%
  mutate(across(phenotype_ck:phenotype_cd4, ~ ifelse(grepl("\\+", .x), 1, 0))) %>%
  relocate(patientImage_id, .before = 1)

#Step 2
cell_types = colnames(spatialTIME_spatial_df) %>% grep("phenotype", ., value = TRUE)

#Step 3
spatial_list = split(spatialTIME_spatial_df, spatial$patientImage_id)

#Step 4
summary_data = lapply(spatial_list, function(df){
  df %>%
    #group by sample ID to maintain ID column
    group_by(patient_id, patientImage_id) %>%
    #find number of positive
    reframe(across(!!cell_types, ~ sum(.x)),
              `Total Cells` = n()) %>%
    #calculate proportion
    mutate(across(!!cell_types, ~.x/`Total Cells` * 100, .names = "{.col} %"))
}) %>%
  #bind the rows together
  do.call(bind_rows, .)

## -----------------------------------------------------------------------------

library(spatialTIME)

#make mif
mif = create_mif(clinical_data = clinical,
                 sample_data = summary_data,
                 spatial_list = spatial_list[1:50],
                 patient_id = "patient_id",
                 sample_id = "patientImage_id")

## -----------------------------------------------------------------------------
mif = NN_G(mif, mnames = cell_types[c(2, 6)],
           r_range = 0:100, num_permutations = 10, 
           edge_correction = "rs", keep_perm_dis = FALSE,
           workers = 1, overwrite = TRUE, xloc = "x", yloc = "y")

## -----------------------------------------------------------------------------
mif$derived$univariate_NN %>%
    ggplot() +
    geom_line(aes(x = r, y = `Degree of Clustering Permutation`, color = patientImage_id), alpha = 0.4) +
    facet_grid(~Marker) +
  theme(legend.position = "none")

## -----------------------------------------------------------------------------
uni_g = mif$derived$univariate_NN %>%
  filter(grepl("cd8", Marker))

## -----------------------------------------------------------------------------
#make mxFDA object 
mxFDA_spatialTIME = make_mxfda(metadata = clinical,
                               spatial = NULL,
                               subject_key = "patient_id",
                               sample_key = "patientImage_id")
#add summary data
mxFDA_spatialTIME = add_summary_function(mxFDAobject,
                                         summary_function_data = uni_g,
                                         metric = "uni g")

## -----------------------------------------------------------------------------
plot(mxFDA_spatialTIME, y = "Degree of Clustering Permutation", what = "uni g")


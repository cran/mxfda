## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>",
  fig.width = 8
)
invisible(suppressPackageStartupMessages(library(tidyverse)))

## ----setup--------------------------------------------------------------------
library(mxfda)
library(tidyverse)
library(ggpubr)
library(broom)

## ----load_data----------------------------------------------------------------
data("ovarian_FDA")

## -----------------------------------------------------------------------------
plot(ovarian_FDA, y = "fundiff", what = "uni g") +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  theme_minimal() +
  ggtitle("Nearest neighbor G-functions for immune cells")
  

## -----------------------------------------------------------------------------
ovarian_FDA <- run_fpca(ovarian_FDA, metric = "uni g", r = "r", value = "fundiff",
                        lightweight = TRUE,
                        pve = .99)
ovarian_FDA

## ----fig.height = 8-----------------------------------------------------------
Gdf_fpc = extract_fpca_scores(ovarian_FDA, 'uni g fpca')

p1 = Gdf_fpc %>%
  mutate(event = factor(event, levels = 0:1, labels = c("censored", "died"))) %>%
  ggplot(aes(fpc1, survival_time, color = event)) +
  geom_point() + 
  labs(y = "survival time (days)", title = "fpc1") +
  theme(legend.position = c(.5, .7))

p2 = Gdf_fpc %>%
  mutate(event = factor(event, levels = 0:1, labels = c("censored", "died"))) %>%
  ggplot(aes(fpc2, survival_time, color = event)) +
  geom_point() + 
  labs(y = "survival time (days)", title = "fpc2") +
  theme(legend.position = "none")


ggarrange(p1, p2, nrow = 1, ncol = 2)

## -----------------------------------------------------------------------------
library(survival)
phmod_fpc = coxph(Surv(survival_time, event) ~ fpc1 + fpc2 + fpc3 + fpc4 + age, 
              data = Gdf_fpc)

## -----------------------------------------------------------------------------
tidy(phmod_fpc, exp = TRUE, conf.int = TRUE) %>%
  mutate(p.value = format.pval(p.value, digits = 1)) %>%
  select(term, hazard_ratio = estimate, conf.low, conf.high, p = p.value) %>%
  knitr::kable(digits = 2)

## -----------------------------------------------------------------------------
ovarian_FDA = run_fcm(ovarian_FDA, model_name = "fit_lfcm",
                      formula = survival_time ~ age, event = "event",
                      metric = "uni g", r = "r", value = "fundiff",
                      afcm = FALSE)

## -----------------------------------------------------------------------------
class(extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_lfcm'))

## -----------------------------------------------------------------------------
summary(extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_lfcm'))

## -----------------------------------------------------------------------------
lfcm_surface = extract_surface(ovarian_FDA, metric = "uni g", model = "fit_lfcm", analysis_vars = c("age"))

plot(lfcm_surface) + ylim(0, 10)

## -----------------------------------------------------------------------------
ovarian_FDA <- run_fcm(ovarian_FDA, model_name = "fit_afcm", 
                       formula = survival_time ~ age, event = "event",
                       metric = "uni g", r = "r", value = "fundiff",
                       afcm = TRUE)

## -----------------------------------------------------------------------------
class(extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_afcm'))

## -----------------------------------------------------------------------------
summary(extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_afcm'))

## ----extract_estimates, fig.with = 12-----------------------------------------
afcm_surface = extract_surface(ovarian_FDA, metric = "uni g", model = "fit_afcm", analysis_vars = c("age"), p = 0.05)

plot(afcm_surface)

## -----------------------------------------------------------------------------
fit_afcm = extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_afcm')
fit_lfcm = extract_model(ovarian_FDA, 'uni g', 'cox', 'fit_lfcm')

c_index = c(
  phmod_fpc$concordance[["concordance"]], 
  extract_c(fit_lfcm, Gdf_fpc$survival_time, Gdf_fpc$event),
  extract_c(fit_afcm, Gdf_fpc$survival_time, Gdf_fpc$event)
)

tibble(model = c("fpc", "lfcm", "afcm"), c_index) %>%
  knitr::kable(digits = 2)


## -----------------------------------------------------------------------------
ovarian_FDA <- run_sofr(ovarian_FDA, 
                        model_name = "fit_sofr_age", 
                        formula = age ~ 1, 
                        metric = "uni g", r = "r", value = "fundiff")

## -----------------------------------------------------------------------------
model = extract_model(ovarian_FDA, 'uni g', type = 'sofr', model_name = 'fit_sofr_age')
plot(model, ylab=expression(paste(beta(t))), xlab="t")

## -----------------------------------------------------------------------------
ovarian_FDA <- run_sofr(ovarian_FDA, 
                        model_name = "fit_sofr_stage", 
                        formula = stage ~ age, 
                        family = "binomial",
                        metric = "uni g", r = "r", value = "fundiff")

## -----------------------------------------------------------------------------
model = extract_model(ovarian_FDA, 'uni g', type = 'sofr', model_name = 'fit_sofr_stage')
plot(model, ylab=expression(paste(beta(t))), xlab="t")


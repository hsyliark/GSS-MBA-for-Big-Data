## Note the MASS package masks select()!
library(tidyverse)
## https://github.com/tidyverse/tibble/issues/395
options(crayon.enabled = FALSE)
## Used for the dataset.
library(survival)
## Used for visualizaiton.
library(survminer)
## Load the Ovarian Cancer Survival Data
data(ovarian)
## Turn into a data_frame
ovarian <- as_data_frame(ovarian)
## Plot
ggsurvplot(survfit(Surv(futime, fustat) ~ 1,
                   data = ovarian),
           risk.table = TRUE,
           break.time.by = 180)

## Fit a Cox model
coxph1 <- coxph(formula = Surv(futime, fustat) ~ pspline(age, df = 4) + factor(resid.ds) +
                  factor(rx) + factor(ecog.ps),
                data    = ovarian)
## Obtain the linear predictor
ovarian$lp <- predict(coxph1, type = "lp")
ovarian

library(survivalROC)
## Define a helper functio nto evaluate at various t
survivalROC_helper <- function(t) {
  survivalROC(Stime        = ovarian$futime,
              status       = ovarian$fustat,
              marker       = ovarian$lp,
              predict.time = t,
              method       = "NNE",
              span = 0.25 * nrow(ovarian)^(-0.20))
}
## Evaluate every 180 days
survivalROC_data <- data_frame(t = 180 * c(1,2,3,4,5,6)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
## Plot
survivalROC_data %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
             mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
  facet_wrap( ~ t) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())

library(risksetROC)
## Define a helper functio nto evaluate at various t
risksetROC_helper <- function(t) {
  risksetROC(Stime        = ovarian$futime,
             status       = ovarian$fustat,
             marker       = ovarian$lp,
             predict.time = t,
             method       = "Cox",
             plot         = FALSE)
}
## Evaluate every 180 days
risksetROC_data <- data_frame(t = 180 * c(1,2,3,4,5,6)) %>%
  mutate(risksetROC = map(t, risksetROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(risksetROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_risksetROC = map(risksetROC, function(obj) {
           ## marker column is too short!
           marker <- c(-Inf, obj[["marker"]], Inf)
           bind_cols(data_frame(marker = marker),
                     as_data_frame(obj[c("TP","FP")]))
         })) %>%
  dplyr::select(-risksetROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
## Plot
risksetROC_data %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_label(data = risksetROC_data %>% dplyr::select(t,auc) %>% unique,
             mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
  facet_wrap( ~ t) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())
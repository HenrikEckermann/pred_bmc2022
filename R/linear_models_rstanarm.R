### README:
# This script was written to run on 
# the HPC Cluster of the Donders Institute
# for Brain, Cognition and Behavior


#####################################################################
##########       packages, settings and data               ##########
#####################################################################

library(tidyverse)
library(projpred)
library(rstanarm)
library(posterior)
library(mice)
library(glmnet)
library(caret)

ncores <- parallel::detectCores(logical = FALSE)
options(mc.cores = ncores)
options(repr.plot.width = 10, repr.plot.height = 10)
load(here::here("rdata/data.Rds"))
source(here::here("R/helper.R"))


#####################################################################
##########               projpred per imp                  ##########
#####################################################################

if (!file.exists(here::here("rdata/pp_allimp_hs.Rds"))) {
  allimp <- map(1:50, function(m) {
    if (!file.exists(glue::glue(here::here("rdata/vs_pp{m}_hs.Rds")))) {
      d <- mutate(dw_imp_all[[m]],
        cortisol_pre_2 = NA,
        cortisol_pre_6 = cortisol_2,
        cortisol_pre_12 = cortisol_6,
        across(contains("season"), function(x) as.numeric(x))
        ) %>%
        # some variables need renaming for succesful transform 
        rename(
          depression_2 = depression_w2,
          depression_6 = depression_w6,
          depression_12 = depression_w12,
          season_2 = season_w2,
          season_6 = season_w6,
          season_12 = season_w12
        ) %>%
        select(-contains("cortisol_diff")) %>%
        # transform to correct format 
        pivot_longer(
            c(
              contains("Collectiontime"),
              contains("postnatalweek"),
              contains("interval_awake"),
              contains("cortisol"),
              contains("sharedactivities_caregiving"),
              contains("season"),
              contains("depression")
            ),
            names_to = c("variable", "week"),
            names_pattern = "(.+)_(\\d+)",
            values_to = "values") %>%
            pivot_wider(
            names_from = variable,
            values_from = values
            ) %>%
            mutate(
            week = factor(week, levels = c("2", "6", "12")),
            season = as.factor(season)) %>%
            # weeks 2 will not be predicted in the longitudinal model 
            filter(week != 2) %>%
            mutate(across(where(is.numeric), function(x) scale(x)[, 1])) 

      formula <-
        cortisol ~
        cortisol_pre +
        postnatalweek +
        interval_awake +
        collectiontime +
        sharedactivities_caregiving +
        father_sens_w6 +
        mother_sens_w6 +
        apl +
        stai +
        iri_pt +
        iri_fs +
        iri_ec +
        iri_pd +
        attach_cat_pren +
        ctq_tot_pren +
        ace_tot_pren +
        eg_tot +
        septi_tot +
        cc_start +
        type_cc_w12 +
        auc_cry +
        season +
        depression +
        sex +
        age +
        edu +
        parity

      # horseshoe prior 
      # Number of regression coefficients:
      D <- 27
      # Prior guess for the number of relevant (i.e., non-zero) regression
      # coefficients:
      p0 <- 5
      # Number of observations:
      N <- nrow(d)
      # Hyperprior scale for tau, the global shrinkage parameter (note that for 
      # Gaussian family, 'rstanarm' will automatically scale this by the residual
      # standard deviation):
      tau0 <- p0 / (D - p0) * 1 / sqrt(N)

      fit <- stan_glm(
        formula = formula,
        family = gaussian(),
        data = d,
        prior = hs(global_scale = tau0),
        chains = 4,
        iter = 4000
      )
      vs <- cv_varsel(
        fit,
        nterms_max = 20,
        seed = 411183
      )
        save(vs, fit, file = glue::glue(here::here("rdata/vs_pp{m}_hs.Rds")))
      } else {
      load(glue::glue(here::here("rdata/vs_pp{m}_hs.Rds")))
    }

    list(vs = vs, refmodel = fit)
  })
  save(allimp, file = here::here("rdata/pp_allimp_hs.Rds"))

 } else {
  load(here::here("rdata/pp_allimp_hs.Rds"))
}

# to inspect variable selection plots 
plots <- map(allimp, function(vs) {
  plot(vs$vs, stats = "rmse")
})

plots 

# Which predictors do most frequently 
# occur as important and how do the models vary in their 
# optimal predictor count?


# less convervative approach (best accuracy score ignoring SD):
solterms1 <- map(allimp, function(vs) {
  result <- summary(vs$vs)$selection
  minsize <- result[result$elpd == max(result$elpd), "size"]
  solterms1 <- filter(result, size <= minsize) %>%
                .$solution_terms %>%
                na.omit()

  solterms1
})


# count predictor frequency
predfreq <- unlist(solterms1) %>% 
  table() %>% 
  as.data.frame() %>%
  rename(Variable = ".") %>%
  mutate(Freq = Freq/50) %>%
  arrange(desc(Freq))
predfreq

# obtain posterior distribution for all predictors with non-zero freq 
prj_pst <- map2_dfr(allimp, solterms, function(ai, st) {
  if (length(st) > 0) {
    prj <- project(ai$refmodel, solution_terms = as.character(st)) 
    prj_mat <- as.matrix(prj)
    prj_tb <- as_tibble(prj_mat)
    prj_tb
  }
})
# summarize posterior 
prj_drws <- as_draws_matrix(prj_pst)
post <- as.data.frame(
  summarize_draws(
    prj_drws, 
    function(x) median(x, na.rm = TRUE), function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  )) %>%
  as_tibble() 
post <- select(post, 
  Variable = variable, 
  Median = "function(x) median(x, na.rm = TRUE)", 
  everything()) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  filter(!Variable %in% c("(Intercept)", "sigma")) 


# for each variable in post$Variable, find out if any subpart of it is in
# predfreq$Variable and if so, grab the Freq 
post$Freq <- rep(0, dim(post)[1])
for (i in seq_along(post$Variable)) {
  for (j in seq_along(predfreq$Variable)) {
    if (str_detect(post$Variable[i], as.character(predfreq$Variable[j]))){
      post[post$Variable == post$Variable[i], "Freq"] <- predfreq[predfreq$Variable == predfreq$Variable[j], "Freq"]
    } 
  }
}

# create table for publication 
projpred_tbl <- arrange(post, desc(Freq)) %>%
  mutate(`95% CI` = glue::glue("{`2.5%`};{`97.5%`}")) %>%
  select(-contains("."), feature = Variable, Median, `95% CI`, Freq) %>%
  rename_vars()
colnames(projpred_tbl) <- str_to_title(colnames(projpred_tbl))
save(projpred_tbl, file = here::here("rdata/projpred_tbl.Rds"))


# how much variation is there in the pred selection between imputed datasets  
pcounts <- map_dbl(solterms1, ~length(.x))
min(pcounts)
max(pcounts)
# although most models suggest a size between 3 and 4 predictors,
# there unfortunately a lot of variation and all model sizes occur 
# leaving us with uncertainty. I think this is yet another reflection 
# of the lack of signal in the data and noisy measurements.
qplot(x = pcounts, geom = "histogram", bins = max(pcounts), color = "white")


# more conservative approach of feature selection 
solterms2<- map(allimp, function(vs) {
  result <- summary(vs$vs)$selection
  minsize <- result[result$elpd == max(result$elpd), "size"]
  if (vs$vs$suggested_size == 0) {
    solterms2 <- NA
    } else {
    solterms2 <- solution_terms(vs$vs)[1:vs$vs$suggested_size]
  }
  solterms2
})

# the same predictors remain, albeit with low frequency. This does not change 
# our conclusion; there seems to be a signal in the data but it is very weak 
unlist(solterms2) %>% 
  table() %>% 
  as.data.frame() %>%
  rename(Variable = ".") %>%
  mutate(Freq = Freq/50) %>%
  arrange(desc(Freq))



# would conclusion change if we did perform the analysis without imputation? 
# father sensitivy is the main variable causing missingness so we must leave 
# this one out:

dnonimp <-  mutate(dw,
        cortisol_pre_2 = NA,
        cortisol_pre_6 = cortisol_2,
        cortisol_pre_12 = cortisol_6,
        across(contains("season"), function(x) as.numeric(x))
        ) %>%
        rename(
          depression_2 = depression_w2,
          depression_6 = depression_w6,
          depression_12 = depression_w12,
          season_2 = season_w2,
          season_6 = season_w6,
          season_12 = season_w12
        ) %>%
        select(-contains("cortisol_diff")) %>%
        pivot_longer(
            c(
              contains("Collectiontime"),
              contains("postnatalweek"),
              contains("interval_awake"),
              contains("cortisol"),
              contains("sharedactivities_caregiving"),
              contains("season"),
              contains("depression")
            ),
            names_to = c("variable", "week"),
            names_pattern = "(.+)_(\\d+)",
            values_to = "values") %>%
            pivot_wider(
            names_from = variable,
            values_from = values
            ) %>%
            mutate(
            week = factor(week, levels = c("2", "6", "12")),
            season = as.factor(season)) %>%
            filter(week != 2) %>%
            mutate(across(where(is.numeric), function(x) scale(x)[, 1])) 


# deselect father sensititvy
formula <-
  cortisol ~
  cortisol_pre +
  postnatalweek +
  interval_awake +
  collectiontime +
  sharedactivities_caregiving +
  # father_sens_w6 +
  mother_sens_w6 +
  apl +
  stai +
  iri_pt +
  iri_fs +
  iri_ec +
  iri_pd +
  attach_cat_pren +
  ctq_tot_pren +
  ace_tot_pren +
  eg_tot +
  septi_tot +
  cc_start +
  type_cc_w12 +
  auc_cry +
  season +
  depression +
  sex +
  age +
  edu +
  parity

# horseshoe prior 
# Number of regression coefficients:
D <- 26
# Prior guess for the number of relevant (i.e., non-zero) regression
# coefficients:
p0 <- 5
# Number of observations:
N <- nrow(d)
# Hyperprior scale for tau, the global shrinkage parameter (note that for 
# Gaussian family, 'rstanarm' will automatically scale this by the residual
# standard deviation):
tau0 <- p0 / (D - p0) * 1 / sqrt(N)

fit <- stan_glm(
  formula = formula,
  family = gaussian(),
  data = select(dnonimp, -father_sens_w6) %>% na.omit(),
  prior = hs(global_scale = tau0),
  chains = 4,
  iter = 4000
)
summary(fit)
vs <- cv_varsel(
  fit,
  nterms_max = 20,
  seed = 411183
)




# given the high uncertainty and the novelty of projpred, I would like to 
# validate these findings using LASSO regression with equal model 
# specification. If then ace does appear, that might be some evidence
set.seed(1)
# lets first determine accuracy of the the algorithm for this task alongside 
# finding our lambda to ultimately obtain feature coefficients accross 
# datasets 


lambda_rsq <- map_dfr(1:50, function(m) {
  # data transformation steps copied from above: 
  d <- mutate(dw_imp_all[[m]],
    cortisol_pre_2 = NA,
    cortisol_pre_6 = cortisol_2,
    cortisol_pre_12 = cortisol_6,
    across(contains("season"), function(x) as.numeric(x))
    ) %>%
    rename(
      depression_2 = depression_w2,
      depression_6 = depression_w6,
      depression_12 = depression_w12,
      season_2 = season_w2,
      season_6 = season_w6,
      season_12 = season_w12
    ) %>%
    select(-contains("cortisol_diff")) %>%
    pivot_longer(
        c(
          contains("Collectiontime"),
          contains("postnatalweek"),
          contains("interval_awake"),
          contains("cortisol"),
          contains("sharedactivities_caregiving"),
          contains("season"),
          contains("depression")
        ),
        names_to = c("variable", "week"),
        names_pattern = "(.+)_(\\d+)",
        values_to = "values") %>%
        pivot_wider(
        names_from = variable,
        values_from = values
        ) %>%
        mutate(
        week = factor(week, levels = c("2", "6", "12")),
        season = as.factor(season)) %>%
        filter(week != 2) %>%
        mutate(across(where(is.numeric), function(x) scale(x)[, 1]))
  
  
  # LASSO 
  y <- d$cortisol
  x <- data.matrix(select(d, -cortisol, -week, -id))
  
  # we only tune lambda 
  grid <- expand.grid(.alpha = 1, .lambda = seq(0, 1, 0.01))
  train_control <- trainControl(method = "LOOCV")
  # train the model
  model <- train(
    x, y, 
    trControl = train_control, 
    method = "glmnet",
    tuneGrid = grid,
    metric = "Rsquared"
  )
  model_metrics <- as_tibble(model$results) %>%
    filter(RMSE == min(RMSE)) 
  list(
    lambda = model_metrics$lambda[1], 
    rsq = model_metrics$Rsquared[1], 
    rmse = model_metrics$RMSE[1])
})
median(lambda_rsq$rsq)
rsqs <- map_dbl(lambda_rsq, ~.x$rsq[1])
rsq <- median(rsqs)
lambda <- median(lambda_rsq$lambda)


# lets get the coefficients using the same lambda across imputed datasets
lasso_allimp <- map(1:50, function(m) {
  d <- mutate(dw_imp_all[[m]],
    cortisol_pre_2 = NA,
    cortisol_pre_6 = cortisol_2,
    cortisol_pre_12 = cortisol_6,
    across(contains("season"), function(x) as.numeric(x))
    ) %>%
    rename(
      depression_2 = depression_w2,
      depression_6 = depression_w6,
      depression_12 = depression_w12,
      season_2 = season_w2,
      season_6 = season_w6,
      season_12 = season_w12
    ) %>%
    select(-contains("cortisol_diff")) %>%
    pivot_longer(
        c(
          contains("Collectiontime"),
          contains("postnatalweek"),
          contains("interval_awake"),
          contains("cortisol"),
          contains("sharedactivities_caregiving"),
          contains("season"),
          contains("depression")
        ),
        names_to = c("variable", "week"),
        names_pattern = "(.+)_(\\d+)",
        values_to = "values") %>%
        pivot_wider(
        names_from = variable,
        values_from = values
        ) %>%
        mutate(
        week = factor(week, levels = c("2", "6", "12")),
        season = as.factor(season)) %>%
        filter(week != 2) %>%
        mutate(across(where(is.numeric), function(x) scale(x)[, 1]))

  y <- d$cortisol
  x <- data.matrix(select(d, -cortisol, -week, -id))
  fit <- glmnet(x, y, alpha = 1, lambda = lambda)
  betas <- as.data.frame(as.matrix(coef(fit)))
  
  list(
    betas = betas
  )
})



# lets count how often each variable has non-zero coefficients 
betas <- map(lasso_allimp, function(listobj) {
  betas <- listobj$betas %>% 
                  rownames_to_column("variable") %>%
                  filter(s0 != 0, variable != "(Intercept)") %>% 
                  .$variable
})
# predictor frequency
nonzerofreq <- unlist(betas) %>% 
  table() %>% 
  as.data.frame() %>%
  rename(variable = ".") %>%
  mutate(Freq = Freq/50) %>%
  arrange(desc(Freq))

# how does optimal predictor number vary for LASSO? 
pcounts <- map_dbl(betas, ~length(.x))
lbins <- max(pcounts)
min(pcounts)
qplot(x = pcounts, geom = "histogram", bins = lbins, color = "white")


# lets calculate the average coefficient size 
betas <- map(lasso_allimp, function(listobj) {
  betas <- listobj$betas %>% 
                  rownames_to_column("variable") %>%
                  filter(variable != "(Intercept)") 
})

# table for publication:
# after rounding, I will filter out those that did not occur at least 50% of 
# the time as non-zero
lasso_tbl <- map_dfr(betas[[1]]$variable, function(var) {
  b <- map_dbl(betas, function(beta) {
    beta[beta$variable == var, "s0"]
  }) %>% mean()
  tibble(variable = var, b = b)}) %>%
  mutate(b = round(b, 3)) %>%
  left_join(nonzerofreq, by = "variable") %>%
  mutate(Freq = ifelse(is.na(Freq), 0, Freq)) %>%
  arrange(desc(Freq), desc(abs(b))) %>%
  select(feature = variable, Beta = b, Freq) %>%
  rename_vars() %>%
  select(Feature = feature, Beta, Freq) %>%
  filter(Freq != 0)
save(lasso_tbl, projpred_tbl, file = here::here("rdata/tbls.Rds"))

lasso_tbl
#####################################################################
##########       packages, data import, settings           ##########
#####################################################################
options(repr.plot.width = 7.5, repr.plot.height = 7.5)
library(tidyverse)


# load ML helper functions
source("https://raw.githubusercontent.com/HenrikEckermann/in_use/master/ml_helper.R")
# helper function for renaming vars specific to this project 
source(here::here("R/helper.R"))
load(here::here("rdata/data.Rds"))


#####################################################################
##########               Random Forrest                    ##########
#####################################################################

# first test code on single imputed dataset     
# create time point specific data sets in wide format

# y = week 2
d0 <- select(
  dw_imp,
  -cortisol_12,
  -cortisol_6,
  -matches("postnatalweek_6|12"),
  -matches("collectiontime_6|12"),
  -matches("interval_awake_6|12"),
  -matches("sharedactivities_caregiving_6|12"),
  -matches("cortisol_diff\\d+"),
  -matches("season_w6|12"),
  -matches("depression_w6|12")
  ) %>%
  rename(cortisol = cortisol_2)


# y = week 6
d1 <- select(
  dw_imp,
  -cortisol_12,
  -matches("postnatalweek_[212]"),
  -matches("collectiontime_[212]"),
  -matches("interval_awake_[212]"),
  -matches("sharedactivities_caregiving_[212]"),
  -matches("cortisol_diff\\d+"),
  -matches("season_w[212]"),
  -matches("depression_w[212]")
  ) %>%
  rename(cortisol_pre = cortisol_2, cortisol = cortisol_6)


#x = week 6, y = week 12
d2 <- select(
  dw_imp,
  -cortisol_2,
  -matches("postnatalweek_[26]"),
  -matches("collectiontime_[26]"),
  -matches("interval_awake_[26]"),
  -matches("sharedactivities_caregiving_[26]"),
  -matches("cortisol_diff\\d+"),
  -matches("season_w[26]"),
  -matches("depression_w[26]")
  ) %>%
  rename(cortisol_pre = cortisol_6, cortisol = cortisol_12)

# x = week 2, y = week 12
d3 <- select(
  dw_imp,
  -cortisol_6,
  -matches("postnatalweek_[26]"),
  -matches("collectiontime_[26]"),
  -matches("interval_awake_[26]"),
  -matches("sharedactivities_caregiving_[26]"),
  -matches("cortisol_diff\\d+"),
  -matches("season_w[26]"),
  -matches("depression_w[26]")
  ) %>%
  rename(cortisol_pre = cortisol_2, cortisol = cortisol_12)



# RF pipeline
if (!file.exists(here::here("rdata/rf2_test.Rds"))) {
  # fit models per timepoints and per imputation method
  result <- map(list(d0, d1, d2, d3), function(d) {


      y <- "cortisol"

      # first base model with only the known predictors
      X_null <- select(
        d,
        contains("cortisol_pre"),
        contains("postnatal"),
        contains("collection"),
        contains("interval_awake")) %>%
        colnames()

      # tune RF hyperparameters
      pars_null <- tune_rf(
        d,
        X_null,
        y,
        regression = TRUE,
        iters = 70,
        iters.warmup = 30,
        ntree = 5000,
        parameters = list(
          replace = FALSE,
          respect.unordered.factors = "order"
        ),
        tune.parameters = c(
          "mtry",
          "min.node.size",
          "sample.fraction"
        ),
        show.info = getOption("mlrMBO.show.info", TRUE)
        )

      # fit model using above hyperparameters
      model_null <- ranger(
        x = select(d, all_of(X_null)),
        y = d[[y]],
        importance = "permutation",
        num.tree = 5000,
        mtry = pars_null$recommended.pars$mtry,
        min.node.size = pars_null$recommended.pars$min.node.size,
        sample.fraction = ifelse(
          pars_null$recommended.pars$sample.fraction < 0.25,
          0.25, pars_null$recommended.pars$sample.fraction)
      )

      # now the full models
      X <- select(d, -id, -cortisol) %>%
        colnames()

      # tune RF hyperparameters
      pars <- tune_rf(
        d,
        X,
        y,
        regression = TRUE,
        iters = 70,
        iters.warmup = 30,
        ntree = 5000,
        parameters = list(
          replace = FALSE,
          respect.unordered.factors = "order"
        ),
        tune.parameters = c(
          "mtry",
          "min.node.size",
          "sample.fraction"
        ),
        show.info = getOption("mlrMBO.show.info", TRUE)
        )

      # fit model using above hyperparameters
      model <- ranger(
        x = select(d, all_of(X)),
        y = d[[y]],
        importance = "permutation",
        num.tree = 5000,
        mtry = pars$recommended.pars$mtry,
        min.node.size = pars$recommended.pars$min.node.size,
        sample.fraction = ifelse(
          pars$recommended.pars$sample.fraction < 0.25,
          0.25, pars$recommended.pars$sample.fraction)
      )


      list(
        model_null = model_null,
        pars_null = pars_null$recommended.pars,
        plot_null = plot_importance(model_null),
        model = model,
        pars = pars$recommended.pars,
        plot = plot_importance(model)
      )
    })
  save(result, file = here::here("rdata/rf2_test.Rds"))
 } else {
  load(file = here::here("rdata/rf2_test.Rds"))
}

result



# the results indicate overfitting at least because of some of our candidate
# predictors. The base models predict better except for y = week 2.
# Nevertheless, we can compute variable importances to evaluate which
# predictors help and which  carry no signal in them. But first lets see how it
# varies between imputations.

impvariation <- map2(1:50, dw_imp_all, function(m, d_imp) {
  if (!file.exists(here::here(glue::glue("rdata/rf2_impvar{m}_out.Rds")))) {

    # create time point specific data sets in wide format
    # predicting at week 2
    d0 <- select(
      d_imp,
      -cortisol_12,
      -cortisol_6,
      -matches("postnatalweek_6|12"),
      -matches("collectiontime_6|12"),
      -matches("interval_awake_6|12"),
      -matches("sharedactivities_caregiving_6|12"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w6|12"),
      -matches("depression_w6|12")
      ) %>%
      rename(cortisol = cortisol_2)

    # change 2 -> 6
    d1 <- select(
      d_imp,
      -cortisol_12,
      -matches("postnatalweek_[212]"),
      -matches("collectiontime_[212]"),
      -matches("interval_awake_[212]"),
      -matches("sharedactivities_caregiving_[212]"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w[212]"),
      -matches("depression_w[212]")
      ) %>%
      rename(cortisol_pre = cortisol_2, cortisol = cortisol_6)


    # change 6 -> 12
    d2 <- select(
      d_imp,
      -cortisol_2,
      -matches("postnatalweek_[26]"),
      -matches("collectiontime_[26]"),
      -matches("interval_awake_[26]"),
      -matches("sharedactivities_caregiving_[26]"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w[26]"),
      -matches("depression_w[26]")
      ) %>%
      rename(cortisol_pre = cortisol_6, cortisol = cortisol_12)

    # change 2 -> 12
    d3 <- select(
      d_imp,
      -cortisol_6,
      -matches("postnatalweek_[26]"),
      -matches("collectiontime_[26]"),
      -matches("interval_awake_[26]"),
      -matches("sharedactivities_caregiving_[26]"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w[26]"),
      -matches("depression_w[26]")
      ) %>%
      rename(cortisol_pre = cortisol_2, cortisol = cortisol_12)


      # fit models per timepoints
      result <- map(list(d0, d1, d2, d3), function(d) {

          y <- "cortisol"
          # first model with only the known predictors
          X_null <- select(
            d,
            contains("cortisol_pre"),
            contains("postnatal"),
            contains("collection"),
            contains("interval_awake")) %>%
            colnames()

          # tune RF hyperparameters
          pars_null <- tune_rf(
            d,
            X_null,
            y,
            regression = TRUE,
            iters = 70,
            iters.warmup = 30,
            ntree = 5000,
            parameters = list(
              replace = FALSE,
              respect.unordered.factors = "order"
            ),
            tune.parameters = c(
              "mtry",
              "min.node.size",
              "sample.fraction"
            ),
            show.info = getOption("mlrMBO.show.info", TRUE)
            )

          # fit model using above hyperparameters
          model_null <- ranger(
            x = select(d, all_of(X_null)),
            y = d[[y]],
            importance = "permutation",
            num.tree = 5000,
            mtry = pars_null$recommended.pars$mtry,
            min.node.size = pars_null$recommended.pars$min.node.size,
            sample.fraction = ifelse(
              pars_null$recommended.pars$sample.fraction < 0.25,
              0.25, pars_null$recommended.pars$sample.fraction)
          )

          # now the full models
          X <- select(d, -id, -cortisol) %>%
            colnames()

          # tune RF hyperparameters
          pars <- tune_rf(
            d,
            X,
            y,
            regression = TRUE,
            iters = 70,
            iters.warmup = 30,
            ntree = 5000,
            parameters = list(
              replace = FALSE,
              respect.unordered.factors = "order"
            ),
            tune.parameters = c(
              "mtry",
              "min.node.size",
              "sample.fraction"
            ),
            show.info = getOption("mlrMBO.show.info", TRUE)
            )

          # fit model using above hyperparameters
          model <- ranger(
            x = select(d, all_of(X)),
            y = d[[y]],
            importance = "permutation",
            num.tree = 5000,
            mtry = pars$recommended.pars$mtry,
            min.node.size = pars$recommended.pars$min.node.size,
            sample.fraction = ifelse(
              pars$recommended.pars$sample.fraction < 0.25,
              0.25, pars$recommended.pars$sample.fraction)
          )


          list(
            model_null = model_null,
            pars_null = pars_null$recommended.pars,
            plot_null = plot_importance(model_null),
            model = model,
            pars = pars$recommended.pars,
            plot = plot_importance(model)
          )
        })
      save(result, file = here::here(glue::glue("rdata/rf2_impvar{m}_out.Rds")))
     } else {
      load(file = here::here(glue::glue("rdata/rf2_impvar{m}_out.Rds")))
    }
    result
})

# summarise accuracy scores 
rsqmedians <- map_dfr(impvariation, function(listobj) {
  map2_dfr(1:4, listobj, function(num, models) {
    rsqmediannull <- models$model_null$r.squared
    rsqmedianfull <- models$model$r.squared
    tibble(num = num, null = rsqmediannull, full = rsqmedianfull)
  })}) %>%
  group_by(num) %>%
  summarise(
    mnull = median(null),
    sdnull = sd(null),
    mfull = median(full),
    sdfull = sd(full)
  )
rsqmedians

# calculate difference between base and full models 
rsqdiff <- map_dfr(impvariation, function(listobj) {
  map2_dfr(1:4, listobj, function(num, models) {
    rsqdiff <- models$model$r.squared - models$model_null$r.squared
    tibble(num = num, rsqdiff = rsqdiff)
  })
})
rsqdiff %>% group_by(num) %>%
  summarise(rsq = median(rsqdiff), sd = sd(rsqdiff))

# accross all the imputed datasets, the rsq difference indicates that the model
# that only has the known covariates fits the data slightly better.
# I think that this is because there are too many noisy variables together
# with a low sample size --> overfitting. Some of the variables might be
# predictive. Therefore, I stick with calculating pvalues of the models after
# which I will calculate variable importances.




#####################################################################
##########      calculate p values for RF models           ##########
#####################################################################
nperms <- 100
nulldist <- map_dfr(1:nperms, function(nperm) {
    map(1:50, function(m) {
      if (!file.exists(here::here(glue::glue("rdata/nulldist2_{nperm}_{m}.Rds")))) {
        d <- dw_imp_all[[m]]

        # create time point specific data sets in wide format
        # predicting at week 2
        d0 <- select(
          d,
          -cortisol_12,
          -cortisol_6,
          -matches("postnatalweek_6|12"),
          -matches("collectiontime_6|12"),
          -matches("interval_awake_6|12"),
          -matches("sharedactivities_caregiving_6|12"),
          -matches("cortisol_diff\\d+"),
          -matches("season_w6|12"),
          -matches("depression_w6|12")
          ) %>%
          rename(cortisol = cortisol_2)

        # change 2 -> 6
        d1 <- select(
          d,
          -cortisol_12,
          -matches("postnatalweek_[212]"),
          -matches("collectiontime_[212]"),
          -matches("interval_awake_[212]"),
          -matches("sharedactivities_caregiving_[212]"),
          -matches("cortisol_diff\\d+"),
          -matches("season_w[212]"),
          -matches("depression_w[212]")
          ) %>%
          rename(cortisol_pre = cortisol_2, cortisol = cortisol_6)


        # change 6 -> 12
        d2 <- select(
          d,
          -cortisol_2,
          -matches("postnatalweek_[26]"),
          -matches("collectiontime_[26]"),
          -matches("interval_awake_[26]"),
          -matches("sharedactivities_caregiving_[26]"),
          -matches("cortisol_diff\\d+"),
          -matches("season_w[26]"),
          -matches("depression_w[26]")
          ) %>%
          rename(cortisol_pre = cortisol_6, cortisol = cortisol_12)

        # change 2 -> 12
        d3 <- select(
          d,
          -cortisol_6,
          -matches("postnatalweek_[26]"),
          -matches("collectiontime_[26]"),
          -matches("interval_awake_[26]"),
          -matches("sharedactivities_caregiving_[26]"),
          -matches("cortisol_diff\\d+"),
          -matches("season_w[26]"),
          -matches("depression_w[26]")
          ) %>%
          rename(cortisol_pre = cortisol_2, cortisol = cortisol_12)


          # fit models per timepoints
          result <- map2_dfr(list(d0, d1, d2, d3), 1:4, function(df, num) {
              d <- df
              y <- "cortisol"
              X <- select(d, -id, -cortisol) %>%
                colnames()

              d[[y]] <- sample(df[[y]], replace = FALSE)

              # tune RF hyperparameters
              pars <- tune_rf(
                d,
                X,
                y,
                regression = TRUE,
                iters = 70,
                iters.warmup = 30,
                ntree = 5000,
                parameters = list(
                  replace = FALSE,
                  respect.unordered.factors = "order"
                ),
                tune.parameters = c(
                  "mtry",
                  "min.node.size",
                  "sample.fraction"
                ),
                show.info = getOption("mlrMBO.show.info", TRUE)
                )

              # fit model using above hyperparameters
              model <- ranger(
                x = select(d, all_of(X)),
                y = d[[y]],
                importance = "permutation",
                num.tree = 5000,
                mtry = pars$recommended.pars$mtry,
                min.node.size = pars$recommended.pars$min.node.size,
                sample.fraction = ifelse(
                  pars$recommended.pars$sample.fraction < 0.25,
                  0.25, pars$recommended.pars$sample.fraction)
              )

              rsq <- model$r.squared


              list(
                num = num,
                rsq = rsq
              )
            })
            save(result, file = here::here(glue::glue("rdata/nulldist2_{nperm}.Rds")))
           } else {
            load(file = here::here(glue::glue("rdata/nulldist2_{nperm}_{m}.Rds")))
          }
          result
    })
})
nulldist_nested <- nulldist %>% group_by(num) %>% nest()
# base models 
map2(rsqmedians$mnull, nulldist_nested$data, function(rsq, dist) {
  mean(rsq <= dist$rsq)
})
map2(rsqmedians$mfull, nulldist_nested$data, function(rsq, dist) {
  mean(rsq <= dist$rsq)
})


#####################################################################
##########            Caculate pvalues for features        ##########
#####################################################################


featpm <- map2(1:50, dw_imp_all, function(m, d_imp) {
  if (!file.exists(here::here(glue::glue("rdata/rf2_altmann_m{m}_out.Rds")))) {

    # create time point specific data sets in wide format
    # predicting at week 2
    d0 <- select(
      d_imp,
      -cortisol_12,
      -cortisol_6,
      -matches("postnatalweek_6|12"),
      -matches("collectiontime_6|12"),
      -matches("interval_awake_6|12"),
      -matches("sharedactivities_caregiving_6|12"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w6|12"),
      -matches("depression_w6|12")
      ) %>%
      rename(cortisol = cortisol_2)

    # change 2 -> 6
    d1 <- select(
      d_imp,
      -cortisol_12,
      -matches("postnatalweek_[212]"),
      -matches("collectiontime_[212]"),
      -matches("interval_awake_[212]"),
      -matches("sharedactivities_caregiving_[212]"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w[212]"),
      -matches("depression_w[212]")
      ) %>%
      rename(cortisol_pre = cortisol_2, cortisol = cortisol_6)


    # change 6 -> 12
    d2 <- select(
      d_imp,
      -cortisol_2,
      -matches("postnatalweek_[26]"),
      -matches("collectiontime_[26]"),
      -matches("interval_awake_[26]"),
      -matches("sharedactivities_caregiving_[26]"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w[26]"),
      -matches("depression_w[26]")
      ) %>%
      rename(cortisol_pre = cortisol_6, cortisol = cortisol_12)

    # change 2 -> 12
    d3 <- select(
      d_imp,
      -cortisol_6,
      -matches("postnatalweek_[26]"),
      -matches("collectiontime_[26]"),
      -matches("interval_awake_[26]"),
      -matches("sharedactivities_caregiving_[26]"),
      -matches("cortisol_diff\\d+"),
      -matches("season_w[26]"),
      -matches("depression_w[26]")
      ) %>%
      rename(cortisol_pre = cortisol_2, cortisol = cortisol_12)


      # fit models per timepoints and per imputation 
      featp <- map2(list(d0, d1, d2, d3), 1:4, function(d, num) {

          model <- impvariation[[m]][[num]]$model
          pars <- impvariation[[m]][[num]]$pars

          pimp <- importance_pvalues(
            model,
            method = "altmann",
            num.permutations = 1000,
            data = select(d, -id),
            formula = cortisol ~ .,
          )
        pimp
        })
      save(featp, file = here::here(glue::glue("rdata/rf2_altmann_m{m}_out.Rds")))
     } else {
      load(file = here::here(glue::glue("rdata/rf2_altmann_m{m}_out.Rds")))
    }
    featp
})


# find median importance value and corresponding p value; 
# create table for paper 
featp_average <- map2_dfr(1:50, featpm, function(m, x) {
  map2_dfr(1:4, x, function(num, d) {
  as.data.frame(d) %>%
  rownames_to_column("feature") %>%
  mutate(num = num, m = m)
  })
})

nums <- 1:3
# rename the variables
table_renamed <- map(nums, function(timepoint) {
  tbl <- featp_average %>% group_by(feature, num) %>% 
  # since median of 50 value will be averaged I find the value closest to that 
  # average value to get the correct pvalue, if there are several, i just pick 
  # randomly 
    mutate(
      median = median(importance),
      mediandist = abs(median - importance),
      mdist = abs(m - runif(1, 50, n = 1))
    ) %>%
    filter(mediandist == min(mediandist)) %>%
    filter(mdist == min(mdist)) %>% # pick one out of duplicates randomly
    ungroup() %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    arrange(num, pvalue, desc(importance)) %>%
    select(feature, importance, pvalue, num) %>%
    filter(num == timepoint)  %>%
    rename_vars() 
  colnames(tbl) <- str_to_title(colnames(tbl))
  tbl
})

save(table_renamed, file = here::here("rdata/tables.Rds"))


#####################################################################
##########             README, packages, settings          ##########
#####################################################################
# In this script, data will be imported, renamed, described, imputed and put in
# formats suitable for the main analyses.

options(repr.plot.width = 15, repr.plot.height = 15)
library(tidyverse)
library(mice)
library(patchwork)



#####################################################################
##########                     import                      ##########
#####################################################################
list.files(here::here("data"))

# the csv import works best for the time vars
temp <- read_csv(here::here("data/Breast milk_paper 2021_varofinterest.csv"))
colnames(temp)
# rename, create time var, reformat
d <- rename(
  temp,
  id = ID,
  sex = Baby_Sex,
  age = MaternalAge,
  edu = EducationalLevel,
  parity = Parity,
  season_w2 = Season2w,
  season_w6 = Season6w,
  season_w12 = Season12w,
  date_w2 = Date_w2,
  date_w6 = Date_w6,
  date_w12 = Date_w12,
  depression_w2 = Depression2w_post,
  depression_w6 = Depression6w_post,
  depression_w12 = Depression12w_post,
  mother_sens_w6 = Mother_sensitivity_6w,
  apl_freq = APLpostFrequency,
  apl = APLPostTotalscore,
  stai = STAIPostTotalscore,
  iri_pt = IRI_PT,
  iri_fs = IRI_FS,
  iri_ec = IRI_EC,
  iri_pd = IRI_PD,
  attach_preoc_pren = Attachment_preoccupied_prenatal,
  attach_dism_pren = Attachment_dismissing_prenatal,
  attach_unres_pren = Attachment_unresolved_prenatal,
  attach_sec_pren = Attachment4_secure_prenatal,
  attach_cat_pren = Attachment_category_prenatal,
  ctq_pn_pren = CTQ_PN_prenatal,
  ctq_en_pren = CTQ_EN_prenatal,
  ctq_pa_pren = CTQ_PA_prenatal,
  ctq_sa_pren = CTQ_SA_prenatal,
  ctq_ea_pren = CTQ_EA_prenatal,
  ctq_tot_pren = CTQ_total_prenatal,
  ace_tot_pren = ACE_tot_prenatal,
  eg_tot = EG_total,
  septi_tot = SEPTI_TOT,
  nurt_septi = Nurturance_SEPTI,
  disc_septi = Discipline_SEPTI,
  play_sept = Play_SEPTI,
  rout_sept = Routine_SEPTI,
  cc_start = Childcare_start,
  cc_start_cat = Childcare_start_dychotonomous,
  type_cc_w12 = Type_daycare_12w,
  auc_cry = AUCg_cry,
  father_sens_w6 = Father_sensitivity_6w,
  auc_cry_m1_m2 = AUCg_cry_m1_m2,
  auc_cry_m2_m3 = AUCg_cry_m2_m3,
  cry_w2 = CryDuration2weeks,
  cry_w6 = CryDuration6weeks,
  cry_w12 = CryDuration12weeks) %>%
  mutate(
    across(contains("Collectiontime"), function(x) {
      parse_time(as.character(x), format = "%H:%M:%S", na = c("", "NA"))
    }),
    across(contains("nterval_awake"), function(x) {
      parse_time(as.character(x), format = "%H:%M:%S", na = c("", "NA"))
    }),
    across(contains("nterval_awake"), function(x) {
      (as.numeric(x)/60)/60  # convert to numeric hours
    }),
    across(contains("Collectiontime"), function(x) {
      (as.numeric(x)/60)/60  # convert to numeric hours
    }),
    across(contains("sharedactivities_caregiving"), function(x) {
      x - 3 # scale so that 0 means equal caregiving
    }),
    across(contains("ortisol"), function(x) {
      log(x) # log cortisol will be needed for regression
    }),
    # the current categories might not suitable for our research question
    cc_start_cat_2 = ifelse(cc_start <= 1, 1, 0),
    cc_start_cat_6 = ifelse(cc_start <= 2, 1, 0),
    cc_start_cat_12 = ifelse(cc_start <= 3, 1, 0),
    attach_cat_pren = as.factor(attach_cat_pren),
    type_cc_w12 = as.factor(type_cc_w12)
  ) %>%
    pivot_longer(
      c(
        contains("Cortisol"),
        contains("Postnatalweek"),
        contains("Collectiontime"),
        matches("^[Ii]nterval_awake_w\\d+$"),
        contains("Sharedactivities"),
        matches("cry_w\\d+"),
        matches("cc_start_cat_\\d")
      ),
      names_to = c("variable", "week"),
      names_pattern = "(.+).*_w*(\\d+)"
    ) %>%
    mutate(
      variable = str_to_lower(variable)
    ) %>%
    select(-cc_start_cat)




# the following ids have issues for the time parsing:
select(temp, ID, contains("time"), contains("nterval")) %>%
   filter(ID %in% c("018m", "096m", "032m"))

# those 3 time values are impossible. I asked Stefania about this and it turned
# out that those are in the raw data wrong (impossible) entries. Therefore, I
# leave them as NA 


# wide format data
dw <- pivot_wider(
  d,
  names_from = c(variable, week),
  names_glue = "{variable}_{week}",
  values_from = value)
  
  


#####################################################################
##########      explore data        ##########
#####################################################################

count(dw, season_w12)
count(dw, parity)
str(dw)


# correlation
p <- mutate(dw, across(where(is.factor), as.integer)) %>%
  select(-id) %>%
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete", method = "pearson") %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r")


ggplot(p, aes(var1, var2, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
   midpoint = 0, limit = c(-1,1), space = "Lab",
   name="Pearson\nCorrelation") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))


# for imputation: change ids names to integer; same for categorical vars
idswitch <- tibble(id = unique(dw$id), new = 1:length(id))

dw <- dw %>% left_join(idswitch, by = "id") %>%
  select(-id,) %>% select(id = new, everything()) %>%
  mutate(
    attach_cat_pren =
      ifelse(attach_cat_pren == "A", 1,
        ifelse(attach_cat_pren == "B", 2,
          ifelse(attach_cat_pren == "C", 3,
            ifelse(attach_cat_pren == "D", 4, NA)))),
    across(c(id, sex, parity, contains("season"), contains("cat")), function(x) as.factor(x))
  )

# visualize each variable for inspection
# first categorical:
select_if(dw, is.factor) %>%
  colnames() %>%
  map(function(x) {
    ggplot(dw, aes_string(x)) +
      geom_bar()
  })


# then continuous:
select_if(dw, is.numeric) %>%
  colnames() %>%
  map(function(var) {
    ggplot(dw, aes_string(var)) +
      geom_histogram() +
      ggtitle(var)
})

# descriptive statistics and missingness
mlr::summarizeColumns(dw)


#####################################################################
##########                  Imputation                     ##########
#####################################################################

# deselect subscales (we explored those as well without different results)
dw <- select(
  dw,
  -contains("date"),
  -attach_preoc_pren,
  -attach_dism_pren,
  -attach_unres_pren,
  -attach_sec_pren,
  -ctq_pn_pren,
  -ctq_en_pren,
  -ctq_pa_pren,
  -ctq_sa_pren,
  -ctq_ea_pren,
  -nurt_septi,
  -disc_septi,
  -play_sept,
  -rout_sept,
  -auc_cry_m1_m2,
  -auc_cry_m2_m3,
  -matches("cry_\\d"),
  -contains("cc_start_cat"),
  -apl_freq
)
colnames(dw)
# pmm imputation
if(!file.exists(here::here("rdata/dimpall.Rds"))) {
  dimp <- mice(dw, m = 50)
  save(dimp, file = here::here("rdata/dimpall.Rds"))
 } else {
  load(here::here("rdata/dimpall.Rds"))
}

# inspect imputations and logged events 
dimp$imp


# to run all analyses first with 1 complete dataset 
dw_imp <- complete(dimp)
# add change in cortisol over time for doublechecking models
dw$cortisol_diff62 <- dw$cortisol_6 - dw$cortisol_2
dw$cortisol_diff126 <- dw$cortisol_12 - dw$cortisol_6
dw_imp$cortisol_diff62 <- dw_imp$cortisol_6 - dw_imp$cortisol_2
dw_imp$cortisol_diff126 <- dw_imp$cortisol_12 - dw_imp$cortisol_6
dw_imp$cortisol_diff122 <- dw_imp$cortisol_12 - dw_imp$cortisol_2

dlpmm <- complete(dimp) %>%
  mutate(across(contains("cc_start_cat"), function(x) as.numeric(x))) %>%
  pivot_longer(
    c(
      contains("Collectiontime"),
      contains("postnatalweek"),
      contains("interval_awake"),
      contains("cortisol"),
      contains("sharedactivities_caregiving"),
      matches("cry_\\d")
    ),
    names_to = c("variable", "week"),
    names_pattern = "(.+)_(\\d+)",
    values_to = "values") %>%
  pivot_wider(
    names_from = variable,
    values_from = values
  ) %>%
  mutate(week = factor(week, levels = c("2", "6", "12")))

# after running this with one imputed dataset, I will want to evaluate 
# using all imputed datasets. For that I create a list
# here with all the datasets
dw_imp_all <- map(1:50, function(m) {
  dw_imp <- complete(dimp, action = m)
  # # add change in cortisol over time for RF models
  # dw$cortisol_diff62 <- dw$cortisol_6 - dw$cortisol_2
  # dw$cortisol_diff126 <- dw$cortisol_12 - dw$cortisol_6
  # dw_imp$cortisol_diff62 <- dw_imp$cortisol_6 - dw_imp$cortisol_2
  # dw_imp$cortisol_diff126 <- dw_imp$cortisol_12 - dw_imp$cortisol_6
  # dw_imp$cortisol_diff122 <- dw_imp$cortisol_12 - dw_imp$cortisol_2
  # dw_imp
})

save(dw, dw_imp, dimp, dw_imp_all, file = here::here("rdata/data.Rds"))



# create plot that shows change over time for all individuals. I create groups
# to avoid cluttering:

dwp <- dw %>%
  mutate(groupmv = ifelse(cortisol_diff62 <=0 & cortisol_diff126 <= 0, 1, ifelse(
    cortisol_diff62 >=0 & cortisol_diff126 <= 0, 2, ifelse(
      cortisol_diff62 <=0 & cortisol_diff126 >= 0, 3, ifelse(
        cortisol_diff62 >=0 & cortisol_diff126 >= 0, 4, NA
)))))

p1 <- dwp %>%
  filter(groupmv == 1) %>%
  pivot_longer(matches("cortisol_\\d+"), names_to = "time", values_to = "cortisol") %>%
  arrange(cortisol_diff62, cortisol_diff126) %>%
  mutate(
    time = as.numeric(str_extract(time, "\\d+")),
    id = fct_reorder(id, cortisol_diff62, min)
  ) %>%
  ggplot(aes(time, cortisol, group = id)) +
    geom_path(position = position_jitter(width = 0.4, seed = 3), color = "lightgrey", size = 2) + 
    geom_point(position = position_jitter(width = 0.4, seed = 3), size = 6, alpha = 0.75) +
  scale_x_continuous(breaks = c(2, 6, 12)) +
  ylim(0.3, 4) +
  xlab("") + ylab("") +
  theme_bw(base_size = 30)
    
p2 <- dwp %>%
  filter(groupmv == 2) %>%
  pivot_longer(matches("cortisol_\\d+"), names_to = "time", values_to = "cortisol") %>%
  arrange(cortisol_diff62, cortisol_diff126) %>%
  mutate(
    time = as.numeric(str_extract(time, "\\d+")),
    id = fct_reorder(id, cortisol_diff62, min)
  ) %>%
  ggplot(aes(time, cortisol, group = id)) +
    geom_path(position = position_jitter(width = 0.2, seed = 1), color = "lightgrey", size = 2) + 
    geom_point(position = position_jitter(width = 0.2, seed = 1), size = 6, alpha = 0.75) +
  scale_x_continuous(breaks = c(2, 6, 12)) +
  ylim(0.3, 4) +
  xlab("") + ylab("") +
  theme_bw(base_size = 30)

p3 <- dwp %>%
  filter(groupmv == 3) %>%
  pivot_longer(matches("cortisol_\\d+"), names_to = "time", values_to = "cortisol") %>%
  arrange(cortisol_diff62, cortisol_diff126) %>%
  mutate(
    time = as.numeric(str_extract(time, "\\d+")),
    id = fct_reorder(id, cortisol_diff62, min)
  ) %>%
  ggplot(aes(time, cortisol, group = id)) +
  geom_path(position = position_jitter(width = 0.2, seed = 1), color = "lightgrey", size = 2) + 
  geom_point(position = position_jitter(width = 0.2, seed = 1), size = 6, alpha = 0.75) +
  scale_x_continuous(breaks = c(2, 6, 12)) +
  ylim(0.3, 4) +
  ylab("Breast Milk Cortisol") + xlab("Sampling Time Point") +
  theme_bw(base_size = 30)


p4 <- dwp %>%
  filter(groupmv == 4) %>%
  pivot_longer(matches("cortisol_\\d+"), names_to = "time", values_to = "cortisol") %>%
  arrange(cortisol_diff62, cortisol_diff126) %>%
  mutate(
    time = as.numeric(str_extract(time, "\\d+")),
    id = fct_reorder(id, cortisol_diff62, min)
  ) %>%
  ggplot(aes(time, cortisol, group = id)) +
  geom_path(position = position_jitter(width = 0.2, seed = 1), color = "lightgrey", size = 2) + 
  geom_point(position = position_jitter(width = 0.2, seed = 1), size = 6, alpha = 0.75) +
  scale_x_continuous(breaks = c(2, 6, 12)) +
  ylim(0.3, 4) +
  ylab("") + xlab("") +
  theme_bw(base_size = 30)


p4


pathall <- (p1 + ggtitle("A") + p2 + ggtitle("B")) / (p3 + ggtitle("C") + p4 + ggtitle("D"))
save(pathall, file = here::here("rdata/pathall.Rds"))
ggsave(pathall, dpi = 300, width = 20, height = 20, filename = here::here("fig/pathall.png"))
pathall


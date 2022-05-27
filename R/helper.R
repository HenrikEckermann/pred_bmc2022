# renaming function for publication table

library(tidyverse)

rename_vars <- function(tb) {
  mutate(tb,
    feature = ifelse(feature == "iri_fs", "Empathy (IRI FS)", ifelse(
    feature == "mother_sens_w6", "Mother Sensitivity", ifelse(
    str_detect(feature, "collectiontime(_\\d)*"), "Collection Time", ifelse(
    feature == "stai", "State-Trait-Anxiety", ifelse(
    feature == "eg_tot", " Ego Resilience", ifelse(
    str_detect(feature, "interval_awake(_\\d)*"), "Time since Awakening", ifelse(
    feature == "iri_pt", "Empathy (IRI PT)", ifelse(
    feature == "attach_cat_pren2", "Attachment (dismissing)", ifelse(
    feature == "attach_cat_pren3", "Attachment (unresolved)", ifelse(
    feature == "attach_cat_pren4", "Attachment (secure)", ifelse(
    feature == "attach_cat_pren", "Attachment", ifelse(
    str_detect(feature, "depression"), "Depression", ifelse(
    feature == "iri_ec", "Empathy (IRI EC)", ifelse(
    feature == "iri_pd", "Empathy (IRI PD)", ifelse(
    feature == "cc_start", "Start of Childcare", ifelse(
    feature == "type_cc_w12", "Type of Childcare", ifelse(
    feature == "type_cc_w121", "Childcare (not centre-based)", ifelse(
    feature == "type_cc_w122", "Childcare (centre-based)", ifelse(
    feature == "age", "Maternal Age", ifelse(
    str_detect(feature, "season_w\\d+"), "Season", ifelse(
    feature == "season2", "Season Summer", ifelse(
    feature == "season3", "Season Fall", ifelse(
    feature == "season4", "Season Winter", ifelse(
    str_detect(feature, "postnatalweek_\\d"), "Postnatal Week", ifelse(
    feature == "ctq_tot_pren", "Child Trauma Questionnaire", ifelse(
    feature == "edu", "Maternal Education", ifelse(
    feature == "ace_tot_pren", "Adverse Childhood Experience", ifelse(
    feature == "septi_tot", "Self-Efficacy for Parenting", ifelse(
    str_detect(feature, "sharedact"), "Shared Caregiving", ifelse(
    feature == "father_sens_w6", "Father Sensitivity", ifelse(
    feature == "cortisol_pre", "Previous Cortisol", ifelse( 
    feature == "apl", "Daily Hassles", ifelse(
    feature == "auc_cry", "Crying (AUC)", str_to_title(feature)))       
  ))))))))))))))))))))))))))))))))
}


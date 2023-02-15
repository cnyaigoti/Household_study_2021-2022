#Analyse symptoms data from the household study
#Last updated 05-Jan-2023
rm(list=ls())

library(tidyverse); library(readxl); library(lubridate); library(scales);library(janitor); 
library(Epi);library(finalfit);library(RColorBrewer); library(artyfarty);  library(sjlabelled)


#Import package that will give visit IDs
library(splitstackshape); ls(package:splitstackshape)

#1. Load epi data from a stata file

ls(package:sjlabelled)

my_dta<- read_data("~/Dropbox/COVID-19/HHSTUDY/epidata/31Oct2022/household_study_updated_311022.dta")%>%
  as_label()

#Notes
#as_label() converts (replaces) values of a variable (also of factors or character vectors) with their associated value labels.

#Using sapply in a simple function to return a variable list as in Stata's Variable Window:
makeVlist <- function(dta) { 
  labels <- sapply(dta, function(x) attr(x, "label"))
  tibble(name = names(labels),
         label = labels)
}
variables <-as_tibble(makeVlist(my_dta))%>%
  apply(2,as.character)

getwd()

#write.csv(variables, file="variables.csv", na="", row.names = F)

#4. Renaming variables, and creating new ones
hh_dta <- my_dta%>%
  select(SampleName, everything())%>%
  rename(hh_no=Household_number, 
         collectiondate=date_collection, 
         subject_id=study_no,
         test_result=Sample_positive)%>%
  mutate(sample_id=str_replace_all(SampleName, "HH", ""),
         sample_id=paste("P",sample_id, sep=""), 
         hh_study="yes", 
         hh_id=paste("HH", hh_no, sep=""), 
         subject_id=str_replace_all(subject_id, "HH",""))%>%
  mutate(test_result=ifelse(test_result=="", "Negative",test_result))%>%
  mutate(test_result=factor(test_result, levels=c("Negative", "Positive")))%>%
  arrange(hh_no, subject_id, collectiondate)%>%
  mutate(hh_id=factor(hh_id, unique(hh_id)), subject_id=factor(subject_id, unique(subject_id)))%>%
  select(hh_no, subject_id, sample_id, collectiondate, test_result, everything())%>%
  getanID("subject_id")%>%rename(visit_no=".id")%>%
  mutate(sample_status=ifelse(test_result=="Positive", 1,0))%>%
  group_by(subject_id)%>%
  mutate(subject_status=factor(ifelse(sum(sample_status)<1, "uninfected", "infected")))%>%
  ungroup()%>%
  group_by(hh_id)%>%
  mutate(hh_status=factor(ifelse(sum(sample_status)<1, "uninfected", "infected")))%>%
  ungroup()%>%
  mutate(agegrp = cut(age_years, 
                      breaks = c(0.0,5,10,20,40,65,140), 
                      right = F,
                      labels=c("0-4", "5-9", "10-19", "20-39", "40-64", "65+")))%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", ""
                                                 )))%>%
  mutate(other_symptom_diary=case_when(other_symtom_diary=="" ~ "None", 
                                       other_symtom_diary=="N/a" ~ "None",
                                       other_symtom_diary=="N/A" ~ "None",
                                       other_symtom_diary=="none" ~ "None",
                                       other_symtom_diary=="No" ~ "None",
                                       other_symtom_diary=="Non" ~ "None",
                                       TRUE ~ other_symtom_diary))


symptom_dta <- hh_dta%>%
  select(c(hh_id, subject_id, subject_status, collectiondate, test_result, contains("Diar", ignore.case = T), -"Diarrhoea"))%>%
  mutate(Diaries_symptom=as.character (Diaries_symptom))%>%
  mutate(symptom_status=case_when(Diaries_symptom=="Yes" ~ 1,
                                  other_symptom_diary!="None" ~ 1,
                                  TRUE ~ 0))%>%
  group_by(subject_id)%>%
  mutate(symptom_status=factor(case_when(sum(symptom_status)<1 ~"Asymptomatic", 
                                        sum(symptom_status)>0 ~"Symptomatic")))%>%
  select(hh_id, subject_id, symptom_status, everything())%>%
  arrange(hh_id, subject_id)%>%
  ungroup()

# Check if symptoms were frequently observed in those 
symptom_dta%>%
  distinct(subject_id, .keep_all=T)%>%
  tabyl(symptom_status, subject_status)%>%
  adorn_totals(where =c("col"))#%>%
  #adorn_percentages(denominator=c("col"))%>%
  #adorn_pct_formatting(digits=1)#%>%
  #chisq.test()


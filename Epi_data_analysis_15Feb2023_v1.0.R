# Analysis of the Kilifi household study data
# Date started; 30 Jun 2021
# Last updated: 01-Dec-2022

#1. Load packages and set working directory
rm(list=ls())

library(tidyverse); library(readxl); library(lubridate); library(scales); library(janitor); library(writexl)
library(phylotools);library(sjlabelled);library(artyfarty); library(data.table); library(patchwork);library(Epi)

Dark2 <-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666", "#FAD0C9FF", "#CE4A7EFF", "#990011FF")

min <- as.Date("2020-03-12")
max <- as.Date("2022-07-31")

setwd("~/Dropbox/COVID-19/HHstudy/")

#2. Load the genomics data: (a) the sequences, (b) Pango lineages information
sequences_dta <- read.csv("~/Dropbox/COVID-19/variants/all_genomes/previous/16Oct2022/allcombined_16Oct2022_info_details.csv")%>%
  mutate(seq_id=str_replace_all(seq_id,c(" "="_")))%>%
  mutate(name_length=str_length(seq_id))%>%
  distinct(seq_id, .keep_all = T)%>%
  tibble()

lineages_dta <-read.csv("~/Dropbox/COVID-19/variants/all_genomes/previous/16Oct2022/lineage_report.csv")%>%
  mutate(taxon=str_replace_all(taxon,c(" "="_")))%>%
  rename(seq_id =taxon)%>%
  mutate(namelength=str_length(seq_id))%>%
  distinct(seq_id, .keep_all = T)%>%
  tibble()

# Merge sequence data to the Pango lineages data
lineage_seq_dta <- full_join(lineages_dta, sequences_dta, by=c("seq_id"), na="", rownames=F)%>%
  separate(seq_id, c("Run","id"), sep="\\|", fill="right",extra = "warn", remove=F)%>%
  separate(id, c("id", "assembler", "reference"), sep="([\\/])", fill="right",extra = "warn")%>%
  mutate(id=ifelse(str_detect(id, "Run|run"), str_split_fixed(id, "_",2)[,1], id))%>%
  select(id,everything())%>%
  rename(sample_id=id)

#3. Load epi data and check the variables and lables

missed_pos=c("92187","92188","92199","95439","95446","95447","95448","96308","95450","96297","96304","96302")

hh_data<- read_data("~/Dropbox/COVID-19/HHSTUDY/epidata/31Oct2022/household_study_updated_311022.dta")%>%
  select(SampleName, everything())%>%
  select(-c(lineage,clade))%>%
  as_label()%>%
  mutate(Sample_positive=ifelse(SampleName%in%missed_pos, "Positive", Sample_positive))

# check location co-ordinates, to plot spatial distribution of the households
hh_location <- hh_data%>%
  select(Household_number, longitude_dd, latitude_dd)%>%
  distinct(Household_number,.keep_all=T)


#Check total number of samples collected, samples positive
tabyl(hh_data, Sample_positive)%>%
  adorn_totals()%>%
  mutate(percent_formated=percent*100)

# determine the range of dates included in this analysis
hh_data%>%
  filter(!is.na(date_collection))%>%
  pull(date_collection)%>%
  range()

# Check the data completeness
#Using sapply in a simple function to return a variable list as in Stata's Variable Window:
makeVlist <- function(dta) { 
  labels <- sapply(dta, function(x) attr(x, "label"))
  tibble(name = names(labels),
         label = labels)
}

#Check variables in the data and their labels
var_list <-as_tibble(makeVlist(hh_data))%>%
  apply(2,as.character)

#write.csv(missing_glimpse(hh_data)%>%filter(missing_n!=0), file="missing_variables.csv", row.names = F, na="")
tabyl(hh_data, Headache)%>%
  adorn_totals()

hh_epi <- hh_data%>%
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
  arrange(hh_no, subject_id)%>%
  mutate(hh_id=factor(hh_id, unique(hh_id)), subject_id=factor(subject_id, unique(subject_id)))%>%
  select(hh_no, subject_id, sample_id, collectiondate, test_result, everything())%>%
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
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))

########################....................Table 1........................#################################
##Table 1 (1.) ##count number of unique  household ids

#RT-PCR results
hh_epi%>%
  tabyl(test_result)%>%
  adorn_totals()


#number of HH
length(unique(hh_epi$hh_id))

##count number of unique subject ids
length(unique(hh_epi$subject_id))

##count number of unique household that were infected
hh_epi%>%
  distinct(hh_id, .keep_all=T)%>%
  tabyl(hh_status)%>%
  adorn_totals()


##count number of unique subjects that were not infected
hh_epi%>%
  distinct(subject_id, .keep_all=T)%>%
  tabyl(subject_status)%>%
  adorn_totals(where=c("row"))

##sex unique subjects that were infected
hh_epi%>%distinct(subject_id, .keep_all=T)%>%
  tabyl(sex,subject_status, show_missing_levels = T)%>%
  adorn_totals(where=c("row", "col"))%>%
  adorn_percentages("col")%>%
  adorn_pct_formatting(digits = 1)%>%
  adorn_ns()

hh_epi%>%distinct(subject_id, .keep_all=T)%>%
  tabyl(sex,subject_status, show_missing_levels = T)%>%
  chisq.test()

## Median distribution
summary(hh_epi$age_years)
tapply(hh_epi$age_years, hh_epi$subject_status, summary)

## Age category distributions
hh_epi%>%distinct(subject_id, .keep_all=T)%>%
  tabyl(agegrp,subject_status, show_missing_levels = T)%>%
  adorn_totals(where=c("row", "col"))%>%
  adorn_percentages("col")%>%
  adorn_pct_formatting(digits = 1)%>%
  adorn_ns()

hh_epi%>%distinct(subject_id, .keep_all=T)%>%
  tabyl(agegrp,subject_status, show_missing_levels = T)%>%
  filter(!is.na(agegrp))%>%
  chisq.test()

###########################.
# Additional epi analysis
hh_epi_2 <- hh_epi%>%
  group_by(hh_id)%>%
  mutate(hh_earliest=min(collectiondate))%>%
  ungroup()%>%
  arrange(hh_earliest)%>%
  select(hh_id, subject_id, collectiondate, hh_earliest, everything())

# Household study period
hh_data%>%
  filter(!is.na(date_collection))%>%
  pull(date_collection)%>%
  range()

#4. combine household dta and sequence data
updated_dta <- full_join(hh_epi, lineage_seq_dta, by=c("sample_id"), all=TRUE )%>%
  filter(hh_study=="yes")%>%
  distinct(sample_id, collectiondate, .keep_all = T)

#Confirm lineages in the data
sequenced_samples <-updated_dta%>%
  filter(actualseq!="")%>%dim()

#Confirm variants in the data
updated_dta%>%
  filter(!is.na(scorpio_call))%>%
  tabyl(scorpio_call, sample_status)%>%
  adorn_totals()%>%
  adorn_percentages(denominator = "col")%>%
  mutate(total =cumsum(`1`)*100)

#view lineages in updated_dta
updated_dta%>%filter(!is.na(lineage))%>%
  tabyl(lineage)%>%
  adorn_totals()

# Generate variables with mean Ct value and variant status of the sequenced samples
master_dta <- updated_dta %>%
  mutate(taxa_id=paste(lineage, hh_id, subject_id, collectiondate, sep="|"))%>%
  mutate(lineage=ifelse(test_result=="Positive" & is.na(lineage), "Not sequenced", lineage))%>%
  mutate(lineage=factor(lineage), participant_id=as.numeric(subject_id))%>%
  pivot_longer(cols=c(30:42), names_to="kit", values_to="Ct")%>%
  mutate(Ct=as.numeric(recode(Ct, "Undetermined"="")))%>%
  pivot_wider(names_from=kit, values_from=Ct)%>%
  mutate(mean_Ct=rowMeans(.[,73:85], na.rm=TRUE))%>%
  #select(test_result, 73:83)%>%filter(!is.na(mean_Ct))%>%
  mutate(Ct_cat=ifelse(mean_Ct<33.0, "low", "high"))%>%
  mutate(Ct_cat=ifelse(test_result=="Positive", Ct_cat, "Neg"))%>%
  mutate(scorpio_call=if_else(lineage=="B.1.617.2", "Delta (B.1.617.2-like)", scorpio_call))%>%
  mutate(voc=fct_collapse(scorpio_call, 
                          Alpha="Alpha (B.1.1.7-like)",
                          Beta="Beta (B.1.351-like)",
                          Delta=c("Delta (B.1.617.2-like)"),
                          Omicron=c("Omicron (BA.1-like)", "Omicron (BA.2-like)", "Omicron (BA.4-like)", "Omicron (BA.5-like)", "Omicron (Unassigned)", "Probable Omicron (BA.1-like)", "Probable Omicron (Unassigned)"),
                          Ancestral=""))%>%
  mutate(voc=factor(voc, levels=c("Ancestral","Alpha", "Beta", "Delta", "Omicron")))%>%
  mutate(count=1)%>%
  mutate(weekly=as.Date(cut(collectiondate,breaks = "1 weeks")))%>%
  mutate(sequencing_status=factor(case_when(is.na(lineage) ~ "Negative",
                                            lineage=="Not sequenced" ~ "Failed \n sequencing", 
                                            TRUE ~ "Sequenced")))


master_dta%>%filter(!is.na(scorpio_call))%>%tabyl(scorpio_call)%>%adorn_totals(where = c("row", "col"))


# Number of households and participants who provided sequence data
master_dta%>%filter(!is.na(scorpio_call))%>%distinct(participant_id, .keep_all = T)%>%dim()
master_dta%>%filter(!is.na(scorpio_call))%>%distinct(participant_id, .keep_all = T)%>%dim()


# Distribution of sequenced sample
master_dta%>%
  filter(!is.na(voc))%>%
  tabyl(voc)%>%
  adorn_totals()

sequence_hh <-master_dta%>%
  dplyr::filter(sequencing_status=="Sequenced")

sequence_hh %>%distinct(hh_id)%>% dim()

master_dta%>%
  dplyr::filter(sequencing_status=="Sequenced")%>%
  distinct(subject_id, .keep_all = T)%>%dim()

range(sequence_hh$collectiondate)

######################......................Supplementary Figure 1....................####################
##Plot the temporal distribution of all the positives and negatives identified during the study
infection_plot<-hh_epi%>%
  ggplot(aes(x=collectiondate, y=subject_id))+
  geom_point(aes(color=test_result), shape=1, size=1)+
  scale_color_manual(values = c("black", "red"))+
  labs(y="Participants", x="")+
  scale_x_date(breaks ="7 day", labels = date_format("%d%b"))+#, limits = c(min, max))+
  scale_y_discrete(expand = c(0, 2))+
  facet_wrap(~hh_id, scales = "free",ncol = 22)+
  theme_scientific()+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text.x = element_text(size = 10, angle=90, vjust = 0.5),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face="bold"),
        legend.position="bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.box.background = element_blank())+
  guides(color=guide_legend(ncol=2, title = "RT-PCR result", title.position = "left"), size=T)
pdf("~/Dropbox/COVID-19/HHstudy/Figures/S2_Figure.pdf", width =22, height =14, family = "Helvetica")
print(infection_plot)
dev.off()

#infection_plot

######################...................... Figure 1....................####################

# Plot the COVID-19 trends from the national data, panel A
national_dta <-read.csv("~/Dropbox/COVID-19/variants/metadata/global/owid-covid-data.csv")%>%
  select(location,date,total_cases,new_cases,new_cases_smoothed, stringency_index)%>%
  filter(location=="Kenya", date<as.Date("2022-07-31"),date>as.Date("2020-03-12"))%>%
  mutate(date=as.Date(date, format="%Y-%m-%d"))%>%
  arrange(date)

national_trend <- ggplot(national_dta, aes(x=date, y=new_cases_smoothed))+
  geom_area(fill="black")+
  labs(y="Nationwide +ves ", x="", tag="A")+
  #geom_line(aes(x=date, y=stringency_index*30), color="#D95F02", linetype="longdash", size=0.6)+
  theme_scientific()+
  scale_y_continuous(limits = c(0,3000), minor_breaks=seq(0, 3000, 250), breaks=seq(0, 3000,500))+#,
                     #sec.axis=sec_axis(~./30, name="Stringency Index", breaks = seq(0, 100, 20)))+
  scale_x_date(breaks ="2 month", date_minor_breaks="1 month", labels = date_format("%b-%y"), limits = c(min, max))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text.x = element_text(size = 10, angle =90, vjust=0.5),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10, face="bold"),
        legend.position = c(0.65, 0.90),
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=4, title = ""), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/Fig 1A.pdf", width = 5, height = 2, family = "Helvetica")
print(national_trend)
dev.off()
#national_trend


# plot weekly test results from the household samples, Panel B
tests<-master_dta%>%
  group_by(weekly)%>%
  ggplot(aes(x=weekly, y=count, fill=test_result))+
  geom_col()+
  labs(y="Total tested ", x="Date", tag = "B")+
  theme_scientific()+
  scale_fill_manual(values = pal("flat")[c(4,3)])+
  scale_y_continuous(limits = c(0,220), minor_breaks=seq(0, 220, 20), breaks=seq(0, 220, 40))+
  scale_x_date(breaks ="2 month", date_minor_breaks="1 month", labels = date_format("%b-%y"), limits = c(min, max))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text.x = element_text(size = 10, angle =90, vjust=0.5),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10, face="bold"),
        legend.position = c(0.47, 0.9),
        #legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=2, title = "Test result", title.position = "left"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/Fig 1B.pdf", width = 5, height = 2,family = "Helvetica")
print(tests)
dev.off()
#tests

# Plot the COVID-19 Oxford Stringency Index, panel C

stringency <- ggplot(national_dta, aes(x=date, y=stringency_index))+
  geom_line(fill="black")+
  labs(y="Stringency index ", x="", tag = "C")+
  #geom_line(aes(x=date, y=stringency_index*30), color="#D95F02", linetype="longdash", size=0.6)+
  theme_scientific()+
  scale_y_continuous(limits = c(0,100), minor_breaks=seq(0, 100, 10), breaks=seq(0, 100,20))+#,
  #sec.axis=sec_axis(~./30, name="Stringency Index", breaks = seq(0, 100, 20)))+
  scale_x_date(breaks ="2 month", date_minor_breaks="1 month", labels = date_format("%b-%y"), limits = c(min, max))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text.x = element_text(size = 10, angle =90, vjust=0.5),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10, face="bold"),
        legend.position = c(0.65, 0.90),
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=4, title = ""), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/Fig 1C.pdf", width = 5, height = 2, family = "Helvetica")
print(stringency)
dev.off()
#stringency

# plot sequences and variants observed during the study period, Panel D
sequenced <-master_dta%>%
  group_by(weekly)%>%
  filter(test_result=="Positive")%>%
  filter(!is.na(sequencing_status))%>%
  mutate(voc=fct_recode(voc, "failed-sequencing"="Not sequenced"))%>%
  filter(voc!="failed-sequencing")%>%
  ggplot(aes(x=weekly, y=count,fill=voc))+
  geom_col(width = 5)+
  labs(y="Sequenced", x="Date", tag = "D")+
  theme_scientific()+
  scale_fill_manual(values = c("red", pal("five38")[c(5,4,3, 1)], "red", "purple"))+
  scale_y_continuous(limits = c(0,30), minor_breaks=seq(0, 30, 5), breaks=seq(0, 30, 5))+
  scale_x_date(breaks ="2 month", date_minor_breaks="1 month", labels = date_format("%b-%y"), limits = c(min, max))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text.x = element_text(size = 10, angle =90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10, face="bold"),
        legend.position = c(0.52, 0.8),
        #legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(nrow = 2, title = "Variant", title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/Fig 1D.pdf",width = 5, height = 2.0, family = "Helvetica")
print(sequenced)
dev.off()
#sequenced

####################................................Create the panel Figure 1.......................#################

Fig_1 <-(national_trend/tests)|(stringency/sequenced)
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/Fig 1.pdf",width = 8.5, height = 7.01, family = "Helvetica")
print(Fig_1)
dev.off()

#Fig_1

######################...................... Supplementary Figure 2....................####################
#(a) Relationship of Ct value and sequencing status
Ct_plot <-master_dta%>%
  filter(test_result=="Positive")%>%
  filter(is.finite(mean_Ct))%>%
  ggplot(aes(x=sequencing_status, y=mean_Ct))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color=sequencing_status))+
  scale_color_manual(values = pal("flat")[c(1,4)])+
  labs(y="Cycle threshold ", x="", tag="A")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12, angle =0),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 11, face="bold"),
        #legend.position = c(0.65, 0.90),
        legend.position="none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 12),
        legend.title =element_text(size = 12),
        legend.box.background = element_blank())+
  guides(color=guide_legend(ncol=1, title = "Key", title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/S2A_Fig.pdf", width = 2.5, height = 2.5, family = "Helvetica")
print(Ct_plot)
dev.off()
#Ct_plot

#(b) temporal distribution  of variants detected 
lineage_plot <-master_dta%>%
  filter(!is.na(lineage), lineage!="Not sequenced")%>%
  group_by(lineage, collectiondate, voc)%>%
  summarize(lineage, Total=n())%>%
  ggplot(aes(x=collectiondate, y=fct_reorder(lineage, collectiondate)))+
  geom_point(aes(fill=voc, size=Total), shape=21)+
  #scale_fill_manual(values = Dark2)+
  scale_fill_manual(values = c("red", pal("five38")[c(5,4,3, 1)], "red", "purple"))+
  scale_x_date(breaks ="2 month", date_minor_breaks="1 month", labels = date_format("%b-%y"))+
  labs(y="", x="Date", tag = "B")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle =90, vjust=0.5),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11, face="bold"),
        #legend.position = c(0.65, 0.90),
        legend.position="right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "Variant", title.position = "top"),
         size=guide_legend(ncol=1, title = "Count", title.position = "top"))
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/S2B_Fig.pdf", width = 6, height = 4, family = "Helvetica")
print(lineage_plot)
dev.off()

#lineage_plot


S2Fig <- Ct_plot+lineage_plot
pdf("~/Dropbox/COVID-19/HHSTUDY/Figures/S2_Fig.pdf", width = 8, height = 7.01, family = "Helvetica")
print(S2Fig)
dev.off()
S2Fig

##############################...................Supplementary table 1................##############################
## Creating the lineage table - S1 Table
#lineage history
lineages_detected <-master_dta%>%pull(lineage)%>%unique()
lineages_detected

lineage_number <-master_dta%>%
  filter(lineage!="Not sequenced")%>%
  rename(Lineage=lineage)%>%mutate(Lineage=factor(Lineage))%>%
  mutate(count=1)%>%group_by(Lineage)%>%
  summarise(Number=sum(count))

lineage_history <-read_xlsx("~/Dropbox/COVID-19/HHstudy/python/data_lineages.xlsx")%>%
  rename(lineage_no="...1")%>%
  tibble()%>%
  filter(Lineage%in%lineages_detected)%>%
  select(Lineage, `Most common countries`, `Earliest date`, Description)%>%
  mutate(`Earliest date`=as.Date(`Earliest date`, "%Y-%m-%d"))%>%
  mutate(VOC_status=case_when(Lineage=="B.1.351" ~ "Beta",
                       Lineage=="B.1.1.7" ~ "Alpha",
                       str_detect(Lineage, "AY|B.1.617.2") ~  "Delta",
                       str_detect(Lineage, "BA.1") ~ "BA.1-like",
                       str_detect(Lineage, "BA.2") ~ "BA.2-like",
                       str_detect(Lineage, "BA.4") ~ "BA.4-like",
                       str_detect(Lineage, "BA.5|BF|BE") ~ "BA.5-like",
                       TRUE ~ "non-VOC"))%>%
  select(Lineage, VOC_status, everything())%>%
  mutate(VOC_status=factor(VOC_status, levels=c("non-VOC", "Alpha","Beta", "Delta","BA.1-like", "BA.2-like", "BA.4-like", "BA.5-like")))

Lineage_table <-full_join(lineage_number, lineage_history, by=c("Lineage"), na="", rownames=F)%>%
  arrange(VOC_status, `Earliest date`)

writexl::write_xlsx(Lineage_table,path="~/Dropbox/COVID-19/HHstudy/manu/Tables/Lineages.xlsx")




#5 Some quick summaries  
hh_summary<-tabyl(master_dta, hh_id, test_result)%>%
  filter(Negative!=0)

# Check out mean Ct values by lineage
stat.table(lineage, list(count(), "Median Ct"=median(mean_Ct)), margins = T, data=master_dta)
lineages <- master_dta%>%filter(!is.na(lineage) & lineage !="Not sequenced")%>%distinct(lineage) %>%as.list(lineage)
lineages

# Plot household infection timelines
house="HH167"
house_dta <- master_dta%>%
  filter(hh_id==house)

p_hh<-house_dta%>%
  ggplot(aes(x=collectiondate, y=subject_id), shape=21)+
  geom_point(aes(x=collectiondate, y=subject_id, fill=test_result),shape=21, color="white",size=4)+
  #facet_wrap(facets=vars(hh_id), ncol=4, scales = "free_y")+
  scale_x_date(breaks ="7 days", date_minor_breaks="1 weeks", labels = date_format("%d %b"), 
               limits =c(min(house_dta$collectiondate)-4, max(house_dta$collectiondate+4)))+
  scale_fill_manual(values = c("grey", "red"))+
  #scale_color_manual(values = c("grey", "red"))+
  labs(x="", y="", title = "")+
  #coord_cartesian(expand = F)
  theme_scientific()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
        #legend.position = c(0.75, 0.750),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 18),
        legend.title =element_text(size = 18),
        strip.background = element_rect(fill="white", color = "white"),
        panel.spacing.x = unit(1.2,"lines"),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=2, title = "Key", title.position = "left"), size=T)
pdf("~/Dropbox/COVID-19/HHstudy/phylogenetics/household/HH167/house167.pdf", width = 3, height = 5)
print(p_hh)
dev.off()


# Find out distribution of the lineages and number of households that yielded the genomic data
master_dta %>%
  filter(sequencing_status=="Sequenced")%>%
  distinct(hh_id, .keep_all=T)%>%
  tabyl(voc)%>%
  adorn_totals(where="row")

master_dta%>%
  arrange(mean_Ct)%>%
  filter(lineage=="Not sequenced")%>%
  mutate(collectiondate=format(collectiondate, "%d-%b-%Y"))%>%
  select(c(sample_id, collectiondate, mean_Ct, 61:71))

master_dta%>%filter(subject_id=="0481")

#Export sequences of the various variants of concern
phylo_dta <- master_dta%>%
  filter(lineage!=""& lineage!="Not sequenced")%>%
  select(taxa_id, subject_id,participant_id, hh_id, lineage, voc,collectiondate, actualseq)%>%
  mutate(Ns=str_count(actualseq, "N"))%>%
  arrange(Ns)%>%
  mutate(duplicated=duplicated(taxa_id))%>%
  select(taxa_id, duplicated, Ns, voc, lineage, hh_id, collectiondate, everything())%>%
  filter(duplicated==FALSE)%>%
  tibble()

phylo_dta%>%tabyl(lineage, voc)%>%adorn_totals(where=c("row", "col"))

Alpha <- phylo_dta %>%
  filter(voc=="Alpha")%>%
  arrange(voc,lineage, hh_id, collectiondate)%>%
  select(taxa_id,actualseq)%>%
  rename(seq.name=taxa_id, seq.text=actualseq)%>%
  dat2fasta(outfile = "~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Alpha_household.fasta")

Beta <- phylo_dta %>%
  filter(voc=="Beta")%>%
  arrange(voc,lineage, hh_id, collectiondate)%>%
  select(taxa_id,actualseq)%>%
  rename(seq.name=taxa_id, seq.text=actualseq)%>%
  dat2fasta(outfile = "~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Beta_household.fasta")

Delta <- phylo_dta %>%
  filter(voc=="Delta")%>%
  arrange(voc,lineage, hh_id, collectiondate)%>%
  select(taxa_id,actualseq)%>%
  rename(seq.name=taxa_id, seq.text=actualseq)%>%
  dat2fasta(outfile = "~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Delta_household.fasta")

Omicron <- phylo_dta %>%
  filter(voc=="Omicron")%>%
  arrange(voc,lineage, hh_id, collectiondate)%>%
  select(taxa_id,actualseq)%>%
  rename(seq.name=taxa_id, seq.text=actualseq)%>%
  dat2fasta(outfile = "~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Omicron_household.fasta")

All_variants <- phylo_dta %>%
  arrange(voc,lineage, hh_id, collectiondate)%>%
  select(taxa_id,actualseq)%>%
  rename(seq.name=taxa_id, seq.text=actualseq)%>%
  dat2fasta(outfile = "~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/All_household.fasta")

hh_89<- phylo_dta %>%
  filter(hh_id=="HH89")%>%
  arrange(voc,lineage, hh_id, collectiondate)%>%
  select(taxa_id,actualseq)%>%
  rename(seq.name=taxa_id, seq.text=actualseq)%>%
  dat2fasta(outfile = "~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH89/hh89.fasta")


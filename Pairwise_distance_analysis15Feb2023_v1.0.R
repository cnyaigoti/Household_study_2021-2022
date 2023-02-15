#Code for plotting pairwise distances
#Written by:  Dr Nyaigoti Agoti
#Last updated 05 Jan 2023

# First import all necessary libraries
rm(list=ls())
library(tidyverse); library(janitor); library(lubridate); library(artyfarty); library(rstatix);library(patchwork);library(wesanderson)



plot_color=c("#000000","#C0C0C0","#696969","#F2D2BD","#800000","#00FF00","#00FFFF","#008000","#FFA500","#9ACD32",
             "#0000FF","#FF00FF","#55ACEE", "#FF0000","#FFFF00","#A4C639","#CCCCFF")

setwd("/Users/cnyaigoti/Dropbox/COVID-19/HHSTUDY/phylogenetics/pairsnp/all/")
hh_dta <- read.csv("./All_household_05Jan2023.csv", stringsAsFactors = F, sep = ",", header=T)

row.names(hh_dta)=colnames(hh_dta)

hh_dta%>%dim()
#hh_dta%>%view()

pairwise_dta <- as.matrix(hh_dta)%>%
  replace_upper_triangle(by=NA, diagonal=TRUE)%>%
  pivot_longer((2:288), names_to="seq_id", values_to="nt_diff")%>%
  filter(!is.na(nt_diff))%>%
  mutate(self=if_else(rowname==seq_id, 0, 1))%>%
  filter(self==1)%>%
  select(rowname, seq_id, nt_diff)%>%
  separate(rowname, into=c("lineage", "hh_id", "subject_id", "date"), remove=FALSE, sep="__")%>%
  mutate(voc=case_when(lineage=="B.1.351" ~ "Beta",
                              lineage=="B.1.1.7" ~ "Alpha",
                              str_detect(lineage, "AY|B.1.617.2") ~  "Delta",
                              str_detect(lineage, "BA.1|BA.2|BA.4|BA.5|BF|BE") ~ "Omicron",
                              TRUE ~ "Ancestral"))%>%
  mutate(voc=factor(voc,levels=c("Ancestral", "Alpha", "Beta", "Delta", "Omicron")))%>%
  separate(seq_id, into=c("lineage2", "hh_id2", "subject_id2", "date2"), remove=F, sep="__")%>%
  mutate(voc2=case_when(lineage2=="B.1.351" ~ "Beta",
                       lineage2=="B.1.1.7" ~ "Alpha",
                       str_detect(lineage2, "AY|B.1.617.2") ~  "Delta",
                       str_detect(lineage2, "BA.1|BA.2|BA.4|BA.5|BF|BE") ~ "Omicron",
                       TRUE ~ "Ancestral"))
  

multiple_intro_hh <-pairwise_dta%>%filter(nt_diff>2, hh_id==hh_id2)%>%distinct(hh_id, .keep_all=T)

summary_hh <-pairwise_dta%>%filter(hh_id==hh_id2, voc==voc2)

summary(summary_hh$nt_diff)


glimpse(pairwise_dta)

nt_compared_alpha <-pairwise_dta%>%
  filter(voc=="Alpha" & voc2=="Alpha")%>%
  mutate(self=ifelse(hh_id==hh_id2, 1, 0))%>%
  filter(self==1)%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))%>%
  mutate(hh_id=as.character(hh_id))%>%
  filter(str_detect(seq_id, hh_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))

within_hh_alpha <-nt_compared_alpha%>%
  group_by(hh_id)%>%
  #filter(!is.na(household_no))%>%
  ggplot(aes(x=hh_id, y=nt_diff))+
  geom_boxplot()+
  geom_jitter(aes(fill=hh_id), shape=21, color="black")+
  #coord_flip()+
  scale_fill_manual(values=plot_color)+
  labs(y="Number of nt differences", x="", title = "Alpha")+
  scale_y_continuous(limits = c(0,4), minor_breaks = seq(0,4,1), breaks = seq(0,4,1))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5, hjust = 1), 
        legend.position = "none",
        #legend.position = c(0.78, 0.80),
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=6, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/pairsnp/Alpha_within_HH2.pdf", width = 3, height = 3.5)
print(within_hh_alpha)
dev.off()

within_hh_alpha

within_hh_alpha

# 2. Beta
nt_compared_beta <-pairwise_dta%>%
  filter(voc=="Beta"&voc2=="Beta")%>%
  mutate(self=ifelse(hh_id==hh_id2, 1, 0))%>%
  filter(self==1)%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
   mutate(hh_id=fct_reorder(hh_id, household_no))%>%
  mutate(hh_id=as.character(hh_id))%>%
  filter(str_detect(seq_id, hh_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))

within_hh_beta <-nt_compared_beta%>%
  group_by(hh_id)%>%
  ggplot(aes(x=hh_id, y=nt_diff))+
  geom_boxplot()+
  geom_jitter(aes(fill=hh_id), color="black", shape=21)+
  #coord_flip()+
  theme_bw()+
  scale_fill_manual(values=plot_color[c(4, 7, 12, 14)])+
  labs(y="Number of nt differences", x="", title = "Beta")+
  scale_y_continuous(limits = c(0,12), minor_breaks = seq(0,12,1), breaks = seq(0,12,2))+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5, hjust = 1), 
        #legend.position = c(0.78, 0.80),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=6, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/pairsnp/Beta/Beta_within_HH2.pdf", width = 1.5, height = 3.5)
print(within_hh_beta)
dev.off()


within_hh_beta

#3. Delta

nt_compared_delta <-pairwise_dta%>%
  filter(voc=="Delta"&voc2=="Delta")%>%
  mutate(self=ifelse(hh_id==hh_id2, 1, 0))%>%
  filter(self==1)%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))%>%
  mutate(hh_id=as.character(hh_id))%>%
  filter(str_detect(seq_id, hh_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))

within_hh_delta <-nt_compared_delta%>%
  group_by(hh_id)%>%
  ggplot(aes(x=hh_id, y=nt_diff))+
  geom_boxplot()+
  geom_jitter(aes(fill=hh_id), color="black", shape=21)+
  theme_bw()+
  scale_y_continuous(limits = c(0,32), minor_breaks = seq(0,32,4), breaks = seq(0,32,4))+
  scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  labs(y="Number of nt differences", x="", title = "Delta")+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5, hjust = 1), 
        #legend.position = c(0.78, 0.80),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=4, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/pairsnp/Delta/Delta_within_HH2.pdf", width = 5, height = 3.5)
print(within_hh_delta)
dev.off()


within_hh_delta



#4. Omicron

nt_compared_omicron <-pairwise_dta%>%
  filter(voc=="Omicron"& voc2=="Omicron")%>%
  mutate(self=ifelse(hh_id==hh_id2, 1, 0))%>%
  filter(self==1)%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))%>%
  mutate(hh_id=as.character(hh_id))%>%
  filter(str_detect(seq_id, hh_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))

within_hh_omicron <-nt_compared_omicron%>%
  group_by(hh_id)%>%
  ggplot(aes(x=hh_id, y=nt_diff))+
  geom_boxplot()+
  geom_jitter(aes(fill=hh_id), color="black", shape=21)+
  theme_bw()+
  scale_y_continuous(limits = c(0,26), minor_breaks = seq(0,26,2), breaks = seq(0,26,4))+
  scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  labs(y="Number of nt differences", x="", title = "Omicron")+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=11, angle=90, vjust=0.5, hjust = 1), 
        #legend.position = c(0.78, 0.80),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=4, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/pairsnp/Omicron/Omicron_within_HH2.pdf", width = 2, height = 3.5)
print(within_hh_omicron)
dev.off()


within_hh_omicron

figure4 <- ((within_hh_alpha|within_hh_beta)/within_hh_delta)/within_hh_omicron
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/pairsnp/Figure 4.pdf", width = 8.5, height = 7.01)
print(figure4)
dev.off()

figure4

print ("I need Madondo!")


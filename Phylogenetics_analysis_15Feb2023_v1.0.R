# Code for variant surveillance
# Code written by Dr Charles N.Agoti
# Last updated 04Jan2023
#1. Clear any data in memory, load necessary packages and set working directory

rm(list=ls())
library(tidyverse); library(stringr); library(lubridate); library(dplyr); library(artyfarty); library(janitor); library(ggtree)
library(readr); library(openxlsx); library(data.table); library(readxl); library(treeio);library(phylotools);
library(wesanderson);library(scales); library(lubridate);library(patchwork)

plot_color=c("#000000","#C0C0C0","#696969","#F2D2BD","#800000","#00FF00","#00FFFF","#008000","#FFA500",
             "#0000FF","#FF00FF","#55ACEE", "#FF0000","#FFFF00","#A4C639","#CCCCFF")
Dark2 <-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666", "#FAD0C9FF", "#CE4A7EFF", "cyan","#990011FF")


setwd("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/")

# 1. Alpha
meta_alpha_dta <- phylotools::read.fasta(file="31May2022/Alpha/Coast_Alpha_31May2022.aligned.fasta")%>%
  separate(seq.name, into=c("lineage", "county", "sample_id", "datecollection"), sep="\\|", remove=FALSE)%>%
  mutate(hh_id=ifelse(str_detect(county, "HH"), county, ""))%>%
  mutate(subject_id=ifelse(str_detect(county, "HH"), sample_id, ""))%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(county=ifelse(str_detect(seq.name, "HH"), "", county), 
         sample_id=ifelse(str_detect(seq.name, "HH"), "", sample_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))

  
tabyl(meta_alpha_dta, county)%>%adorn_totals()
tabyl(meta_alpha_dta, hh_id)

alpha_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Alpha/Coast_Alpha_31May2022.aligned.fasta.treefile")
alpha_p <- ggtree(alpha_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 400)+
  theme(axis.text.x = element_text(size=10,angle=0))
alpha_p

alpha_tre <-  alpha_p%<+% meta_alpha_dta+ 
  geom_tippoint(aes(subset=(hh_id!=""),fill=hh_id,),size=2, stroke=0.2, color="black", shape=21)+
  labs(y="No of sequences", x="Genetic distance", title="Alpha")+
  scale_fill_manual(values=plot_color)+
  #scale_color_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1"))+
  theme_bw()+
  scale_y_continuous(limits=c(0,375), minor_breaks = seq(0 , 375, 50), breaks = seq(0,375, 100))+
  scale_x_continuous(labels = comma, limits = c(0, 0.00075), breaks = seq(0, 0.0009, 0.0003))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text = element_text(size=10), 
        #legend.position = c(0.85, 0.7),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=3, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Alpha/Alpha_coast_2022.pdf", width = 2.2, height = 6.02)
print(alpha_tre)
dev.off()

alpha_tre


meta_beta_dta <- read.fasta(file="~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Beta/Coast_Beta_31May2022.aligned.fasta")%>%
  separate(seq.name, into=c("lineage", "county", "sample_id", "datecollection"), sep="\\|", remove=FALSE)%>%
  mutate(hh_id=ifelse(str_detect(county, "HH"), county, ""))%>%
  mutate(subject_id=ifelse(str_detect(county, "HH"), sample_id, ""))%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(county=ifelse(str_detect(seq.name, "HH"), "", county), 
         sample_id=ifelse(str_detect(seq.name, "HH"), "", sample_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))



tabyl(meta_beta_dta, hh_id)%>%adorn_totals()

beta_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Beta/Coast_Beta_31May2022.aligned.fasta.treefile")
beta_p <- ggtree(beta_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 180)+
  theme(axis.text.x = element_text(size=10,angle=0))
beta_p

beta_tre <-  beta_p%<+% meta_beta_dta+ 
  geom_tippoint(aes(subset=(hh_id!=""),fill=hh_id,),size=2, stroke=0.2, color="black", shape=21)+
  labs(y="No of sequences", x="Genetic distance", title ="Beta")+
  scale_fill_manual(values=plot_color[c(2,4,5, 7, 11, 13, 15)])+
  theme_bw()+
  scale_y_continuous(limits=c(0,180), minor_breaks = seq(0 , 180, 25), breaks = seq(0 , 180, 50))+
  scale_x_continuous(labels = comma, limits = c(0, 0.00075), breaks = seq(0, 0.0009, 0.0003))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text = element_text(size=10), 
        #legend.position = c(0.85, 0.7),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=3, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Beta/Beta_coast_2022.pdf", width = 2.2, height = 6.02)
print(beta_tre)
dev.off()

beta_tre

#3. Delta tree
meta_delta_dta <- phylotools::read.fasta(file="~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Delta/Coast_Delta_31May2022.aligned.fasta")%>%
  separate(seq.name, into=c("lineage", "county", "sample_id", "datecollection"), sep="\\|", remove=FALSE)%>%
  mutate(hh_id=ifelse(str_detect(county, "HH"), county, ""))%>%
  mutate(subject_id=ifelse(str_detect(county, "HH"), sample_id, ""))%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(county=ifelse(str_detect(seq.name, "HH"), "", county), 
         sample_id=ifelse(str_detect(seq.name, "HH"), "", sample_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))



tabyl(meta_delta_dta, hh_id)%>%adorn_totals()

delta_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Delta/Coast_Delta_31May2022.aligned.fasta.treefile")
delta_p <- ggtree(delta_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 550)+
  theme(axis.text.x = element_text(size=10,angle=0))
delta_p

delta_tre <-  delta_p%<+% meta_delta_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),fill=hh_id),size=2, stroke=0.2, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance", title = "Delta")+
  scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  theme_bw()+
  scale_y_continuous(limits=c(0,550), minor_breaks = seq(0 , 550, 50), breaks = seq(0 , 550, 100))+
  scale_x_continuous(labels = comma, limits = c(0, 0.0018), breaks = seq(0, 0.0018, 0.0008))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text = element_text(size=10), 
        #legend.position = c(0.85, 0.6),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=3, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/31May2022/Delta/Delta_coast_2021.pdf", width = 2.2, height = 6.02)
print(delta_tre)
dev.off()

delta_tre

#4. Omicron tree
meta_omicron_dta <- phylotools::read.fasta(file="~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/omicron/ALIGNMENT/Combined_omicron_COAST+HH.aligned.fasta")%>%
  separate(seq.name, into=c("lineage", "county", "sample_id", "datecollection"), sep="\\|", remove=FALSE)%>%
  mutate(hh_id=ifelse(str_detect(county, "HH"), county, ""))%>%
  mutate(subject_id=ifelse(str_detect(county, "HH"), sample_id, ""))%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(county=ifelse(str_detect(seq.name, "HH"), "", county), 
         sample_id=ifelse(str_detect(seq.name, "HH"), "", sample_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))

tabyl(meta_omicron_dta, hh_id)%>%adorn_totals()
tabyl(meta_omicron_dta, lineage)%>%adorn_totals()

omicron_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/omicron/iqtree/Combined_omicron_COAST+HH.aligned.fasta.treefile")
omicron_p <- ggtree(omicron_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 600)+
  theme(axis.text.x = element_text(size=10,angle=0))
omicron_p

omicron_tre <- omicron_p%<+% meta_omicron_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),fill=hh_id),size=2, stroke=0.2, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance", title="Omicron")+
  #scale_fill_manual(values=plot_color[-1])+
  scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  theme_bw()+
  scale_y_continuous(limits=c(0,580), minor_breaks = seq(0 , 580, 50), breaks = seq(0 , 580, 100))+
  scale_x_continuous(labels = comma, limits = c(0, 0.0025), breaks = seq(0, 0.0025, 0.001))+
  theme(axis.title.x = element_text(size = 10, face="bold"),
        axis.title.y = element_text(size = 10, face="bold"),
        axis.text = element_text(size=10), 
        #legend.position = c(0.85, 0.6),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=3, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/omicron/Omicron_coast_2022.pdf", width = 2.2, height = 6.02)
print(omicron_tre)
dev.off()

omicron_tre

figure3 <-(alpha_tre|beta_tre | delta_tre|omicron_tre)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/Figure 3.pdf", height = 7.02, width = 8.8)
print(figure3)
dev.off()


#5. all variants tree

meta_all_dta <- phylotools::read.fasta(file="~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/seq_data/Combined_COAST+HH.fasta")%>%
  separate(seq.name, into=c("lineage", "county", "sample_id", "datecollection"), sep="\\|", remove=FALSE)%>%
  mutate(hh_id=ifelse(str_detect(county, "HH"), county, ""))%>%
  mutate(subject_id=ifelse(str_detect(county, "HH"), sample_id, ""))%>%
  mutate(household_no=as.integer(str_replace_all(hh_id, "HH", "")))%>%
  mutate(county=ifelse(str_detect(seq.name, "HH"), "", county), 
         sample_id=ifelse(str_detect(seq.name, "HH"), "", sample_id))%>%
  mutate(hh_id=fct_reorder(hh_id, household_no))%>%
  mutate(VOC_status=case_when(lineage=="B.1.351" ~ "Beta",
                              lineage=="B.1.1.7" ~ "Alpha",
                              str_detect(lineage, "AY|B.1.617.2") ~  "Delta",
                              str_detect(lineage, "BA.1") ~ "BA.1-like",
                              str_detect(lineage, "BA.2") ~ "BA.2-like",
                              str_detect(lineage, "BA.4") ~ "BA.4-like",
                              str_detect(lineage, "BA.5|BF|BE") ~ "BA.5-like",
                              TRUE ~ "non-VOC"))%>%
  select(lineage, VOC_status, everything())%>%
  mutate(VOC_status=factor(VOC_status, levels=c("non-VOC", "Alpha","Beta", "Delta","BA.1-like", "BA.2-like", "BA.4-like", "BA.5-like")))%>%
  arrange(VOC_status, lineage)%>%
  mutate(lineage=factor(lineage, levels=unique(lineage)))%>%
  mutate(date=decimal_date(as.Date(datecollection)))%>%
  mutate(datecollection=as.Date(datecollection))%>%
  rename(strain=seq.name)%>%
  select(strain, date, everything(), -seq.text)%>%
  arrange(date)%>%
  mutate(participant=paste(hh_id, subject_id, sep="_"))%>%
  mutate(participant=ifelse(participant=="_", "non_HH", participant))%>%
  mutate(Source=ifelse(str_detect(hh_id, "HH"), as.character(hh_id), "non_HH"))

  

glimpse(meta_all_dta)

tabyl(meta_all_dta,participant)%>%mutate(n=1, count=cumsum(n))
tabyl(meta_all_dta,hh_id)%>%mutate(n=1, count=cumsum(n))


tabyl(meta_all_dta,VOC_status)
tabyl(meta_all_dta,lineage)
tabyl(meta_all_dta,Source)




range(meta_all_dta$datecollection)
range(meta_all_dta$date)


write.csv(meta_all_dta, file="~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport/all_variants_03Jan2023.csv", na="", row.names = F)

write.csv(meta_all_dta, file="~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport2/all_variants_04Jan2023.csv", na="", row.names = F)

tabyl(meta_all_dta, hh_id)%>%adorn_totals()

all_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/iqtree/Combined_COAST+HH.aligned.fasta.treefile")
all_p <- ggtree(all_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 2900)+
  theme(axis.text.x = element_text(size=10,angle=0))
all_p

all_tre <-  all_p%<+% meta_all_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),shape=VOC_status, fill=lineage),size=2, stroke=0.2, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance", tag = "A", title="Mutation tree")+
  scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  #scale_fill_manual(values=c("grey", "skyblue", "orange", "maroon", "red"))+
  #scale_fill_manual(values=Dark2)+
  #coord_flip()+
  theme_scientific()+
  scale_y_continuous(limits=c(0,2900), minor_breaks = seq(0 , 2900, 250), breaks = seq(0 , 2900, 500))+
  scale_x_continuous(labels = comma, limits = c(0, 0.0035), breaks = seq(0, 0.0035, 0.001))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size=12), 
        #legend.position = c(0.85, 0.8),
        legend.position = "none",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 12),
        legend.title =element_text(size = 12),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "Variant",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/coast_household_2023.pdf", width = 8.5, height = 5.02)
print(all_tre)
dev.off()

all_tre

#### time tree equivalent of all variants tree
recent_all_tree <-as.Date(meta_all_dta%>%
                            select(date)%>%max()%>%
                            date_decimal()%>%
                            format("%Y-%m-%d"))
recent_all_tree

time_all_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/treetime/2023-01-03_treetime/all_variants_timetree.nwk")
time_all_p <- ggtree(time_all_tree, mrsd=recent_all_tree,as.Date=TRUE,color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 2900)+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "3 month", minor_breaks = "1 month")+
  theme(axis.text.x = element_text(size=10,angle=0))
time_all_p

time_all_tre <-  time_all_p%<+% meta_all_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),shape=voc, fill=lineage),size=2, stroke=0.2, color="black",shape=21)+
  labs(y="No of sequences", x="", tag = "B", title="Time tree")+
  scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  #scale_fill_manual(values=c("grey", "skyblue", "orange", "maroon", "red"))+
  #scale_fill_manual(values=Dark2)+
  #coord_flip()+
  theme_scientific()+
  scale_y_continuous(limits=c(0,2900), minor_breaks = seq(0 , 2900, 250), breaks = seq(0 , 2900, 500))+
  #scale_x_continuous(labels = comma, limits = c(0, 0.00275), breaks = seq(0, 0.00275, 0.0005))+
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size=12), 
        legend.position = c(0.85, 0.5),
        #legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 12),
        legend.title =element_text(size = 12),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "Variant",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/time_coast_household_2023.pdf", width = 8.5, height = 5.02)
print(time_all_tre)
dev.off()

time_all_tre




Supple_Figure_2 <-(all_tre)|(time_all_tre)+plot_layout(guides="collect")
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/Supplementary Figure 2.pdf", height = 7.4, width = 12)
print(Supple_Figure_2)
dev.off()
Supple_Figure_2
######

hh9_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH9/iqtree/HH9.aligned.fasta.treefile")
hh9_p <- ggtree(hh9_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 15)+
  theme(axis.text.x = element_text(size=10,angle=0))
hh9_p

hh9_tre <-  hh9_p%<+% meta_all_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),shape=voc, fill=subject_id),size=3, stroke=0.4, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance")+#, tag = "C", title="Delta")+
  #scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  scale_fill_manual(values=c("grey", "skyblue", "orange", "maroon", "red"))+
  scale_fill_manual(values=c("blue", Dark2[c(4,6,12)]))+
  #coord_flip()+
  theme_bw()+
  scale_y_continuous(limits=c(0,15), minor_breaks = seq(0 , 15, 2), breaks = seq(0 , 15, 4))+
  scale_x_continuous(labels = comma, limits = c(0, 0.0004), breaks = seq(0, 0.0004, 0.00015))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=14), 
        legend.position = c(0.75, 0.4),
        #legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 14),
        legend.title =element_text(size = 14),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "Participant",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH9/iqtree/hh9_tree.pdf", width = 3, height = 4)
print(hh9_tre)
dev.off()

hh9_tre


###
hh53_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH53/iqtree/hh53.aligned.fasta.treefile")
hh53_p <- ggtree(hh53_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 15)+
  theme(axis.text.x = element_text(size=10,angle=0))
hh53_p

hh53_tre <-  hh53_p%<+% meta_all_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),shape=voc, fill=subject_id),size=3, stroke=0.4, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance")+#, tag = "C", title="Delta")+
  #scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  scale_fill_manual(values=c("grey", "skyblue", "orange", "maroon", "red"))+
  scale_fill_manual(values=c("blue", "cyan","red", Dark2[c(6,8,9,12)]))+
  #coord_flip()+
  theme_bw()+
  scale_y_continuous(limits=c(0,12), minor_breaks = seq(0 , 12, 2), breaks = seq(0 , 12, 4))+
  scale_x_continuous(labels = comma, limits = c(0, 0.00015), breaks = seq(0, 0.00015, 0.0001))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=14), 
        legend.position = c(0.75, 0.35),
        #legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 14),
        legend.title =element_text(size = 14),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "Participant",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH53/iqtree/hh53_tree.pdf", width = 3, height = 4)
print(hh53_tre)
dev.off()

hh53_tre

###
hh64_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH64/iqtree/hh64.aligned.fasta.treefile")
hh64_p <- ggtree(hh64_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 14)+
  theme(axis.text.x = element_text(size=10,angle=0))
hh64_p

hh64_tre <-  hh64_p%<+% meta_all_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),shape=voc, fill=subject_id),size=3, stroke=0.4, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance")+#, tag = "C", title="Delta")+
  #scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  scale_fill_manual(values=c("grey", "skyblue", "orange", "maroon", "red"))+
  scale_fill_manual(values=c("blue", "cyan","red", Dark2[c(6,8,9,12)]))+
  #coord_flip()+
  theme_bw()+
  scale_y_continuous(limits=c(0,14), minor_breaks = seq(0 , 14, 2), breaks = seq(0 , 14, 4))+
  scale_x_continuous(labels = comma, limits = c(0, 0.00015), breaks = seq(0, 0.00015, 0.0001))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=14), 
        legend.position = c(0.75, 0.25),
        #legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 14),
        legend.title =element_text(size = 14),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "Participant",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH64/iqtree/hh64_tree.pdf", width = 3, height = 4)
print(hh64_tre)
dev.off()

hh64_tre

###
hh114_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH114/iqtree/hh114.aligned.fasta.treefile")
hh114_p <- ggtree(hh114_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 11)+
  theme(axis.text.x = element_text(size=10,angle=0))
hh114_p

hh114_tre <-  hh114_p%<+% meta_all_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),shape=voc, fill=subject_id),size=3, stroke=0.4, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance")+#, tag = "C", title="Delta")+
  #scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  scale_fill_manual(values=c("grey", "skyblue", "orange", "maroon", "red"))+
  scale_fill_manual(values=c("blue", "cyan","red", Dark2[c(6,8,9,12)]))+
  #coord_flip()+
  theme_bw()+
  scale_y_continuous(limits=c(0,11), minor_breaks = seq(0 , 11, 1), breaks = seq(0 , 11, 2))+
  scale_x_continuous(labels = comma, limits = c(0, 0.0009), breaks = seq(0, 0.0009, 0.0004))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size=14), 
        legend.position = c(0.75, 0.25),
        #legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 14),
        legend.title =element_text(size = 14),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "Participant",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH114/iqtree/hh114_tree.pdf", width = 3, height = 4)
print(hh114_tre)
dev.off()

hh114_tre


hh167_tree <-read.newick("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH167/iqtree/household_167.aligned.fasta.treefile")
hh167_p <- ggtree(hh167_tree, color='grey40',size=0.2)+
  theme_tree2()+
  expand_limits(y = 11)+
  theme(axis.text.x = element_text(size=10,angle=0))
hh167_p

hh167_tre <-  hh167_p%<+% meta_all_dta+ 
  #geom_tippoint(aes(subset=(hh_id=="")),fill='whi',size=1.5,stroke=0.2,color='black',shape=21)+
  geom_tippoint(aes(subset=(hh_id!=""),shape=voc, fill=subject_id),size=4, stroke=0.4, color="black",shape=21)+
  labs(y="No of sequences", x="Genetic distance")+#, tag = "C", title="Delta")+
  #scale_fill_manual(values=c(wes_palette(5, name = "Darjeeling1", type = "discrete")[c(-1,-2)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1",plot_color))+
  scale_fill_manual(values=c("grey", "skyblue", "orange", "maroon", "red"))+
  scale_fill_manual(values=c("blue", "cyan","red", Dark2[c(6,8,9,12)]))+
  #coord_flip()+
  theme_scientific()+
  scale_y_continuous(limits=c(0,6), minor_breaks = seq(0 , 6, 1), breaks = seq(0 , 6, 2))+
  scale_x_continuous(labels = comma, limits = c(0, 0.0022), breaks = seq(0, 0.0022, 0.001))+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size=18), 
        legend.position = c(0.75, 0.75),
        #legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 18),
        legend.title =element_text(size = 18),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=1, title = "ID",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/household/HH167/iqtree/hh167_tree.pdf", width = 3, height = 5)
print(hh167_tre)
dev.off()

hh167_tre




glimpse(meta_all_dta)
tabyl(meta_all_dta, hh_id)%>%
  adorn_totals()


range(meta_all_dta$datecollection)


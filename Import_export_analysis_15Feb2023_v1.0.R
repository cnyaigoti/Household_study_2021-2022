# Code written by Dr Charles Agoti
# Analysis household import export events
#Last update 05 Jan 2023

# Clear the work space and upload required packages

rm(list=ls())
library(tidyverse);library(janitor); library(lubridate); library(ggthemes); library(artyfarty); library(ggalluvial);library(ggrepel); library(patchwork)
plot_color <- c("#FF0000", "#FFFF00", "#00EAFF", "#AA00FF", "#FF7F00", "#BFFF00", "#0095FF","#FF00AA","#FFD400", "#6AFF00", "#0040FF","#EDB9B9", "#B9D7ED", 
                "#E7E9B9", "#DCB9ED", "#B9EDE0", "#8F2323", "#23628F", "#8F6A23", "#6B238F", "#4F8F23", "#000000", "#737373", "#CCCCCC")

# specify working directory and upload data
setwd("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport/2023-01-03_mugration/")
import_dta <- read.csv("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport/2023-01-03_mugration/annottated_tree_events.csv")%>%
  tibble()

glimpse(import_dta)

tabyl(import_dta, Destination)%>%
  filter(Destination!="nonHH")%>%
  mutate(cumtotal=cumsum(n))%>%
  mutate(hh_count=row_number())%>%
  mutate(number_introductions=factor(n))%>%
  tabyl(number_introductions)%>%
  arrange(desc(number_introductions))%>%
  mutate(sumsum=cumsum(percent))%>%
  adorn_totals()

total_imports <- import_dta%>%filter(Origin!="UNKNOWN", str_detect(Origin, "nonHH"))%>%tabyl(Destination)%>%arrange(desc(n))%>%adorn_totals()

##specify missing values are non-household 
hh_dta <-import_dta%>%
  filter(Origin!="UNKNOWN")%>%
  mutate(Origin=if_else(str_detect(Origin, "non"),"non-HH", Origin))%>%
  mutate(Destination=if_else(str_detect(Destination, "non"),"non-HH", Destination))%>%
  mutate(EventType=case_when(str_detect(Destination,"non") ~ "Export",
                             str_detect(Origin,"non") ~ "Import",
                             TRUE ~ "Interhousehold"))


tabyl(hh_dta, EventType)%>% adorn_totals()

hh_dta%>%filter(EventType!="Export")%>%
  tabyl(Destination)%>%
  mutate(count=1)%>%
  group_by(Destination)%>%
  summarise(count2=sum(n))%>%
  ggplot(aes(x=Destination, y=count2))+
               geom_col()

hh_intros <-hh_dta%>%filter(Destination!="non-HH")%>%tabyl(Destination)%>%
  mutate(household_no=as.integer(str_replace_all(Destination, "HH", "")))%>%
  mutate(hh_id =fct_reorder(Destination, household_no))
  
  
hh_intros_plot <-hh_intros%>%
  ggplot()+
  geom_point(aes(x=hh_id, y=n), size=3)+
  geom_segment(aes(x=hh_id, xend=hh_id, y=0, yend=n))+
  labs(x="", y="No. of introductions")+
  #coord_flip()+
  theme_scientific()+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=8,angle=90, vjust=0.5, hjust = 1), 
        legend.position = "none",
        #legend.position = c(0.78, 0.80),
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())#+
  #guides(fill=guide_legend(ncol=6, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport/patterns.pdf", width = 8.5, height = 3.5)
print(hh_intros_plot)
dev.off()

hh_intros_plot

hh_intros_summary <-hh_intros%>%
  mutate(no_intro =factor(n))%>%
  tabyl(no_intro)%>%
  ggplot()+
  geom_col(aes(x=no_intro, y=n))+
  labs(x="No. of introductions", y="Frequency", tag="B")+
  scale_y_continuous(limits=c(0, 60), breaks=seq(0,60,5))+
  #coord_flip()+
  theme_scientific()+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=11),#angle=90, vjust=0.5, hjust = 1), 
        legend.position = "none",
        #legend.position = c(0.78, 0.80),
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())#+
#guides(fill=guide_legend(ncol=6, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport/summary.pdf", width = 4, height = 7)
print(hh_intros_summary)
dev.off()

hh_intros_summary

#Generating the alluvium plot/flow diagram
p1<-hh_dta%>%
  filter(Origin!=Destination)%>%
  ggplot(aes(axis1 = Origin, axis2 = Destination)) +
  geom_alluvium(width = 1/8, aes(fill=EventType),decreasing = TRUE, alpha=0.65) +
  geom_stratum(width = 1/8, fill='white', color='black', decreasing = TRUE) +
  geom_label_repel(label.padding=0.25, size=3,stat = "stratum", min.y = 3,decreasing = TRUE, colour = "black", aes(label = after_stat(stratum)))+ ## Note: Can be uncommented to show country labels
  #scale_fill_manual(values=c(plot_color, wes_palette(2, name = "Darjeeling1", type = "discrete")[c(-4,-5)],"black",wes_palette(5, name = "Cavalcanti1", type = "discrete"),terrain.colors(2),"khaki","gray77",rainbow(10)[c(-1,-2,-4)],as.character(wes_palette(4, name = "Royal1", type = "discrete")),as.character(wes_palette(5, name = "Zissou1", type = "discrete")),"pink","gray32","tomato1"))+
  scale_fill_manual(values=plot_color[c(1, 7,19)])+
  scale_x_discrete(limits = c("Origin", "Destination"), expand = c(.08, .08)) +
  scale_y_continuous(limits = c(0,170), minor_breaks = seq(0,170, 10), breaks = seq(0,170,20))+
  labs(y="Number of Events", tag="A")+
  #theme_minimal()+
  theme_scientific() +
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.position = "bottom",
        #legend.position = c(0.80, 0.80),
        legend.key.size = unit(0.35, "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())+
  guides(fill=guide_legend(ncol=3, title = "Event Type", title.position = "left"), size=T)

pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport/import_export.pdf", width = 4.5, height =7)
print(p1)
dev.off()

p1

figure3 <- p1+hh_intros_summary
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/figure 3.pdf", width = 12, height =7.02)
print(figure3)
dev.off()
figure3

hh_dta%>%filter(Destination!="non-HH")%>%tabyl(Origin)%>%adorn_totals()

### analysis 2


transmission_dta <- read.csv("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport2/2023-01-04_mugration/annottated_tree_events.csv")%>%
  filter(Origin!="UNKNOWN")%>%
  separate(Origin, into=c("from_HH", "from_individual"), sep="_", remove=F)%>%tibble%>%
  separate(Destination, into=c("to_HH", "to_individual"), sep="_", remove=F)%>%
  mutate(EventType=case_when(from_HH==to_HH ~ "intra_household",
                             Origin=="non-HH" & str_detect(Destination, "_") ~ "importation",
                             from_HH!=to_HH &!is.na(from_HH) & !is.na(to_HH) ~"inter_household",
                             TRUE ~ "other"))%>%
  select(Origin, Destination, EventType, everything())%>%
  mutate(to_hh_no=as.integer(str_replace_all(to_HH, "HH", "")))%>%
  mutate(to_HH=fct_reorder(to_HH, to_hh_no))

tabyl(transmission_dta, EventType)%>%adorn_percentages()

glimpse(transmission_dta)



intra_hh <-transmission_dta%>%
  tabyl(to_HH,EventType)%>%
  filter(intra_household>0)%>%
  mutate(rownumber=row_number())

intra_hh

transmission_dta%>%tabyl(EventType)


intra_household <-intra_hh%>%ggplot()+
  geom_point(aes(x=to_HH, y=intra_household), size=3)+
  geom_segment(aes(x=to_HH, xend=to_HH, y=0, yend=intra_household))+
  labs(x="", y="No. of introductions")+
  #coord_flip()+
  theme_scientific()+
  scale_y_continuous(limits=c(0,8), breaks = seq(0, 8,1))+
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size=8,angle=90, vjust=0.5, hjust = 1), 
        legend.position = "none",
        #legend.position = c(0.78, 0.80),
        legend.key.size = unit(0.25, "cm"),
        legend.spacing.x = unit(0.15, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.text = element_text(size = 11),
        legend.title =element_text(size = 11),
        legend.background = element_rect(fill="#FFFFFF", color = NA),
        legend.box.background = element_blank())#+
#guides(fill=guide_legend(ncol=6, title = "Household",title.position = "top"), size=T)
pdf("~/Dropbox/COVID-19/HHSTUDY/phylogenetics/22Dec2022/all/importexport2/intra_household_patterns.pdf", width = 3.5, height = 2.5)
print(intra_household)
dev.off()



glimpse(transmission_dta)

from_HH <- transmission_dta%>%pull(from_HH)%>%unique()
str(from_HH)

from_individuals <- transmission_dta%>%pull(from_individual)%>%unique()
str(from_individuals)

to_HH <-transmission_dta%>%mutate(to_HH=as.character(to_HH))%>%pull(to_HH)%>%unique%>%as.list()
str(to_HH)


HH_analysed <- to_HH+from_HH

to_individuals <-transmission_dta%>%pull(to_individual)%>%unique
str(to_individuals)



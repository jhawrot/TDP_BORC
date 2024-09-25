library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggrepel)
library(shiny)
library(ggiraph)
library(shinyWidgets)
library(scales)
library(GO.db)
library(org.Hs.eg.db)
library(GOSim)
library(biomaRt)
library(readxl)

cellrox_x_i3N_meta_analysis_results <- read.csv("~/Downloads/cellrox_x_i3N_meta_analysis_results.csv")
liperfluo_x_i3N_meta_analysis_results <- read.csv("~/Downloads/liperfluo_x_i3N_meta_analysis_results.csv")

rox_novel <- annotated_cellrox_x_i3N_meta_analysis_results %>% 
  filter(novel_hit==1) %>% 
  dplyr::select(Gene)

lipo_novel <- annotated_liperfluo_x_i3N_meta_analysis_results %>% 
  filter(novel_hit==1) %>% 
  dplyr::select(Gene)

Combined<-cellrox_x_i3N_meta_analysis_results %>% 
  merge(liperfluo_x_i3N_meta_analysis_results, by="Gene",all=T)


Combined_mean<-cellrox_x_i3N_meta_analysis_results %>%
  group_by(Gene) %>% 
  dplyr::summarise(across(c(beta_1:fixed_effect_p_FDR_adjusted), mean,na.rm = TRUE)) %>% 
  merge(liperfluo_x_i3N_meta_analysis_results %>%
          group_by(Gene) %>% 
          dplyr::summarise(across(c(beta_1:fixed_effect_p_FDR_adjusted), mean,na.rm = TRUE)), by="Gene",all=T)
##check number of duplicates
Combined %>% 
  group_by(Gene) %>% 
  filter(n()>1)
Combined_mean %>% 
  group_by(Gene) %>% 
  filter(n()>1)

### add negative z scores
Combined_mean<-Combined_mean %>% 
  mutate(
    z_cell_rox = 
      case_when(
        beta_1.x<0 ~ z_score_1.x* -1,
        T ~ z_score_1.x
      ),
    z_i3N = 
      case_when(
        beta_2.x<0 ~ z_score_2.x* -1,
        T ~ z_score_2.x
      ),
    z_lipo = 
      case_when(
        beta_1.y<0 ~ z_score_1.y* -1,
        T ~ z_score_1.y
      ),
    z_i3N2 = 
      case_when(
        beta_2.y<0 ~ z_score_2.y* -1,
        T ~ z_score_2.y
      )
      
  )
  
  



  
validation_filtered <- read_excel("input/microscopy_validation_filtered.xlsx")

#validation_filtered[is.na(validation_filtered)]<-0
BORC<-toupper(borc)

Screen_1_a<-both %>% dplyr::select(
  Gene=id, LFC=lfc,pvalue=pvalue) %>% 
  mutate(Gene = toupper(Gene)) 
Screen_1_a<-Screen_1_a %>% 
  dplyr::filter(Gene%in%validation_filtered$Gene) %>% 
  mutate(Group="CRISPRi_Screen") %>% 
  merge(validation_filtered %>% 
          dplyr::select(Gene, Category))

x<-Combined_mean$fixed_effect_z.x
y<-Combined_mean$fixed_effect_z.y
cor(x,y, use='complete.obs')
metaZ<-Combined_mean %>% 
  mutate(Gene=toupper(Gene)) %>% 
  ggplot(aes(fixed_effect_z.x,fixed_effect_z.y,label=Gene))+
  geom_point(alpha=0.5, color="grey")+
  geom_point(data=Combined_mean %>% filter(Gene %in% lipo_novel$Gene), color="chartreuse3",alpha=0.8)+
  geom_point(data=Combined_mean %>% filter(Gene %in% rox_novel$Gene), color="chocolate2",alpha=0.8)+
  geom_point(data=Combined_mean %>% 
               mutate(Gene=toupper(Gene)) %>% 
               filter(Gene %in% BORC), color="deepskyblue2",alpha=0.8)+
  #geom_point(data=Combined_mean %>% 
  #             filter(Gene=="TARDBP"), color="red",alpha=0.8)+
  geom_label_repel(data=Combined_mean %>%
                     mutate(Gene=toupper(Gene)) %>% 
                     filter(Gene %in% BORC),
                   max.overlaps = Inf,box.padding = 0.5)+
  labs(title = "Meta analysis comparing Z",
       x="CellRox and TDP fixed Z-score",
       y="Liperfluo and TDP fixed Z-score",
       caption = "Novel hits from liperfluo in green
       Novel hits from cellrox in orange
       BORC hits in blue
       Correlation coefficient = 0.56")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
metaZ
metaZ %>% ggsave("output/metaZ.svg", ., height = 30, width = 30, units = "cm")
metaZ %>% ggsave("output/metaZ.tiff", ., height = 10, width = 10, units = "cm")


metaZ_cellrox<-Combined_mean %>% 
  mutate(Gene=toupper(Gene)) %>% 
  ggplot(aes(z_cell_rox,z_i3N,label=Gene))+
  geom_point(alpha=0.5, color="grey")+
  geom_point(data=Combined_mean %>% filter(Gene %in% lipo_novel$Gene), color="chartreuse3",alpha=0.5)+
  geom_point(data=Combined_mean %>% filter(Gene %in% rox_novel$Gene), color="chocolate2",alpha=0.5)+
  geom_point(data=Combined_mean %>% 
               mutate(Gene=toupper(Gene)) %>% 
               filter(Gene %in% BORC), color="deepskyblue2",alpha=0.5)+
  geom_label_repel(data=Combined_mean %>%
                     mutate(Gene=toupper(Gene)) %>% 
                     filter(Gene %in% BORC),
                   max.overlaps = Inf,box.padding = 0.5)+
  labs(title = "Meta analysis comparing Z",
       x="CellRox z-score",
       y="TDP z-score",
       caption = "Novel hits from liperfluo in green
       Novel hits from cellrox in orange
       BORC hits in blue")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
metaZ_cellrox
metaZ_cellrox %>% ggsave("output/metaZ_cellrox_1.svg", ., height = 30, width = 30, units = "cm")
metaZ_cellrox %>% ggsave("output/metaZ_cellrox_1.tiff", ., height = 10, width = 10)



metaZ_lipo<-Combined_mean %>% 
  mutate(Gene=toupper(Gene)) %>% 
  ggplot(aes(z_lipo,z_i3N,label=Gene))+
  geom_point(alpha=0.5, color="grey")+
  geom_point(data=Combined_mean %>% filter(Gene %in% lipo_novel$Gene), color="chartreuse3",alpha=0.5)+
  geom_point(data=Combined_mean %>% filter(Gene %in% rox_novel$Gene), color="chocolate2",alpha=0.5)+
  geom_point(data=Combined_mean %>% 
               mutate(Gene=toupper(Gene)) %>% 
               filter(Gene %in% BORC), color="deepskyblue2",alpha=0.5)+
  geom_label_repel(data=Combined_mean %>%
                     mutate(Gene=toupper(Gene)) %>% 
                     filter(Gene %in% BORC),
                   max.overlaps = Inf,box.padding = 0.5)+
  labs(title = "Meta analysis comparing Z",
       x="Lipfluo z-score",
       y="TDP z-score",
       caption = "Novel hits from liperfluo in green
       Novel hits from cellrox in orange
       BORC hits in blue")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
metaZ_lipo
metaZ_lipo %>% ggsave("output/metaZ_lipo_1.svg", ., height = 30, width = 30, units = "cm")
metaZ_lipo %>% ggsave("output/metaZ_lipo_1.tiff", ., height = 10, width = 10)



metaZ_lipo_rox<-Combined_mean %>% 
  mutate(Gene=toupper(Gene)) %>% 
  ggplot(aes(z_lipo,z_cell_rox,label=Gene))+
  geom_point(alpha=0.5, color="grey")+
  geom_point(data=Combined_mean %>% filter(Gene %in% lipo_novel$Gene), color="chartreuse3",alpha=0.5)+
  geom_point(data=Combined_mean %>% filter(Gene %in% rox_novel$Gene), color="chocolate2",alpha=0.5)+
  geom_point(data=Combined_mean %>% 
               mutate(Gene=toupper(Gene)) %>% 
               filter(Gene %in% BORC), color="deepskyblue2",alpha=0.5)+
  geom_label_repel(data=Combined_mean %>%
                     mutate(Gene=toupper(Gene)) %>% 
                     filter(Gene %in% BORC),
                   max.overlaps = Inf,box.padding = 0.5)+
  labs(title = "Meta analysis comparing Z",
       x="Lipfluo z-score",
       y="Cellrox z-score",
       caption = "Novel hits from liperfluo in green
       Novel hits from cellrox in orange
       BORC hits in blue")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
metaZ_lipo_rox
metaZ_lipo_rox %>% ggsave("output/metaZ_lipo_rox.svg", ., height = 30, width = 30, units = "cm")
metaZ_lipo_rox %>% ggsave("output/metaZ_lipo_rox.tiff", ., height = 10, width = 10)









Combined_mean %>% 
  ggplot(aes(fixed_effect_z.x,fixed_effect_z.y,label=Gene))+
  geom_point(alpha=0.5)+
  geom_point(data=Combined_mean %>% filter(Gene %in% Screen_1_a$Gene), color="chartreuse3",alpha=0.5)+
  geom_label_repel(data=Combined_mean %>%filter(Gene %in% Screen_1_a$Gene,
                                           abs(fixed_effect_z.x)>4|
                                             abs(fixed_effect_z.y)>4), color="chartreuse3",
                   max.overlaps = Inf)+
  labs(title = "Comparing Z",
       x="LFC (CellRox/TDP)",
       y="LFC (Liperfluo/TDP)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


Combined_mean %>% 
  ggplot(aes(fixed_effect_z.x,fixed_effect_z.y,label=Gene))+
  geom_point(alpha=0.5)+
  geom_point(data=Combined_mean %>% filter(abs(fixed_effect_z.x)>4&
                                             abs(fixed_effect_z.y)>4), color="chartreuse3",alpha=0.5)+
  geom_label_repel(data=Combined_mean %>%filter(fixed_effect_z.x < -4&
                                                  fixed_effect_z.y< -4), color="chartreuse3",
                   max.overlaps = Inf)+
  labs(title = "Comparing Z",
       x="LFC (CellRox/TDP)",
       y="LFC (Liperfluo/TDP)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

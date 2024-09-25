# all dots ----------------------------------------------------------------
# James Hawrot
# Dotplots for each group and ordered alphabetical with TDP43 at top
# 240531


#### all of the dots
library(tidyverse)
library(viridis)
library(readxl)
library(png)
validation_filtered <- read_excel("input/microscopy_validation_filtered.xlsx")

#validation_filtered[is.na(validation_filtered)]<-0
BORC<-toupper(borc)

Screen_1_a<-both %>% select(
  Gene=id, LFC=lfc,pvalue=pvalue) %>% 
  mutate(Gene = toupper(Gene)) 
Screen_1_a<-Screen_1_a %>% 
  filter(Gene%in%validation_filtered$Gene) %>% 
  mutate(Group="CRISPRi_Screen") %>% 
  merge(validation_filtered %>% 
          select(Gene, Category))

Halo_values_a<-validation_filtered %>%
  filter(Gene%in%validation_filtered$Gene) %>% 
  select(Gene,LFC=`log2 halo`,pvalue=`p halo`,Category) %>% 
  mutate(
    Group="Halo")

options(scipen = 999)

IF_values_a<-validation_filtered %>%
  filter(Gene%in%validation_filtered$Gene) %>% 
  select(Gene,LFC=`log2 IF`,pvalue=`p IF`,Category) %>% 
  mutate(
    Group="IF"
  ) 
all_values_a<-Screen_1_a %>% rbind(Halo_values_a) %>% 
  rbind(IF_values_a)


borc_lib <- read.csv("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/input/borc_lib.csv")

all_values_a<-all_values_a %>% merge(
  borc_lib %>% 
    mutate(gene=str_sub(gene,1,6)), by="Gene", all = T)



### change the names of the borc complex genes
all_values_a<- all_values_a %>%
  rowwise() %>% 
  mutate(Gene= ifelse(Gene%in%BORC,gene,Gene))
all_values_a<- all_values_a %>%
  mutate(Gene=as.factor(Gene))



# BORC Individual dot plots ----------------------------------------------------


level_order <- c('BORCS8', 'BORCS7',
                 'BORCS6', 'BORCS5',
                 'BORCS4', 'BORCS3',
                 'BORCS2', 'BORCS1',
                 'TARDBP') 
dotplot_BORC<-all_values_a %>%
  filter(Category=="BORC") %>% 
  ggplot() +
  geom_point(aes(x = Group, y = Gene, size = -log10(pvalue),fill=LFC),
             shape=21)+
  scale_fill_viridis(limits = c(-3, 3), oob = scales::squish,
                     option = "turbo")+
  #scale_fill_gradient(low = "black",high = "white")+
  scale_y_discrete(limits = level_order)+
  scale_size(range=c(2, 10),
             breaks = c(1.3,2,3,4),
             limits = c(0,4))+
  #facet_wrap(~Category)+
  labs(title="BORC dotplot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dotplot_BORC

dotplot_BORC %>% ggsave("dotplott_BORC.svg", ., height = 20, width = 20, units = "cm")
width=20
height=20
dotplot_BORC %>% ggsave("dotplot_BORC.tiff",., width = width, height = height, device='tiff', dpi=700)

# misc Individual dot plots ----------------------------------------------------

all_values_a %>% distinct(Category)
all_values_a %>% filter(Category=="misc") %>% 
  distinct(Gene)

level_order_misc <- c('ZYFVE26', 'UVRAG',
                      'TMEM192', 'RBM27',
                      'OGFOD1', 'OCRL',
                      'ATP6V0A1', 'ALG2',
                      'TARDBP') 
dotplot_misc<-all_values_a %>%
  filter(Category=="misc") %>% 
  ggplot() +
  geom_point(aes(x = Group, y = Gene, size = -log10(pvalue),fill=LFC),
             shape=21)+
  scale_fill_viridis(limits = c(-3, 3), oob = scales::squish,
                     option = "turbo")+
  #scale_fill_gradient(low = "black",high = "white")+
  scale_y_discrete(limits = level_order_misc)+
  scale_size(range=c(2, 10),
             breaks = c(1.3,2,3,4),
             limits = c(0,4))+
  #facet_wrap(~Category)+
  labs(title="misc dotplot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dotplot_misc

dotplot_misc %>% ggsave("dotplot_misc.svg", ., height = height, width = width, units = "cm")
dotplot_misc %>% ggsave("dotplot_misc.tiff",., width = width, height = height, device='tiff', dpi=700)

# m6a Individual dot plots ----------------------------------------------------

all_values_a %>% distinct(Category)
all_values_a %>% filter(Category=="m6A") %>% 
  distinct(Gene)

level_order_m6a <- c('ZC3H18', 'ZC3H13',
                     'YTHDF2', 'YTHDF1',
                     'SFPQ', 'METTL3',
                     'METTL21A', 'METTL14',
                     'METLL4','KIAA1429',
                     'CNOT3','CNOT1',
                     'TARDBP') 
dotplot_m6a<-all_values_a %>%
  filter(Category=="m6A") %>% 
  ggplot() +
  geom_point(aes(x = Group, y = Gene, size = -log10(pvalue),fill=LFC),
             shape=21)+
  scale_fill_viridis(limits = c(-3, 3), oob = scales::squish,
                     option = "turbo")+
  #scale_fill_gradient(low = "black",high = "white")+
  scale_y_discrete(limits = level_order_m6a)+
  scale_size(range=c(2, 10),
             breaks = c(1.3,2,3,4),
             limits = c(0,4))+
  #facet_wrap(~Category)+
  labs(title="m6a dotplot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dotplot_m6a

dotplot_m6a %>% ggsave("dotplot_m6a.svg", ., height = height, width = width, units = "cm")
dotplot_m6a %>% ggsave("dotplot_m6a.tiff",., width = width, height = height, device='tiff', dpi=700)


# mito Individual dot plots ----------------------------------------------------

all_values_a %>% distinct(Category)
all_values_a %>% filter(Category=="mito") %>% 
  distinct(Gene)

level_order_mito <- c('UQCRQ', 'UQCRC2',
                      'RPIA', 'PMPCB',
                      'PITRM1', 'NDUFA8',
                      'MRPL39', 'MFN2',
                      'GRSF1','FASTK',
                      'DTYMK','CPOX',
                      'COX6B1','COX10',
                      'TARDBP') 
dotplot_mito<-all_values_a %>%
  filter(Category=="mito") %>% 
  ggplot() +
  geom_point(aes(x = Group, y = Gene, size = -log10(pvalue),fill=LFC),
             shape=21)+
  scale_fill_viridis(limits = c(-3, 3), oob = scales::squish,
                     option = "turbo")+
  #scale_fill_gradient(low = "black",high = "white")+
  scale_y_discrete(limits = level_order_mito)+
  scale_size(range=c(2, 10),
             breaks = c(1.3,2,3,4),
             limits = c(0,4))+
  #facet_wrap(~Category)+
  labs(title="mito dotplot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dotplot_mito

dotplot_mito %>% ggsave("dotplot_mito.svg", ., height = height, width = width, units = "cm")
dotplot_mito %>% ggsave("dotplot_mito.tiff",., width = width, height = height, device='tiff', dpi=700)

# Ubiquitin Individual dot plots ----------------------------------------------------

all_values_a %>% distinct(Category)
all_values_a %>% filter(Category=="Ubiquitin") %>% 
  distinct(Gene)

level_order_Ubiquitin <- c('WDR26', 'UCHL5',
                           'UBE2I', 'UBE2H',
                           'UBC', 'UBA3',
                           'UBA2', 'SKP2',
                           'SAE1','RMND5A',
                           'RANBP9','NEDD8',
                           'NAE1','MAEA',
                           'GID8','CUL1',
                           'TARDBP') 
dotplot_Ubiquitin<-all_values_a %>%
  filter(Category=="Ubiquitin") %>% 
  ggplot() +
  geom_point(aes(x = Group, y = Gene, size = -log10(pvalue),fill=LFC),
             shape=21)+
  scale_fill_viridis(limits = c(-3, 3), oob = scales::squish,
                     option = "turbo")+
  #scale_fill_gradient(low = "black",high = "white")+
  scale_y_discrete(limits = level_order_Ubiquitin)+
  scale_size(range=c(2, 10),
             breaks = c(1.3,2,3,4),
             limits = c(0,4))+
  #facet_wrap(~Category)+
  labs(title="Ubiquitin dotplot")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dotplot_Ubiquitin

dotplot_Ubiquitin %>% ggsave("dotplot_Ubiquitin.svg", ., height = height, width = width, units = "cm")
dotplot_Ubiquitin %>% ggsave("dotplot_Ubiquitin.tiff",., width = width, height = height, device='tiff', dpi=700)

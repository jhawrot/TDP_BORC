library(tidyverse)
library(ggrepel)
library(svglite)
library(plotly)

rra.i3N.gene_summary <- read.delim("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/input/CP6277_first/rra-i3N.gene_summary.txt")
rra.i3N.sgrna_summary <- read.delim("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/input/CP6277_first/rra-i3N.sgrna_summary.txt")

borc<-c('BLOC1S1',"BLOC1S2","SNAPIN","KXD1","LOH12CR1","C17orf59","C10orf32","MEF2BNB")
##transcript plot - need to add noise modeling too
rra.i3N.sgrna_summary %>% 
  ggplot(aes(control_mean,treat_mean, label=Gene))+
  geom_point(data = rra.i3N.sgrna_summary %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.i3N.sgrna_summary %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.7)+
  geom_point(data = rra.i3N.sgrna_summary %>% 
               filter(Gene=="TARDBP"|
                        Gene%in%borc),
             color="purple",
             alpha=0.7)+
  #geom_point(data = rra.i3N.sgrna_summary %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$id),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = rra.i3N.sgrna_summary %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$id),
  #           color="green",
#           alpha=0.7)+
geom_label_repel(data = rra.i3N.sgrna_summary %>% 
                   filter(Gene=="TARDBP"|
                            Gene%in%borc)
                 ,
                 box.padding = 0.5,
                 segment.color="black",
                 alpha=0.6,
                 max.overlaps = Inf
)+
  scale_x_log10()+
  scale_y_log10()+
  coord_fixed(1)+
  labs(
    title = "Halo TDP-43 Expression Screen",
    x="KD decreases TDP43 (median transcript count)",
    y="KD increases TDP43 (median transcript count)",
    caption = 
      "Grey: NT guides"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


neg<-rra.i3N.gene_summary %>% 
  filter(neg.lfc<=0) %>% 
  mutate(lfc=neg.lfc,
         fdr=neg.fdr,
         score = neg.score,
         pvalue = neg.p.value) 
neg<-neg%>% 
  select(id,num,lfc,neg.rank,fdr,score,pvalue)

pos<-rra.i3N.gene_summary %>% 
  filter(neg.lfc>0) %>% 
  mutate(lfc=pos.lfc,
         fdr=pos.fdr,
         score = pos.score,
         pvalue = pos.p.value) 
pos<-pos%>% 
  select(id,num,lfc,neg.rank,fdr,score,pvalue)

both <- neg %>% rbind(pos)

both %>% 
  ggplot(aes(lfc,-log(pvalue)))+
  geom_point(data = both %>% 
               filter(fdr<0.05),
             color="black",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter(fdr>=0.05),
  #           color="grey",
  #           alpha=0.2)+
  geom_point(data = both %>% 
               filter(id=="TARDBP"|
                        id%in%borc),
             color="purple",
             alpha=0.7)+
  theme_classic()


rra.i3N.sgrna_summary_filtered<-rra.i3N.sgrna_summary %>%
  filter(control_count>100|
     treatment_count>100) %>% 
  arrange(LFC)
rra.i3N.sgrna_summary_filtered<-rra.i3N.sgrna_summary_filtered %>% 
  mutate(
    zscore = scale(LFC),
    row = as.numeric(rownames(rra.i3N.sgrna_summary_filtered))
  )
  


##update gene name
rra.i3N.sgrna_summary_filtered_a<-rra.i3N.sgrna_summary_filtered %>%
  mutate(Gene=case_when(Gene=="BLOC1S1"~"BORCS1",
                     Gene=="BLOC1S2"~"BORCS2",
                     Gene=="SNAPIN"~"BORCS3",
                     Gene=="KXD1"~"BORCS4",
                     Gene=="LOH12CR1"~"BORCS5",
                     Gene=="C17orf59"~"BORCS6",
                     Gene=="C10orf32"~"BORCS7",
                     Gene=="MEF2BNB"~"BORCS8",
                     TRUE~Gene)) %>% 
  rbind(rra.i3N.sgrna_summary_filtered %>% 
          filter(Gene=="negative_control"))

subset(rra.i3N.sgrna_summary_filtered, !(sgrna %in% rra.i3N.sgrna_summary_filtered_a$sgrna))
subset(rra.i3N.sgrna_summary_filtered_a, !(sgrna %in% rra.i3N.sgrna_summary_filtered$sgrna))

## Rank plot i3N filtered ----
rank_i3N_T<-rra.i3N.sgrna_summary_filtered_a %>% 
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.i3N.sgrna_summary_filtered_a %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.i3N.sgrna_summary_filtered_a %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  #geom_point(data = rra.i3N.sgrna_summary_filtered_a %>% 
  #             filter(Gene=="TARDBP"|
  #                      grepl("BORCS",Gene)),
  #           color="purple",
  #           alpha=0.7)+
  geom_point(data = rra.i3N.sgrna_summary_filtered_a %>% 
               filter(Gene=="TARDBP"),
             color="purple",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$Gene),
  #           color="green",
#           alpha=0.7)+
  #geom_label_repel(data = rra.i3N.sgrna_summary_filtered_a %>% 
  #                 filter(grepl("BORCS",Gene)),
  #               fill="purple"
  #               ,
  #               box.padding = 1,
  #               segment.color="black",
  #               alpha=0.6,
  #               max.overlaps = Inf)+
  geom_label_repel(data = rra.i3N.sgrna_summary_filtered_a %>% 
                     filter(Gene=="TARDBP"),
                   fill="purple"
                   ,
                   box.padding = 1,
                   segment.color="black",
                   alpha=0.6,
                   max.overlaps = Inf)+
  labs(
    title = "TDP-43 rank plot filtered at 100",
    caption = 
      "Grey: NT guides
    Purple: TARDBP"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
rank_i3N_T


rank_i3N_T %>% ggsave("rank_i3N_T.svg", ., height = 10, width = 20, units = "cm")

rank_i3N_BORC<-rra.i3N.sgrna_summary_filtered_a %>% 
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.i3N.sgrna_summary_filtered_a %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.i3N.sgrna_summary_filtered_a %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  geom_point(data = rra.i3N.sgrna_summary_filtered_a %>% 
               filter(grepl("BORCS",Gene)),
             color="deepskyblue2",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$Gene),
  #           color="green",
#           alpha=0.7)+
geom_label_repel(data = rra.i3N.sgrna_summary_filtered_a %>% 
                 filter(grepl("BORCS",Gene)),
               fill="deepskyblue2"
               ,
               box.padding = 1,
               segment.color="black",
               alpha=0.6,
               max.overlaps = Inf)+
  labs(
    title = "BORC Complex rank plot filtered at 100",
    caption = 
      "Grey: NT guides
    Blue: BORC Complex"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
rank_i3N_BORC


rank_i3N_BORC %>% ggsave("rank_i3N_BORC.svg", ., height = 10, width = 20, units = "cm")


##### ipscs
# Rank plot iPSC ----

rra.ipsc.gene_summary <- read.delim("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/input/CP6277_first/rra-ipsc.gene_summary.txt")
rra.ipsc.sgrna_summary <- read.delim("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/input/CP6277_first/rra-ipsc.sgrna_summary.txt")



rra.ipsc.sgrna_summary_filtered<-rra.ipsc.sgrna_summary %>%
  filter(control_count>100|
           treatment_count>100) %>% 
  arrange(LFC)
rra.ipsc.sgrna_summary_filtered<-rra.ipsc.sgrna_summary_filtered %>% 
  mutate(
    zscore = scale(LFC),
    row = as.numeric(rownames(rra.ipsc.sgrna_summary_filtered))
  )



##update gene name
rra.ipsc.sgrna_summary_filtered_a<-rra.ipsc.sgrna_summary_filtered %>%
  mutate(Gene=case_when(Gene=="BLOC1S1"~"BORCS1",
                        Gene=="BLOC1S2"~"BORCS2",
                        Gene=="SNAPIN"~"BORCS3",
                        Gene=="KXD1"~"BORCS4",
                        Gene=="LOH12CR1"~"BORCS5",
                        Gene=="C17orf59"~"BORCS6",
                        Gene=="C10orf32"~"BORCS7",
                        Gene=="MEF2BNB"~"BORCS8",
                        TRUE~Gene)) %>% 
  rbind(rra.ipsc.sgrna_summary_filtered %>% 
          filter(Gene=="negative_control"))



rank_ipsc_T<-rra.ipsc.sgrna_summary_filtered_a %>% 
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.ipsc.sgrna_summary_filtered_a %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered_a %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  #geom_point(data = rra.ipsc.sgrna_summary_filtered_a %>% 
  #             filter(Gene=="TARDBP"|
  #                      grepl("BORCS",Gene)),
  #           color="purple",
  #           alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered_a %>% 
               filter(Gene=="TARDBP"),
             color="purple",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$Gene),
  #           color="green",
#           alpha=0.7)+
#geom_label_repel(data = rra.ipsc.sgrna_summary_filtered_a %>% 
#                   filter(grepl("BORCS",Gene)),
#                 fill="purple"
#                 ,
#                 box.padding = 1,
#                 segment.color="black",
#                 alpha=0.6,
#                 max.overlaps = Inf)+
  geom_label_repel(data = rra.ipsc.sgrna_summary_filtered_a %>% 
                     filter(Gene=="TARDBP"),
                   fill="purple"
                   ,
                   box.padding = 1,
                   segment.color="black",
                   alpha=0.6,
                   max.overlaps = Inf)+
  labs(
    title = "TDP-43 rank plot iPSC filtered at 100",
    caption = 
      "Grey: NT guides
    Purple: TARDBP KD"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))
rank_ipsc_T
ggsave(file="test.svg", plot=rank_ipsc_T, width=20, height=16)














borc<-c('BLOC1S1',"BLOC1S2","SNAPIN","KXD1","LOH12CR1","C17orf59","C10orf32","MEF2BNB")
##transcript plot - need to add noise modeling too
rra.ipsc.sgrna_summary %>% 
  ggplot(aes(control_mean,treat_mean, label=Gene))+
  geom_point(data = rra.ipsc.sgrna_summary %>% 
               filter(Gene!="negative_control"),
             color="black")+
  geom_point(data = rra.ipsc.sgrna_summary %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary %>% 
               filter(Gene=="TARDBP"|
                        Gene%in%borc),
             color="purple",
             alpha=0.7)+
  #geom_point(data = rra.ipsc.sgrna_summary %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$id),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = rra.ipsc.sgrna_summary %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$id),
  #           color="green",
#           alpha=0.7)+
geom_label_repel(data = rra.ipsc.sgrna_summary %>% 
                   filter(Gene=="TARDBP"|
                            Gene%in%borc)
                 ,
                 box.padding = 0.5,
                 segment.color="black",
                 alpha=0.6,
                 max.overlaps = Inf
)+
  scale_x_log10()+
  scale_y_log10()+
  coord_fixed(1)+
  labs(
    title = "Halo TDP-43 Expression Screen - iPSCs",
    x="KD decreases TDP43 (median transcript count)",
    y="KD increases TDP43 (median transcript count)",
    caption = 
      "Grey: NT guides"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))


neg<-rra.ipsc.gene_summary %>% 
  filter(neg.lfc<=0) %>% 
  mutate(lfc=neg.lfc,
         fdr=neg.fdr,
         score = neg.score,
         pvalue = neg.p.value) 
neg<-neg%>% 
  select(id,num,lfc,neg.rank,fdr,score,pvalue)

pos<-rra.ipsc.gene_summary %>% 
  filter(neg.lfc>0) %>% 
  mutate(lfc=pos.lfc,
         fdr=pos.fdr,
         score = pos.score,
         pvalue = pos.p.value) 
pos<-pos%>% 
  select(id,num,lfc,neg.rank,fdr,score,pvalue)

both <- neg %>% rbind(pos)

both %>% 
  ggplot(aes(lfc,-log(pvalue)))+
  geom_point(data = both %>% 
               filter(fdr<0.05),
             color="black",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter(fdr>=0.05),
  #           color="grey",
  #           alpha=0.2)+
  geom_point(data = both %>% 
               filter(id=="TARDBP"|
                        id%in%borc),
             color="purple",
             alpha=0.7)+
  theme_classic()


rra.ipsc.sgrna_summary_filtered<-rra.ipsc.sgrna_summary %>%
  filter(control_count>100|
           treatment_count>100) %>% 
  arrange(LFC)
rra.ipsc.sgrna_summary_filtered<-rra.ipsc.sgrna_summary_filtered %>% 
  mutate(
    zscore = scale(LFC),
    row = as.numeric(rownames(rra.ipsc.sgrna_summary_filtered))
  )


## Rank plot iPSC filtered ----

rra.ipsc.sgrna_summary_filtered %>%
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
               filter(Gene=="TARDBP"|
                        Gene%in%borc),
             color="purple",
             alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
               filter(Gene=="TARDBP"),
             color="green",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$Gene),
  #           color="green",
#           alpha=0.7)+
geom_label_repel(data = rra.ipsc.sgrna_summary_filtered %>% 
                   filter(Gene%in%borc),
                 fill="purple"
                 ,
                 box.padding = 1,
                 segment.color="black",
                 alpha=0.6,
                 max.overlaps = Inf)+
  geom_label_repel(data = rra.ipsc.sgrna_summary_filtered %>% 
                     filter(Gene=="TARDBP"),
                   fill="green"
                   ,
                   box.padding = 0.5,
                   segment.color="black",
                   alpha=0.6,
                   max.overlaps = Inf)+
  labs(
    title = "TDP-43 iPSC rank plot filtered at 100 - iPSCs",
    caption = 
      "Grey: NT guides
    Purple: BORC complex subunits"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

rank_iPSC_T<-rra.ipsc.sgrna_summary_filtered %>% 
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  #geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
  #             filter(Gene=="TARDBP"|
  #                      grepl("BORCS",Gene)),
  #           color="purple",
  #           alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered %>% 
               filter(Gene=="TARDBP"),
             color="purple",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$Gene),
  #           color="green",
#           alpha=0.7)+
#geom_label_repel(data = rra.ipsc.sgrna_summary_filtered %>% 
#                 filter(grepl("BORCS",Gene)),
#               fill="purple"
#               ,
#               box.padding = 1,
#               segment.color="black",
#               alpha=0.6,
#               max.overlaps = Inf)+
geom_label_repel(data = rra.ipsc.sgrna_summary_filtered %>% 
                   filter(Gene=="TARDBP"),
                 fill="purple"
                 ,
                 box.padding = 1,
                 segment.color="black",
                 alpha=0.6,
                 max.overlaps = Inf)+
  labs(
    title = "TDP-43 rank plot iPSC filtered at 100",
    caption = 
      "Grey: NT guides
    Purple: TARDBP"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
rank_iPSC_T

rank_iPSC_T %>% ggsave("rank_iPSC_T.svg", ., height = 10, width = 20, units = "cm")


rra.ipsc.sgrna_summary_filtered_a<-rra.ipsc.sgrna_summary_filtered %>%
  mutate(Gene=case_when(Gene=="BLOC1S1"~"BORCS1",
                        Gene=="BLOC1S2"~"BORCS2",
                        Gene=="SNAPIN"~"BORCS3",
                        Gene=="KXD1"~"BORCS4",
                        Gene=="LOH12CR1"~"BORCS5",
                        Gene=="C17orf59"~"BORCS6",
                        Gene=="C10orf32"~"BORCS7",
                        Gene=="MEF2BNB"~"BORCS8",
                        TRUE~Gene)) %>% 
  rbind(rra.ipsc.sgrna_summary_filtered %>% 
          filter(Gene=="negative_control"))

rank_iPSC_BORC<-rra.ipsc.sgrna_summary_filtered_a %>% 
  ggplot(aes(row,zscore, label=Gene))+
  geom_point(data = rra.ipsc.sgrna_summary_filtered_a %>% 
               filter(Gene!="negative_control"),
             color="black",
             alpha=0.7)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered_a %>% 
               filter(Gene=="negative_control"),
             color="grey",
             alpha=0.2)+
  geom_point(data = rra.ipsc.sgrna_summary_filtered_a %>% 
               filter(grepl("BORCS",Gene)),
             color="deepskyblue2",
             alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%transcription$Gene),
  #           color="blue",
  #           alpha=0.7)+
  #geom_point(data = both %>% 
  #             filter((control_mean>=3000|
  #                       treat_mean>=3000),
  #                    Gene%in%translation$Gene),
  #           color="green",
#           alpha=0.7)+
geom_label_repel(data = rra.ipsc.sgrna_summary_filtered_a %>% 
                 filter(grepl("BORCS",Gene)),
               fill="deepskyblue2"
               ,
               box.padding = 1,
               segment.color="black",
               alpha=0.6,
               max.overlaps = Inf)+
  labs(
    title = "BORC Complex rank plot iPSC filtered at 100",
    caption = 
      "Grey: NT guides
    Light Blue: BORC Complex"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
rank_iPSC_BORC

rank_iPSC_BORC %>% ggsave("rank_iPSC_BORC.svg", ., height = 10, width = 20, units = "cm")




# Import meta results -----------------------------------------------------




library(tidyverse)


liperfluo_x_i3N_meta_analysis_results <- read.csv("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/meta/liperfluo_x_i3N_meta_analysis_results.csv")
annotated_liperfluo_x_i3N_meta_analysis_results <- read.csv("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/meta/annotated_liperfluo_x_i3N_meta_analysis_results.csv")
annotated_cellrox_x_i3N_meta_analysis_results <- read.csv("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/meta/annotated_cellrox_x_i3N_meta_analysis_results.csv")

liperfluo_x_i3N_meta_analysis_results %>% 
  summary()
annotated_liperfluo_x_i3N_meta_analysis_results %>% 
  summary()
annotated_liperfluo_x_i3N_meta_analysis_results %>% 
  filter(novel_hit==1) %>% 
  summary()


#BiocManager::install("tools")

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)
library(msigdbr)
library(patchwork)
library(ggnewscale)


novel<- annotated_liperfluo_x_i3N_meta_analysis_results %>% 
  dplyr::distinct(Gene,.keep_all = T)%>%
  filter(novel_hit==1)
#file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
#                 "testdata/rra.gene_summary.txt")
#gdata = ReadRRA(file1)

annotated_liperfluo_x_i3N_meta_analysis_results %>% 
  filter(novel_hit==1)

gdata_lipo_up<-annotated_liperfluo_x_i3N_meta_analysis_results %>% 
  dplyr::distinct(Gene,.keep_all = T)%>%
  filter(novel_hit==1,
         fixed_effect_z> 0)
genelist_lipo_up = gdata_lipo_up$fixed_effect_z
names(genelist_lipo_up) = gdata_lipo_up$Gene
genelist_lipo_up[1:5]
names(genelist_lipo_up)


gdata_lipo_down<-annotated_liperfluo_x_i3N_meta_analysis_results %>% 
  dplyr::distinct(Gene,.keep_all = T)%>%
  filter(novel_hit==1,
         fixed_effect_z< 0)
genelist_lipo_down = gdata_lipo_down$fixed_effect_z
names(genelist_lipo_down) = gdata_lipo_down$Gene
genelist_lipo_down[1:5]
names(genelist_lipo_down)

###Hypergeometric test
# Alternative functions EnrichAnalyzer and enrich.HGT.
hgtRes1 = EnrichAnalyzer(genelist_lipo_down, method = "HGT")
head(hgtRes1@result)
# hgtRes2 = enrich.HGT(genelist[genelist< -z])


hgtRes1_up = EnrichAnalyzer(genelist_lipo_up, method = "HGT")
head(hgtRes1_up@result)
# hgtRes2 = enrich.HGT(genelist[genelist< -z])

###over-representation test
# Alternative functions EnrichAnalyzer and enrich.ORT.
ortRes1 = EnrichAnalyzer(genelist_lipo_down, method = "ORT")
head(ortRes1@result)
# ortRes2 = enrich.ORT(genelist[genelist< -z])


###other plots
hgtRes1_up@result$geneID = hgtRes1_up@result$geneName
hgtRes1@result$geneID = hgtRes1@result$geneName

###Remove redundant results using EnrichedFilter.
enrich1 = EnrichAnalyzer(genelist_lipo_down, type = "GOMF+GOBP+GOCC")
metaL_go_all_down<-EnrichedView(enrich1, bottom = 10)+
  scale_y_discrete(limits=rev)+
  labs(title = "Lipofluo Meta: TDP43 negative GOALL")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

enrich2 = EnrichAnalyzer(genelist_lipo_up, type = "GOMF+GOBP+GOCC")
metaL_go_all_up<-EnrichedView(enrich2, top = 10)+
  scale_y_discrete(limits=rev)+
  labs(title = "Lipofluo Meta: TDP43 positive GOALL")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

width=20
height=20
metaL_go_all_down %>% ggsave("output/metaL_go_all_down.svg", ., height = 20, width = 30, units = "cm")
metaL_go_all_down %>% ggsave("output/metaL_go_all_down.tiff",., width = 30, height = 20, device='tiff')

metaL_go_all_up %>% ggsave("output/metaL_go_all_up.svg", ., height = 20, width = 30, units = "cm")
metaL_go_all_up %>% ggsave("output/metaL_go_all_up.tiff",., width = 30, height = 20, device='tiff')






### patchwork
a1<-barplot(hgtRes1, showCategory = 10, title = "TDP43 negative modulating genes enrichplot")
a2<-barplot(hgtRes1_up, showCategory = 10, title = "TDP43 positive modulating genes enrichplot")
b1<-dotplot(hgtRes1, showCategory = 10)+
  labs(title = "TDP43 negative modulating HGT analyis")
b2<-dotplot(hgtRes1_up, showCategory = 10)+
  labs(title = "TDP43 positive modulating HGT analyis")

c1<-cnetplot(hgtRes1, 5, 
             foldChange= genelist, 
             circular = TRUE, colorEdge = TRUE,
             max.overlaps = Inf)+
  labs(title = "TDP43 negative modulating HGT analyis")
c2<-cnetplot(hgtRes1_up, 5, foldChange= genelist, circular = TRUE, colorEdge = TRUE)+
  labs(title = "TDP43 positive modulating HGT analyis")
?cnetplot

d1<-EnrichedView(enrichPro, bottom = 10)+
  labs(title = "TDP43 negative modulating CORUM analyis")
d2<-EnrichedView(enrichPro2, top = 10)+
  labs(title = "TDP43 positive modulating genes CORUM analyis")

p<-a1+a2+b1+b2+
  plot_layout(ncol = 2)

?patchwork
p
a1
a2
b1
b1 %>% ggsave("output/HGT_metaL_negative.svg", ., height = width, width = height, units = "cm")
b1 %>% ggsave("output/HGT_metaL_negative.tiff",., width = 20, height = 10, device='tiff')


b2
b2 %>% ggsave("output/HGT_metaL_positive.svg", ., height = width, width = height, units = "cm")
b2 %>% ggsave("output/HGT_metaL_positive.tiff",., width = 20, height = 10, device='tiff')

c1
c1 %>% ggsave("output/HGT_metaL_negative.svg", ., height = width, width = height, units = "cm")
c1 %>% ggsave("output/HGT_metaL_negative.tiff",., width = 20, height = 20, device='tiff')

c2
c2 %>% ggsave("output/HGT_metaL_positive.svg", ., height = width, width = height, units = "cm")
c2 %>% ggsave("output/HGT_metaL_positive.tiff",., width = 20, height = 20, device='tiff')

d1
d2




# Cellrox -----------------------------------------------------------------


annotated_cellrox_x_i3N_meta_analysis_results %>% 
  filter(novel_hit==1)

gdata_cellrox_up<-annotated_cellrox_x_i3N_meta_analysis_results %>% 
  dplyr::distinct(Gene,.keep_all = T)%>%
  filter(novel_hit==1,
         fixed_effect_z> 0)
genelist_cellrox_up = gdata_cellrox_up$fixed_effect_z
names(genelist_cellrox_up) = gdata_cellrox_up$Gene
genelist_cellrox_up[1:5]
names(genelist_cellrox_up)


gdata_cellrox_down<-annotated_cellrox_x_i3N_meta_analysis_results %>% 
  dplyr::distinct(Gene,.keep_all = T)%>%
  filter(novel_hit==1,
         fixed_effect_z< 0)
genelist_cellrox_down = gdata_cellrox_down$fixed_effect_z
names(genelist_cellrox_down) = gdata_cellrox_down$Gene
genelist_cellrox_down[1:5]
names(genelist_cellrox_down)

###Hypergeometric test
# Alternative functions EnrichAnalyzer and enrich.HGT.
hgtRes1 = EnrichAnalyzer(genelist_cellrox_down, method = "HGT")
head(hgtRes1@result)
# hgtRes2 = enrich.HGT(genelist[genelist< -z])


hgtRes1_up = EnrichAnalyzer(genelist_cellrox_up, method = "HGT")
head(hgtRes1_up@result)
# hgtRes2 = enrich.HGT(genelist[genelist< -z])

###over-representation test
# Alternative functions EnrichAnalyzer and enrich.ORT.
ortRes1 = EnrichAnalyzer(genelist_cellrox_down, method = "ORT")
head(ortRes1@result)
# ortRes2 = enrich.ORT(genelist[genelist< -z])


###other plots
hgtRes1_up@result$geneID = hgtRes1_up@result$geneName
hgtRes1@result$geneID = hgtRes1@result$geneName

###Remove redundant results using EnrichedFilter.
enrich1 = EnrichAnalyzer(genelist_cellrox_down, type = "GOMF+GOBP+GOCC")
metaC_go_all_down<-EnrichedView(enrich1, bottom = 10)+
  scale_y_discrete(limits=rev)+
  labs(title = "Cellrox Meta: TDP43 negative GOALL")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

enrich2 = EnrichAnalyzer(genelist_cellrox_up, type = "GOMF+GOBP+GOCC")
metaC_go_all_up<-EnrichedView(enrich2, top = 10)+
  scale_y_discrete(limits=rev)+
  labs(title = "Cellrox Meta: TDP43 positive GOALL")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

metaC_go_all_down %>% ggsave("output/metaC_go_all_down.svg", ., height = 20, width = 30, units = "cm")
metaC_go_all_down %>% ggsave("output/metaC_go_all_down.tiff",., width = 30, height = 20, device='tiff', dpi=700)

metaC_go_all_up %>% ggsave("output/metaC_go_all_up.svg", ., height = 20, width = 30, units = "cm")
metaC_go_all_up %>% ggsave("output/metaC_go_all_up.tiff",., width = 30, height = 20, device='tiff', dpi=700)






### patchwork
a1<-barplot(hgtRes1, showCategory = 10, title = "TDP43 negative modulating genes enrichplot")
a2<-barplot(hgtRes1_up, showCategory = 10, title = "TDP43 positive modulating genes enrichplot")
b1<-dotplot(hgtRes1, showCategory = 10)+
  labs(title = "TDP43 negative modulating HGT analyis")
b2<-dotplot(hgtRes1_up, showCategory = 10)+
  labs(title = "TDP43 positive modulating HGT analyis")

c1<-cnetplot(hgtRes1, 5, 
             foldChange= genelist, 
             circular = TRUE, colorEdge = TRUE,
             max.overlaps = Inf)+
  labs(title = "TDP43 negative modulating HGT analyis")
c2<-cnetplot(hgtRes1_up, 5, foldChange= genelist, circular = TRUE, colorEdge = TRUE)+
  labs(title = "TDP43 positive modulating HGT analyis")
?cnetplot

d1<-EnrichedView(enrichPro, bottom = 10)+
  labs(title = "TDP43 negative modulating CORUM analyis")
d2<-EnrichedView(enrichPro2, top = 10)+
  labs(title = "TDP43 positive modulating genes CORUM analyis")

p<-a1+a2+b1+b2+
  plot_layout(ncol = 2)

?patchwork
p
a1
a2
b1
b1 %>% ggsave("output/HGT_metaC_negative.svg", ., height = width, width = height, units = "cm")
b1 %>% ggsave("output/HGT_metaC_negative.tiff",., width = 20, height = 10, device='tiff')


b2
b2 %>% ggsave("output/HGT_metaC_positive.svg", ., height = width, width = height, units = "cm")
b2 %>% ggsave("output/HGT_metaC_positive.tiff",., width = 20, height = 10, device='tiff')

c1
c1 %>% ggsave("output/HGT_metaC_negative.svg", ., height = width, width = height, units = "cm")
c1 %>% ggsave("output/HGT_metaC_negative.tiff",., width = 20, height = 20, device='tiff')

c2
c2 %>% ggsave("output/HGT_metaC_positive.svg", ., height = width, width = height, units = "cm")
c2 %>% ggsave("output/HGT_metaC_positive.tiff",., width = 20, height = 20, device='tiff')

d1
d2

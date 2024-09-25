#BiocManager::install("tools")

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)
library(msigdbr)
library(patchwork)
library(ggnewscale)

#file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
 #                 "testdata/rra.gene_summary.txt")
#gdata = ReadRRA(file1)


dSILAC <- read.csv("~/Desktop/Projects/Others/Veronica/Halo_TDP43_screens/input/dSILAC_halflife_hits.csv",header = T)
dSILAC<-dSILAC %>% 
  filter(P.value.Day.7.KD.NT<= 0.05)


genelist = dSILAC$Day.7.KD.NT
names(genelist) = dSILAC$Gene.names
genelist[1:5]
names(genelist)



#####go terms and pathways
## KEGG and REACTOME pathways

enrich_up = EnrichAnalyzer(geneList = genelist[genelist>1], type = "KEGG+REACTOME")
EnrichedView(enrich_up, top = 10)

## Only KEGG pathways
enrich_up = EnrichAnalyzer(geneList = genelist[genelist>1], type = "KEGG")
EnrichedView(enrich_up, bottom = 5)


enrich2 = EnrichAnalyzer(genelist[genelist> 1], type = "KEGG")
dSILAC_kegg<-EnrichedView(enrich2, top = 10)+
  labs(title="dSILAC KEGG pathway analysis")+
  theme_classic()

dSILAC_kegg %>% ggsave("output/dSILAC_kegg.svg", ., height = 20, width = 30, units = "cm")
dSILAC_kegg %>% ggsave("output/dSILAC_kegg.tiff",., width = 30, height = 20, device='tiff', dpi=700)


enrich2 = EnrichAnalyzer(genelist[genelist> 1], type = "REACTOME")
dSILAC_reactome<-EnrichedView(enrich2, top = 10)+
  labs(title="dSILAC REACTOME pathway analysis")+
  theme_classic()

dSILAC_reactome %>% ggsave("output/dSILAC_reactome.svg", ., height = 20, width = 30, units = "cm")
dSILAC_reactome %>% ggsave("output/dSILAC_reactome.tiff",., width = 30, height = 20, device='tiff', dpi=700)



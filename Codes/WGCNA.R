library(tidyverse)
library(purrrlyr)
library(gplots)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("GO.db")
# BiocManager::install("impute")
# BiocManager::install("preprocessCore")
# BiocManager::install("RcppArmadillo")

#install.packages("WGCNA")
library(WGCNA)
# BiocManager::install("sva")
library(sva)
#install.packages("flashClust")
library(flashClust)
library(cluster)
library(igraph)
#options(stringsAsFactors = FALSE);
allowWGCNAThreads(10) #enableWGCNAThreads() doesn't work

dir <- "/$HOME/PROJECTDIR"
df_c <- readr::read_delim(paste0(dir, "/Data/Diff_expression_all_comparisons.tsv"), delim="\t") %>%
  rename(Geneid = "Row-names")

df_tpm <- readr::read_delim(paste0(dir, "/Data/Mouse_RNAseq_tpm.tsv"), delim="\t") %>%
  inner_join(df_c %>% select(Geneid, Symbol), by = "Geneid") %>%
  dplyr::select(-TT_amygdala_4_17, Geneid, Symbol) %>%
  dplyr::select(contains("amygdala"), Geneid, Symbol)

genelist2 <- read.table(paste0(dir, "/Data/All_gene_lists_GMT_47samples_intersect_amygdala.txt"), sep="\t", head=F) %>%
  rename(number = V1, comp = V2, Geneid = V3, symbol = V4)

df_exp <- df_tpm
#### WGCNA ####

# Remove confounding variable effects
df_meta <- as.data.frame(colnames(df_tpm)[1:15]) %>%
  rename(sample = 'colnames(df_tpm)[1:15]') %>%
  mutate(genotype =  str_split(sample, "_")  %>% map_chr(., 1) %>% as.factor(),
         age =  str_split(sample, "_")  %>% map_chr(., 3) %>% paste(., "month") %>% as.factor(),
         batch = 1) %>%
  column_to_rownames(var = "sample")
annotation <- df_meta
rownames(annotation) <- df_meta$sample

data <- df_tpm %>% dplyr::select(!c(Geneid, Symbol)) %>% as.matrix()
df_exp = ComBat(data, batch = annotation$age) %>% as.data.frame() %>%
  bind_cols(data.frame(Geneid = df_tpm$Geneid,
                       Symbol = df_tpm$Symbol))

# select 1000 most variable genes by cv (coefficient of variation)
N <- 1000
df_cv <- df_exp %>% 
  dplyr::select(!c(Geneid, Symbol)) %>%
  mutate(mean = apply(., 1, mean),
         sd = apply(., 1, sd),
         cv = sd / mean)

rownames(df_cv) <- df_exp$Geneid
df_cv2 <- df_cv %>% arrange(desc(cv))
df_1000genes <- rownames(df_cv2)[1:N] %>% as.data.frame()
colnames(df_1000genes) <- "Geneid"

datExpr_t <- df_exp %>% 
  inner_join(df_1000genes, by = "Geneid") %>% 
  dplyr::select(!c(Geneid, Symbol))

dft <- df_exp %>% inner_join(df_1000genes, by = "Geneid") %>% as.data.frame()
write.table(dft, paste0(dir, "/Data/WGCNA/WGCNA_most_variable_1000genes.tsv"), sep="\t", row.names = F, quote = F)

SubGeneNames=dft[1:1000,"Geneid"]

datExpr <- t(log2(datExpr_t+1))

powers = c(1:50)
sft=pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers,corFnc = cor, corOptions = list(use = 'p'), networkType = "signed")
#sft=pickSoftThreshold(datExpr, powerVector=powers, networkType = "unsigned", RsquaredCut=0.8, verbose=5)
#sft$powerEstimate

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#softPower = 30 #30 for all tissues; 31 is the best but can not specify more than 30.
softPower = 21 #for only amygdala

#calclute the adjacency matrix
adj= adjacency(datExpr, type = "signed", power = softPower)

#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr, networkType = "signed", TOMType = "signed", power = softPower)

colnames(TOM) = rownames(TOM) = SubGeneNames
dissTOM=1-TOM

#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize = 20;

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#set the diagonal of the dissimilarity to NA 
diag(dissTOM) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
# sizeGrWindow(7,7)
# TOMplot(dissTOM^4, geneTree, as.character(dynamicColors))
# 
# module_colors= setdiff(unique(dynamicColors), "grey")
# for (color in module_colors){
#   module=SubGeneNames[which(dynamicColors==color)]
#   write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
# }

module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
sample_amy <- c("AA_amygdala_10_12", "AA_amygdala_10_13", "AA_amygdala_4_20", "AA_amygdala_4_24", 
                "TT_amygdala_10_10", "TT_amygdala_10_14", "TT_amygdala_4_21", 
                "TI_amygdala_10_15", "TI_amygdala_4_19", "TI_amygdala_4_23", "TI_amygdala_10_26",
                "II_amygdala_10_16", "II_amygdala_4_18", "II_amygdala_4_22", "II_amygdala_10_25")
# m2 <- m[sample_amy,]
# heatmap(t(m2),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])

G.expression = datExpr[sample_amy,]
G.expressionColor = numbers2colors(G.expression, signed = F, colors = colorpanel(100, low = "#594255", high = "#dccb18", mid = "#e6eae3"))#gray.colors(100)) #colorpanel(100,low = "#839b5c", high = "#74325c", mid = "#E7E6D5")
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

# dynamicColors2 <- data.frame(dynamicColors, G.expressionColor[1,], G.expressionColor[2,], 
#                              G.expressionColor[3,], G.expressionColor[4,], G.expressionColor[5,], 
#                              G.expressionColor[6,], G.expressionColor[7,], G.expressionColor[8,], 
#                              G.expressionColor[9,], G.expressionColor[10,], G.expressionColor[11,], 
#                              G.expressionColor[12,], G.expressionColor[13,], G.expressionColor[14,], 
#                              G.expressionColor[15,])
# pdf("./result_47samples/WGCNA/WGCNA_tree_all.pdf", w=7, h=5)
# plotDendroAndColors(geneTree, dynamicColors2, groupLabels = c("Dynamic Tree Cut", sample_amy), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
# dev.off()

#### Correlation with behavioral phenotype ####

df_Test <- df_tpm %>% 
  dplyr::select(matches("amygdala"))
rownames(df_Test) <- df_tpm$Geneid

Zt <- read.delim(paste0(dir, "/Data/Behavior_Zscore.tsv"), h=T) %>% 
  mutate(
    ID = case_when(
      ID == "Vmat1_090419-4d" ~ "AA_amygdala_10_12",
      ID == "Vmat1_091919-4c" ~ "AA_amygdala_10_13",
      ID == "Vmat1_032420-4a" ~ "AA_amygdala_4_20",
      ID == "Vmat1_032420-4b" ~ "AA_amygdala_4_24",
      ID == "Vmat1_091019-1b" ~ "TT_amygdala_10_10",
      ID == "Vmat1_091419-1e" ~ "TT_amygdala_10_14",
      ID == "Vmat1_032420-1a" ~ "TT_amygdala_4_17",
      ID == "Vmat1_032420-1b" ~ "TT_amygdala_4_21",
      ID == "Vmat1_091619-3a" ~ "TI_amygdala_10_15",
      ID == "Vmat1_032420-3a" ~ "TI_amygdala_4_19",
      ID == "Vmat1_032420-3b" ~ "TI_amygdala_4_23",
      ID == "Vmat1_083119-3d" ~ "TI_amygdala_10_26",
      ID == "Vmat1_091619-2a" ~ "II_amygdala_10_16",
      ID == "Vmat1_032420-2a" ~ "II_amygdala_4_18",
      ID == "Vmat1_032420-2b" ~ "II_amygdala_4_22",
      ID == "Vmat1_090119-2d" ~ "II_amygdala_10_25")) %>% 
  #  filter(str_detect(ID, "amygdala"))
  na.omit()

df_z <- Zt %>% dplyr::select(Z1N, Z2N)
rownames(df_z) <- Zt$ID

df_Test2 <- t(log2(df_Test+1))
GZ1.expression <- as.numeric(cor(df_z[,1], df_Test2[rownames(df_z),1:13608], method="spearman")[,SubGeneNames])
GZ2.expression <- as.numeric(cor(df_z[,2], df_Test2[rownames(df_z),1:13608], method="spearman")[,SubGeneNames])
GZ.expression <- cbind(GZ1.expression, GZ2.expression)
GZ.expressionColor <- numbers2colors(GZ.expression, signed = T, colors = colorpanel(100,low = "#223a70", high = "#bb5548", mid = "#ede4cd"))

dynamicColors2 <- data.frame(dynamicColors, GZ.expressionColor[,1], GZ.expressionColor[,2], 
                             G.expressionColor[1,], G.expressionColor[2,], G.expressionColor[3,], 
                             G.expressionColor[4,], G.expressionColor[5,], G.expressionColor[6,], 
                             G.expressionColor[7,], G.expressionColor[8,], G.expressionColor[9,], 
                             G.expressionColor[10,], G.expressionColor[11,], G.expressionColor[12,], 
                             G.expressionColor[13,], G.expressionColor[14,], G.expressionColor[15,])
pdf(paste0(dir, "/Data/WGCNA/Figure4c.pdf"), w=7, h=5)
plotDendroAndColors(geneTree, dynamicColors2, groupLabels = c("Dynamic Tree Cut", "Correlation with Activity", "Correlation with Anxiety", sample_amy), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, autoColorHeight = T, main = "Gene dendrogram and module colors")
dev.off()


#install.packages("heatmap3")
library("heatmap3")

#show color bar
#pdf("../Figure_WGCNA_colorlegend.pdf", w=7, h=5)
#colByValue(GZ.expression, col = colorpanel(100,low = "#223a70", high = "#bb5548", mid = "#e6eae3"))
#colByValue(G.expression, col = colorpanel(100, low = "#594255", high = "#dccb18", mid = "#e6eae3"))
#dev.off()

#barplot(matrix(rep(1,15),15),horiz=T,col=G.expressionColor[, dft$Geneid=="ENSMUSG00000022144"])


#### Module genes ####
#### Module turquoise ####
tt<- dft[rownames(dynamicColors2[dynamicColors=="turquoise",]),] %>%
  select(matches("amygdala")) %>% as.matrix()

tt2 <- log2(tt+1) %>%
  cbind(dft[rownames(dynamicColors2[dynamicColors=="turquoise",]),] %>% select(Geneid, Symbol))

nrow(tt)
tt2 %>% inner_join(genelist2 %>% filter(comp=="AA-II"), by="Geneid") %>% nrow() #38 / 189
dft %>% inner_join(genelist2 %>% filter(comp=="AA-II"), by="Geneid") %>% nrow() #53 / 1000
fisher.test(matrix(c(1000, 53, 189, 38), nrow=2)) #turquoise enrichment: p-value = 2.876-08
fisher.test(matrix(c(947, 53, 151, 38), nrow=2)) #turquoise enrichment: p-value = 2.876-08

write.table(tt2, paste0(dir, "/Data/WGCNA/Module_amy_turquoise.tsv"), sep="\t", row.names = F, quote = F)

#### Module purple ####
tt<- dft[rownames(dynamicColors2[dynamicColors=="purple",]),] %>%
  select(matches("amygdala")) %>% as.matrix()

tt2 <- log2(tt+1) %>%
  cbind(dft[rownames(dynamicColors2[dynamicColors=="purple",]),] %>% select(Geneid, Symbol))

nrow(tt2)
tt2 %>% inner_join(genelist2, by="Geneid") %>% nrow() #0 / 35
dft %>% inner_join(genelist2, by="Geneid") %>% nrow() #72 / 1000
fisher.test(matrix(c(1000, 72, 35, 0), nrow=2)) #purple enrichment: p-value = 0.3567

#write.table(tt2, "../Data/Module_amy_purple.tsv", sep="\t", row.names = F, quote = F)

# #### Using all the regions for analysis ####
# #### Module blue ####
# tt<- dft[rownames(dynamicColors2[dynamicColors=="blue",]),] %>%
#   select(matches("amygdala")) %>% as.matrix()
# 
# tt2 <- log2(tt+1) %>%
#   cbind(dft[rownames(dynamicColors2[dynamicColors=="blue",]),] %>% select(Geneid, Symbol))
# dft
# tt2 %>% inner_join(genelist2 %>% filter(comp=="AA-II"), by="Geneid") %>% nrow()
# dft %>% inner_join(genelist2 %>% filter(comp=="AA-II"), by="Geneid") %>% nrow()
# dft %>% inner_join(genelist2, by="Geneid") %>% nrow()
# #nrow(genelist2 %>% filter(comp=="AA-II"))
# nrow(genelist2) #all DEGs
# nrow(tt2) #genes belonging to blue module
# write.table(tt2, "../Data/Module_blue.tsv", sep="\t", row.names = F, quote = F)
# fisher.test(matrix(c(1000, 174, 44, 19), nrow=2)) #blue
# 
# #### Module greenyellow ####
# tt<- dft[rownames(dynamicColors2[dynamicColors=="greenyellow",]),] %>%
#   select(matches("amygdala")) %>% as.matrix()
# 
# tt2 <- log2(tt+1) %>%
#   cbind(dft[rownames(dynamicColors2[dynamicColors=="greenyellow",]),] %>% select(Geneid, Symbol))
# fisher.test(matrix(c(1000, 24, 44, 10), nrow=2)) #greenyellow
# 
# write.table(tt2, "../Data/Module_greenyellow.tsv", sep="\t", row.names = F, quote = F)

#### Eigengenes ####
datME = moduleEigengenes(datExpr,dynamicColors)$eigengenes
datME
#as.numeric(cor(df_z[,1], datME[rownames(df_z),1:12], method="spearman"))
#as.numeric(cor(df_z[,2], datME[rownames(df_z),1:12], method="spearman"))
as.numeric(cor(df_z[,1], datME[rownames(df_z),1:12], method="spearman"))
as.numeric(cor(df_z[,2], datME[rownames(df_z),1:12], method="spearman"))

#### Create networks using the genes above by STRING ####
#https://string-db.org

#### Network ####
#install.packages("GGally")
library(GGally)
#install.packages("sna")
library(sna)
#install.packages("network")
library(network)
library(ggrepel)

# tt<- dft[rownames(dynamicColors2[dynamicColors=="greenyellow" | dynamicColors=="blue",]),] %>%
#   select(matches("amygdala")) %>% as.matrix()
# 
# tt2 <- log2(tt+1) %>%
#   cbind(dft[rownames(dynamicColors2[dynamicColors=="greenyellow"| dynamicColors=="blue",]),] %>% select(Geneid, Symbol))

tt<- dft[rownames(dynamicColors2[dynamicColors=="turquoise",]),] %>%
  select(matches("amygdala")) %>% as.matrix()

tt2 <- log2(tt+1) %>%
  cbind(dft[rownames(dynamicColors2[dynamicColors=="turquoise",]),] %>% select(Geneid, Symbol))


tt3 <- tt2 %>% left_join(genelist2, by="Geneid") %>%
  mutate(lab.col = if_else(is.na(number), "#383c3c", "#a22041"),
         Symbol = if_else(Geneid == "ENSMUSG00000118401", "Gpr52", Symbol),
         Symbol = if_else(Symbol == "Drd1", "Drd1a", Symbol),
         #         Symbol = if_else(Symbol == "Stk26", "Mst4", Symbol),
         Symbol = if_else(Symbol == "Ccn2", "Ctgf", Symbol),
         lab.size = 3)

tt4 <- tt3 %>%  
  left_join(read.delim(paste0(dir, "/Data/Cytoscape/CytoScape_result_format.tsv"), head=T), by="Symbol") %>%
  mutate(col = case_when(
    is.na(Cluster) ~ "#383c3c",
    Cluster == 1 ~ "#9079ad", 
    Cluster == 2 ~ "#dccb18",
    Symbol == "Mme" ~ "#b94047",
    Symbol == "Nt5e" ~ "#b94047",
    Symbol == "Gdnf" ~ "#b94047",
    Symbol == "Foxp1" ~ "#b94047",
    Symbol == "Foxp2" ~ "#b94047",
    Symbol == "Grb7" ~ "#b94047",
    Symbol == "Six3" ~ "#b94047")
  ) %>%
  distinct(Symbol, .keep_all = TRUE)

S#df_net <- read.delim("../Data/string_interactions.tsv", head=T, stringsAsFactors=FALSE)
df_net <- read.delim(paste0(dir, "/Data/STRING/string_interactions_amy.tsv"), head=T, stringsAsFactors=FALSE)
df.net <- network(df_net[,1:2], directed = TRUE)

rownames(tt4) <- tt4$Symbol
df.net %v% "col" <- tt4[network.vertex.names(df.net),]$col
df.net %v% "lab.col" <- tt4[network.vertex.names(df.net),]$lab.col
df.net %v% "lab.size" <- tt4[network.vertex.names(df.net),]$lab.size
network::set.edge.attribute(df.net, "weights", df_net$combined_score*1.2)

# set.seed(15)
# set.seed(4)
set.seed(20)

g1 <- ggnet2(df.net, color = "col", label.color = "lab.col", label.size = "lab.size",
             edge.color = "#17184b", edge.size = "weights", edge.alpha = .4, node.size = 4,
             node.alpha = 1, node.label = FALSE)  + 
  geom_text_repel(size = tt4[network.vertex.names(df.net),]$lab.size,
                  color = tt4[network.vertex.names(df.net),]$lab.col,
                  label = network.vertex.names(df.net),
                  #                  fontface = data_amy[network.vertex.names(amy.net),]$font, 
                  #                  box.padding = .05, segment.size = .1,
                  nudge_x       = -.05,
                  nudge_y       = -.0001,
                  max.iter      = 10000, 
                  segment.size  = .2,
                  segment.color = "grey50",
                  segment.alpha = .8,
                  box.padding   = .05,
                  direction     = "y",
                  hjust         = .5)
g1

ggsave(paste0(dir, "/Figures/Figure4d.pdf"), g1, w=4, h=4)


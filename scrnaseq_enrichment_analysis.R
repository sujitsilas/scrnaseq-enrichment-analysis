
###### Mouse Integrated Dataset ########
# Marco Olfr2+ 
load("MPC_Reclustered.Robj")

## remove cDC1 cDC2 and Mature DCs
x <- MPC_Reclustered@meta.data
select_cells <- rownames(x)[!(x$Annotation_Final %in% c("cDC1", "cDC2", "Mature_DCs"))]
MPC_Reclustered <- subset(MPC_Reclustered,  cells = select_cells)
DimPlot(MPC_Reclustered)
dev.off()

# Volcano Plot
library(EnhancedVolcano)

to_show <- c("Olr1", "Vegf3",  "Il1b", "Il6", "Il12a", "Il23a", "Ccl22", "Ccl5", "Cx3cl1", "Cxcl9" ,"Vegf3", "Il1b", "Il23a", "Ccl22", "Ccl5", "Cx3cl1", "Cxcl9") %>% unique()
to_remove <- c("Tox4", "Unc119", "Tip2", "Stra6l", "Tamalin", "DCn", "Brat1", "F10", "Cdc45", "Selp")

df_olfr2pos_vs_olfr2neg <- read.csv("olfr2_Exp_Positive_vs_Negative.csv")

t <- df_olfr2pos_vs_olfr2neg %>% mutate(Type="Olfr2+ vs Olfr2-", Directionality=ifelse(log2FoldChange>0 & padj<0.05, yes="Up in Olfr2+", no=NA),
                                   Directionality=ifelse(log2FoldChange<0 & padj<0.05, yes="Up in Olfr2-", no=Directionality))

t <- t[grep(pattern =paste(to_remove,collapse =  "|"), t$gene_symbol, invert = T),]


keyvals = ifelse((t$log2FoldChange< 0 & t$padj <0.05), 'blue',
                 ifelse((t$log2FoldChange >0 & t$padj <0.05), 'red', NA))


names(keyvals)[keyvals == 'black'] <- "NS"
names(keyvals)[keyvals == 'red'] <- unique(t$Directionality)[2]
names(keyvals)[keyvals == 'blue'] <- unique(t$Directionality)[3]

pointSize <- 14
lineWidth <- 1 / 2.835

na_values_x <- function(log2FC, limit){
  a <- as.vector(log2FC)
  b <- ifelse(a < -limit, yes=-limit, no=a)
  c <- ifelse(b > limit, yes=limit, no=b)
  return(c)
  
}

na_values_y <- function(padj, limit){
  a <- as.vector(-log10(padj))
  b <- ifelse(a > limit, yes=limit, no=a)
  return(b)
}


t$label <- ifelse(t$log2FoldChange>10 & t$padj<0.05| t$log2FoldChange< -10 & t$padj<0.05, yes ="lable",no="no")
t$gene_symbol[t$label=="lable"]

e = EnhancedVolcano(t, lab =t$gene_symbol, x = "log2FoldChange", y = "padj",
                    legendPosition ="top",
                    subtitle = paste0('Adj. P-Value Cutoff = 0.05'),
                    selectLab = ,
                    xlab = bquote("average" ~log[2]~ "fold change"),
                    ylab = bquote('-log10 (adj. P-value)'),
                    labSize =5,
                    xlim = c(-15,15),
                    ylim=c(0,16),
                    boxedLabels = FALSE, 
                    labFace = 'bold', 
                    pCutoff = 0.05,
                    FCcutoff =0,
                    axisLabSize = 25,
                    pointSize = 2,
                    drawConnectors = F,
                    widthConnectors =0.1,
                    arrowheads = F,
                    colCustom = keyvals,
                    colAlpha = 0.8,
                    title= unique(t$Type), gridlines.major = F)+ scale_y_continuous(breaks = scales::pretty_breaks(n = 5),limits =c(0,16), na.value =) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5),limits =c(-14,14), na.value =) + theme_light() +
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(size = pointSize, colour = "black",face = "bold"),
    axis.title  = element_text(size =20, colour = "black"),
    axis.text.x  = element_text(size = 18, colour = "black",vjust = 0.5, hjust=1),
    axis.text.y  = element_text(size = 18, colour = "black", hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize * 0.8, colour = "black", face = "bold"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+ guides(colour = guide_legend(override.aes = list(size=4)))


png(filename = paste0("Olfr2+_vs_Olfr2-.png"), height=10, width=11, units="in", res=300)
plot(e)
dev.off()

## Pseudobulk
DefaultAssay(MPC_Reclustered) <- "RNA"
df_pseudo <- FetchData(MPC_Reclustered, vars = c("Annotation_Final", rownames(MPC_Reclustered)), slot = "counts")
df_pseudo %>% group_by(Annotation_Final) %>% summarize_all(sum) %>% as.data.frame() -> df_pseudo
rownames(df_pseudo) <- df_pseudo$Annotation_Final
df_pseudo$Annotation_Final <- NULL
df_pseudo <- df_pseudo %>% t() %>% as.data.frame()


add_noise_to_data <- function(x){
  for(i in names(x)){
    print(i)
    x[paste0(i,"_S1")] <- x[i] + 2
    x[paste0(i,"_S2")] <- x[i] + 1
    x[paste0(i,"_S3")] <- x[i] + 3
  }
  return(x[,order(names(x))])
}

names(df_pseudo)
df_pseudo <- add_noise_to_data(df_pseudo)
clipr::write_clip(unique(MPC_Reclustered$Annotation_Final), breaks = "\t")

## vst transform
vst_MPC <- DESeq2::varianceStabilizingTransformation(as.matrix(df_pseudo))
write.csv(vst_MPC, "Zernecke_PseudoBulk.csv")

## ssGSEA
library(escape)
library(dittoSeq)
library(Seurat)

up500 <- df_olfr2pos_vs_olfr2neg %>% filter(padj<0.05) %>% arrange(-log2FoldChange) %>% top_n(n = 500, wt = log2FoldChange) %>% select(gene_symbol)
down500 <- df_olfr2pos_vs_olfr2neg %>% filter(padj<0.05) %>% arrange(log2FoldChange) %>% top_n(n = 500, wt = -log2FoldChange) %>% select(gene_symbol)

gene_sets <- list(
  `Olfr2_Pos_Macs` = up500$gene_symbol,
  `Olfr2_Neg_Macs` = down500$gene_symbol
)

saveRDS(gene_sets, "Olfr2_Gene_Sets.robj")

length(colnames(MPC_Reclustered))
gene_sets <- readRDS("Olfr2_Gene_Sets.robj")
ES.seurat <- enrichIt(obj = MPC_Reclustered, 
                      gene.sets = gene_sets, cores = 8,
                      groups = 2000, 
                      min.size = 5)
MPC_Reclustered<- AddMetaData(MPC_Reclustered, metadata = ES.seurat)

names(MPC_Reclustered@meta.data) <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Protocol","Species","Williams_C57_ID","S.Score","G2M.Score","Phase",
                                      "IEG_Score1", "Annotation_Simplified","Annotation_Final", "Olfr2+ Macs" ,"Olfr2- Macs")

x <- MPC_Reclustered@meta.data
png("Olfr2+ and Olfr2- scEnrichment.png", width=17, height = 6,units = "in",bg = "white",res = 400)
multi_dittoPlot(MPC_Reclustered, vars = c("Olfr2+ Macs", 
                              "Olfr2- Macs"),
                group.by = "Annotation_Final", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", min=0, max=1500,
                theme = theme_classic() + theme(text = element_text(size = pointSize, colour = "black"),
                                                rect = element_blank(),
                                                line = element_line(size = lineWidth, colour = "black"),
                                                plot.title  = element_text(size = pointSize, colour = "black",face = "bold"),
                                                axis.title  = element_text(size =20, colour = "black"),
                                                axis.text.x  = element_text(size = 16, colour = "black",vjust = 0.5, hjust=1),
                                                axis.text.y  = element_text(size = 18, colour = "black", hjust = 1)))
dev.off()



# GSEA Enrichment Dot Plot
files <- grep(list.files(path="Zerneke_GSEA",recursive = T, full.names = T), pattern="gsea_report_for", value = T)

es_reports <- list()
for(i in files[files %like% ".tsv"]){
  df<- read_tsv(i) %>% as.tibble()
  name <- str_split(i, "/", simplify = T)[,8]
  df$Comparison <- str_split(name, "\\.", simplify = T)[,1]
  if(nrow(df)<1){
    next
  }
  es_reports[[i]] <- df %>% as.data.frame()
}


write.csv(df_es_reports, "GSEA_NES_Score_Reports.csv")
png("ES_Socre_Plot.png", width =5, height = 8, units = "in",bg = "white", res = 400)
ggplot(df_es_reports, aes(x=NAME, y=Comparison)) + geom_point(aes(color=NES), size=10) +scale_color_gradient2() + xlab("")+ylab("")+theme(text = element_text(size = pointSize, colour = "black"),
                                                                                                                                          line = element_line(size = lineWidth, colour = "black"),
                                                                                                                                          plot.title  = element_text(size = pointSize, colour = "black",face = "bold"),
                                                                                                                                          axis.title  = element_text(size =10, colour = "black"),
                                                                                                                                          axis.text.x  = element_text(size = 13, colour = "black",vjust = 0.5, hjust=1, angle=90),
                                                                                                                                          axis.text.y  = element_text(size = 13, colour = "black", hjust = 1))
dev.off()


# Feature plots
Idents(MPC_Reclustered)
png("Olfr2+ Enrichment FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(MPC_Reclustered, features ="Olfr2+ Macs", min.cutoff = 800, max.cutoff =  max(MPC_Reclustered$`Olfr2+ Macs`), pt.size =1, order = T, label = T, label.size =8,repel = T) +scale_color_gradientn(colours = c(alpha("grey80",0.5), alpha("white",0.1),alpha("red",1),"red","red"),
                                                                                                                                                                        limits=c(800,max(MPC_Reclustered$`Olfr2+ Macs`)),
                                                                                                                                                                        labels=c("","  "," maximum  \n enrichment"),
                                                                                                                                                                        breaks = c(800,1000, max(MPC_Reclustered$`Olfr2+ Macs`))) +
  theme(text = element_text(size =12, colour = "black"),
         line = element_line(size = lineWidth, colour = "black"),
         plot.title  = element_text(size = 24, colour = "black",face = "bold"),
         axis.title  = element_text(size =18, colour = "black", face="bold"),
         axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
         axis.text.y  = element_text(size = 20, colour = "black", hjust = 1))

dev.off()
getwd()
setwd("/Users/sujitsilas/Desktop/Ley Lab/Marco_GSEA_Olfr2/")
png("Olfr2- Enrichment FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(MPC_Reclustered, features ="Olfr2- Macs", min.cutoff = 800, max.cutoff =  max(MPC_Reclustered$`Olfr2- Macs`), pt.size =1, order = T, label = T, label.size = 8, repel = T) +scale_color_gradientn(colours = c(alpha("grey85",0.5), alpha("white",0.5),alpha("blue",0.8),"blue","blue"),
                                                                                                                                                                                                     limits=c(800,max(MPC_Reclustered$`Olfr2- Macs`)),
                                                                                                                                                                                                     labels=c("","  "," maximum  \n enrichment"),
                                                                                                                                                                                                     breaks = c(800,1000, max(MPC_Reclustered$`Olfr2- Macs`))) +
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1))

dev.off()
?ggplot2


png("Ccr2_FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(MPC_Reclustered, features ="Ccr2", pt.size =1, min.cutoff = 0,order = T, label = T, label.size = 8, repel = T) +
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1), legend.text = element_text(size=15))
dev.off()




png("Trem2_FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(MPC_Reclustered, features ="Trem2", pt.size =1, min.cutoff = 0,order = T, label = T, label.size = 8, repel = T) +
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1), legend.text = element_text(size=15))
dev.off()



png("Itgax_FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(MPC_Reclustered, features ="Itgax", pt.size =1, min.cutoff = 0,order = T, label = T, label.size = 8, repel = T) +
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1), legend.text = element_text(size=15))
dev.off()



png("Lyve1_FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(MPC_Reclustered, features ="Lyve1", pt.size =1, min.cutoff = 0,order = T, label = T, label.size = 8, repel = T) +
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1), legend.text = element_text(size=15))
dev.off()


###### Human Dataset ########
load("/Users/sujitsilas/Downloads/Figure3_Human_MPC.Robj")
DimPlot(Figure3_Human_MPC)
dev.off()

Figure3_Human_MPC

## Pseudobulk
DefaultAssay(Figure3_Human_MPC) <- "RNA"
x <- Figure3_Human_MPC@meta.data
df_pseudo <- FetchData(Figure3_Human_MPC, vars = c("Annotation_hCCA", rownames(Figure3_Human_MPC)), slot = "counts")
df_pseudo %>% group_by(Annotation_hCCA) %>% summarize_all(sum) %>% as.data.frame() -> df_pseudo
rownames(df_pseudo) <- df_pseudo$Annotation_hCCA
df_pseudo$Annotation_hCCA <- NULL
df_pseudo <- df_pseudo %>% t() %>% as.data.frame()

?ggplot2
add_noise_to_data <- function(x){
  for(i in names(x)){
    print(i)
    x[paste0(i,"_S1")] <- x[i] + 2
    x[paste0(i,"_S2")] <- x[i] + 1
    x[paste0(i,"_S3")] <- x[i] + 3
  }
  return(x[,order(names(x))])
}

names(df_pseudo)
df_pseudo <- add_noise_to_data(df_pseudo)
clipr::write_clip(unique(MPC_Reclustered$Annotation_Final), breaks = "\t")
df_pseudo <- round(df_pseudo, digits = 0)

## vst transform
vst_hcca <- DESeq2::varianceStabilizingTransformation(as.matrix(df_pseudo))
write.csv(vst_hcca, "Zernecke_Human_PseudoBulk.csv")

## ssGSEA
library(dplyr)
library(stringr)

up <- df_olfr2pos_vs_olfr2neg %>% filter(padj<0.05 & log2FoldChange>0)%>% arrange(-log2FoldChange)
down <- df_olfr2pos_vs_olfr2neg %>% filter(padj<0.05 & log2FoldChange<0) %>% arrange(log2FoldChange)

gene_sets_all <- list(
  `Olfr2_Pos_Macs` = str_to_upper(up500$gene_symbol),
  `Olfr2_Neg_Macs` = str_to_upper(down500$gene_symbol)
)


saveRDS(gene_sets_all, "Olfr2_Gene_Sets_AllSig.robj")


length(colnames(Figure3_Human_MPC))

ES.seurat <- enrichIt(obj = Figure3_Human_MPC, 
                      gene.sets = gene_sets_all, cores = 8,
                      groups = 500, 
                      min.size = 5)

Figure3_Human_MPC@meta.data$`OLFR2+ Macs` <- NULL
Figure3_Human_MPC@meta.data$`OLFR2- Macs` <- NULL
Figure3_Human_MPC@meta.data$`NA` <- NULL
Figure3_Human_MPC@meta.data$`NA` <- NULL

Figure3_Human_MPC<- AddMetaData(Figure3_Human_MPC, metadata = ES.seurat)

names(Figure3_Human_MPC@meta.data) <- c("nCount_RNA","nFeature_RNA","percent.mt" ,"Patient","Protocol","Species","Vessel","Annotation_hCCA","OLFR2+ Macs","OLFR2- Macs" )

x <- Figure3_Human_MPC@meta.data

# Enrichmet Plot
png("OLFR2+ and OLFR2- scEnrichment.png", width=17, height = 6,units = "in",bg = "white",res = 400)
multi_dittoPlot(Figure3_Human_MPC, vars = c("OLFR2+ Macs", 
                                          "OLFR2- Macs"),
                group.by = "Annotation_hCCA", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", min=-300, max=2000,
                theme = theme_classic() + theme(text = element_text(size = pointSize, colour = "black"),
                                                rect = element_blank(),
                                                line = element_line(size = lineWidth, colour = "black"),
                                                plot.title  = element_text(size = pointSize, colour = "black",face = "bold"),
                                                axis.title  = element_text(size =20, colour = "black"),
                                                axis.text.x  = element_text(size = 16, colour = "black",vjust = 0.5, hjust=1),
                                                axis.text.y  = element_text(size = 18, colour = "black", hjust = 1)))
dev.off()

# Feature plots
Idents(Figure3_Human_MPC)
png("Olfr2+ Enrichment FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(Figure3_Human_MPC, features ="OLFR2+ Macs", min.cutoff = 800, max.cutoff =  max(Figure3_Human_MPC$`OLFR2+ Macs`), pt.size =1, order = T, label = T, label.size = 6, ) +scale_color_gradientn(colours = c(alpha("grey80",0.5), alpha("white",0.1),alpha("red",1),"red","red"),
                                                                                                                                                                                                     limits=c(800,max(Figure3_Human_MPC$`OLFR2+ Macs`)),
                                                                                                                                                                                                     labels=c("","  "," maximum  \n enrichment"),
                                                                                                                                                                                                     breaks = c(800,1000, max(Figure3_Human_MPC$`OLFR2+ Macs`))) +
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1))

dev.off()
?harmony

png("Olfr2- Enrichment FeaturePlot.png", width = 11,height = 8, units = "in",res = 400,bg = "white")
FeaturePlot(MPC_Reclustered, features ="Olfr2- Macs", min.cutoff = 800, max.cutoff =  max(MPC_Reclustered$`Olfr2- Macs`), pt.size =1, order = T, label = T, label.size = 6, ) +scale_color_gradientn(colours = c(alpha("grey85",0.5), alpha("white",0.5),alpha("red",0.8),"red","red"),
                                                                                                                                                                                                     limits=c(800,max(MPC_Reclustered$`Olfr2- Macs`)),
                                                                                                                                                                                                     labels=c("","  "," maximum  \n enrichment"),
                                                                                                                                                                                                     breaks = c(800,1000, max(MPC_Reclustered$`Olfr2- Macs`))) +
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1))

dev.off()

## SS GSEA
suppressMessages({
  library(splatter)
  library(Seurat)
  library(msigdbr)
  library(singleseqgset)
  library(heatmap3)
})



gene_sets_new <- gdata::cbindX(as.data.frame(gene_sets[[1]][gene_sets[[1]] %in%  rownames(MPC_Reclustered@assays$RNA)]), as.data.frame(gene_sets[[2]][gene_sets[[2]] %in%  rownames(MPC_Reclustered@assays$RNA)]))
names(gene_sets_new) <- c("Olfr2+ Macs", "Olfr2- Macs")
write.csv(gene_sets_new, "Gene_Sig_New.csv")


logfc.data <- logFC(cluster.ids=MPC_Reclustered$Annotation_Final,expr.mat=MPC_Reclustered@assays$RNA@data)
names(logfc.data)
logfc.data$cluster.cells
logfc.data$log.fc.cluster
gse.res <- wmw_gsea(expr.mat=MPC_Reclustered@assays$RNA@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=gene_sets)
names(gse.res)

res.stats <- gse.res[["GSEA_statistics"]]
res.pvals <- gse.res[["GSEA_p_values"]]

res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
res.stats[order(res.stats[,1],decreasing=TRUE),] #Top gene sets enriched by z scores


write.csv(as.data.frame(cbind(t(res.stats),t(res.pvals))),"GSEA_Enrichment.csv")
df <- read.csv("GSEA_Enrichment.csv")
df$gene_set <- ifelse(df$gene_set=="Olfr2_Neg_Macs", yes="Olfr2- Macs", no = "Olfr2+ Macs")


png("ES_Socre_Plot.png", width =4, height = 6, units = "in",bg = "white", res = 400)
ggplot(df, aes(x=gene_set, y=Cluster)) + geom_point(aes(color=ES,size=adj.P.values)) +scale_colour_gradientn(colours = c("blue","red"), na.value = NA)+ xlab("")+ylab("")+theme(text = element_text(size = pointSize, colour = "black"),
                                                                                                                                          line = element_line(size = lineWidth, colour = "black"),
                                                                                                                                          plot.title  = element_text(size = pointSize, colour = "black",face = "bold"),
                                                                                                                                          axis.title  = element_text(size =10, colour = "black"),
                                                                                                                                          axis.text.x  = element_text(size = 13, colour = "black",vjust = 0.5, hjust=1, angle=90),
                                                                                                                                          axis.text.y  = element_text(size = 13, colour = "black", hjust = 1)) +
  scale_size(name="adj.P.values", breaks = c(0.001,0.01,0.05),
             range=c(5,1), limits = c(2.03e-15, 0.05), labels = c("p \u2264 0.001", "p \u2264 0.01", "p \u2264 0.05"))
dev.off()
length(unique(seurat_1$nFeature_RNA))
length((seurat_1$nFeature_RNA))

###### Integrated Analysis #####
load("mmu_hsa_Integrated_MPC.Robj")
row.names(mmu_hsa_Integrated_MPC@assays$RNA@counts)

DimPlot(mmu_hsa_Integrated_MPC, group.by = "Species")
dev.off()
library(escape)
load("Olfr2_Gene_Sets_AllSig.robj")
ES.seurat <- escape::enrichIt(obj = mmu_hsa_Integrated_MPC, 
                      gene.sets = gene_sets_all, cores = 8,
                      groups = 2000, 
                      min.size = 5)

mmu_hsa_Integrated_MPC<- AddMetaData(mmu_hsa_Integrated_MPC, metadata = ES.seurat)

names(mmu_hsa_Integrated_MPC@meta.data) <- c("nCount_RNA","nFeature_RNA" ,"percent.mt" , "Sample" ,"Patient","Protocol" , "Species", "Vessel" ,"Annotation_hCCA" ,"nCount_integrated" , "nFeature_integrated","Integrated_Annotation", "OLFR2+ Macs" ,"OLFR2- Macs" )

mmu_hsa_Integrated_MPC <- subset(mmu_hsa_Integrated_MPC, cell=rownames(mmu_hsa_Integrated_MPC@meta.data)[!(mmu_hsa_Integrated_MPC$Integrated_Annotation %like% "T/B/NK|Non_Leuko|MoDC/cDC2|FSN1/CCR7_DC")])

x <- mmu_hsa_Integrated_MPC@meta.data
png("OLFR2+ and OLFR2- scEnrichment Integrated.png", width=17, height = 6,units = "in",bg = "white",res = 400)
multi_dittoPlot(mmu_hsa_Integrated_MPC, vars = c("OLFR2+ Macs", 
                                          "OLFR2- Macs"),
                group.by = "Integrated_Annotation", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", min=-550, max=1500,
                theme = theme_classic() + theme(text = element_text(size = pointSize, colour = "black"),
                                                rect = element_blank(),
                                                line = element_line(size = lineWidth, colour = "black"),
                                                plot.title  = element_text(size = pointSize, colour = "black",face = "bold"),
                                                axis.title  = element_text(size =20, colour = "black"),
                                                axis.text.x  = element_text(size = 16, colour = "black",vjust = 0.5, hjust=1),
                                                axis.text.y  = element_text(size = 18, colour = "black", hjust = 1)))
dev.off()


png("OLFR2+ Enrichment FeaturePlot Integrated.png", width = 16,height = 7, units = "in",res = 400,bg = "white")
FeaturePlot(mmu_hsa_Integrated_MPC, features ="OLFR2+ Macs", min.cutoff =-100, max.cutoff =  max(mmu_hsa_Integrated_MPC$`OLFR2+ Macs`), pt.size =1, order = T, label = T, label.size = 6, split.by = "Species")&scale_color_gradientn(colours = c(alpha("grey80",0.5), alpha("white",0.1),alpha("red",1),"red","red"),
                                                                                                                                                                                                         limits=c(-100,max(mmu_hsa_Integrated_MPC$`OLFR2+ Macs`)),
                                                                                                                                                                                                         labels=c("","  "," maximum  \n enrichment"),
                                                                                                                                                                                                         breaks = c(-100,1000, max(mmu_hsa_Integrated_MPC$`OLFR2+ Macs`)))&
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1))

dev.off()



png("OLFR2- Enrichment FeaturePlot Integrated.png", width = 16,height = 7, units = "in",res = 400,bg = "white")
FeaturePlot(mmu_hsa_Integrated_MPC, features ="OLFR2- Macs", min.cutoff = -100, max.cutoff =  max(mmu_hsa_Integrated_MPC$`OLFR2+ Macs`), pt.size =1, order = T, label = T, label.size = 6, split.by = "Species")&scale_color_gradientn(colours = c(alpha("grey80",0.5), alpha("white",0.1),alpha("red",1),"red","red"),
                                                                                                                                                                                                                                      limits=c(-100,max(mmu_hsa_Integrated_MPC$`OLFR2- Macs`)),
                                                                                                                                                                                                                                      labels=c("","  "," maximum  \n enrichment"),
                                                                                                                                                                                                                                      breaks = c(-100,1000, max(mmu_hsa_Integrated_MPC$`OLFR2- Macs`)))&
  theme(text = element_text(size =12, colour = "black"),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(size = 24, colour = "black",face = "bold"),
        axis.title  = element_text(size =18, colour = "black", face="bold"),
        axis.text.x  = element_text(size = 20, colour = "black",vjust = 0.5,),
        axis.text.y  = element_text(size = 20, colour = "black", hjust = 1))

dev.off()


unique(mmu_hsa_Integrated_MPC$Integrated_Annotation)
unique(mmu_hsa_Integrated_MPC$Annotation_hCCA)

logfc.data <- logFC(cluster.ids=mmu_hsa_Integrated_MPC$Integrated_Annotation,expr.mat=mmu_hsa_Integrated_MPC@assays$integrated@data)
names(logfc.data)
logfc.data$cluster.cells
logfc.data$log.fc.cluster
gse.res <- wmw_gsea(expr.mat=mmu_hsa_Integrated_MPC@assays$integrated@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=gene_sets_all)
names(gse.res)


res.stats <- gse.res[["GSEA_statistics"]]
res.pvals <- gse.res[["GSEA_p_values"]]

res.pvals <- apply(res.pvals,2,p.adjust,method="fdr") #Correct for multiple comparisons
res.stats[order(res.stats[,1],decreasing=TRUE),] #Top gene sets enriched by z scores


write.csv(as.data.frame(cbind(t(res.stats),t(res.pvals))),"GSEA_Enrichment_Human_Integrated.csv")
df <- read.csv("GSEA_Enrichment_Human_Integrated.csv")
df$gene_set <- ifelse(df$gene_set=="Olfr2_Neg_Macs", yes="OLFR2- Macs", no = "OLFR2+ Macs")
df <- filter(df, !(df$ES==0))
df <- df[ !(df$Cluster %in% c("XCR1_DC", "FSCN1/CCR7_DC")),]
setwd("/Users/sujitsilas/Desktop/")
png("ES_Socre_Plot_Human_Integrated.png", width =5, height = 6, units = "in",bg = "white", res = 400)
ggplot(df, aes(x=gene_set, y=Cluster)) + geom_point(aes(color=ES,size=adj.P.values)) +scale_colour_gradientn(colours = c("darkgrey","red"), na.value = NA)+ xlab("")+ylab("")+theme(text = element_text(size = pointSize, colour = "black"),
                                                                                                                                                                                line = element_line(size = lineWidth, colour = "black"),
                                                                                                                                                                                plot.title  = element_text(size = pointSize, colour = "black",face = "bold"),
                                                                                                                                                                                axis.title  = element_text(size =10, colour = "black"),
                                                                                                                                                                                axis.text.x  = element_text(size = 13, colour = "black",vjust = 0.5, hjust=1, angle=90),
                                                                                                                                                                                axis.text.y  = element_text(size = 13, colour = "black", hjust = 1)) +
  scale_size(name="adj.P.values", breaks = c(0.001,0.01,0.05),
             range=c(5,1), limits = c(2.03e-15, 0.05), labels = c("p \u2264 0.001", "p \u2264 0.01", "p \u2264 0.05"))
dev.off()


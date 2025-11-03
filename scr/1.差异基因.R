# /home/yuchen/miniconda3/envs/R4.0/bin/R
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(reshape2)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(ggsci)
library(biomaRt)
##############################################################################################################################################################################################################################################################################
data_integretion = function(filepath){
  files = list.files(filepath, pattern = ".count$")
  data_list = list()
  file_num = length(files)
  for(i in 1:file_num){
    data_list[[i]] = read.table(paste(filepath,files[i],sep = "/"), sep = "\t", col.names = c("gene",files[i]))
  }
  RawCount = reduce(data_list, inner_join, by = "gene")
  return(RawCount)
}

## pheatmap使用的是grid图形体系，不是ggplot2体系
save_pheatmap_pdf = function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

##############################################################################################################################################################################################################################################################################
## 读取数据
RawCount = data_integretion("/Volumes/GaoxyData/data/RNA-seq/counts")
colnames(RawCount)[2:ncol(RawCount)] = c("CON_1", "CON_2", "CON_3",
                                         "LPS_1","LPS_2","LPS_3",
                                         "LPS+NC_1","LPS+NC_2","LPS+NC_3",
                                         "LPS+SiR_1", "LPS+SiR_2", "LPS+SiR_3")

## gene ID转换
# 选择数据库与数据集（人类）
mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# 批量转换
mapping = getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = RawCount$gene,
  mart = mart
)
colnames(mapping) = c('gene', 'symbol')
RawCount = merge(RawCount, mapping, by = "gene")

RawCount = RawCount[,c("symbol", "CON_1", "CON_2", "CON_3",
                                         "LPS_1","LPS_2","LPS_3",
                                         "LPS+NC_1","LPS+NC_2","LPS+NC_3",
                                         "LPS+SiR_1", "LPS+SiR_2", "LPS+SiR_3")]

group = factor(c(rep("CON",3), rep("LPS",3), rep("LPS+NC",3), rep('LPS+SiR',3)))
RawCountFilter = RawCount[rowSums(RawCount[,2:ncol(RawCount)]) > 180,]
rangelist.filted = DESeqDataSetFromMatrix(RawCountFilter[,2:ncol(RawCountFilter)], DataFrame(group), design= ~ group )
rangelist.res = DESeq(rangelist.filted)
rld = rlog(rangelist.res, blind=FALSE)

## 绘制PCA
pca = plotPCA(rld, intgroup=c("group"), returnData=T)
head(pca)
percentVar = round(100 * attr(pca, "percentVar"))

## plot
pdf("/Volumes/GaoxyData/data/RNA-seq/All_pca.pdf",width = 10, height = 8)
ggplot(pca,aes(PC1,PC2))+
  geom_point(aes(colour = group), size = 5)+
  # scale_shape_manual(values = c(15, 4,  16, 17, 18))+
  scale_color_manual(values = c("#FFC636", "#00ADA9", "#FF6444", "#CEF09D")) +
  xlab(paste0("PCA1: ", percentVar[1], "% variance"))+
  ylab(paste0("PCA2: ", percentVar[2], "% variance"))+
  theme_bw()+
  theme(legend.background = element_blank(),
        legend.text=element_text(size=16),
        legend.position="right",
        legend.title=element_blank(),
        legend.key.size=unit(0.3,'cm'),
        panel.grid = element_blank(),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_text(face = "bold", colour = "black", size = 16),
        axis.title.y = element_text(face = "bold", colour = "black", size = 16)
  )
dev.off()


## 差异基因
## statistic
DEgeneRes_a_b_results = function(rangelist.res,Aname,Bname){
  DEgeneRes_a_b_tmp = results(rangelist.res, contrast=c("group",Aname,Bname))
  DEgeneRes_a_b_tmp$symbol = RawCountFilter[row.names(DEgeneRes_a_b_tmp),]$symbol
  return(DEgeneRes_a_b_tmp)
}

DEgeneRes_LPS_CON = DEgeneRes_a_b_results(rangelist.res, "LPS", "CON")
DEgeneRes_LPSSiR_LPS = DEgeneRes_a_b_results(rangelist.res, "LPS+SiR", "LPS")

## read counts normalized
Normalized_counts = counts(rangelist.res, normalized = T)
Normalized_counts = data.frame(Normalized_counts)
Normalized_counts$symbol = RawCountFilter[row.names(Normalized_counts),]$symbol

## Functions
pairwise_Compare = function(DEgenesRes,FC){
  df = as.data.frame(DEgenesRes)
  df$express = "NO"
  df[which(df$padj < 0.05 & df$log2FoldChange<=-FC),]$express = "DOWN"
  df[which(df$padj < 0.05 & df$log2FoldChange>=FC),]$express = "UP"
  DOWN_UP = table(df$express)
  return(df)
}

DE_data_LPS_CON = pairwise_Compare(DEgeneRes_LPS_CON, 1)
DE_data_LPSSiR_LPS = pairwise_Compare(DEgeneRes_LPSSiR_LPS, 1)

saveResults = function(DE_data, filepath, filename){
  write.csv(DE_data, file = paste(filepath, paste0(filename,".csv"), sep = "/"), row.names = F)
  write.table(DE_data[which(DE_data$express=="UP"),"symbol"], 
              file = paste(filepath, paste0(filename, "_UP.txt"), sep = "/"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(DE_data[which(DE_data$express=="DOWN"),"symbol"],
              file = paste(filepath, paste0(filename, "_DOWN.txt"), sep = "/"),
              sep = "\t", quote = F, row.names = F, col.names = F)
}

saveResults(DE_data = DE_data_LPS_CON, filepath = "/Volumes/GaoxyData/data/RNA-seq", filename = "LPS.vs.CON.DEgene")
saveResults(DE_data = DE_data_LPSSiR_LPS, filepath = "/Volumes/GaoxyData/data/RNA-seq", filename = "LPSSiR.vs.LPS.DEgene")

##############################################################################################################################################################################################################################################################################
## Protein coding gene的表达分别
LPS_CON_data = merge(DE_data_LPS_CON, Normalized_counts[,c("CON_1", "CON_2", "CON_3", "LPS_1", "LPS_2", "LPS_3", "symbol")], by = "symbol")
LPSSiR_LPS_data = merge(DE_data_LPSSiR_LPS, Normalized_counts[,c("LPS_1", "LPS_2", "LPS_3", "LPS.SiR_1", "LPS.SiR_2", "LPS.SiR_3", "symbol")], by = "symbol")
LPS_CON_data$WTcounts = rowSums(LPS_CON_data[,c("CON_1", "CON_2", "CON_3")])/3
LPS_CON_data$KOcounts = rowSums(LPS_CON_data[,c("LPS_1", "LPS_2", "LPS_3")])/3
LPSSiR_LPS_data$WTcounts = rowSums(LPSSiR_LPS_data[,c("LPS_1", "LPS_2", "LPS_3")])/3
LPSSiR_LPS_data$KOcounts = rowSums(LPSSiR_LPS_data[,c("LPS.SiR_1", "LPS.SiR_2", "LPS.SiR_3")])/3

pdf("/Volumes/GaoxyData/data/RNA-seq/Scatter_LPS_LPSSiR.pdf", width = 5, height = 5)
ggplot(LPSSiR_LPS_data, aes(log2(WTcounts),log2(KOcounts))) +
  geom_point(aes(color=express),size=0.5) +
  geom_abline(intercept = 0, slope = 1, lty=2) +
  scale_color_manual(values = c("#d01c8b","grey70","#4dac26")) +
  scale_x_continuous(limits = c(0,15), labels = c(0,5,10,15)) +
  scale_y_continuous(limits = c(0,15), labels = c(0,5,10,15)) +
  labs(x="log2(LPS Normalized counts)", y="log2(LPS+SiR Normalized counts)") +
  annotate(geom = "text", label = "UP:472", x = 1, y = 15, fontface = "bold") +
  annotate(geom = "text", label = "DOWN:380", x = 13.5, y = 0, fontface = "bold") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        )
dev.off()


##############################################################################################################################################################################################################################################################################
##heatmap
tsv = read.delim('/Volumes/GaoxyData/data/RNA-seq/HALLMARK_INFLAMMATORY_RESPONSE.v2025.1.Mm.tsv', header = T, sep = '\t')
gene_symbols = unlist(strsplit(tsv[17,2], ","))  # 按逗号分割
gene_symbols = trimws(gene_symbols)  # 去除首尾空格
gene_symbols = unique(gene_symbols)  # 去重

LPS.vs.CON_heatmap = LPS_CON_data[which(LPS_CON_data[which(LPS_CON_data$express != 'NO'), 'symbol'] %in% gene_symbols), c("CON_1", "CON_2", "CON_3", "LPS_1", "LPS_2", "LPS_3", "symbol")]
LPSSiR.vs.LPS_heatmap = LPSSiR_LPS_data[which(LPSSiR_LPS_data[which(LPSSiR_LPS_data$express != 'NO'), 'symbol'] %in% gene_symbols), c("LPS_1", "LPS_2", "LPS_3", "LPS.SiR_1", "LPS.SiR_2", "LPS.SiR_3", "symbol")]
rownames(LPS.vs.CON_heatmap) = LPS.vs.CON_heatmap$symbol
rownames(LPSSiR.vs.LPS_heatmap) = LPSSiR.vs.LPS_heatmap$symbol




# 先准备一个包含 express 的小表
direction_info <- LPSSiR_LPS_data[, c("symbol", "express")]

# 合并到热图数据
LPS.vs.CON_heatmap2 <- merge(
  LPS.vs.CON_heatmap,
  direction_info,
  by = "symbol",
  all.x = TRUE
)
# merge 后列顺序会变一下：symbol, CON_1, CON_2, ..., express
# 我们重新设置行名
rownames(LPS.vs.CON_heatmap2) <- LPS.vs.CON_heatmap2$symbol
# 明确一个排序优先级：UP first, then DOWN
LPS.vs.CON_heatmap2$express <- factor(
  LPS.vs.CON_heatmap2$express,
  levels = c("UP", "DOWN")
)

# 得到按 express 排序后的基因顺序
row_order <- order(LPS.vs.CON_heatmap2$express)

# 按这个顺序重排矩阵（只保留表达列）
mat_plot <- as.matrix(LPS.vs.CON_heatmap2[row_order, c("CON_1", "CON_2", "CON_3", "LPS_1", "LPS_2", "LPS_3")])


p = pheatmap(mat_plot,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         scale = "row", 
         cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100))

save_pheatmap_pdf(p, filename = '/Volumes/GaoxyData/data/RNA-seq/LPS.vs.CON炎症基因热图.pdf', height = 12, width = 4)

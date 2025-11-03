## GO enrichment
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(forcats)

DE_LPS_CON = read.csv("/Volumes/GaoxyData/data/RNA-seq/LPS.vs.CON.DEgene.csv")
DE_LPS_LPSSiR = read.csv("/Volumes/GaoxyData/data/RNA-seq/LPSSiR.vs.LPS.DEgene.csv")

LPS_UP = DE_LPS_CON[DE_LPS_CON$express == "UP", "symbol"]
LPS_DOWN = DE_LPS_CON[DE_LPS_CON$express == "DOWN", "symbol"]

LPSSiR_UP = DE_LPS_LPSSiR[DE_LPS_LPSSiR$express == "UP", "symbol"]
LPSSiR_DOWN = DE_LPS_LPSSiR[DE_LPS_LPSSiR$express == "DOWN", "symbol"]

LPS_up_entrez = bitr(LPS_UP, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
LPS_down_entrez = bitr(LPS_DOWN, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
LPSSiR_up_entrez = bitr(LPSSiR_UP, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
LPSSiR_down_entrez = bitr(LPSSiR_DOWN, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)

LPS_up_kegg = enrichKEGG(
  gene         = LPS_up_entrez$ENTREZID,
  organism     = "mmu",          # 小鼠
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)
LPS_down_kegg = enrichKEGG(
  gene         = LPS_down_entrez$ENTREZID,
  organism     = "mmu",          # 小鼠
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)
LPSSiR_up_kegg = enrichKEGG(
  gene         = LPSSiR_up_entrez$ENTREZID,
  organism     = "mmu",          # 小鼠
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)
LPSSiR_down_kegg = enrichKEGG(
  gene         = LPSSiR_down_entrez$ENTREZID,
  organism     = "mmu",          # 小鼠
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

LPS_UP_KEGG_df = LPS_up_kegg@result
LPS_DOWN_KEGG_df = LPS_down_kegg@result
LPSSiR_UP_KEGG_df = LPSSiR_up_kegg@result
LPSSiR_DOWN_KEGG_df = LPSSiR_down_kegg@result

pathway = c('mmu04062', 'mmu04620', 'mmu04621', 'mmu04622', 'mmu04623', 'mmu04625', 'mmu04657', 'mmu04668', 'mmu04610', 'mmu04612', 'mmu04658', 'mmu04659', 'mmu04662', 'mmu04664', 'mmu04666', 'mmu04670')

LPS_UP_s_df = LPS_UP_KEGG_df[which(LPS_UP_KEGG_df$ID %in% pathway), ]
LPS_DOWN_s_df = LPS_DOWN_KEGG_df[which(LPS_DOWN_KEGG_df$ID %in% pathway), ]
LPSSiR_UP_s_df = LPSSiR_UP_KEGG_df[which(LPSSiR_UP_KEGG_df$ID %in% pathway), ]
LPSSiR_DOWN_s_df = LPSSiR_DOWN_KEGG_df[which(LPSSiR_DOWN_KEGG_df$ID %in% pathway), ]

LPS_UP_s_df$lgpadj = -log10(LPS_UP_s_df$p.adjust)
LPS_DOWN_s_df$lgpadj = -log10(LPS_DOWN_s_df$p.adjust)
LPSSiR_UP_s_df$lgpadj = -log10(LPSSiR_UP_s_df$p.adjust)
LPSSiR_DOWN_s_df$lgpadj = -log10(LPSSiR_DOWN_s_df$p.adjust)

LPS_UP_s_df$cluster = 'LPS.vs.CON UP'
LPS_DOWN_s_df$cluster = 'LPS.vs.CON DOWN'
LPSSiR_UP_s_df$cluster = 'LPSSiR.vs.LPS UP'
LPSSiR_DOWN_s_df$cluster = 'LPSSiR.vs.LPS DOWN'

df = rbind(LPS_UP_s_df[, c('Description', 'cluster', 'Count', 'GeneRatio', 'FoldEnrichment', 'p.adjust', 'lgpadj')], 
           LPS_DOWN_s_df[, c('Description', 'cluster', 'Count', 'GeneRatio', 'FoldEnrichment', 'p.adjust', 'lgpadj')], 
           LPSSiR_UP_s_df[, c('Description', 'cluster', 'Count', 'GeneRatio', 'FoldEnrichment', 'p.adjust', 'lgpadj')], 
           LPSSiR_DOWN_s_df[, c('Description', 'cluster', 'Count', 'GeneRatio', 'FoldEnrichment', 'p.adjust', 'lgpadj')])

df$GeneRatio_num = sapply(strsplit(df$GeneRatio, "/"), function(x) (as.numeric(x[1]) / as.numeric(x[2]))*100)

p = ggplot(df) +
        geom_segment(aes(x = Description, xend = Description, y = 0, yend = GeneRatio_num),linetype = "solid", color = "#d3d3d3", linetype = "solid") +
        geom_point(aes(x = Description, y = GeneRatio_num, size = lgpadj, color = cluster)) +
        scale_color_manual(values=c("LPS.vs.CON UP"="#a8dcd4", "LPS.vs.CON DOWN"="#f29e95", "LPSSiR.vs.LPS UP"="#f6c48a", "LPSSiR.vs.LPS DOWN"="#9fc3dd"))+
        #scale_size_continuous(range = c(0, 4),breaks = c(1, 2, 3, 4),labels = c(1, 2, 3, 4)) +
        scale_y_continuous(expand = c(0,0),limits = c(0, 5.5)) +
        #scale_x_discrete(labels = scales::label_wrap(width = 23)) +
        labs(x = '', y = 'Gene Ratio', title = 'Inflammatory KEGG pathways')+
        guides(size = guide_legend(title = '-log10(P adj)'),color = guide_legend(title = '')) +
        coord_flip()+
        theme_classic() +
        theme(plot.title.position = "plot",
              plot.title = element_text(size = 22, hjust=0.5, face = 'bold'),
              legend.text = element_text(size = 12),
              axis.text.x = element_text(size=18),
              axis.text.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18))

pdf('/Volumes/GaoxyData/data/RNA-seq/KEGG_MacthStick.pdf', width = 14, height = 16)
print(p)
dev.off()




## GSEA
LPS.vs.CON = DE_LPS_CON[which(DE_LPS_CON$express %in% c("UP", "DOWN")), c("symbol", "log2FoldChange")]
LPSSiR.vs.LPS = DE_LPS_LPSSiR[which(DE_LPS_LPSSiR$express %in% c("UP", "DOWN")), c("symbol", "log2FoldChange")]

LPS.vs.CON_tbl_clean <- DE_LPS_CON[, c("symbol", "log2FoldChange")] %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%          # 高表达在前
  distinct(symbol, .keep_all = TRUE)  # 如果重复基因取一次

LPS.vs.CON_tbl_clean_gene_map <- bitr(LPS.vs.CON_tbl_clean$symbol,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Mm.eg.db)

# 合并，把 ENTREZID 加回表
LPS.vs.CON_tbl_clean_ens <- merge(LPS.vs.CON_tbl_clean,
                     LPS.vs.CON_tbl_clean_gene_map,
                     by.x = "symbol",
                     by.y = "SYMBOL")

geneList <- LPS.vs.CON_tbl_clean_ens$log2FoldChange
names(geneList) <- LPS.vs.CON_tbl_clean_ens$ENTREZID

# 排序：必须是降序
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)

gsea_kegg <- gseKEGG(
  geneList     = geneList,
  organism     = "mmu",      # 小鼠用 mmu，人用 hsa
  minGSSize    = 10,         # 最小基因集大小
  maxGSSize    = 2000,        # 最大基因集大小
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# 看结果
head(as.data.frame(gsea_kegg)[, c("ID","Description","NES","p.adjust")])

pdf('/Volumes/GaoxyData/data/RNA-seq/LPS.vs.CON下调GSEA.pdf', width = 8, height = 8)
gseaplot2(gsea_kegg,
          geneSetID = "mmu04110",    # 比如 IL-17 signaling pathway
          title = "Cell cycle")
dev.off()



## GSEA

LPS.vs.LPSSiR_tbl_clean <- DE_LPS_LPSSiR[, c("symbol", "log2FoldChange")] %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%          # 高表达在前
  distinct(symbol, .keep_all = TRUE)  # 如果重复基因取一次

LPS.vs.LPSSiR_tbl_clean_gene_map <- bitr(LPS.vs.LPSSiR_tbl_clean$symbol,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Mm.eg.db)

# 合并，把 ENTREZID 加回表
LPS.vs.LPSSiR_tbl_clean_ens <- merge(LPS.vs.LPSSiR_tbl_clean,
                     LPS.vs.LPSSiR_tbl_clean_gene_map,
                     by.x = "symbol",
                     by.y = "SYMBOL")

geneList <- LPS.vs.LPSSiR_tbl_clean_ens$log2FoldChange
names(geneList) <- LPS.vs.LPSSiR_tbl_clean_ens$ENTREZID

# 排序：必须是降序
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)

gsea_kegg2 <- gseKEGG(
  geneList     = geneList,
  organism     = "mmu",      # 小鼠用 mmu，人用 hsa
  minGSSize    = 10,         # 最小基因集大小
  maxGSSize    = 2000,        # 最大基因集大小
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

# 看结果
head(as.data.frame(gsea_kegg2)[, c("ID","Description","NES","p.adjust")])

pdf('/Volumes/GaoxyData/data/RNA-seq/LPSSiR.vs.LPS下调GSEA.pdf', width = 8, height = 8)
gseaplot(gsea_kegg2,
          geneSetID = "mmu04060",    # 比如 IL-17 signaling pathway
          title = "Cytokine-cytokine receptor interaction")
dev.off()



library(GseaVis)

gsekegg1 = DOSE::setReadable(gsea_kegg,
                             OrgDb = "org.Mm.eg.db",
                             keyType = "ENTREZID")
gsekegg2 = DOSE::setReadable(gsea_kegg2,
                             OrgDb = "org.Mm.eg.db",
                             keyType = "ENTREZID")
gseakegg_df1 = as.data.frame(gsekegg1)
gseakegg_df2 = as.data.frame(gsekegg2)

geneset1 = strsplit(gseakegg_df1[which(gseakegg_df1$ID == 'mmu04657'), 'core_enrichment'], split = '/') %>% unlist()
geneset2 = strsplit(gseakegg_df2[which(gseakegg_df2$ID == 'mmu04060'), 'core_enrichment'], split = '/') %>% unlist()

IL17_geneset = c('Cebpb', 'Ccl2', 'Ccl7', 'Cxcl2', 'Mmp9', 'S100a8', 'Lcn2')


LPS.vs.CON_mat = LPS_CON_data[which(LPS_CON_data$symbol %in% geneset1), c("symbol", "CON_1", "CON_2", "CON_3", "LPS_1", "LPS_2", "LPS_3")]
rownames(LPS.vs.CON_mat) = LPS.vs.CON_mat$symbol
colnames(LPS.vs.CON_mat)[1] = 'gene_name'

LPSSiR.vs.LPS_mat = LPSSiR_LPS_data[which(LPSSiR_LPS_data$symbol %in% geneset2), c("symbol", "LPS_1", "LPS_2", "LPS_3", "LPS.SiR_1", "LPS.SiR_2", "LPS.SiR_3")]
rownames(LPSSiR.vs.LPS_mat) = LPSSiR.vs.LPS_mat$symbol
colnames(LPSSiR.vs.LPS_mat)[1] = 'gene_name'


pdf('/Volumes/GaoxyData/data/RNA-seq/IL17_GSEA.pdf', height = 8, width = 10)
gseaNb(object = gsekegg1,
       geneSetID = 'mmu04657',
       newGsea = T,
       add.geneExpHt = T,
       exp.col = c("#3299FF", "white", "#FF2A00"),
       exp = LPS.vs.CON_mat,
       kegg = T,
       addPval = T,
       ght.geneText.size = 12,
       ght.sampleText.size = 12,
       rmSegment = T)
dev.off()


pdf('/Volumes/GaoxyData/data/RNA-seq/Cytokine_GSEA.pdf', height = 8, width = 10)
gseaNb(object = gsekegg2,
       geneSetID = 'mmu04060',
       newGsea = T,
       add.geneExpHt = T,
       exp.col = c("#3299FF", "white", "#FF2A00"),
       exp = LPSSiR.vs.LPS_mat,
       kegg = T,
       addPval = T,
       ght.geneText.size = 12,
       ght.sampleText.size = 12,
       rmSegment = T)
dev.off()

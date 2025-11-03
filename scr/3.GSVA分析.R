library(GSVA)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)


NFKB = c("Tnf","Il1b","Nfkb1","Rela","Ikbkb","Ccl2","Ccl5","Cxcl2","Nfkbia","Tnfaip3")
#IL17 = 

LPSSiR_LPS_data_df = LPSSiR_LPS_data[, c("LPS_1", "LPS_2", "LPS_3", "LPS.SiR_1", "LPS.SiR_2", "LPS.SiR_3")]
rownames(LPSSiR_LPS_data_df) = LPSSiR_LPS_data$symbol

LPS_CON_data_df = LPS_CON_data[, c("CON_1", "CON_2", "CON_3", "LPS_1", "LPS_2", "LPS_3")]
rownames(LPS_CON_data_df) = LPS_CON_data$symbol

expr_mat.LPS.vs.CON = as.matrix(LPS_CON_data_df)
expr_mat.LPSSiR.vs.LPS = as.matrix(LPSSiR_LPS_data_df)

gene_sets <- list(NFKB = NFKB)

scores <- gsva(expr_mat, gene_sets, method="ssgsea", kcdf="Gaussian")


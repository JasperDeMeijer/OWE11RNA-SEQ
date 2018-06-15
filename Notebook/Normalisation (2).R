################################################
# Normaliseren van de count data doormiddel van EdgeR en de
# significante genen van de experimenten ophalen.
################################################

counts <- read.delim(header=TRUE, file= file.choose(),row.names = 1, skip =1)
setNormalisation(counts)

setNormalisation <- function(RNAseq){
  library(edgeR)
  exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib", "NC8.glc", "NC8.glc", "NC8.rib", "NC8.rib")
  group <- factor(exp)
  # Inlezen van de data en de definitie van de kolommen
  dge <- DGEList(counts= as.matrix(RNAseq[,1:8]), group=group)
  
  keep.genes <- rowSums(cpm(dge) > 2) >= 2
  dge <- dge[keep.genes,]
  dge$samples$lib.size <- colSums(dge$counts)
  
  # Normaliseren op basis van de totale RNA concentratie. Het berekent een normalisatie factor doormiddel van TMM
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Beschrijven hoe de replicates zijn gedefinieerd.
  design <- model.matrix(~0+group, data=dge$samples)
  colnames(design) <- levels(dge$samples$group)
  
  # Schatten van de dispersie. De breedte van de verdeling wordt bepaald.
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge,design, method="power")
  dge <- estimateGLMTagwiseDisp(dge,design)
  
  #pdf("Results/EdgeR.pdf")
  plotMDS(dge)
  plotBCV(dge)
  dev.off()

  fit <- glmFit(dge, design)
  
  # Bepalen van de Fold-Changes voor de WCFS1 stam
  mc_WCFS1 <- makeContrasts(exp.r=WCFS1.rib-WCFS1.glc, levels = design)
  fit_WCFS1 <- glmLRT(fit, contrast = mc_WCFS1)
  fit_WCFS1 <- as.data.frame(fit_WCFS1$table)
  WCFS1_FC <<- fit_WCFS1
  
  # Bepalen van de Fold-Change voor de NC8 stam
  mc_NC8 <- makeContrasts(exp.r=NC8.rib-NC8.glc, levels = design)
  fit_NC8 <- glmLRT(fit, contrast = mc_NC8)
  fit_NC8 <- as.data.frame(fit_NC8$table)
  
  # Bepalen van de Fold-Change voor de glucose conditie
  mc_Glc <- makeContrasts(exp.r=NC8.glc-WCFS1.glc, levels = design)
  fit_Glc <- glmLRT(fit, contrast = mc_Glc)
  fit_Glc <- as.data.frame(fit_Glc$table)
  
  # Bepalen van de Fold-Change voor de ribose conditie
  mc_Rib <- makeContrasts(exp.r=NC8.rib-WCFS1.rib, levels = design)
  fit_Rib <- glmLRT(fit, contrast = mc_Rib)
  fit_Rib <- as.data.frame(fit_Rib$table)
  
  # Aanroepen functie getSignificantGenes.
  WCFS1_sign <<- getSignificantGenes(fit_WCFS1)
  NC8_sign <<- getSignificantGenes(fit_NC8)
  Glc_sign <<- getSignificantGenes(fit_Glc)
  Rib_sign <<- getSignificantGenes(fit_Rib)
  
  sigGenenWCFS <<- rownames(WCFS1_sign)
  
  # Samenvoegen van de NC8 stam en WCFS1 stam in een dataframe.
  WCFS1.NC8_FC_all_genes <- cbind(fit_WCFS1$logFC, fit_NC8$logFC)
  rownames(WCFS1.NC8_FC_all_genes) <- rownames(fit_WCFS1)
  colnames(WCFS1.NC8_FC_all_genes) <- c("WFCS1", "NC8")
  
  # QQ-plot maken voor de normalisatie en na de normalisatie
  #source("QQplot.R")
  #makeQQPlot(RNAseq, WCFS1.NC8_FC_all_genes)
  
  # Samenvoegen van de NC8 stam en WCFS1 stam in een dataframe.
  WCFS1.NC8_FC <- merge(WCFS1_sign, NC8_sign, by = "row.names", all=TRUE)
  rownames(WCFS1.NC8_FC) <- WCFS1.NC8_FC[,1]
  WCFS1.NC8_FC <<- WCFS1.NC8_FC[,-1]
  colnames(WCFS1.NC8_FC) <<- c("WFCS1", "NC8")
  
  # Samenvoegen van de glucose conditie en ribose conditie in een dataframe.
  glc.rib_FC <- merge(Glc_sign, Rib_sign, by = "row.names", all = TRUE)
  rownames(glc.rib_FC) <- glc.rib_FC[,1]
  glc.rib_FC <<- glc.rib_FC[,-1]
  colnames(glc.rib_FC) <<- c("Glucose", "Ribose")
  write.csv(glc.rib_FC, file = 'test.csv')
}

################################################
# Significante genen bepalen (absolute FC > 1 EN de FDR is kleiner dan 0.05).
################################################
getSignificantGenes <- function(stem){
  # Berekenen van de FDR dmv de p-value met de BH methode.
  stem$FDR <- p.adjust(stem$PValue, method = "BH")
  sign_logFC <- subset(stem, abs(logFC)>1)
  sign_pVal <- subset(sign_logFC, FDR < 0.05)
  
  # Alleen de Fold-Changes overhouden
  significantGenes <-data.frame(sign_pVal$logFC)
  rownames(significantGenes) <- rownames(sign_pVal)

  return(significantGenes)
}

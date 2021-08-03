
getSet1DEGS = function(x, foldChange = 2, padj = 0.01) {

  ax4_set1_raw <- utils::read.table(x)

  ax4_set1_veg_raw <- ax4_set1_raw[,c(1:3, 7:9)]

  colData <- data.frame(sample = colnames(ax4_set1_veg_raw), condition = c(rep("ax4_veg", 3), rep("set1_veg", 3)))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(ax4_set1_veg_raw), colData = colData, design = ~ condition)

  dds <- DESeq2::DESeq(dds)

  ax4_set1_counts <- DESeq2::counts(dds, norm = T)

  ax4_set1_res = DESeq2::results(dds)
  ax4_set1_res = as.data.frame(ax4_set1_res)
  ax4_set1_res_omit <- na.omit(ax4_set1_res)

  set1_up <- ax4_set1_res_omit[ax4_set1_res_omit$log2FoldChange > log2(foldChange) & ax4_set1_res_omit$padj < padj, ]
  set1_down <- ax4_set1_res_omit[ax4_set1_res_omit$log2FoldChange < -(log2(foldChange)) & ax4_set1_res_omit$padj < padj, ]

  set1DegList = list(set1_up = set1_up,
                     set1_down = set1_down)

  return(set1DegList)

}


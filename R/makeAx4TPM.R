
#' makeAx4TPM
#'
#' TPM (Transcript per Million) normalisation of Ax4 RNA-Seq data. Returns Genes and norm expression for those genes with above zero readcount.
#' TPM normalises gene expression with gene length included.
#' See StatQuest RPKM, FPKM and TPM, Clearly Explained!!! on YouTube for good explanation.
#'
#' This could maybe be split into two functions one which gets the width of genes from the dicty TxDb object
#' another which actually does the TPM calculation
#'
#' Extra info: https://support.bioconductor.org/p/91218/ Mike Love guy who created DESEq2. All columns now = 1e6 as he says it should
#'
#' @param x should be a 6 column un-normalised ax4 set1 readcounts file. Can be found in the dicty_chip_seq Rproj or in RNA-Seq folder in Dropbox
#'
#' @return TPM matrix
#' @export
makeAx4TPM = function(x) {

  ax4_set1_raw <- read.table(x)
  ax4_veg_raw <- ax4_set1_raw[,c(1:3)]

  dicty_genes = GenomicFeatures::genes(makeDictyGrangesfromTxDb())
  dictyGenesWidth = GenomicRanges::width(dicty_genes)
  names(dictyGenesWidth) = dicty_genes@elementMetadata$gene_id

  ax4VegRawInDictyGenesWidth = ax4_veg_raw[row.names(ax4_veg_raw) %in% names(dictyGenesWidth), ]
  dictyGenesWidthInAx4VegRaw = dictyGenesWidth[names(dictyGenesWidth) %in% row.names(ax4VegRawInDictyGenesWidth)]
  ax4VegRawInDictyGenesWidth = ax4VegRawInDictyGenesWidth[order(row.names(ax4VegRawInDictyGenesWidth)), ]
  dictyGenesWidthInAx4VegRaw = dictyGenesWidthInAx4VegRaw[order(names(dictyGenesWidthInAx4VegRaw))]

  ax4VegRawInDictyGenesWidthT = t(ax4VegRawInDictyGenesWidth)

  for (i in 1:length(dictyGenesWidthInAx4VegRaw)) {

    if (colnames(ax4VegRawInDictyGenesWidthT)[i] == names(dictyGenesWidthInAx4VegRaw)[i]) {

      ax4VegRawInDictyGenesWidthT[ , i] = ax4VegRawInDictyGenesWidthT[ , i] / dictyGenesWidthInAx4VegRaw[i]

    }
  }

  tpm.mat <- t(ax4VegRawInDictyGenesWidthT * (1e6 / colSums(t(ax4VegRawInDictyGenesWidthT))))
  tpm.mat = as.data.frame(tpm.mat)
  tpm.mat$mean = apply(tpm.mat, 1, mean)

  tpm.mat = tpm.mat[order(tpm.mat$mean), ]
  dim(tpm.mat[tpm.mat$mean == 0, ]) # 2667

  tpmMatNoZero = tpm.mat
  tpmMatNoZero$geneId = row.names(tpmMatNoZero)
  tpmMatNoZero = tpmMatNoZero[tpmMatNoZero$mean > 0 , 4:5] # could maybe turn into regex call for anything with Ax4, ax4 or WT


  return (tpmMatNoZero)

}

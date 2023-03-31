#file <- "~/Downloads/LungSquamous.ie.sample.DEG_output_IE_Excluded_vs_Hot.csv"
#outpath <- '~/Parthenon/test/GSEA'
#GSEAnalysis <- function(
#    file,
#    outpath
#    ) {
#  if (missing(file)) {
#    stop("File not provided. Please provide a valid input.")
#  }
#  if (missing(outpath)) {
#    stop("Output path not provided. Please provide a valid input.")
#  }
#  if (file.exists(outpath)==FALSE){
#    dir.create(outpath)
#  }
#}

#' Different Expression Genes Data Prepare
#'
#' @title DEGsData
#' @description Load differentially expressed genes and rank the genes by log2FC.
#' @param file the file path and name, it should be a csv file.
#' @param EntrezID If using Symbol name or Ensembl ID, please convert them first.
#' @param log2FC log2 Fold Change
#'
#' @importFrom  readr read_csv
#' @importFrom stats na.omit
#' @return a gene list ranked by log2FC
#' @export

DEGsData <- function(file,
                     EntrezID=NULL,
                     log2FC=NULL){
  if (missing(file)) {
    stop("File not provided. Please provide a valid input.")
  }
  dfg <- readr::read_csv(file)
  ranks <- dfg[,c(EntrezID, log2FC)]
  ranks <- stats::na.omit(ranks)
  rnk <- ranks$log2FoldChange
  names(rnk) <- ranks$EntrezID
  rnk.list <- stats::na.omit(rnk)
  rnk.list <- sort(rnk.list,decreasing = T)
  return(rnk.list)
}






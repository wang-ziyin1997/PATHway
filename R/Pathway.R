#' GSEA Analysis on log2FC ranked DEGs
#'
#' @title GSEAnalysis
#' @description
#' @param data
#' @param outpath
#' @importFrom clusterProfiler GSEA
#' @import enrichplot
#' @import magrittr
#'
#' @return
#' @export
GSEAnalysis <- function(
    data,
    outpath
) {
  if (missing(data)) {
    stop("Data not provided. Please provide a valid input.")
  }
  if (missing(outpath)) {
    stop("Output path not provided. Please provide a valid input.")
  }

  if (file.exists(outpath)==FALSE){
    dir.create(outpath)
  }
  gseapath <- paste(outpath,'GSEA',sep='/')
  if (file.exists(gseapath)==FALSE){
    dir.create(gseapath)}

  data(c2)
  data(go.bp)
  data(go.cc)
  data(go.mf)
  data_lists <- list(c2, go.bp, go.cc, go.mf)
  result_list <- list()
  Class <- c('C2','GO.BP','GO.CC','GO.MF')
  for (i in 1:length(data_lists)) {
    gsea.result <- clusterProfiler::GSEA(data,TERM2GENE =data_lists[[i]], pvalueCutoff = 0.5)
    gsea.result <- clusterProfiler::setReadable(gsea.result, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType="ENTREZID")
    df.gsea <- gsea.result@result
    df.count <- df.gsea %>% dplyr::group_by(ID) %>% dplyr::summarise(count = sum(stringr::str_count(core_enrichment, "/")) + 1)
    df.gsea <- dplyr::left_join(df.gsea, df.count, by = "ID") %>% dplyr::mutate(GeneRatio = count/setSize)
    df.gsea <- df.gsea %>% dplyr::filter(pvalue <= 0.05)
    df.gsea$classification <- Class[i]
    df.gsea <- df.gsea %>%
      dplyr::mutate(up_down = ifelse(enrichmentScore > 0, "up-regulated", "down-regulated"))
    result_list[[i]] <- df.gsea

    up_regulated <- df.gsea %>%
      dplyr::filter(up_down == "up-regulated") %>%
      dplyr::arrange(desc(enrichmentScore)) %>%
      dplyr::slice(1:min(10, nrow(.)))
    down_regulated <- df.gsea %>%
      dplyr::filter(up_down == "down-regulated") %>%
      dplyr::arrange(enrichmentScore) %>%
      dplyr::slice(1:min(10, nrow(.)))
    sub.gsea <- dplyr::bind_rows(up_regulated, down_regulated)

    up_regulated <- df.gsea %>%
      dplyr::filter(up_down == "up-regulated") %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice(1:min(10, nrow(.)))
    down_regulated <- df.gsea %>%
      dplyr::filter(up_down == "down-regulated") %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice(1:min(10, nrow(.)))
    sub.gsea.bbp <- dplyr::bind_rows(up_regulated, down_regulated)

    pw.list <- unique(c(sub.gsea.bbp$ID,sub.gsea$ID))

    count <- 0
    len <- as.character(length(pw.list))
    for (n in pw.list){
      temp.fig <- enrichplot::gseaplot2(gsea.result, n, pvalue_table = T, base_size  =  9)
      name <- paste(gseapath,'/', Class[i],'.',n,".png", sep = "")
      png(
        filename = name,
        width = 3152,
        height = 2000,
        units = "px",
        bg = "white",
        res = 300)
      print(temp.fig)
      dev.off()
      count <- count +1
      print(paste(as.character(count),'/',len))}
  }
  res <- dplyr::bind_rows(result_list)
  return(res)
}


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

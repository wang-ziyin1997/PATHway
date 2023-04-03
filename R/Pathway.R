#' GSEA Analysis on log2FC ranked DEGs
#'
#' @title GSEAnalysis
#' @description Generates GSEA results and GSEA plots for the top ten pathways
#'   sorted by p-value and Enrichment Score, respectively.
#' @param data a named gene list with EntrezID format sorted by log2FC.
#' @param outpath output pathway
#' @importFrom clusterProfiler GSEA setReadable
#' @importFrom dplyr group_by summarise filter mutate arrange slice left_join bind_rows
#' @importFrom enrichplot gseaplot2
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom stringr str_count
#' @importFrom grDevices png dev.off
#' @importFrom readr write_csv
#'
#' @return a set of GSEA plots and a data table
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
      dplyr::arrange(dplyr::desc(enrichmentScore)) %>%
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
  readr::write_csv(res, paste(outpath,"GSEA_results.csv",sep='/'))
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

#' GSEA Histogram Plotting
#'
#' This function generates a GSEA histogram plot using the input data
#' and saves it to the specified output path.
#'
#' @param data A data frame containing the following columns: 'ID',
#' 'enrichmentScore', 'up_down', and 'classification'.
#' @param outpath A string representing the output path where
#' the histogram plot will be saved.
#'
#' @importFrom dplyr mutate if_else group_by arrange case_when desc slice_head ungroup
#' @importFrom ggplot2 ggplot aes geom_bar facet_grid scale_fill_manual theme
#' element_blank ylab xlab scale_y_discrete ggsave
#' @importFrom stringr str_sub str_wrap
#' @importFrom magrittr %>%
#' @importFrom stats reorder
#'
#' @export
GSEA_histogram <- function(data, outpath){
  if (missing(data)) {
    stop("Data not provided. Please provide a valid input.")
  }
  required_columns <- c('ID', 'enrichmentScore', 'up_down', 'classification')
  if (!all(required_columns %in% colnames(data))) {
    stop("Input data frame must contain the following columns: ",
         paste(required_columns, collapse = ", "), ".")
  }
  if (missing(outpath)) {
    stop("Output path not provided. Please provide a valid input.")
  }

  if (file.exists(outpath)==FALSE){
    dir.create(outpath)
  }
  bar.data <- data %>%
    dplyr::mutate(
      Path = gsub("_", " ", ID),
      label = stringr::str_sub(Path, 1, 45)
    ) %>%
    dplyr::mutate(
      label = dplyr::if_else(nchar(label) >= 40, paste0(label, "..."), label)
    )

  result <- bar.data %>%
    dplyr::group_by(classification, up_down) %>%
    dplyr::arrange(dplyr::case_when(up_down == "up" ~ dplyr::desc(enrichmentScore),
                                    up_down == "down" ~ enrichmentScore)) %>%
    dplyr::slice_head(n = min(20, nrow(.))) %>%
    dplyr::ungroup()

  barplot <- ggplot2::ggplot(result, ggplot2::aes(x = enrichmentScore, y = stats::reorder(label, enrichmentScore), fill = up_down)) +
    ggplot2::geom_bar(stat = "identity", position = "identity", colour = "black", linewidth = 0.1) +
    ggplot2::facet_grid(classification~., scales = "free")+
    ggplot2::scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = "none") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('Pathway')+ ggplot2::xlab('Enrichment Score') +
    ggplot2::scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=50))
  ggplot2::ggsave(barplot,filename = paste(outpath,"bar_plot.png",sep='/'), width = 10, height = 18)
}

#' GSEA Bubble Plot
#'
#' This function creates a GSEA bubble plot for each specified classification
#' ('C2', 'GO.BP', 'GO.CC', 'GO.MF') and saves it to the provided output path.
#'
#' @param data A data frame containing the required columns: 'ID',
#' 'enrichmentScore', 'up_down', 'classification', 'GeneRatio', 'count',
#' and 'pvalue'.
#' @param outpath The output path where the bubble plot images will be saved.
#'
#' @return Saves the bubble plot images to the specified output path.
#' @export
#'
#' @importFrom dplyr filter mutate group_by arrange case_when slice_head ungroup
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradient xlab ylab scale_y_discrete theme_void ggsave
#' @importFrom stringr str_wrap
#' @importFrom gridExtra grid.arrange
GSEA_bubble <- function(data, outpath){
  if (missing(data)) {
    stop("Data not provided. Please provide a valid input.")
  }
  required_columns <- c('ID', 'enrichmentScore', 'up_down', 'classification','GeneRatio','count','pvalue')
  if (!all(required_columns %in% colnames(data))) {
    stop("Input data frame must contain the following columns: ",
         paste(required_columns, collapse = ", "), ".")
  }
  if (missing(outpath)) {
    stop("Output path not provided. Please provide a valid input.")
  }

  if (file.exists(outpath)==FALSE){
    dir.create(outpath)
  }

  for(i in c('C2','GO.BP','GO.CC','GO.MF')){
    subdata <- data %>%
      dplyr::filter(classification == i) %>%
      dplyr::mutate(Path = gsub("_", " ", ID)) %>%
      dplyr::group_by(up_down) %>%
      dplyr::arrange(dplyr::case_when(up_down == "up" ~ pvalue,
                                      up_down == "down" ~ pvalue)) %>%
      dplyr::slice_head(n = min(10, nrow(.))) %>%
      dplyr::ungroup()

    empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()
    up_plot <- empty_plot
    down_plot <- empty_plot
    if (any(subdata$up_down == "up-regulated")) {
      up_plot <- ggplot2::ggplot(subdata[subdata$up_down == "up-regulated", ], ggplot2::aes(x = GeneRatio, y = reorder(Path, GeneRatio))) +
        ggplot2::geom_point(aes(size = count, color = `pvalue`)) +
        ggplot2::scale_colour_gradient(low = "red", high = "blue") +
        ggplot2::xlab("Gene Ratio: numbers of overlap gene/numbers of genes in the pathway") +
        ggplot2::ylab("Pathways (Up)") +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
    }

    if (any(subdata$up_down == "down-regulated")) {
      down_plot <- ggplot2::ggplot(subdata[subdata$up_down == "down-regulated", ], ggplot2::aes(x = GeneRatio, y = reorder(Path, GeneRatio))) +
        ggplot2::geom_point(aes(size = count, color = `pvalue`)) +
        ggplot2::scale_colour_gradient(low = "red", high = "blue") +
        ggplot2::xlab("Gene Ratio: numbers of overlap gene/numbers of genes in the pathway") +
        ggplot2::ylab("Pathways (Down)") +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))
    }

    bbp <- gridExtra::grid.arrange(up_plot, down_plot, ncol = 1)

    ggplot2::ggsave(bbp,filename = paste(outpath,'/bubble.plot.',i,'.png',sep=''), width = 7.5, height = 10.5)
    print(paste(i,'bubble plot saved.'))
  }
}


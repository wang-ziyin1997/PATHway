
# PATHway

`PATHway` is an R package for performing Gene Set Enrichment Analysis (GSEA) and visualizing the results.

## Installation

To install the latest version of `PATHway` directly from GitHub, you can use the `devtools` package: [GitHub](https://github.com/wang-ziyin1997/PATHway) with:

``` r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("wang-ziyin1997/PATHway")
```  
  
Package clusterProfiler (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), org.Hs.eg.db (https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html), enrichplot(https://bioconductor.org/packages/release/bioc/html/enrichplot.html) need to be download by BiocManager.  

## Usage  
After installing the package, you can load it using `library(PATHway)`.  

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PATHway)
DEGsData <- PATHway::DEGsData(Yourfilepath, EntrezID=EntrezID, log2FC=log2FoldChange)
GSEAres <- PATHway::GSEAnalysis(DEGsData,Youroutpath)
PATHway::GSEA_histogram(GSEAres,Youroutpath)
PATHway::GSEA_bubble(GSEAres,Youroutpath)
```

Use test csv to test the code.
``` r
library(PATHway)
test_file <- system.file("extdata", "test_data.csv", package = "PATHway")
DEGs <- DEGsData(file = test_file, EntrezID = "EntrezID", log2FC = "log2FoldChange")
GSEAres <- PATHway::GSEAnalysis(DEGsData,Youroutpath)
PATHway::GSEA_histogram(GSEAres,Youroutpath)
PATHway::GSEA_bubble(GSEAres,Youroutpath)
```


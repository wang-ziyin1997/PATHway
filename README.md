
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
  
## Usage  
After installing the package, you can load it using `library(PATHway)`.  

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PATHway)
DEGsData <- DEGsData(Yourfilepath, EntrezID=EntrezID, log2FC=log2FoldChange)
GSEAres <- GSEAnalysis(DEGsData,Youroutpath)
GSEA_histogram(GSEAres,Youroutpath)
GSEA_bubble(GSEAres,Youroutpath)
```


library(testthat)
# test for DEGsData
context("DEGsData tests")

test_that("DEGsData requires a file as input", {
  expect_error(DEGsData(), "File not provided. Please provide a valid input.")
})

test_that("DEGsData returns a list of ranked genes", {
  test_file <- system.file("extdata", "test_data.csv", package = "PATHway")
  result <- DEGsData(file = test_file, EntrezID = "EntrezID", log2FC = "log2FoldChange")
  expect_is(result, "numeric")
  expect_named(result)
})

# test for GSEAnalysis
context("GSEAnalysis tests")

test_that("GSEAnalysis requires data and output path as input", {
  expect_error(GSEAnalysis(), "Data not provided. Please provide a valid input.")
  expect_error(GSEAnalysis(data = data.frame()), "Output path not provided. Please provide a valid input.")
})

test_that("GSEAnalysis creates output files and returns a data frame", {
  test_data <- data.frame() # Replace with appropriate test data
  test_outpath <- tempdir()
  result <- GSEAnalysis(data = test_genelist, outpath = test_outpath)
  expect_is(result, "data.frame")
  expect_file_exists(paste0(test_outpath, "/GSEA_results.csv"))
})

# test for GSEA_histogram
context("GSEA_histogram tests")

test_that("GSEA_histogram requires data and output path as input", {
  expect_error(GSEA_histogram(), "Data not provided. Please provide a valid input.")
  expect_error(GSEA_histogram(data = data.frame()), "Output path not provided. Please provide a valid input.")
})

test_that("GSEA_histogram creates output files", {
  test_data <- data.frame() # Replace with appropriate test data
  test_outpath <- tempdir()
  GSEA_histogram(data = test_data, outpath = test_outpath)
  expect_file_exists(paste0(test_outpath, "/bar_plot.png"))
})

# test for GSEA_bubble
context("GSEA_bubble tests")

test_that("GSEA_bubble requires data and output path as input", {
  expect_error(GSEA_bubble(), "Data not provided. Please provide a valid input.")
  expect_error(GSEA_bubble(data = data.frame()), "Output path not provided. Please provide a valid input.")
})

test_that("GSEA_bubble creates output files", {
  test_data <- data.frame() # Replace with appropriate test data
  test_outpath <- tempdir()
  GSEA_bubble(data = test_data, outpath = test_outpath)
  for (i in c('C2', 'GO.BP', 'GO.CC', 'GO.MF')) {
    expect_file_exists(paste0(test_outpath, "/bubble.plot.", i, ".png"))
  }
})

test_that("GSEA_bubble function handles different classifications", {
  test_data <- DEGsData(file = "test_data1.csv")
  gsea_results <- GSEAAnalysis(data = test_data, outpath = "test_output1")

  # Create a custom data frame to include all classifications
  custom_gsea_results <- gsea_results %>%
    rbind(gsea_results %>%
            dplyr::mutate(classification = "GO.BP"),
          gsea_results %>%
            dplyr::mutate(classification = "GO.CC"),
          gsea_results %>%
            dplyr::mutate(classification = "GO.MF"))

  # Test GSEA_bubble function with custom data frame
  expect_silent(GSEA_bubble(data = custom_gsea_results, outpath = "test_output1"))
})

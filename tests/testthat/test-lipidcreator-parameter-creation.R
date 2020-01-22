context("test-lipidcreator-parameter-creation.R")
library(dplyr)

test_that("Creation of LipidCreator parameters works", {
  paramsFile <- system.file("testdata","12-HETE-d8-_M-H_1--qex_fip-raw-parameters.tsv",package="flipr")
  params <- readr::read_tsv(paramsFile)
  params$combinationId <- as.factor(params$combinationId)
  lipidCreatorParams <- createLipidCreatorParameters(params)
  expect(ncol(lipidCreatorParams)==9, "Expected 9 columns")
  expectedColumnNames <- c("instrument", "class", "precursorAdduct", "fragment", "adduct", "ppmMassRange", "group", "ParKey", "ParValue")
  expect(
    setequal(expectedColumnNames, colnames(lipidCreatorParams)),
    paste("Expected the following columns: ", paste(sort(expectedColumnNames), collapse=","),
      "Got:", paste(sort(colnames(lipidCreatorParams)),collapse=","), sep=" "))
  expect(!is.null(lipidCreatorParams), sprintf("lipidCreatorParams should not be null!"))
  # print(lipidCreatorParams)

})

context("test-lipidcreator-parameter-creation.R")
library(dplyr)

test_that("Creation of LipidCreator parameters works", {
  paramsFile <- system.file("testdata","PE-PE_17_0-17_0-_M-H_1-_fip-1-of-2-NEGATIVE-raw-parameters.tsv",package="flipr")
  params <- readr::read_tsv(paramsFile)
  params$combinationId <- as.factor(params$combinationId)
  lipidCreatorParams <- createLipidCreatorParameters(params)
  expect(ncol(lipidCreatorParams)==8)
  expect(!is.null(lipidCreatorParams))
  print(lipidCreatorParams)

})
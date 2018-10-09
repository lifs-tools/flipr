context("test-flipr.R")

test_that("processing of Transition Extractor data works", {
  testDir <- tempdir()
  testFile <- "PE_17_0-17_0_qexactive-fip.tsv"
  qexFip <- system.file("testdata",testFile,package="flipr")
  file.copy(qexFip, file.path(testDir, testFile), overwrite=TRUE, copy.mode=TRUE)
  minPrecursorCollisionEnergy <- 0
  start_lower=c(
    meanlog = -10,
    sdlog = 0.01,
    scale = 0,
    shift = -minPrecursorCollisionEnergy+1
  )
  start_upper = c(
    meanlog = 10,
    sdlog = 10,
    scale = 1,
    shift = 200
  )
  lower = c(
    meanlog = -20,
    sdlog = 0.0001,
    scale = 0.000001,
    shift = -minPrecursorCollisionEnergy+1
  )
  upper = c(
    meanlog = 20,
    sdlog = 20,
    scale = 5,
    shift = 1000
  )
  results <- flip(projectDir = testDir,
    plotFormat = "png",
    filePattern = testFile,
    dataPlots = TRUE,
    minPrecursorCollisionEnergy = minPrecursorCollisionEnergy, 
    start_lower = start_lower,
    start_upper = start_upper,
    lower = lower,
    upper = upper,
    trainModel = TRUE,
    minDataPoints = 10,
    max_iter = 1000
    )
  expect(!is.null(results), "Expected non-null results!")
  print(results)
  #expect(results$name=="PE_17_0-17_0_qexactive-fip", "Expected result name to be 'PE_17_0-17_0_qexactive-fip'")
#  expect(length(names(results))==2,"Expected two named results for list")
  #expect(!is.null(results$fits),"Expected non-null results for list member 'fits'")
  #expect(!is.null(results$fits$fits),"Expected non-null results for fits list member 'fits'")
  #expect(!is.null(results$fits$params),"Expected non-null results for fits list member 'params'")
  #expect(!is.null(results$fits$CI),"Expected non-null results for fits list member 'CI'")
  #expect(!is.null(results$fits$preds),"Expected non-null results for fits list member 'preds'")
  #expect(!is.null(results$fits$nls.tibble.unfiltered),"Expected non-null results for fits list member 'nls.tibble.unfiltered'")
  #expect(!is.null(results$fits$nls.tibble),"Expected non-null results for fits list member 'nls.tibble'")
  #expect(!is.null(results$fits$preds_from_data),"Expected non-null results for fits list member 'preds_from_data'")
  #expect(!is.null(results$fits$fitinfo),"Expected non-null results for fits list member 'fitinfo'")
})

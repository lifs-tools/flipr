context("test-flipr.R")

test_that("processing of Transition Extractor data works", {
  testDir <-  file.path(file.path(tempdir(),".."), "flipr", "tests", "test-flipr-processing", strftime(Sys.time(), format="%Y-%b-%d_%H-%M-%S"))
  if(!dir.exists(testDir)) {
    dir.create(testDir, recursive = TRUE)
  }
  testFileName <- "12-HETE-d8-_M-H_1--qex_fip.tsv"
  qexFip <- system.file("testdata",testFileName,package="flipr")
  fipFile <- file.path(testDir, testFileName)
  file.copy(qexFip, fipFile, overwrite=TRUE, copy.mode=TRUE)
  start_lower=c(
    meanlog = 0,
    sdlog = 0.00001,
    scale = 0,
    shift = -50
  )
  start_upper = c(
    meanlog = 10,
    sdlog = 2,
    scale = 100,
    shift = 250
  )
  lower = c(
    meanlog = 0,
    sdlog = 0.0001,
    scale = 0,
    shift = -50
  )
  upper = c(
    meanlog = 10,
    sdlog = 2,
    scale = 100,
    shift = 250
  )

    expect(file.exists(qexFip), sprintf("Source file '%s' for test does not exist!",qexFip))
    expect(file.exists(fipFile), sprintf("Source file '%s' in test directory does not exist!", fipFile))
    results <- flip(projectDir = testDir,
      plotFormat = "png",
      filePattern = testFileName,
      dataPlots = TRUE,
      minPrecursorCollisionEnergy = 0,
      start_lower = start_lower,
      start_upper = start_upper,
      lower = lower,
      upper = upper,
      trainModel = TRUE,
      minDataPoints = 10,
      max_iter = 5000
      )
    expect(!is.null(results), "Expected non-null results!")
    # print(results)
  #expect(results$name=="PE_17_0-17_0_qexactive-fip", "Expected result name to be 'PE_17_0-17_0_qexactive-fip'")
  expect(length(names(results))==1, paste0("Expected one named results for list, but got ",length(names(results))))
  # expect(!is.null(results$fits),"Expected non-null results for list member 'fits'")
  #expect(!is.null(results$fits$fits),"Expected non-null results for fits list member 'fits'")
  #expect(!is.null(results$fits$params),"Expected non-null results for fits list member 'params'")
  #expect(!is.null(results$fits$CI),"Expected non-null results for fits list member 'CI'")
  #expect(!is.null(results$fits$preds),"Expected non-null results for fits list member 'preds'")
  #expect(!is.null(results$fits$nls.tibble.unfiltered),"Expected non-null results for fits list member 'nls.tibble.unfiltered'")
  #expect(!is.null(results$fits$nls.tibble),"Expected non-null results for fits list member 'nls.tibble'")
  #expect(!is.null(results$fits$preds_from_data),"Expected non-null results for fits list member 'preds_from_data'")
  #expect(!is.null(results$fits$fitinfo),"Expected non-null results for fits list member 'fitinfo'")
})

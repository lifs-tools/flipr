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
    results <- flipr::flip(projectDir = testDir,
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
    #print(results)
  expect(results[[1]]$name=="12-HETE-d8-_M-H_1--qex_fip", "Expected result name to be '12-HETE-d8-_M-H_1--qex_fip'")
  expect(length(results)==1, paste0("Expected one result in list, but got ",length(results)))
  expect(!is.null(results[[1]]$fits),"Expected non-null results for list member 'fits'")
  expect(!is.null(results[[1]]$fits[["QExHF03_NM_0001295.mzML"]]),"Expected non-null results for fits list raw data file 'QExHF03_NM_0001295.mzML'")
  expect(!is.null(results[[1]]$fits[["QExHF03_NM_0001295.mzML"]]$fitResults),"Expected non-null fitResults for fits list raw data file 'QExHF03_NM_0001295.mzML'")
  expect(!is.null(results[[1]]$fits[["QExHF03_NM_0001295.mzML"]]$dataPlots),"Expected non-null dataPlots for fits list raw data file 'QExHF03_NM_0001295.mzML'")
})

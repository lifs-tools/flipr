#' Calculate non linear regression models.
#'
#' Calculates non linear regression models by performing an iterative grid search within the coordinate bounds provided
#' by \code{lower} and \code{upper} vectors, starting at \code{start_lower} and \code{start_upper}.
#' Intermediate models are scored by the AIC value until convergence has been achieved as to the defaults of nls.multstart.
#'
#' @param projectDir the path to the project directory containing the '_fip.tsv' files as input.
#' @param plotFormat the plot format, as supported by \code{ggplot2::ggsave}.
#' @param filePattern the file pattern for 'fip' files.
#' @param dataPlots whether data plots (diagnostics) should be created.
#' @param minPrecursorCollisionEnergy the minimum precursor collision energy to consider for model training.
#' @param start_lower the lower bound to start the parameter grid search, argument is passed to nls_multstart.
#' @param start_upper the upper bound to start the parameter grid search, argument is passed to nls_multstart.
#' @param lower the lower bound for the parameter estimates, argument is passed to nlsLM.
#' @param upper the upper bound for the parameter estimates, argument is passed to nlsLM.
#' @param trainModel whether the model should be calculated.
#' @param minDataPoints the minimum number of data points required per fragment / adduct / ppm combination to be considered for model calculation.
#' @param max_iter the number of combinations for grid expansion starting parameters.
#' @return The list of fit tables which includes in named element \code{name} the input file name for the model data, in named element \code{fits} a list with the named elements \code{fits} (nonlinear fits), \code{params} (parameters), \code{CI} (confidence intervals for params), \code{preds} (immediate predictions), \code{nls.tibble.unfiltered} (unfiltered input data), \code{nls.tibble} (data for calculation of fits), and \code{preds_from_data} (predictions from equidistantly resampled x-value range).
#' @importFrom magrittr %>%
#' @export
flip <-
  function(projectDir = getwd(),
           plotFormat = "png",
           filePattern = "*_fip.tsv",
           dataPlots = TRUE,
           minPrecursorCollisionEnergy = 0,
           start_lower,
           start_upper,
           lower,
           upper,
           trainModel = FALSE,
           minDataPoints = 0,
           max_iter = 500) {
    message(paste("Project dir is", projectDir))
    oldwd <- getwd()
    fip_fits <- tryCatch({
      setwd(projectDir)
      fip_files <-
        list.files(path = projectDir,
                   pattern = filePattern,
                   full.names = TRUE)
      fip_fits <- purrr::map(fip_files, function(fip_file) {
        message(paste("Creating plots for file", fip_file, "\n", sep = " "))
        colspec <- readr::cols(
          .default = readr::col_double(),
          instrument = readr::col_character(),
          localDateTimeCreated = readr::col_datetime(format = ""),
          origin = readr::col_character(),
          scanNumber = readr::col_integer(),
          polarity = readr::col_character(),
          basePeakMz = readr::col_double(),
          basePeakIntensity = readr::col_double(),
          id = readr::col_character(),
          scanDefinition = readr::col_character(),
          msLevel = readr::col_integer(),
          precursorActivationType = readr::col_character(),
          precursorCollisionEnergyUnit = readr::col_character(),
          `ionInjectionTime[0]` = readr::col_character(),
          `precursorCharge[0]` = readr::col_integer(),
          msFunction = readr::col_character(),
          spectrumType = readr::col_character(),
          group = readr::col_character(),
          `foundMassRange[ppm]` = readr::col_integer(),
          species = readr::col_character(),
          precursorAdduct = readr::col_character(),
          fragment = readr::col_character(),
          adduct = readr::col_character()
        )

        originalData <-
          readr::read_tsv(fip_file, col_types = colspec)
        message(
          paste(
            "Filtering",
            nrow(originalData),
            "rows with precursorCollisionEnergy >=",
            minPrecursorCollisionEnergy
          )
        )
        originalData <-
          originalData %>% dplyr::filter(precursorCollisionEnergy >= minPrecursorCollisionEnergy)
        message(paste(nrow(originalData), "rows retained after filter!"))
        baseFileName <-
          gsub(pattern = "\\.tsv$", "", basename(fip_file))
        a4 <- list(width = 8.27, height = 11.69)
        a4r <- list(width = 11.69, height = 8.27)
        nlsFitOutputList <- list()
        if (length(unique(originalData$instrument)) != 1) {
          warning(paste(
            "instrument column must have one distinct value only, has:",
            unique(originalData$instrument)
          ))
          stopifnot(length(unique(originalData$instrument) == 1))
        }
        instrumentId <- unique(originalData$instrument)[[1]]
        skipGroupOutput <- TRUE
        originalData$originId <-
          as.integer(as.factor(originalData$origin))
        originIndexMap <-
          data.frame("origin" = originalData$origin, "id" = originalData$originId)
        readr::write_tsv(unique(originIndexMap), path = file.path(paste(
          baseFileName, "-origin-index-map.tsv", sep = ""
        )))
        # split by origin to allow better comparability / training
        splitOriginalData <-
          split(originalData, originalData$origin)
        message(paste(
          "Split data into",
          length(splitOriginalData),
          "partitions with levels",
          paste(levels(as.factor(
            originalData$origin
          )), collapse = ",")
        ))
        sodLists <-
          purrr::map(splitOriginalData, function(splitOriginalData,
                                                 lengthOfIndex) {
            sodList <- list()
            subSetData <- splitOriginalData %>%
              dplyr::mutate(
                precursorCollisionEnergyUnitLabel = dplyr::case_when(
                  precursorCollisionEnergyUnit == "electronvolt" ~ "Collision Energy [eV]",
                  precursorCollisionEnergyUnit == "normalized" ~ "Normalized Collision Energy",
                  TRUE ~ as.character(precursorCollisionEnergyUnit)
                )
              ) # captures anything else as character
            dataIndex <- unique(subSetData$originId)[[1]]
            message(paste0("Processing data from ", unique(subSetData$origin)))
            subSetData$fragadd <-
              paste(subSetData$fragment, subSetData$adduct, sep = " ")
            subSetData$fragadd <-
              factor(subSetData$fragadd, levels = unique(subSetData[order(subSetData$calculatedMass),]$fragadd))
            message(paste0(
              "Using Fragment+Adduct levels ",
              paste0(levels(subSetData$fragadd), collapse = " | ")
            ))
            color_scale <-
              ggplot2::scale_colour_hue(
                name = "Fragment",
                limits = as.character(levels(subSetData$fragadd)),
                aesthetics = c("colour", "fill")
              )
            fileName <-
              paste0(baseFileName, "-", dataIndex, "-of-", lengthOfIndex)
            message(paste0("Writing to file ", fileName))
            # browser()
            if (nrow(subSetData) > 0) {
              if (dataPlots) {
                dataPlotsList <- flipr::createDataPlots(
                  fileName = fileName,
                  subSetData = subSetData,
                  plotFormat = plotFormat,
                  plotDimensions = a4r,
                  color_scale = color_scale
                )
                sodList[["dataPlots"]] <- dataPlotsList
              } else {
                message(
                  paste(
                    "Skipping creation of data plots. Set argument 'dataPlots=TRUE' to create!"
                  )
                )
                sodList[["dataPlots"]] <- list()
              }

              if (trainModel) {
                sodList[["fitResults"]] <-
                  flipr::createFits(
                    fileName,
                    subSetData,
                    instrumentId,
                    skipGroupOutput,
                    start_lower,
                    start_upper,
                    lower,
                    upper,
                    minDataPoints,
                    max_iter,
                    plotFormat,
                    plotDimensions = a4r,
                    color_scale = color_scale
                  )
              } else {
                message(
                  paste(
                    "Skipping training of models. Set argument 'trainModel=TRUE' to create!"
                  )
                )
                sodList[["fitResults"]] <- list()
              }
            } else {
              warning(paste("No data available for", fip_file, sep = " "))
            }
            return(sodList)
          }, lengthOfIndex = length(splitOriginalData))
        return(list(name = baseFileName, fits = sodLists))
      })
    }, error = function(e) {
      print(e)
      return(list())
    }, finally = {
      setwd(oldwd)
    })
    fip_fits
  }

#' Create the fragment model fits and fit plots for an FIP table.
#'
#' @param fileName the base filename of the FIP table, without the file extension.
#' @param subsetData the data subset.
#' @param instrumentId the instrumentId of this data.
#' @param skipGroupOutput whether to skip writing of grouping information (for debugging).
#' @param start_lower the lower bound to start the parameter grid search, argument is passed to nls_multstart.
#' @param start_upper the upper bound to start the parameter grid search, argument is passed to nls_multstart.
#' @param lower the lower bound for the parameter estimates, argument is passed to nlsLM.
#' @param upper the upper bound for the parameter estimates, argument is passed to nlsLM.
#' @param minDataPoints the minimum number of data points required per fragment / adduct / ppm combination to be considered for model calculation.
#' @param max_iter the maximum number of iterations of the model to calculate.
#' @param plotFormat the plot format, as supported by \code{ggplot2::ggsave}.
#' @param plotDimensions the plot dimensions to use, when printing the plot object to the output device.
#' @param color_scale the ggplot2 color scale to use for fragments.
#' @return a list with a named entry of 'fits', containing the model fitting results, and a named entry 'fitPlots', containing ggplot grob object for the corresponding plot.
#' @export
createFits <- function(fileName,
                       subSetData,
                       instrumentId,
                       skipGroupOutput,
                       start_lower,
                       start_upper,
                       lower,
                       upper,
                       minDataPoints,
                       max_iter,
                       plotFormat = "PNG",
                       plotDimensions = list(width = 11.69,
                                             height = 8.27),
                       color_scale = ggplot2::scale_colour_hue()) {
  dfl <- split(subSetData, subSetData$polarity)
  purrr::map(dfl, function(splitData,
                       splitFileName,
                       plotFormat,
                       plotDimensions,
                       color_scale,
                       sodList) {
    polarity <- unique(splitData$polarity)
    sodList <- list()
    splitFileName <-
      paste(splitFileName, polarity, sep = "-")
    stopifnot(length(unique(splitData$group)) == 1)
    sodList[["fits"]] <- flipr::fits(
      tibble = splitData,
      outputPrefix = splitFileName,
      group = unique(splitData$group)[[1]],
      instrumentId = instrumentId,
      skipGroupOutput = skipGroupOutput,
      start_lower = start_lower,
      start_upper = start_upper,
      lower = lower,
      upper = upper,
      minDataPoints = minDataPoints,
      max_iter = max_iter
    )
    sodList[["fitPlots"]] <- flipr::plotFits(
      sodList[["fits"]],
      splitFileName,
      plotFormat = plotFormat,
      plotDimensions = plotDimensions,
      color_scale = color_scale
    )
    return(sodList)
  }, splitFileName = fileName, plotFormat = plotFormat, plotDimensions = plotDimensions, color_scale = color_scale, sodList = sodList)
}

#' Create the data qc plots for an FIP table.
#'
#' @param fileName the base filename of the FIP table, without the file extension.
#' @param subsetData the data subset.
#' @param plotFormat the plot format, as supported by \code{ggplot2::ggsave}.
#' @param plotDimensions the plot dimensions to use, when printing the plot object to the output device.
#' @param color_scale the ggplot2 color scale to use for fragments.
#' @return a list with a named entry of a ggplot grob object for the corresponding plot.
#' @export
createDataPlots <-
  function(fileName,
           subSetData,
           plotFormat = "PNG",
           plotDimensions = list(width = 11.69,
                                 height = 8.27),
           color_scale = ggplot2::scale_colour_hue()) {
    # split dataframe based on polarity
    dfl <- split(subSetData, subSetData$polarity)

    purrr::map(dfl, function(splitData) {
      dataPlots <- list()
      customTheme <-
        ggplot2::theme_get() + ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5))
      ggplot2::theme_set(customTheme)
      polarity <- unique(splitData$polarity)
      fileName <- paste(fileName, polarity, sep = "-")

      message("fragment-ppm-boxplot")
      ppms <- sort(unique(splitData$`foundMassRange[ppm]`))
      splitData$`foundMassRange[ppm]` <-
        paste(splitData$`foundMassRange[ppm]`, "[ppm]", sep = " ")
      splitData$`foundMassRange[ppm]` <-
        factor(splitData$`foundMassRange[ppm]`,
               levels = paste(ppms, "[ppm]", sep = " "))

      splitData$fragadd <-
        paste(splitData$fragment, splitData$adduct, sep = " ")
      splitData$fragadd <-
        factor(splitData$fragadd,
               levels = unique(splitData[order(splitData$calculatedMass), ]$fragadd))

      splitData$fragment <-
        factor(splitData$fragment,
               levels = unique(splitData[order(splitData$calculatedMass), ]$fragment))
      readr::write_tsv(splitData, path = file.path(paste(
        fileName, "-raw-plots-input.tsv", sep = ""
      )))
      dataPlots[["rawTicVsTotalIonCurrent"]] = flipr::plotRawTicVsTotalIonCurrent(
        data = splitData,
        basename = fileName,
        plotFormat = plotFormat,
        plotDimensions = plotDimensions
      )
      dataPlots[["fragmentPpmBoxplot"]] <-
        flipr::plotFragmentPpmBoxplot(splitData,
                                      fileName,
                                      plotFormat = plotFormat,
                                      plotDimensions = plotDimensions)

      splitData.NAs <-
        subset(splitData, subset = is.na(splitData$foundMass))
      if (nrow(splitData.NAs) > 0) {
        readr::write_tsv(splitData.NAs, path = file.path(paste(
          fileName, "-no-ions-found.tsv", sep = ""
        )))
      }
      splitData.noNAs <-
        splitData[!is.na(splitData$foundMass),]

      readr::write_tsv(splitData, path = file.path(paste(
        fileName, "-data-plots-input.tsv", sep = ""
      )))
      dataPlots[["precCollEnergyVsFoundIntensity"]] <-
        flipr::plotPrecCollEnergyVsFoundIntensity(
          splitData.noNAs,
          fileName,
          plotFormat = plotFormat,
          plotDimensions = plotDimensions,
          color_scale = color_scale
        )
      dataPlots[["precCollEnergyVsScanRelativeIntensityNormalized"]] <-
        flipr::plotPrecCollEnergyVsScanRelativeIntensityNormalized(
          splitData.noNAs,
          fileName,
          plotFormat = plotFormat,
          plotDimensions = plotDimensions,
          color_scale = color_scale
        )
      dataPlots[["precCollEnergyVsScanRelativeIntensityOverlay"]] <-
        flipr::plotPrecCollEnergyVsScanRelativeIntensityOverlay(
          splitData.noNAs,
          fileName,
          plotFormat = plotFormat,
          plotDimensions = plotDimensions,
          color_scale = color_scale
        )
      dataPlots[["precCollEnergyVsMassErrorPpm"]] <-
        flipr::plotPrecCollEnergyVsMassErrorPpm(
          splitData.noNAs,
          fileName,
          plotFormat = plotFormat,
          plotDimensions = plotDimensions,
          color_scale = color_scale
        )
      dataPlots[["massDensityDistribution"]] <-
        flipr::plotMassDensityDistribution(
          splitData.noNAs,
          fileName,
          plotFormat = plotFormat,
          plotDimensions = plotDimensions,
          color_scale = color_scale
        )
      dataPlots[["mzVsMerrPpm"]] <- flipr::plotMzVsMerrPpm(
        splitData.noNAs,
        fileName,
        plotFormat = plotFormat,
        plotDimensions = plotDimensions,
        color_scale = color_scale
      )
      dataPlots[["scanRelativeIntenstityHistogram"]] <-
        flipr::plotScanRelativeIntensityHistogram(
          splitData.noNAs,
          fileName,
          plotFormat = plotFormat,
          plotDimensions = plotDimensions,
          color_scale = color_scale
        )
      return(dataPlots)
    })
  }
#' @importFrom magrittr %>%
#' @export
flip <- function(projectDir=getwd(), plotFormat="png", filePattern="*_fip.tsv", dataPlots=TRUE, minPrecursorCollisionEnergy=0) {
  setwd(projectDir)
  fip_files <-
    list.files(path = projectDir,
               pattern = filePattern,
               full.names = TRUE)
  multiDataMap <- data.frame()
  fip_fits <- sapply(fip_files, function(fip_file) {
    message(paste("Creating plots for file", fip_file, "\n", sep = " "))
    colspec <- readr::cols(
      .default = readr::col_double(),
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
      `foundMassRange[ppm]` = readr::col_integer(),
      species = readr::col_character(),
      fragment = readr::col_character(),
      adduct = readr::col_character()
    )
    originalData <- readr::read_tsv(fip_file, col_types = colspec)
    message(paste("Retaining rows with precursorCollisionEnergy >=",minPrecursorCollisionEnergy))
    originalData <- originalData %>% filter_all(precursorCollisionEnergy>=minPrecursorCollisionEnergy)
    message(paste(nrow(originalData)," rows retained after filter!"))
    baseFileName <- gsub(pattern = "\\.tsv$", "", basename(fip_file))
    a4 <- list(width=8.27, height=11.69)
    a4r <- list(width=11.69, height=8.27)
    nlsFitOutputList <- list()
    # TODO: split by origin to allow better comparability / training
    splitOriginalData <- split(originalData, originalData$origin)
    splitOriginalDataIndex <- seq_along(splitOriginalData)
    multiDataMap <- rbind(multiDataMap, data.frame(name=baseFileName, index=splitOriginalDataIndex, origin=splitOriginalData))
    lapply(splitOriginalData, function(data, dataIndex, lengthOfIndex) {
      message(paste0("Processing data from ", unique(data$origin)))

      fileName <- paste0(baseFileName, "-", dataIndex,"-of-",lengthOfIndex)
      if (nrow(data) > 0) {
        nlsFitOutputList <- flipr::fits(data, fileName)
        flipr::plotFits(nlsFitOutputList, fileName, format=plotFormat)
        if(dataPlots) {
          # split dataframe based on polarity
          dfl <- split(data, data$polarity)

          lapply(dfl, function(splitData) {
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

            splitData$fragadd <- paste(splitData$fragment, splitData$adduct, sep = " ")
            splitData$fragadd <-
              factor(splitData$fragadd, levels = unique(splitData[order(splitData$calculatedMass),]$fragadd))

            splitData$fragment <-
              factor(splitData$fragment, levels = unique(splitData[order(splitData$calculatedMass),]$fragment))

            flipr::plotRawTicVsTotalIonCurrent(data=splitData, basename = fileName, plotFormat=plotFormat, plotDimensions=a4r)
            flipr::plotFragmentPpmBoxplot(splitData, fileName, plotFormat=plotFormat, plotDimensions=a4r)

            splitData.NAs <- subset(splitData, subset = is.na(splitData$foundMass))
            if (nrow(splitData.NAs) > 0) {
              readr::write_tsv(splitData.NAs, path = file.path(paste(
                fileName, "-no-ions-found.tsv", sep = ""
              )))
            }
            splitData.noNAs <- splitData[!is.na(splitData$foundMass), ]

            flipr::plotPrecCollEnergyVsFoundIntensity(splitData.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
            flipr::plotPrecCollEnergyVsScanRelativeIntensityNormalized(splitData.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
            flipr::plotPrecCollEnergyVsScanRelativeIntensityOverlay(splitData.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
            flipr::plotPrecCollEnergyVsMassErrorPpm(splitData.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
            flipr::plotMassDensityDistribution(splitData.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
            flipr::plotMzVsMerrPpm(splitData.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          })
        } else {
          message(paste("Skipping creation of data plots. Set argument 'dataPlots=TRUE' to create!"))
        }
        # charge and adduct dependency, mass error influence on intensity ?

        # subset by fragment
      } else {
        warning(paste("No data available for", fip_file, sep = " "))
      }
    }, dataIndex=splitOriginalDataIndex, lengthOfIndex=length(splitOriginalDataIndex))
    list(name=fileName, fits=nlsFitOutputList)
  }, simplify = FALSE, USE.NAMES = TRUE)
  readr::write_tsv(multiDataMap, path=file.path("origin-index-data-map.tsv"))
  fip_fits
}

#' @importFrom magrittr %>%
#' @export
flip <- function(projectDir=getwd(), plotFormat="png", filePattern="*_fip.tsv", dataPlots=TRUE, minPrecursorCollisionEnergy=0,
start_lower, start_upper, lower, upper, trainModel=FALSE) {
  setwd(projectDir)
  fip_files <-
    list.files(path = projectDir,
               pattern = filePattern,
               full.names = TRUE)
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
      instrument = readr::col_character(),
      group = readr::col_character(),
      `foundMassRange[ppm]` = readr::col_integer(),
      species = readr::col_character(),
      fragment = readr::col_character(),
      adduct = readr::col_character()
    )
    originalData <- readr::read_tsv(fip_file, col_types = colspec)
    message(paste("Filtering",nrow(originalData), "rows with precursorCollisionEnergy >=",minPrecursorCollisionEnergy))
    originalData <- originalData %>% dplyr::filter(precursorCollisionEnergy>=minPrecursorCollisionEnergy)
    message(paste(nrow(originalData),"rows retained after filter!"))
    baseFileName <- gsub(pattern = "\\.tsv$", "", basename(fip_file))
    a4 <- list(width=8.27, height=11.69)
    a4r <- list(width=11.69, height=8.27)
    nlsFitOutputList <- list()
    if(length(unique(originalData$instrument))!=1) {
      warning(paste("instrument column must have one distinct value only, has:",unique(originalData$instrument)))
      stopifnot(length(unique(originalData$instrument)==1))
    }
    instrumentId <- unique(originalData$instrument)[[1]]
    skipGroupOutput <- TRUE
    originalData$originId <- as.integer(as.factor(originalData$origin))
    originIndexMap <- data.frame("origin"=originalData$origin, "id"=originalData$originId)
    readr::write_tsv(unique(originIndexMap), path=file.path(paste(baseFileName, "-origin-index-map.tsv", sep ="")))
    # TODO: split by origin to allow better comparability / training
    splitOriginalData <- split(originalData, originalData$origin)
    #print(splitOriginalData)
    message(paste("Split data into",length(splitOriginalData),"partitions with levels", paste(levels(as.factor(originalData$origin)),collapse=",")))
    lapply(splitOriginalData, function(splitOriginalData, lengthOfIndex) {
      subSetData <- splitOriginalData
      dataIndex <- unique(subSetData$originId)[[1]]
      print(dataIndex)
      message(paste0("Processing data from ", unique(subSetData$origin)))

      fileName <- paste0(baseFileName, "-", dataIndex,"-of-",lengthOfIndex)
      message(paste0("Writing to file ",fileName))
      if (nrow(subSetData) > 0) {

        if(dataPlots) {
          # split dataframe based on polarity
          dfl <- split(subSetData, subSetData$polarity)

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
            flipr::plotScanRelativeIntensityHistogram(splitData.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          })
        } else {
          message(paste("Skipping creation of data plots. Set argument 'dataPlots=TRUE' to create!"))
        }

        if(trainModel) {
          dfl <- split(subSetData, subSetData$polarity)
          lapply(dfl, function(splitData, splitFileName, plotFormat){
            polarity <- unique(splitData$polarity)
            splitFileName <- paste(splitFileName, polarity, sep="-")
            stopifnot(length(unique(splitData$group))==1)
            nlsFitOutputList <- flipr::fits(tibble=splitData,
                                            outputPrefix=splitFileName,
                                            group=unique(splitData$group)[[1]],
                                            instrumentId=instrumentId,
                                            skipGroupOutput=skipGroupOutput,
                                            start_lower=start_lower,
                                            start_upper=start_upper,
                                            lower=lower,
                                            upper=upper)
            flipr::plotFits(nlsFitOutputList, splitFileName, plotFormat=plotFormat)
          },fileName, plotFormat)
        }
        # charge and adduct dependency, mass error influence on intensity ?

        # subset by fragment
      } else {
        warning(paste("No data available for", fip_file, sep = " "))
      }
    }, lengthOfIndex=length(splitOriginalData))
    list(name=baseFileName, fits=nlsFitOutputList)
  }, simplify = FALSE, USE.NAMES = TRUE)
  fip_fits
}

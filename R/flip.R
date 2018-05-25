#' @importFrom magrittr %>%
#' @export
flip <- function(projectDir=getwd(), plotFormat="png", filePattern="*_fip.tsv", dataPlots=TRUE) {
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
      `foundMassRange[ppm]` = readr::col_integer(),
      species = readr::col_character(),
      fragment = readr::col_character(),
      adduct = readr::col_character()
    )
    data <- readr::read_tsv(fip_file, col_types = colspec)
    fileName <- gsub(pattern = "\\.tsv$", "", basename(fip_file))
    a4 <- list(width=8.27, height=11.69)
    a4r <- list(width=11.69, height=8.27)
    nlsFitOutputList <- list()
    if (nrow(data) > 0) {
      nlsFitOutputList <- flipr::fits(data, fileName)
      flipr::plotFits(nlsFitOutputList, fileName, format=plotFormat)
      if(dataPlots) {
        # split dataframe based on polarity
        dfl <- split(data, data$polarity)

        lapply(dfl, function(dt) {
          customTheme <-
            ggplot2::theme_get() + ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5))
          ggplot2::theme_set(customTheme)
          polarity <- unique(dt$polarity)
          fileName <- paste(fileName, polarity, sep = "-")

          message("fragment-ppm-boxplot")
          ppms <- sort(unique(dt$`foundMassRange[ppm]`))
          dt$`foundMassRange[ppm]` <-
            paste(dt$`foundMassRange[ppm]`, "[ppm]", sep = " ")
          dt$`foundMassRange[ppm]` <-
            factor(dt$`foundMassRange[ppm]`,
                   levels = paste(ppms, "[ppm]", sep = " "))

          dt$fragadd <- paste(dt$fragment, dt$adduct, sep = " ")
          dt$fragadd <-
            factor(dt$fragadd, levels = unique(dt[order(dt$calculatedMass),]$fragadd))

          dt$fragment <-
            factor(dt$fragment, levels = unique(dt[order(dt$calculatedMass),]$fragment))

          flipr::plotRawTicVsTotalIonCurrent(dt, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          flipr::plotFragmentPpmBoxplot(dt, fileName, plotFormat=plotFormat, plotDimensions=a4r)

          dt.NAs <- subset(dt, subset = is.na(dt$foundMass))
          if (nrow(dt.NAs) > 0) {
            readr::write_tsv(dt.NAs, path = file.path(paste(
              fileName, "-no-ions-found.tsv", sep = ""
            )))
          }
          dt.noNAs <- dt[!is.na(dt$foundMass), ]

          flipr::plotPrecCollEnergyVsFoundIntensity(dt.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          flipr::plotPrecCollEnergyVsScanRelativeIntensityNormalized(dt.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          flipr::plotPrecCollEnergyVsScanRelativeIntensityOverlay(dt.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          flipr::plotPrecCollEnergyVsMassErrorPpm(dt.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          flipr::plotMassDensityDistribution(dt.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
          flipr::plotMzVsMerrPpm(dt.noNAs, fileName, plotFormat=plotFormat, plotDimensions=a4r)
        })
      } else {
        message(paste("Skipping creation of data plots. Set argument 'dataPlots=TRUE' to create!"))
      }
      # charge and adduct dependency, mass error influence on intensity ?

      # subset by fragment
    } else {
      warning(paste("No data available for", fip_file, sep = " "))
    }
    list(name=fileName, fits=nlsFitOutputList)
  }, simplify = FALSE, USE.NAMES = TRUE)
  fip_fits
}

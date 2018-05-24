#library("ggplot2")
#library("readr")
#library("plyr")
#library("dplyr")
#source("fitDistributions.R")

#' @export
flip <- function(projectDir=getwd(), plotFormat="png", filePattern="*_fip.tsv", dataPlots=TRUE) {
  setwd(projectDir)
  fip_files <-
    list.files(path = projectDir,
               pattern = filePattern,
               full.names = TRUE)
  sapply(fip_files, function(fip_file) {
    print(paste("Creating plots for file", fip_file, "\n", sep = " "))
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
    basename <- gsub(pattern = "\\.tsv$", "", fip_file)
    a4 <- list(width=8.27, height=11.69)
    a4r <- list(width=11.69, height=8.27)
    if (nrow(data) > 0) {
      nlsFitOutputList <- flipr::fits(data, basename)
      flipr::plotFits(nlsFitOutputList, basename, format=plotFormat)
      if(dataPlots) {
        # split dataframe based on polarity
        dfl <- split(data, data$polarity)

        lapply(dfl, function(dt) {
          customTheme <-
            ggplot2::theme_get() + ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5))
          ggplot2::theme_set(customTheme)
          polarity <- unique(dt$polarity)
          basename <- paste(basename, polarity, sep = "-")

          print("fragment-ppm-boxplot")
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

          plotRawTicVsTotalIonCurrent(dt, basename, plotFormat=plotFormat, plotDimensions=a4r)
          plotFragmentPpmBoxplot(dt, basename, plotFormat=plotFormat, plotDimensions=a4r)

          dt.NAs <- subset(dt, subset = is.na(dt$foundMass))
          if (nrow(dt.NAs) > 0) {
            write_tsv(dt.NAs, path = file.path(paste(
              basename, "-no-ions-found.tsv", sep = ""
            )))
          }
          dt.noNAs <- dt[!is.na(dt$foundMass), ]

          plotPrecCollEnergyVsFoundIntensity(dt.noNAs, basename, plotFormat=plotFormat, plotDimensions=a4r)
          plotPrecCollEnergyVsScanRelativeIntensityNormalized(dt.noNAs, basename, plotFormat=plotFormat, plotDimensions=a4r)
          plotPrecCollEnergyVsScanRelativeIntensityOverlay(dt.noNAs, basename, plotFormat=plotFormat, plotDimensions=a4r)
          plotPrecCollEnergyVsMassErrorPpm(dt.noNAs, basename, plotFormat=plotFormat, plotDimensions=a4r)
          plotMassDensityDistribution(dt.noNAs, basename, plotFormat=plotFormat, plotDimensions=a4r)
          plotMzVsMerrPpm(dt.noNAs, basename, plotFormat=plotFormat, plotDimensions=a4r)
        })
      } else {
        print(paste("Skipping creation of data plots. Set argument 'dataPlots=TRUE' to create!"))
      }
      # charge and adduct dependency, mass error influence on intensity ?

      # subset by fragment
    } else {
      print(paste("No data available for", basename, sep = " "))
    }

  })
}

plotRawTicVsTotalIonCurrent <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  rawTicVsTotalIonCurrentPlot <-
    ggplot2::qplot(
      x = rawTic,
      y = totalIonCurrent,
      data = data,
      color = fragadd,
      geom = "point",
      main = paste(unique(dt$species), unique(data$polarity), sep = " ")
    )  + ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap( ~ data$`foundMassRange[ppm]`, ncol =6) +
    ggplot2::ylab("Total Ion Current [a.u.]") + ggplot2::xlab("Raw Total Ion Current [a.u.]")
  ggplot2::ggsave(
    rawTicVsTotalIonCurrentPlot,
    filename = paste0(basename, "-rawTic-vs-Tic.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
}

plotFragmentPpmBoxplot <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  fpbPlot <-
    ggplot2::qplot(
      x = fragadd,
      y = (calculatedMass - foundMass) ,
      data = data,
      color = fragadd,
      geom = c("violin"),
      main = paste(unique(data$species), unique(data$polarity), sep = " ")
    ) + ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap( ~ data$`foundMassRange[ppm]`, ncol =
                  6) + ggplot2::ylab(expression(paste(
                    Delta, " of calculated and found m/z", sep = " "
                  ))) + ggplot2::xlab("Fragment") + ggplot2::geom_jitter(ggplot2::aes(
                    x = fragadd,
                    y = (calculatedMass - foundMass),
                    color = fragadd
                  ),
                  data = data,
                  alpha = 0.1)
  ggplot2::ggsave(
    fpbPlot,
    filename = paste0(basename, "-frag-ppm-box.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  fpbPlot
}

plotPrecCollEnergyVsFoundIntensity <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  print("precursorCollisionEnergy-vs-foundIntensity")

  pceFiPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = foundIntensity,
      data = data,
      color = fragadd,
      geom = c("point"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Absolute Intensity [a.u.]") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6)
  ggplot2::ggsave(
    pceFiPlot,
    filename = paste0(basename, "-precCE-vs-I.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
}

plotPrecCollEnergyVsScanRelativeIntensityNormalized <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  print("precursorCollisionEnergy-vs-foundIntensity-scan-relative-normalized")
  precCeVsIsrnPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = scanRelativeIntensity,
      data = data,
      color = fragadd,
      geom = c("point"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=round_any(max(dt.noNAs$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Scan Relative Intensity") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6)
  ggplot2::ggsave(
    precCeVsIsrnPlot,
    filename = paste0(basename, "-precCE-vs-I-srn.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  precCeVsIsrnPlot
}

plotPrecCollEnergyVsScanRelativeIntensityOverlay <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {

  print(
    "precursorCollisionEnergy-vs-foundIntensity-scan-relative-normalized-overlay"
  )
  precCeVsIsrnOverlayPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = scanRelativeIntensity,
      data = data,
      color = fragadd,
      fill = fragadd,
      # geom = c("boxplot"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) +
    ggplot2::geom_area(position = "identity", alpha = 0.3) +
    ggplot2::scale_x_continuous(breaks = seq(from=-10, to=round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Scan Relative Intensity") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment', fill = 'Fragment') +
    ggplot2::facet_wrap(. ~ data$`foundMassRange[ppm]`, ncol = 2)
  ggplot2::ggsave(
    precCeVsIsrnOverlayPlot,
    filename = paste0(basename, "-precCE-vs-I-srn-overlay.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  precCeVsIsrnOverlayPlot
}

plotPrecCollEnergyVsMassErrorPpm <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  print("precursorCollisionEnergy-vs-mass-error-ppm")
  precCeVsMerrPpmmPlot <-
    ggplot2::qplot(
      x = precursorCollisionEnergy,
      y = `foundMassError[ppm]`,
      data = data,
      color = fragadd,
      geom = c("point"),
      main = paste(
        unique(data$species),
        unique(data$polarity),
        sep = " "
      )
    ) + ggplot2::geom_smooth(
      method = "loess",
      colour = "blue",
      se = TRUE,
      level = 0.95
    ) + ggplot2::scale_x_continuous(breaks = seq(from=-10, to=round_any(max(data$precursorCollisionEnergy), 10, f = ceiling)+10, by=10)) +
    ggplot2::ylab("Mass Error [ppm]") +
    ggplot2::xlab("Precursor Collision Energy [eV]") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(fragadd ~ data$`foundMassRange[ppm]`, ncol = 6)
  ggplot2::ggsave(
    precCeVsMerrPpmmPlot,
    filename = paste0(basename, "-precCE-vs-merr-ppm.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  precCeVsMerrPpmmPlot
}

plotMassDensityDistribution <- function(data, basename, plotFormat="png", nppmLevels, plotDimensions=list(width=11.69, height=8.27)) {
  print("m/z density distribution")
  stopifnot(!missing(data))
  stopifnot(!missing(basename))
  nppms <- length(unique(data$`foundMassRange[ppm]`))
  #plot m/z density distributions
  mddPlot <-
    ggplot2::ggplot(ggplot2::aes(
      x = foundMass,
      color = fragadd,
      fill = fragadd
    ), data = data) +
    ggplot2::ggtitle(paste(
      unique(data$species),
      unique(data$polarity),
      sep = " "
    )) + ggplot2::geom_density(ggplot2::aes(y = ..count.. / sum(..count..)), data = data) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = calculatedMass),
               data = data,
               linetype = 3) +
    ggplot2::geom_text(
      x = data$calculatedMass,
      y = 0,
      label = round(data$calculatedMass, digits = 4),
      check_overlap = TRUE,
      colour = "blue",
      hjust = "left",
      vjust = "top",
      size = 3
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = foundMassLowerBound),
      colour = "darkblue",
      data = data,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = foundMassUpperBound),
      colour = "darkblue",
      data = data,
      linetype = 2
    ) +
    ggplot2::facet_wrap( ~ `foundMassRange[ppm]`, nrow = nppms) +
    ggplot2::xlab("Fragment m/z") +
    ggplot2::ylab("Normalized Count") +
    ggplot2::labs(colour = 'Fragment', fill = 'Fragment')
  ggplot2::ggsave(
    mddPlot,
    filename = paste0(basename, "-mass-dens-distr.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  mddPlot
}

plotMzVsMerrPpm <- function(data, basename, plotFormat="png", plotDimensions=list(width=11.69, height=8.27)) {
  # plot m/z vs mass error
  print("m/z-vs-mass-error-ppm")
  mzVsMerrPpmPlot <- ggplot2::qplot(
    x = foundMass,
    y = `foundMassError[ppm]`,
    data = data,
    color = fragadd,
    geom = c("violin"),
    main = paste(
      unique(data$species),
      unique(data$polarity),
      sep = " "
    )
  ) + ggplot2::ylab("Mass Error [ppm]") +
    ggplot2::xlab("Fragment m/z") +
    ggplot2::labs(colour = 'Fragment') +
    ggplot2::facet_wrap(~ data$`foundMassRange[ppm]`, ncol = 6) +
    ggplot2::geom_jitter(
      ggplot2::aes(x = foundMass, y = `foundMassError[ppm]`, color = fragadd),
      data = data,
      alpha = 0.1
    )
  ggplot2::ggsave(
    mzVsMerrPpmPlot,
    filename = paste0(basename, "-mz-vs-merr-ppm.", plotFormat),
    width = plotDimensions$width,
    height = plotDimensions$height
  )
  plotMzVsMerrPpm
}

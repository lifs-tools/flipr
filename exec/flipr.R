library("optparse")
library("flipr")
message(paste0("Using flipr ", packageVersion("flipr")))
options(show.error.locations = TRUE)
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
  make_option("--projectDir", default=getwd(),
              help = "The projectDir containing the *_fip.tsv files [default \"%default\"]"),
  make_option("--plotFormat", default="png",
              help = "The file format used for plots (png, pdf, svg) [default \"%default\"]"),
  make_option("--filePattern", default="*_fip.tsv$", help = "The filePattern of files to process [default \"%default\"]"),
  make_option("--dataPlots", default=TRUE, help = "If true, diagnostic data plots will be generated [default \"%default\"]"),
  make_option("--minPrecursorCollisionEnergy", default="0", help = "Minimum precursor collision energy to use for training [default \"%default\"]"),
  make_option("--config", help = "If provided (as an .R file), read settings for non-linear regression search bounds [default \"%default\"]"),
  make_option("--trainModel", default=FALSE, help = "If true, train model based on provided data [default \"%default\"]"),
  make_option("--minDataPoints", default=0, help = "Minimum number of data points required for model calculation [default \"%default\"]"),
  make_option("--maxCombinations", default=100, help = "Number of combinations for grid expansion starting parameters. [default \"%default\"]")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- optparse::parse_args(OptionParser(option_list=option_list))

minPrecursorCollisionEnergy <- as.numeric(opt$minPrecursorCollisionEnergy)


if(!is.null(opt$config)) {
  source(opt$config, local = TRUE)
  message(paste("Loaded config file", system.file(opt$config)))
} else {
  message(paste("Using default settings for parameter optimization range!"))
  # defined default search bounds for nls

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

}

message("Using parameter optimization values:")
message(paste(
  "start_lower=",
  paste(
    names(start_lower),
    start_lower,
    sep = "=",
    collapse = ", "
  )
))
message(paste(
  "start_upper=",
  paste(
    names(start_upper),
    start_upper,
    sep = "=",
    collapse = ", "
  )
))
message(paste(
  "lower=",
  paste(
    names(lower), lower, sep = "=", collapse = ", "
  )))
message(paste(
  "upper=",
  paste(
    names(upper), upper, sep = "=", collapse = ", "
  )
))

flipFits <- flipr::flip(projectDir = opt$projectDir,
                        plotFormat = opt$plotFormat,
                        filePattern = opt$filePattern,
                        dataPlots = opt$dataPlots,
                        minPrecursorCollisionEnergy=minPrecursorCollisionEnergy,
                        start_lower=start_lower,
                        start_upper=start_upper,
                        lower=lower,
                        upper=upper,
                        trainModel=opt$trainModel,
                        minDataPoints=opt$minDataPoints,
                        max_iter=opt$maxCombinations)

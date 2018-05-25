library("optparse")
library("flipr")
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
  make_option("--dataPlots", default=TRUE, help = "If true, diagnostic data plots will be generated [default \"%default\"]")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))
flipFits <- flip(projectDir = opt$projectDir, plotFormat = opt$plotFormat, filePattern = opt$filePattern, dataPlots = opt$dataPlots)

# FlipR <img src="man/figures/flipr-site.png" align="right" /> #

[![Build Status](https://travis-ci.org/lifs-tools/flipr.svg?branch=master)](https://travis-ci.org/lifs-tools/flipr)

FlipR is a package to help in optimal collision energy selection for molecules based on MS2 data that was gathered with increasing collision energies. It therefor fits a parameterized log-normal distribution to the relative intensity of one particular fragment ion (incl. adduct) over the range of collision energies. We use the nls.multstart package to select a model that minimizes the AIC value within a given range of parameters that are used in a grid search, starting in an initially user-determined range. The log-normal distribution is able to adapt quite well to different fragments and their distribution profiles. Some cases though seemingly show different fragmentation behaviour that is usually visible through a less steep decline (flatter curve) of the scan relative intensity over the range of fragmentation energies.

It has been tested on Thermo QExactive HF and Waters QTof instruments.

## Installation ##

Install devtools:

For Ubuntu, please install the following libraries:

```
  sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
```

Your will require R>=3.4.4 for flipr. Start administrator's R session as follows, to install the packages into the site-wide (all users) library:

```
  sudo -i R
```

Install devtools and wait until installation has finished:

```R
  install.packages("devtools")
```

Load the devtools library:

```R
  library(devtools)
```
  
Run

```R
  install_git("https://github.com/lifs-tools/flipr.git")
```

Done!

## Usage ##

### Creating a transition list with LipidCreator



### Measuring a standard sample on an MS platform

### Converting the raw MS data

We use msConvert to convert the raw MS data for different vendor platforms into mzML.


### Extracting fragment ion traces from the MS file 

Transition extractor 

### Fitting a model for multiple fragment ion traces over a CE sequence

```R
  library("flipr")
  xics_file <- system.file("extdata", "12-HETE-d8-_M-H_1--qex_fip.tsv", package = "flipr", mustWork = TRUE)
  fits <- flip(projectDir = dirname(xics_file), plotFormat = "pdf", filePattern = basename(xics_file), dataPlots = TRUE)
```



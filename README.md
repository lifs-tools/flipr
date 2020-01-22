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

## Building the Docker image ##

The flipr package can also be packaged as a Docker container. To build it for your local use, use the following command under Linux:

```
./build-docker.sh -d=docker/flipr/ -v=MYTAG 
```

`MYTAG` would usually be the package version, e.g. `1.0.6`. This command will build the image and will tag it under the version provided with `-v`.

If you want to change the default registry prefix used for tagging, edit the `config.sh` file. Also, the tagging base (the first part after the registry) can 
be defined in that file.

### Running the Docker image ###

To run flipr using the docker image from your current directory (under Linux), run it as follows (assuming you have not changed `config.sh`):

```
docker run -v "$PWD":/home/data/ do1-aps-feris.isas.de:5000/lifs/flipr:1.0.6 --projectDir=/home/data --trainModel=TRUE --dataPlots=TRUE
```

This will mount the contents of the current directory (PWD) into the `/home/data` folder within the container. Then, Docker will run the image tagged 
`do1-aps-feris.isas.de:5000/lifs/flipr:1.0.6` on your local machine. The commands following are passed to the `flipr.R` script, which is located in the `exec` folder of the package.
To see, which options are available, add `--help` to the command. It is advisable to add a custom `config.R` file to define the boundaries for the optimization. You can use it as follows:


```
docker run -v "$PWD":/home/data/ do1-aps-feris.isas.de:5000/lifs/flipr:1.0.6 --projectDir=/home/data --trainModel=TRUE --dataPlots=TRUE --config=config.R
```
This assumes, that `config.R` is available within the `projectDir` directory. All file and directory arguments apply to the file system **within** the container!


Example files are available from the `inst/testdata` folder of this project.

1. Create a new folder `flipr-test`
2. Copy the files `inst/testdata/qex-HF-config.R` and `inst/testdata/12-HETE-d8-_M-H_1-qex_fip.tsv` into `flipr-test`
3. Change to the `flipr-test` folder
4. Run the flipr Docker container: `docker run --rm -v "$PWD":/home/data/ do1-aps-feris.isas.de:5000/lifs/flipr:1.0.6 --projectDir=/home/data --trainModel=TRUE --dataPlots=TRUE --config=/home/data/qex-HF-config.R`
5. Examine the execution of the container and open the results in `flipr-test`

You will see multiple diagnostic plots of the extracted signals, as well as plots generated to display the model fit and various metrics for the fit quality. Each plot will contain display information for each of the fragments defined for one precursor molecule and adduct combination. 


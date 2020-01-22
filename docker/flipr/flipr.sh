#!/bin/bash

FLIPR_PATH=$(R -e "system.file('exec','flipr.R',package='flipr');" | egrep -oe "(/[A-Za-z0-9_.-]+)+")
if [ -z $FLIPR_PATH ]; then
  echo "Please make sure that the flipr package has been installed on your system!"
  echo "Consult https://gitlab.isas.de/hoffmann/flipr for installation instructions!"
  exit 1
fi

RSCRIPT=$(command -v Rscript)

if [ -z $RSCRIPT ]; then
  echo "Rscript is not available! Please install an R distribution on your system."
  echo "R is available free of charge from https://www.r-project.org/"
  exit 1
fi

echo "Using flipr script at $FLIPR_PATH"
$RSCRIPT $FLIPR_PATH "$@"
exit $?

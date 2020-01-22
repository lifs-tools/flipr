FROM rocker/tidyverse:3.6.2

MAINTAINER nils.hoffmann@isas.de

ADD install.R /tmp/
# invalidates cache every 24 hours
ADD http://master.bioconductor.org/todays-date /tmp/
#need to install the packages globally
RUN R -f /tmp/install.R

ADD install-custom.R /tmp/
RUN R -f /tmp/install-custom.R

# Biocontainers labels
LABEL base.image="rocker/tidyverse"
LABEL version="${APP_VERSION}"
LABEL software="flipr"
LABEL software.version="${APP_VERSION}"
LABEL description="A command line interface for the flipr library for CE model training"
LABEL website="https://github.com/lifs-tools/flipr"
LABEL documentation="https://github.com/lifs-tools/flipr"
LABEL license="https://github.com/lifs-tools/flipr/blob/master/LICENSE"
LABEL tags="Metabolomics Lipidomics"

RUN mkdir /home/data && chmod 755 /home/data

ADD flipr.sh /tmp/

RUN cp /tmp/flipr.sh /bin/flipr && chmod 755 /bin/flipr
ENTRYPOINT ["bash","/bin/flipr"]
LABEL application.version="1.0.6-testing"
LABEL application.build.timestamp="2020-01-22T17:40:22+01:00"
LABEL application.name="flipr"
LABEL application.authors="Nils Hoffmann"
LABEL application.maintainer=""
LABEL application.license="MIT"
LABEL application.from.image=""
LABEL application.scm=""

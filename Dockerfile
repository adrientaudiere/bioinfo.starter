# syntax=docker/dockerfile:1

FROM rocker/rstudio:4.3.2
LABEL maintainer="Adrien Taudière <adrientaudiere@zaclys.net>"

## Install system dependencies ----
RUN apt-get update \
    && apt-get install -y libz-dev \
      libglpk-dev \
      libudunits2-dev \
      libgdal-dev \
      libgeos-dev \
      libproj-dev \
      libbz2-dev \
      curl \
    && apt-get clean \
    && apt purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


## Install vsearch
RUN wget https://github.com/torognes/vsearch/releases/download/v2.25.0/vsearch-2.25.0-linux-x86_64.tar.gz \
   && tar xzf vsearch-2.25.0-linux-x86_64.tar.gz
ENV PATH=/vsearch-2.25.0-linux-x86_64/bin:$PATH

## Install MUMU
RUN git clone https://github.com/frederic-mahe/mumu.git \
  && cd ./mumu/ \
  && make && make check && make install \
  && mv mumu ~/mumu \
  && cd .. \
  && rm -r mumu
ENV PATH=~/mumu:$PATH

## Install miniconda to install with conda
ENV CONDA_PATH="/opt/miniforge3"

RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  && bash Miniforge3.sh -b -p "/opt/conda" \
  && rm -rf Miniforge3 \
  && source "opt/conda/etc/profile.d/conda.sh"


#RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh \
#  && bash /opt/miniconda.sh -b -p ${CONDA_PATH} \
#  && rm -rf /opt/miniconda3/miniconda.sh \
#  && ${CONDA_PATH}/bin/conda init bash
#ENV PATH=${CONDA_PATH}/bin:$PATH

## Parametrize conda and install softwares
RUN conda config --add channels bioconda \
  && conda config --add channels conda-forge  \
  && conda install -y -c bioconda itsxpress \
  && conda install -y -c bioconda blast \
  && conda install -y -c bioconda multiqc \
  && conda install -y -c bioconda swarm \
  && conda install -y -c bioconda fastqc \
  && conda install -y -c bioconda falco \
  && conda create -n cutadaptenv cutadapt \
  && conda init bash


RUN R -e "install.packages('pak', repos = c(CRAN = 'https://cloud.r-project.org'))" \
  && R -e  "pak::pkg_install(c('DT', 'emojifont', 'flexdashboard', 'forcats','formattable','ggimage', 'ggplot2', 'ggVennDiagram','grid', 'gridExtra','here', 'kableExtra', 'knitr', 'magick', 'networkD3', 'optparse', 'pbapply', 'plotly', 'qs', 'quarto', 'rmarkdown', 'sessioninfo', 'stringr','targets', 'tarchetypes', 'this.path', 'tinytex', 'treemap','vegan','viridis','visNetwork'))" \
  && R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

ENV FOLDER="/home/rstudio/"

## Set working directory ----
WORKDIR $FOLDER
ENV PROJECT="Bioinformatic analysis"

# Allowing only uncached the git clone starter_folder and update MiscMetabar
# using docker build --build-arg CACHE_DATE="$(date)"
# explanation : https://github.com/moby/moby/issues/22832
ARG CACHE_DATE=2016-01-01

# Install MiscMetabar
RUN R -e  "pak::pkg_install('adrientaudiere/MiscMetabar')" \
  && rm -rf /tmp/downloaded_packages/

## Clone the starter folder from gitlab ----
RUN git clone https://github.com/adrientaudiere/bioinfo.starter ${FOLDER}/$PROJECT

WORKDIR ${FOLDER}/$PROJECT

RUN chown -R rstudio:rstudio $FOLDER

VOLUME /app/data

EXPOSE 8787


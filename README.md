# bioinfo.starter

<!-- badges: start -->

<!-- badges: end -->

The goal of bioinfo.starter is to ...

## Installation

You can install the development version of bioinfo.starter like so:

``` r
devtools::install_github("adrientaudiere/bioinfo.starter")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bioinfo.starter)
## basic example code
```

## Colophon

bioinfo.starter is inspired by:
- {[rrtools]()} developed by Ben Marwick et al.
- {[rcompendium](https://frbcesab.github.io/rcompendium/)} developed by Nicolas Casajus.
- The [Turing way](https://the-turing-way.netlify.app/) recommendation

<!-- 

## Docker recipe

### For each new build

```sh
version_build=0.1.7
docker build -t adrienta/mycea_starter:$version_build -t adrienta/mycea_starter:latest --build-arg CACHE_DATE="$(date)" .

docker push adrienta/mycea_starter:$version_build --all-tags
```


#### Exemple avec projet XXXX

```sh
docker run --rm --env PROJECT="XXXX" -p 8787:8787 -e PASSWORD=221310 -e ROOT=TRUE -ti -v /media/adrien/homeMX3/kDrive/BIO_INFO/data_raw_mycea:/home/rstudio/data/data_raw:ro -v /home/adrien/Nextcloud/IdEst/Projets/BIOINFORMATIQUE/Mycea/Rendus_clients_RetD/ANALYSES_JANVIER_2024:/home/rstudio/data/  adrienta/mycea_starter:latest
```

```sh
mv Dossier_analyse ${PROJECT}
cd ${PROJECT}/

bash ${FOLDER}${PROJECT}/init.sh --sam_data ${FOLDER}/data/sam_data_pont_a_mousson.csv --name_folder $PROJECT  --raw_data ${FOLDER}/data/data_raw --min_reads_samp 0 --path_folder ${FOLDER} --path_params ${FOLDER}data/pont_a_mousson.yaml

Rscript run.R
```
-->


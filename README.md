# bioinfo.starter

<!-- badges: start -->

<!-- badges: end -->

The goal of bioinfo.starter is to [...]

## Installation

### Clone from github

```sh
git clone git@github.com:adrientaudiere/bioinfo.starter.git name_analyse
cd name_analyse
git checkout -b name_analyse
```

### Adapt to your pipeline

- Replace `data/data_raw/metadata/sam_data.csv` with good metadata file. 
  - Must be a true comma separated csv. If you prefer tabulation or ; you may want to add parameter to function `sam_data_matching_names()`.
  - The name of the file and the the column names indicating the samples are defined at the start of the `_targets.R` file
- Copy fastq files in `data/data_raw/rawseq`
- Add references database in `data/data_raw/refseq`
- Modify the `_targets.R` files (at least modify primers sequences and name of the reference database)
- Modify (if necessary) params `pattern_remove_sam_data` and `pattern_remove_fastq_files` to make matching fastq files and sample names in metadata

### Install R packages and run targets pipeline

```r
renv::install()
targets::tar_make()
```

## Colophon

bioinfo.starter is inspired by:
- {[rrtools]()} developed by Ben Marwick et al.
- {[lumo](https://github.com/holtzy/lumo)} developed by Yan Holtz.
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


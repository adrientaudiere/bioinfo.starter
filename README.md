# bioinfo.starter

<!-- badges: start -->

<!-- badges: end -->

The goal of bioinfo.starter is to [...]

```sh
git clone https://github.com/adrientaudiere/bioinfo.starter.git
git checkout -b name_analyse
```

- Replace `data/data_raw/metadata/sam_data.csv` with good metadata file
- Copy fastq files in `data/data_raw/rawseq`
- Add references database in `data/data_raw/refseq`
- Modify the `_targets.R` files (at least modify primers sequences and name of the reference database) 

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


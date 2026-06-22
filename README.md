# bioinfo.starter

## Install R packages and run targets pipeline

```r
install.packages("pak")
pak::local_install_deps(dependencies = TRUE)
fastqcr::fastqc_install()
```

```r
targets::tar_make()
```

## Multiple marker pipelines

This project includes three independent `targets` pipelines managed via `_targets.yaml`:

| Project | Script | Marker | Primers | Store |
|---------|--------|--------|---------|-------|
| `main`  | `_targets.R` | Archaea 16S (Arch349F / Arch806R) | `CCCTACGGGGTGCASCAG` / `GGACTACVSGGGTATCTAAT` | `_targets/` |
| `mcrA`  | `_targets_mcrA.R` | mcrA — methane production | `GGTGGTGTMGGDTTCACMCARTA` / `CGTTCATBGCGTAGTTVGGRTAGT` | `_targets_mcrA/` |
| `pmoA`  | `_targets_pmoA.R` | pmoA — methane oxidation | `GGNGACTGGGACTTCTGG` / `CCGGMGCAACGTCYTTACC` | `_targets_pmoA/` |

The mcrA and pmoA pipelines have no taxonomic assignment step (no reference database available for these markers).

To run a specific pipeline, set the `TAR_PROJECT` environment variable:

```r
Sys.setenv(TAR_PROJECT = "mcrA")
targets::tar_make()
```

Or without side effects:

```r
withr::with_envvar(c(TAR_PROJECT = "mcrA"), targets::tar_make())
```

All three pipelines use the same raw FASTQ files (`data/data_raw/rawseq/`) — the markers were pooled before sequencing. Cutadapt extracts marker-specific reads using the primer sequences. Intermediate and final outputs are written to marker-prefixed folders (e.g., `data/data_intermediate/mcrA_*/`, `data/data_final/quality_fastqc/pmoA_*/`) to avoid collisions between pipelines.


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


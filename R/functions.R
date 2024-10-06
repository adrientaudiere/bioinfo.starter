# Add some functions with documentation here

#fastq_dir <- system.file("extdata", package = "dada2")
#fastqc_agg(fastq_dir, qc.dir = "~/essai_fastqc")

fastqc_agg <- function(fq.dir = getwd(),
                       qc.dir = NULL,
                       threads = 4,
                       fastqc.path = "~/bin/FastQC/fastqc",
                       aggregate = TRUE) {
  fastqcr::fastqc(
    fq.dir = fq.dir,
    qc.dir = qc.dir,
    threads = threads,
    fastqc.path = fastqc.path
  )

  if (aggregate) {
    res <-  fastqcr::qc_aggregate(qc.dir)
  }

  message(paste("Fastqc reports are available in the folder", qc.dir, "."))
  return(res)
}



#fastq_dir <- system.file("extdata", package = "dada2")
#fastqc_agg(fastq_dir, qc.dir = "~/essai_fastqc")
#fastqc_plot("~/essai_fastqc")

fastqc_plot <- function(qc.dir,
                        modules = c(
                          "Basic Statistics",
                          "Per base sequence quality",
                          "Per base sequence content",
                          "Sequence Length Distribution",
                          "Adapter Content"
                        )) {
  p <- fastqcr::qc_plot_collection(qc_read_collection(
    list.files(qc.dir, pattern = "*.zip", full.names = TRUE),
    list.files(qc.dir, pattern = "*.zip"),
    modules = modules
  ),
  modules = modules)

  return(p)
}

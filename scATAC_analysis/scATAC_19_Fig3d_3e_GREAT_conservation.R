# Authors: Jack Huey

library(rGREAT)
library(rtracklayer)

ont=c("GO Molecular Function","GO Biological Process","GO Cellular Component",
      "Mouse Phenotype","Mouse Phenotype Single KO","Human Phenotype","Ensembl Genes")

do.ontology = function(filepath) {
  tmp = import.bed(filepath)
  if (length(tmp) == 0) {
    return()
  }
  if (file.exists(paste0(gsub(".bed", "", filepath), ".", gsub(" ", "_", ont[1]), ".csv"))) {
    return()
  }
  job = submitGreatJob(tmp,species='mm10')
  tb1 = getEnrichmentTables(job,ontology = ont)
  for (info in ont){
    write.csv(data.frame(tb1[info]), paste0(gsub(".bed", "", filepath), ".", gsub(" ", "_", info), ".csv"))
  }
}

mclapply(list.files("cluster_conservation_euar_top", pattern = ".bed", full.names = T), do.ontology, mc.cores = 2)

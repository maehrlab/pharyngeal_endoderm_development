# author: Jack Huey
library(ArchR)
library(ggrepel)
library(universalmotif)
library(TFBSTools)

## Read motifs

files <- list.files(path="Mus_musculus_2021_05_06_12_37_pm/pwms_all_motifs", full.names=TRUE)

motifs.raw = lapply(files, function(f) {
  id = gsub(".txt", "", basename(f))
  tryCatch({
    motif = universalmotif::read_cisbp(f)
    motif@name = id
    convert_motifs(motif, "TFBSTools-PWMatrix")
  }, error=function(e){
    #print(paste0(id, " ", name))
  })
})
motifs.raw <- motifs.raw[!sapply(motifs.raw,is.null)]
names(motifs.raw) = sapply(motifs.raw, function(f) { f@name })

all.motifs = list()
con = file("Mus_musculus_2021_05_06_12_37_pm/TF_Information.txt", "r")
readLines(con, n = 1)
while ( TRUE ) {
  line = readLines(con, n = 1)
  if (length(line) == 0) {
    break
  }
  cols = strsplit(line, "\t")
  id = cols[[1]][4]
  if (id == ".") {
    next
  }
  name = cols[[1]][7]
  motif = motifs.raw[[id]]
  if (is.null(motif)) {
    print(name)
    next
  }
  idx = length(grep(name, names(all.motifs))) + 1
  if (idx > 1) {
    name = paste0(name, "_", idx)
  }
  motif@name = name
  all.motifs[name] = motif
}
close(con)

#all.motifs["Foxn1_2"] = convert_motifs(read_meme("MA1684.1.meme"), "TFBSTools-PWMatrix") this is to add a custom foxn1 motif from meme, didnt use in E11_E12
motifsList = do.call(PWMatrixList, all.motifs)


## E11_E12
proj_e11_e12_cisbp2 = loadArchRProject("ML_0414_ABSOLUTELY_FINAL_object_whole_analysis")
proj_e11_e12_cisbp2 <- addMotifAnnotations(ArchRProj = proj_e11_e12_cisbp2, motifPWMs = motifsList, name = "MotifCisbp2_only")
proj_e11_e12_cisbp2 <- addDeviationsMatrix(ArchRProj = proj_e11_e12_cisbp2, peakAnnotation = "MotifCisbp2_only")
saveArchRProject(proj_e11_e12_cisbp2, "E11_E12_cisbp2_deviations")




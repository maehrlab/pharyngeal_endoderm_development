# author: Eric Kernfeld
results_path = Sys.getenv( "FREEZR_DESTINATION" )

# Getting these packages installed SUCKED! Here's what ultimately worked.
# install.packages("~/Dropbox/scRNA_data_analysis/pharynx_GRN/renv/library/R-3.5/x86_64-apple-darwin15.6.0/tgutil", repos = NULL, type = "source")
# install.packages("~/Dropbox/scRNA_data_analysis/pharynx_GRN/renv/library/R-3.5/x86_64-apple-darwin15.6.0/tgstat", repos = NULL, type = "source")
# install.packages(c("Rgraphviz", "pdist", "doMC", "RSvgDevice", "dbscan", "entropy"))
# install.packages("~/Downloads/metacell-445a346401873ba679ad7c19b632bcdd275e1996", repos = NULL, type = "source")
unloadNamespace("metacell")
library("metacell")

#' Partition cells into metacells. Don't call this directly; use record_metacell_assignments.
#'
#'
#' @param dge A Matrix storing expression values.
ComputeMetaCellAssignments = function( dge ){
  
  # Make temp directory
  host_folder = file.path(results_path, "temp" ,"metacell")
  if(!dir.exists(host_folder)) dir.create(host_folder, recursive = T, showWarnings = F)
  scdb_init(host_folder, force_reinit=T)
  scfigs_init(host_folder)
  # Import data
  dge %>% Write10X(data.dir = file.path( host_folder, "10x_v2"))
  mcell_import_scmat_10x(force = T, base_dir = file.path( host_folder, "10x_v2"), mat_nm = "temp")
  
  #ignore mito genes and low-depth cells
  mat = scdb_mat( "temp")
  nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
  rm(mat); gc()
  bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T)))
  mcell_mat_ignore_genes(new_mat_id="temp_filtered", mat_id="temp", bad_genes, reverse=F) 
  
  # Feature selection
  mcell_add_gene_stat(       gstat_id = "temp_gstats", mat_id="temp_filtered", force=T)
  mcell_gset_filter_varmean( gset_id  = "temp_features", gstat_id = "temp_gstats", T_vm=0.08, force_new=T)
  mcell_gset_filter_cov(     gset_id  = "temp_features", gstat_id = "temp_gstats", T_tot=100, T_top3=2)
  mcell_plot_gstats(         gset_id  = "temp_features", gstat_id = "temp_gstats"  )
  
  # KNN and metacell construction
  mcell_add_cgraph_from_mat_bknn(mat_id="temp_filtered", 
                                 gset_id = "temp_features", 
                                 graph_id="temp_graph",
                                 K=100,
                                 dsamp=T)
  mcell_coclust_from_graph_resamp(
    coc_id="temp_coc500", 
    graph_id="temp_graph",
    min_mc_size=20, 
    p_resamp=0.75, n_resamp=500)
  mcell_mc_from_coclust_balanced(
    coc_id="temp_coc500", 
    mat_id= "temp",
    mc_id= "temp_mc", 
    K=30, min_mc_size=30, alpha= 2)
  
  #report
  mc = scdb_mc("temp_mc")
  X = data.frame( METACELL_ID = mc@mc, CELL_ID = names(mc@mc), stringsAsFactors = F )
  unlink(list.files(host_folder, full.names = T), recursive = T)
  X
}
  
#' Aggregate expression using labels from the 'metacells_compute.R'
#'
#'
#'
ComputeMetaCellCounts = function( dge_mat, metadata, skip_none = F ){
  assertthat::assert_that( "METACELL_ID" %in% colnames(metadata) )
  assertthat::assert_that(all(!is.na(metadata[["METACELL_ID"]])))
  metacell_indicators = model.matrix(data = metadata, object = ~METACELL_ID + 0)
  colnames(metacell_indicators) %<>% gsub("^METACELL_ID", "", .)
  dge_mat_metacell = dge_mat %*% metacell_indicators
  assertthat::are_equal( ncol( dge_mat_metacell ), length( unique( metadata[["METACELL_ID"]] ) ) )
  assertthat::are_equal( nrow( dge_mat_metacell ), nrow( dge_mat_metacell ) )
  assertthat::are_equal( sum(dge_mat_metacell), sum(dge_mat))
  if(skip_none){
    idx_skip = colnames(dge_mat_metacell)=="NONE"
    if(any(idx_skip)){
      warning("Omitting a column NONE.\n")
      dge_mat_metacell = dge_mat_metacell[, !idx_skip]
    }
  }
  return(dge_mat_metacell)
}

#' Aggregate metadata using labels from the 'metacells_compute.R'
#'
#'
#'
ComputeMetaCellMetadata = function( metadata ){
  # Name, cell count
  metacell_metadata = 
    metadata[["METACELL_ID"]] %>% 
    table %>% 
    as.data.frame() %>% 
    set_colnames(c("METACELL_ID", "CELL_COUNT"))
  
  # for categorical variables, record the highest and its proportions.
  # for quantitative variables, take mean. (For nUMI, use sum as well.)
  get_prop_max = function(x) max(table(x) / sum(table(x)))
  get_mode = function(x) names(which.max(table(x)))
  aggregate_by_metacell = function(field, FUN){
    aggregate(x = metadata[[field]], 
              by = metadata["METACELL_ID"], 
              FUN = FUN)$x
  }
  metacell_metadata[["nUMI_total"]] = aggregate_by_metacell("nUMI", sum)
  is_quantitative = sapply(metadata, function(x) is.numeric(x) & !is.factor(x) )
  for( v in colnames(metadata)[is_quantitative] ){
    metacell_metadata[[v %>% paste0("_mean")]] = 
      aggregate_by_metacell(v, mean)
  }
  for( v in  colnames(metadata)[!is_quantitative] ){
    metacell_metadata[[v %>% paste0("_highest")]] = 
      aggregate_by_metacell(v, get_mode)
    metacell_metadata[[v %>% paste0("_prop_highest")]] = 
      aggregate_by_metacell(v, get_prop_max)
  }
  
  assertthat::assert_that( "METACELL_ID" %in% colnames(metacell_metadata) )
  assertthat::assert_that(all(!is.na(metadata[["METACELL_ID"]])))
  return(metacell_metadata)
}


dge = inventory_get( tag = "initial_exploration" ) %>% readRDS %>% add_cluster_names
metacell_assignments = list()
for( overview_cluster in dge@ident %>% levels %>% rev ){
  cat(overview_cluster, "\n")
  metacell_assignments[[overview_cluster]] = ComputeMetaCellAssignments(SubsetData(dge, ident.use = overview_cluster))
  metacell_assignments$METACELL_ID %<>% paste0(overview_cluster, "__", .)
}
metacell_assignments %<>% Reduce(f = rbind)
dge@meta.data[["METACELL_ID"]] = "NONE"
cu = metacell_assignments$CELL_ID
dge@meta.data[cu,"METACELL_ID"] = paste0( "pharynx_cluster_", dge@ident[cu], "__metacell_", metacell_assignments[cu, "METACELL_ID"] )
mcx = ComputeMetaCellCounts(dge@raw.data, dge@meta.data)
mcm = dge@meta.data %>% cbind(dge@dr$umap@cell.embeddings) %>% ComputeMetaCellMetadata
mcm %<>% set_rownames(mcm$METACELL_ID)
# All metacells are uniform in metacell ID (duh) and cluster (because of the loop above)
assertthat::assert_that( all(1==mcm$METACELL_ID_prop_highest))
assertthat::assert_that( all(1==mcm$res.1_prop_highest))
dim(mcm)
dim(mcx)
summary(mcm)
metacells_object = CreateSeuratObject(mcx) %>% AddMetaData(mcm)
metacells_object %<>% NormalizeData
save_feature_plots(metacells_object, 
                   results_path,
                   get_pharynx_genes(), 
                   axes = UMAP_AXES %>% paste0("_mean"),
                   axes_description = "metacell_mean_umap_coords")
inventory_save_and_add(metacells_object, "metacells_object", extra = "Seurat v2 object containing metacells (summed raw counts and mean or mode of metadata). Metacells were generated within each overview cluster.")
fp_10x = file.path(results_path, "metacells_object_10x")
metacells_object@raw.data %<>% as.matrix()
metacells_object@raw.data %<>% as("Matrix")
thymusatlastools2::Write10X(metacells_object, fp_10x, slot = "raw.data")
inventory_add(
  tag = fp_10x %>% basename, filename = fp_10x,
  extra = "10x triplet format containing metacells (summed raw counts and mean or mode of metadata). Metacells were generated within each overview cluster.")
)
# author: Eric Kernfeld

### One-time setup
options(repos = rev(BiocManager::repositories()))
packrat::set_opts(local.repos = c("~/Desktop/software_projects"))
packrat::init("~/Desktop/scRNA_data_analysis/pharynx_analysis3")
packrat::restore()
packrat::install_local("Seurat")
packrat::install_local("thymusatlastools2")
packrat::install_local("freezr")
packrat::snapshot()
packrat::set_opts(vcs.ignore.src = TRUE)
packrat::status()

### Daily setup
{
  proj_dir = path.expand("~/Desktop/scRNA_data_analysis/pharynx_analysis3")
  setwd( proj_dir )
  packrat::on()
  library(magrittr)
  library(freezr)
  library(Seurat)
  library(Matrix)
  library(thymusatlastools2)
  
  #' This is specific to my folder layout. 
  #' When run without args, this function runs only setup.Rmd and saves to the "interactive" folder.
  #' Otherwise, it runs both setup.Rmd and your analysis script, and it 
  #' names the results subfolder after the source code subfolder. This way, the same
  #' organizational structure is used for the analysis scripts and the results folders.
  #' 
  flash_freeze = configure_flash_freeze( project_directory = proj_dir, 
                                         setup_scripts = "setup.Rmd",
                                         repos_to_track =
                                           c( proj_dir,
                                              "~/Desktop/software_projects/thymusatlastools2",
                                              "~/Desktop/software_projects/freezr" ) )
  check_results = inventory_check(inv_location = file.path(proj_dir, "results"))
  check_results %>% subset(!exists, select = "tag", drop = T) 
  flash_freeze()
}


# 10x pharynx data
flash_freeze( analyses_to_run = "atlas/assemble.Rmd" )
flash_freeze( analyses_to_run = "atlas/subset_explore.Rmd" )
flash_freeze( analyses_to_run = "atlas/clean.Rmd" )
flash_freeze( analyses_to_run = "atlas/explore.Rmd" )
flash_freeze( analyses_to_run = "atlas/annotate_clusters.Rmd" )
flash_freeze( analyses_to_run = "atlas/subset_explore.Rmd" )
flash_freeze( analyses_to_run = "atlas/markers.Rmd" )
flash_freeze( analyses_to_run = "atlas/markers_terminal.Rmd" )

# Interface with Theis group's code
flash_freeze( analyses_to_run = "interface_with_pseudotime_pouch_subset_analysis/theis_group_export/mark_pouch_subsets.Rmd" )
flash_freeze( analyses_to_run = "interface_with_pseudotime_pouch_subset_analysis/theis_group_export/export_for_theis_group.Rmd" )
flash_freeze( analyses_to_run = "interface_with_pseudotime_pouch_subset_analysis/theis_group_export/hana_pouch3_gene_selection.Rmd" )
flash_freeze( analyses_to_run = "interface_with_pseudotime_pouch_subset_analysis/theis_group_export/plot_hana_results.Rmd" )

# GRN-related
flash_freeze( analyses_to_run = "atlas/export_for_GRN.R" )
flash_freeze( analyses_to_run = "atlas/make_metacells.R" )

# Foxn1 KO
flash_freeze( analyses_to_run = c(foxn1_tools, "foxn1_ko/foxn1_ko_initial.Rmd" ))
flash_freeze( analyses_to_run = c(foxn1_tools, "foxn1_ko/foxn1_ko_clean.Rmd" ))
flash_freeze( analyses_to_run = c(foxn1_tools, "foxn1_ko/foxn1_ko_visualize.Rmd" ))
flash_freeze( analyses_to_run = c(foxn1_tools, "foxn1_ko/foxn1_ko_isolate_via_classifier.Rmd" ))
flash_freeze( analyses_to_run = c(foxn1_tools, "foxn1_ko/foxn1_ko_isolate_via_clustering.Rmd" ))
flash_freeze( analyses_to_run = c(foxn1_tools, "foxn1_ko/foxn1_ko_isolate_with_pth.Rmd" ))
flash_freeze( analyses_to_run = c(foxn1_tools, "foxn1_ko/foxn1_ko_atlas_counterparts.Rmd" ))


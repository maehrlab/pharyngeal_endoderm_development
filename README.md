# pharyngeal_endoderm_development

This repository contains scripts for our manuscript titled "Integrated single-cell transcriptomic and chromatin accessibility approaches reveal regulatory programs driving pharyngeal organ development"

You can find the data at GEO: 

The repository is broadly organized as follows. Each directory requires a different set of packages as documented.
- scRNA_analysis

    a) atlas, foxn1_ko and interface_with_pseudotime_pouch_subset_analysis
       Contains code for Fig 1, 5, 6b and 6c
    
        Package requirements:
            - Seurat version 2 (from the Satija lab) 
            - thymusatlastools2 available at https://github.com/maehrlab/thymusatlastools2
            - https://github.com/ekernf01/freezr  to save code and session info and to track processed data for use in downstream scripts.
            - packrat (for dependency management)
            - magrittr
            - Matrix
    
       How to run:
       - in main.R edit the packrat(optional, can be deleted)
       - in main.R set the proj_dir to point to this repo
       - edit the flashfreeze paths to point to the location of the freezr and thymusatlastools2 packages
       - in setup.Rmd, edit the PATH_TO_DATA to point to the  CellRanger output directory (from GEO)
       
    b) pouch_object_pseudotime
    uses output of interface_with_pseudotime_pouch_subset_analysis
    Fig6a and FigS6; 
    output used in 6b and 6c
    
       Package requirements:
        - scanpy (https://github.com/theislab/scanpy)
    
        How to run:
        - self contained jupyter notebook
        
    c)  The standalone script FigS2_dotplot.R  is used to make the dotplot in supplementary figure 2. 
        
        Package requirements: Seurat version3, ggplot2 and cowplot
    
    d) notebooks "Fig6e_EnrichR_interface" and "Fig6e_dotplot_from_curated_terms"  are used to make the dotplots in Fig 6 
        
        Package requirements : enrichR, dplyr, magrittr, ggplot2, viridis, scales
        
- scATAC_qc_analysis
Contains code for Fig 2, 3, 5

        Package requirements:
        - ArchR (ArchR_1.0.1 in R 3.6.3)
        - MACS2 (part of ArchR suggested packages)
        - Cis-BP version 2 database (http://cisbp.ccbr.utoronto.ca)
        - Seurat - version 3
        - scanpy (https://github.com/theislab/scanpy)
        - rGREAT (https://www.bioconductor.org/packages/release/bioc/html/rGREAT.html)
        - universalmotif (https://bioconductor.org/packages/release/bioc/html/universalmotif.html)
        - TFBSTools (https://bioconductor.org/packages/release/bioc/html/TFBSTools.html)
        - ChromVar (part of ArchR suggested packages)
        - lineup (https://github.com/kbroman/lineup)
        - Phastcon score download (http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phastCons60way/mm10.60way.phastCons60wayEuarchontoGlire.bw)
        - enrichR (https://github.com/wjawaid/enrichR)
        - jupyter
        
        How to run:
        - series of self contained jupyter notebooks and R scripts

- gene_regulatory_networks
Contains code for Fig 4, Fig 5a and Fig 5e
NOTE: Code for creating the metacells is in the scRNA_analysis directory

    a) GENIE3 network (run from SCENIC package)
    
        Package requirements:
        - R version (R 3.6.3)
        - metacell_0.3.6
        - Seurat (v3)
        - Matrix
        - magrittr
        - parallel
        - RcisTarget
        - GENIE3
        - SCENIC_1.2.4 
        - igraph
        - network
        - sna
        - visNetwork
        - leiden
        - ggplot2
        NOTE: Heatmap in Fig 4e was made with scanpy, panda and anndata
      
      How to run:
          - series of self contained jupyter notebooks which also contain further details on packages used

        
    b) CellOracle GRN
        
        Package requirements:
        - CellOracle 0.6.6 and all the prescribed dependencies
        - scvelo 0.2.2.dev51+ga7de78a and all the prescribed dependencies
        NOTE: Fig 5e was made using R3.6.3, ggplot2 and Seurat
    
        How to run:
        - series of self contained jupyter notebooks which also contain further details on packages used


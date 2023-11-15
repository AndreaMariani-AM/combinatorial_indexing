# Function to load and convert ensembl genes to gene Symbol
load_split_seq <- function(fryDir, outputFormat, mixedExpr=FALSE){
  # This function load a alevinFry experiment and map ens2symbol genes
  # moreover makes it 1:1 mapping by filtering out genes not present in either
  # one of the annotation.
  require(fishpond)
  require(magrittr)
  require(stringr)
  library(SummarizedExperiment)

  #check fryDir is legit
  mtx_file <- file.path(fryDir, "alevin", "quants_mat.mtx")
  if(!file.exists(mtx_file)){
    stop("The 'fryDir' doesn't look like a directory generate by Alevin-fry:\n",
         sprintf("Quantification file is missing: %s", mtx_file)
         )
  }

  #Check annotation file exists
  annot_file <- file.path(fryDir, "bc_ex_mapping", "ensembl2symbol.txt")
  if(!file.exists(annot_file)){
    stop("The file 'ensembl2symbol' doesn't exist or it's not in the right directory:\n",
         sprintf("annotation file is missing: %s", annot_file))
  }

  # Load in experiment and clean up gene names. I can have more files in the counts directory.
  sce <- loadFry(fryDir = fryDir, outputFormat = outputFormat)
  rownames(sce) <- rownames(sce) %>%
    str_replace(pattern = "\\.\\d*", replacement = "") #remove transcript info

  # Clean up gene names
  rowData(sce)$gene_ids <- rowData(sce)$gene_ids %>%
    str_replace(pattern = "\\.\\d*", replacement = "")

  # Get the mapping table and filter it to be = to rownames(sce). This mapping file
  # will be supplied with the pipeline in the bc_ex_mapping dir
  gene_annotation_table <- data.table::fread(annot_file) %>%
    mutate(Geneid = str_replace(Geneid, pattern = "\\.\\d*", replacement = ""))


  gene_annotation_table <- gene_annotation_table %>%
    filter(Geneid %in% rownames(sce))

  # filter out genes that aren't present in reduced gene annotation table
  sce <- sce[rownames(sce) %in% gene_annotation_table$Geneid,]

  if(mixedExpr == FALSE){
    # Trim Mouse-## in front of the Gene Symbol
    # Map back Symbols as rownames and colData
    geneNames <- gene_annotation_table$GeneSymbol[match(rownames(sce),gene_annotation_table$Geneid)]
    geneNames <- stringr::str_remove(geneNames, "Mouse-")
    rownames(sce) <- geneNames
    rowData(sce)$SYMBOL <- geneNames

  } else {
    # Do not trim Gene Symbol name
    # Map back Symbols as rownames and colData
    geneNames <- gene_annotation_table$GeneSymbol[match(rownames(sce),gene_annotation_table$Geneid)]
    rownames(sce) <- geneNames
    rowData(sce)$SYMBOL <- geneNames
  }

  return(sce)
}


add_sample_info <- function(sce, fryDir, outDir, filename){
  # This function takes in the output of `load_split_seq` function, adds the samples info
  # and saves the output to a location
  require(fishpond)
  require(zellkonverter)

  #Check bc_sample_mapping file exists
  bc_sample_map <- file.path(fryDir, "bc_ex_mapping", "bc_sample_mapping.txt")
  if(!file.exists(bc_sample_map)){
    stop("The file 'bc_sample_mapping' doesn't exist or has not been generated yet:\n",
         sprintf("bc_sample_mapping doesn't exist: %s", bc_sample_map))
  }

  # read in sample info to add the information
  sample_df <- read.table(bc_sample_map, strip.white = TRUE, header = TRUE,
                          col.names = c("bc", "third_bc", "sample_info"),
                          sep = "\t")

  #check that colData(sce) and sample_df order is the same
  if(identical(sample_df$bc, colData(sce)$barcodes) != TRUE){
    stop("The order of barcodes in sce object and in the 'bc_sample_mapping' doesn't match")
  }

  colData(sce)$sample_info <- sample_df$sample_info


  # Deduplicate rows
  to_drop <- duplicated(rownames(sce))
  sce <- sce[!to_drop, ]

  # Save obj as h5ad that can be read in R or Python
  zellkonverter::writeH5AD(sce, paste0(outDir, "/", filename, ".h5ad"))
}



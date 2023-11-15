# Function to load and convert ensembl genes to gene Symbol
load_split_seq <- function(fryDir, outputFormat, expType){
  # This function load a alevinFry experiment and map ens2symbol genes
  # moreover makes it 1:1 mapping by filtering out genes not present in either
  # one of the annotation.
  require(fishpond)
  require(magrittr)
  require(stringr)
  require(SummarizedExperiment)

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
    stringr::str_replace(pattern = "\\.\\d*", replacement = "") #remove transcript info

  # Clean up gene names
  rowData(sce)$gene_ids <- rowData(sce)$gene_ids %>%
    stringr::str_replace(pattern = "\\.\\d*", replacement = "")

  # Get the mapping table and filter it to be = to rownames(sce). This mapping file
  # will be supplied with the pipeline in the bc_ex_mapping dir
  gene_annotation_table <- data.table::fread(annot_file) %>%
    dplyr::mutate(Geneid = stringr::str_replace(Geneid, pattern = "\\.\\d*", replacement = ""))


  gene_annotation_table <- gene_annotation_table %>%
    dplyr::filter(Geneid %in% rownames(sce))

  # filter out genes that aren't present in reduced gene annotation table
  sce <- sce[rownames(sce) %in% gene_annotation_table$Geneid,]

  # Check experiment type, valid values are: "Mixed, Mouse, Human"
  valid_list <- c("Mixed", "Mouse", "Human")
  if(!expType %in% valid_list){
    stop("Experiment type is not valid. Accepted values: Mixed, Mouse, Human")
  }

  if(expType == "Mixed"){
    # Trim Mouse-## in front of the Gene Symbol
    # Map back Symbols as rownames and colData, keep both annotations
    geneNames <- gene_annotation_table$GeneSymbol[match(rownames(sce),gene_annotation_table$Geneid)]
    geneNames <- stringr::str_remove(geneNames, "Mouse-")
    rownames(sce) <- geneNames
    rowData(sce)$SYMBOL <- geneNames

  } else if (expType == "Mouse") {
    # Select mouse
    # trim Gene Symbol name
    # Map back Symbols as rownames and colData
    geneNames <- gene_annotation_table$GeneSymbol[match(rownames(sce),gene_annotation_table$Geneid)]
    rownames(sce) <- geneNames
    geneNames <- stringr::str_subset(geneNames, "Mouse-")
    sce <- sce[geneNames,]
    geneNames_trimmed <- stringr::str_remove(geneNames, "Mouse-")
    rownames(sce) <- geneNames_trimmed
    rowData(sce)$SYMBOL <- geneNames_trimmed

  } else if (expType == "Human") {
    # Select Human
    # Map back Symbols as rownames and colData
    geneNames <- gene_annotation_table$GeneSymbol[match(rownames(sce),gene_annotation_table$Geneid)]
    rownames(sce) <- geneNames
    geneNames <- stringr::str_subset(geneNames, "Mouse-", negate = TRUE)
    sce <- sce[geneNames,]
    rowData(sce)$SYMBOL <- geneNames

  }

  print(sce)
  return(sce)
}


add_sample_info <- function(sce, fryDir, bc_filter){
  # This function takes in the output of `load_split_seq` function,
  # adds the samples info
  require(SummarizedExperiment)


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


  # Filter out others
  if(bc_filter == TRUE){
    #this means that bc annotated with "other" have been removed
    # need to align colData(sce)$barcodes, i'll operate on colnames(sce)
    drop_other <- colnames(sce) %in% sample_df$bc
    # this effectively drops BCs not in the sample_df_bc
    sce <- sce[,drop_other]

  }
  # Otherwise everything stays the same, as every BCs is kept
  # check that the order is the same
  if(identical(sample_df$bc, colData(sce)$barcodes) != TRUE){
    stop("The order of barcodes in sce object and in the 'bc_sample_mapping' doesn't match:\n",
         sprintf("have you filtered out barcodes annotated with 'other'?"))
  }

  colData(sce)$sample_info <- sample_df$sample_info


  # Deduplicate rows
  to_drop <- duplicated(rownames(sce))
  sce <- sce[!to_drop, ]

  return(sce)
}



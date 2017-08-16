library(ImmuneSpaceR)
library(Rlabkey)
library(data.table)
library(Biobase) # why are we loading library AND doing @importFrom?
library(tools) # why are we loading library AND doing @importFrom?

# read the job info
jobInfo <- read.table("${pipeline, taskInfo}",
                      col.names=c("name", "value", "type"),
                      header=FALSE,
                      check.names=FALSE,
                      stringsAsFactors=FALSE,
                      sep="\t",
                      quote="",
                      fill=TRUE,
                      na.strings="")

## Handle ALL inputs (useful for working locally on the code)
labkey.url.base <- jobInfo$value[jobInfo$name == "baseUrl"]
labkey.url.path <- jobInfo$value[jobInfo$name == "containerPath"]
pipeline.root   <- jobInfo$value[jobInfo$name == "pipeRoot"]
inputFiles      <- jobInfo$value[ grep("input\\.", jobInfo$name)]

selectedBiosamples <- selectedGEOs <- NULL
# From LABKEY.Pipeline.startAnalysis in views/CreateMatrix.html
selectedBiosamples <- "${selected-biosamples}"
selectedGEOs <- "${selected-GEOs}"

################################################
###             HELPER FN                    ###
################################################

setGEO <- function(GEOs, gef, con, inputFiles){
  # TEMPORARY list of studies where GEO should be avoided
  noGEO <- c("SDY224")
  
  # Decide whether we use GEO or the files
  if( all(GEOs %in% gef$geo_accession) & !con$study %in% noGEO ){
    return(TRUE)
  } else if( all(basename(inputFiles) %in% gef$file_info_name) ){
    return(FALSE)
  } else{
    stop("Could not decide between GEO and files. Check selected rows to make sure that all selected rows have a file or they all have a GEO accession number.")
  }
}

#' @importFrom preprocessCore normalize.quantiles
logNorm <- function(norm_exprs){
  # Norm process removes col / rownames
  rnames <- norm_exprs[ , feature_id ]
  norm_exprs[ , feature_id := NULL]
  cnames <- colnames(tmp)
  norm_exprs <- preprocessCore::normalize.quantiles(as.matrix(norm_exprs))
  colnames(norm_exprs) <- cnames
  
  # pmax ensures that no values are below 1, which would generate neg numbers when log2
  norm_exprs <- pmax(norm_exprs, 1)
  
  # if NOT rna-seq, then log2, handle NAs? Yes, for now.
  if( max(norm_exprs, na.rm = T) > 100 ){ norm_exprs <- log2(norm_exprs) }

  # need to do this last b/c need probes as col, not rownames and want dt output
  norm_exprs <- data.table(norm_exprs)
  norm_exprs[, feature_id := rnames ]
  
  return(norm_exprs)
}

.clean_colnames <- function(table){
  setnames(table, tolower(chartr(" ", "_", names(table))))
}

# Map colnames from expsample_accession to biosample_accession
.es2bs <- function(con, table){
  ess <- grep("^ES", colnames(table), value = TRUE)
  if(length(ess) > 0){
    esFilter <- makeFilter(c("expsample_accession", "IN", paste0(ess, collapse = ";")))
    bs2es <- data.table(labkey.selectRows(baseUrl = con$config$labkey.url.base,
                                          folderPath = con$config$labkey.url.path,
                                          schemaName = "immport", 
                                          queryName = "expsample_2_biosample",
                                          colFilter = esFilter,
                                          colNameOpt = "rname"))
    bss <- bs2es[match(ess, bs2es$expsample_accession), biosample_accession]
    setnames(table, ess, bss)
  }
  return(table)
}

################################################
###        EXT-SPECIFIC PROCESSING FN        ###
################################################
# Process GEO accession
# Returns a data.table with a feature_id column and one column per expsample_accession
#' @importFrom Biobase exprs sampleNames
.process_GEO <- function(gef){
  gef <- gef[geo_accession != ""] # Make sure we only have rows with GEO.
  gsm <- getGEO(gef[1, geo_accession])
  gse <- gsm@header$series_id[1] # Assumes seriesID is the first one (SuperSeries with higher ID)
  es <- getGEO(gse)[[1]]
  es <- es[, sampleNames(es) %in% gef$geo_accession]
  if( all(gef$geo_accession %in% sampleNames(es)) ){ # All selected samples are in the series
    #sampleNames are GEO accession
    sampleNames(es) <- gef[match(sampleNames(es), geo_accession), expsample_accession] 
    exprs <- data.table(exprs(es))
    exprs <- exprs[, feature_id := featureNames(es)]
    cnames <- colnames(exprs)[ colnames(exprs) != "feature_id"] # b/c dt is by ref!
    setcolorder(exprs, c("feature_id", cnames))
  } else{
    stop(paste0("Some of the selected gene_expression_files rows are not part of ",
                gse, ". Add code!"))
  }
  return(exprs)
}

# Process CEL files
# @param gef A \code{data.table} the gene_expression_files table or a subset of
#  it.
# @param inputFiles A \code{character}. The filenames.
# @return A \code{matrix} with biosample_accession as cols and feature_id as rownames
#' @importFrom affy ReadAffy rma
.process_CEL <- function(con, gef, inputFiles){
  affybatch <- ReadAffy(filenames = inputFiles)
  eset <- rma(affybatch)
  exprs <- exprs(eset)
  if (all(file_ext(colnames(exprs)) == "CEL")) { #If filenames used as samplenames
    colnames(exprs) <- gef[match(colnames(exprs), gef$file_info_name), biosample_accession]
  }
  return(exprs)
}

# This will work for files that follow the standards from immport
# Eventually, all tsv files should be rewritten to follow this standard.
# @return A matrix with biosample_accession as cols and feature_id as rownames
#' @importFrom reshape2 acast
.process_TSV <- function(gef, inputFiles){
  exprs <- fread(inputFiles, header = TRUE)
  if(min(exprs[, lapply(.SD, min), .SDcols = grep("^BS", names(exprs))]) == 0 &
     max(exprs[, lapply(.SD, max), .SDcols = grep("^BS", names(exprs))]) > 100){
    #RNA-Seq data. Assume it is already count base.
  } else{
    exprs <- .clean_colnames(exprs)
    if(!all(c("target_id", "raw_signal") %in% colnames(exprs))){
      stop("The file does not follow HIPC standards.")
    }
    try(setnames(exprs, "experiment_sample_accession", "expsample_accession"))
    try(setnames(exprs, "target_id", "feature_id"))
    if("expsample_accession" %in% colnames(exprs)){
      EorB <- "expsample_accession"
    } else if("biosample_accession" %in% colnames(exprs)){
      EorB <- "biosample_accession"
    } else{
      stop("The input file must contain either biosample_accession or expsample_accession")
    }
    exprs <- acast(exprs, formula = paste("feature_id ~", EorB), value.var = "raw_signal")
    exprs <- exprs[, c(colnames(exprs) %in% gef[[EorB]])]
  }
  return(exprs)
}

# biosample_accession as colnames instead of experiment sample
# Works for SDY212 & 162
.process_others <- function(gef, inputFiles){
  
  if( length(inputFiles) > 1 ){
    lf <- lapply(inputFiles, fread)
    names(lf) <- basename(inputFiles)
    exprs <- rbindlist(lf, idcol = TRUE)
    exprs <- .clean_colnames(exprs)
    try(setnames(exprs, "target_id", "probe_id"))
    try(setnames(exprs, "raw_signal", "expression_value"))
    file2sample <- unique(gef[, list(file_info_name, expsample_accession)])
    exprs[, sample := file2sample[match(.id, file_info_name), expsample_accession]]
    exprs <- exprs[, list(sample, probe_id, gene_symbol, expression_value)]
    
    exprs <- dcast.data.table(exprs, formula = "probe_id ~ sample", value.var="expression_value")
    rnames <- exprs[, probe_id]
    exprs <- exprs[, probe_id := NULL]
  } else {
    exprs <- fread(inputFiles, header = TRUE)
    sigcols <- grep("Signal", colnames(exprs), value = TRUE)
    rnames <- exprs[, PROBE_ID]
    if(length(sigcols) > 0){
      exprs <- exprs[, sigcols, with = FALSE]
      setnames(exprs, gsub(".AVG.*$", "", colnames(exprs)))
    } else if( all(gef$biosample_accession %in% colnames(exprs)) ){
      exprs <- exprs[, gef$biosample_accession, with = FALSE]
    } else{
      stop("Unknown format: check data and add code if needed.")
    }
  }
  norm_exprs <- norm_exprs[, c(colnames(norm_exprs) %in% gef$expsample_accession)]
  setnames(exprs, rnames, "feature_id")
  
  return(exprs)
}

################################################
###           HIGHER LEVEL FN                ###
################################################
#' makeMatrix
#'
#' Create a standard expression matrix from a given set of files or GEO
#' accession numbers.
#'
#' @param labkey.url.base string of test or production server url
#' @param labkey.url.path string of file path to sdy level
#' @param inputFiles vector of filepaths from sdy level to rawdata
#' @param selectedBiosamples vector of biosample ids to process
#' @param selectedGEOs vector of GEO accessions to process
#'
#' @name makeMatrix
#' @importFrom tools file_ext
#' @importFrom GEOquery getGEO
#' @export
#'
makeMatrix <- function(gef, isGEO, con){
  
  if( length( unique(gef$arm_name) ) > 1){
    message("There are more than one cohort selected in this HIPCMatrix run")
  }
  
  # Create expression matrix based on file extension or GEO
  if(isGEO){
    library(GEOquery)
    exprs <- .process_GEO(gef)
  }else{
    inputFiles <- unique(gef$file_info_name) # why not use the global via arg?
    ext <- unique(file_ext(inputFiles))
    inputFiles <- file.path("/share/files/",
                            con$config$labkey.url.path,
                            "@files/rawdata/gene_expression",
                            inputFiles)
    
    # Filetypes
    # After this step exprs is a matrix with features as rownames and expsample as colnames
    if(length(ext) > 1){
      stop(paste("There is more than one file extension:", paste(ext, collapse = ",")))
    } else if(ext == "CEL"){
      library(affy)
      exprs <- .process_CEL(con, gef, inputFiles)
    } else if(ext == "txt" | con$study %in% c("SDY162", "SDY180", "SDY212")){
      exprs <- .process_others(gef, inputFiles)
    } else if(ext == "tsv"){
      exprs <- .process_TSV(gef, inputFiles)
    } else{
      stop("File extension not supported.")
    }
    
    # want to return as data.table with feature_id instead of rownames?
    if( !is(exprs, "data.table") ){
      exprs <- data.table(exprs, keep.rownames = TRUE)
      setnames(exprs, "rn", "feature_id")
    }
  }
  
  # map experiment samplenames to biosample names
  exprs <- .es2bs(con, exprs)
}

summarizeMatrix <- function(norm_exprs, f2g){
  norm_exprs[ , gene_symbol := f2g[match(norm_exprs$feature_id, f2g$featureid), genesymbol] ]
  norm_exprs <- norm_exprs[ !is.na(gene_symbol) & gene_symbol != "NA"]
  norm_exprs <- norm_exprs[ ,lapply(.SD, mean),
                            by="gene_symbol",
                            .SDcols = grep("^BS", colnames(norm_exprs))]
}

writeMatrix <- function(pipeline.root, exprs, bygene){
  # - EM
  setnames(exprs, "feature_id", " ")
  write.table(exprs, 
              file = file.path(pipeline.root, 
                               "analysis/exprs_matrices", 
                               "${output.tsv}"), 
              sep = "\t", 
              quote=FALSE, 
              row.names=FALSE)
  
  # - summary EM
  write.table(bygene, 
              file = file.path(pipeline.root, 
                               "analysis/exprs_matrices", 
                               paste0("${output.tsv}", 
                                      ".summary")), 
              sep = "\t", 
              quote=FALSE, 
              row.names=FALSE)
  
  # - EM used for pipeline (not moved to the right location)
  write.table(exprs, 
              file = "${output.tsv}", 
              sep = "\t", 
              quote=FALSE, 
              row.names=FALSE)
}

#----------------------------------------------------------------
#                        PIPELINE
#----------------------------------------------------------------
co <- labkey.setCurlOptions(ssl.verifyhost = 2, sslversion=1)
FAS_filter <- makeFilter(c("FeatureAnnotationSetId/RowId", "IN", "${assay run property, featureSet}"))
f2g <- data.table(labkey.selectRows(baseUrl = labkey.url.base,
                                    folderPath = "/Studies/",
                                    schemaName = "Microarray",
                                    queryName = "FeatureAnnotation",
                                    colFilter = FAS_filter,
                                    colNameOpt = "rname", 
                                    colSelect = c("featureid", "genesymbol")))

if(nrow(f2g) == 0){ stop("The downloaded feature annotation set has 0 rows.") }

con <- CreateConnection() # sdy specific b/c labkey.url.path

# Subset gef according to inputFiles or GEO
filter <- makeFilter(c("biosample_accession", "IN", gsub(",", ";", selectedBiosamples)))
gef <- con$getDataset("gene_expression_files", 
                      colFilter = filter, 
                      original_view = TRUE, 
                      reload = TRUE)
inputFiles <- inputFiles[file.exists(inputFiles)]
GEOs <- unlist(strsplit(selectedGEOs, ","))
gef <- gef[file_info_name %in% basename(inputFiles) | geo_accession %in% GEOs]

# decide whether to use GEO
isGEO <- setGEO(GEOs, gef, con, inputFiles)

# want probe level to be non-normalized so user can have access to true rawdata
# output (exprs) is a data.table with "feature_id" as col, not rownames
exprs <- makeMatrix(gef, isGEO, con) # > dt output

# normalize (if not RNA-seq) and summary by gene 
norm_exprs <- copy(exprs) # need to keep exprs safe from data.table by-ref changes!
if( unique(file_ext(inputFiles)) != "CEL" ){ norm_exprs <- logNorm(norm_exprs) } # dt > matrix > dt
bygene <- summarizeMatrix(norm_exprs, f2g) # dt expected and output

writeMatrix(pipeline.root, exprs, bygene)

outProps = file(description="${pipeline, taskOutputParams}", open="w")
cat(file=outProps, sep="", "name\tvalue\n")
cat(file=outProps, sep="", "assay run property, cohort\t", unique(gef$cohort), "\n")
flush(con=outProps)
close(con=outProps)
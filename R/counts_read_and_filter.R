# Functions to rad and filter counts
#
# Author: vidal
###############################################################################

utils::globalVariables("expression.raw")

#' Function to read counts as produced by htseq-count
#' @param expr.file \code{character} or \code{list}. String or vector or list of input file paths.
#' Allowed are individual text files with counts for one sample each, with gene IDs in the first and counts in the second column or
#' a single counts matrix file containing read counts for all samples with rows corresponding to genes (genomic features) and columns to samples.
#' Negative values or NAs are not allowed and gene IDs are expected in the first column.
#' If \code{miRSEQ=TRUE} this expects the output from CAP-miRSEQ summary script which is also a counts matrix.
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet.
#' @param miRSEQ \code{logical}. Is the input data the output from CAP-miRSEQ summary script?
#' @export
diff_expr_read_counts <-
  function(expr.file, samp.info, miRSEQ = FALSE)
  {
    # expr.file <- as.character(expr.file)
    if (length(unlist(expr.file)) > 1) {
      ## check sample names
      cat("Checking expression file names...\n")
      if (is.null(names(expr.file))) {
        names(expr.file) <- sub("\\.[a-z]{1,5}$", "", basename(expr.file))
      }
      if (length(grep("\\.txt$|\\.count$", names(expr.file)))) {
        names(expr.file) <- sub("\\.[a-z]{1,5}$", "", names(expr.file))
      }
      expr.file <- expr.file[order(names(expr.file))]
      if (any(names(expr.file) %in% samp.info$SampleNames)) {
        snam <- which(names(expr.file) %in% samp.info$SampleNames)
        if (!identical(names(expr.file)[snam], as.character(samp.info$SampleNames))) {
          stop("Input vector of count file names must be identical to all or a subset of sample names provided in 'samp.info', i.e., it must be the same names in the same order!")
        } else {
          if (nrow(samp.info) < length(expr.file)) {
            cat("  Importing subset of samples...\n")
          }
          expr.file <- expr.file[as.character(samp.info$SampleNames)]
        }
      } else {
        stop("Sample names not found in input files")
      }
      # import individual count files
      cat("Reading count files for samples:\n")
      print(samp.info$SampleNames)
      counts <- edgeR::readDGE(expr.file)$counts
    } else if (miRSEQ) {
      # reading output from CAP-miRSEQ summary script
      cat("Reading output from CAP-miRSEQ summary script...\n")
      counts <- read.table(expr.file, sep="\t", header=T, stringsAsFactors=F)
      # assign unique mature_precursor id to each row
      row.names(counts) <- expression.raw$Mature.miRNA
      # extract samples used in this comparison
      samples <- make.names(samp.info$SampleNames)
      counts <- counts[, samples]
    } else {
      # reading counts matrix with rows corresponding to genes and columns to samples
      cat("Reading counts matrix...\n")
      counts <- read.delim(expr.file, row.names = 1)
      samples <- make.names(samp.info$SampleNames)
      counts <- counts[, samples]
    }
    # in case of untypical gene identifiers change those
    ## for "gene:" in front of the Ensembl Gene ID
    rownames(counts) <- sub("^gene:", "", rownames(counts))
    return(counts)
  }

#' Function to filter counts
#' @param counts Count matrix.
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet.
#' @param strict Logical; only keep a miRNA if there are > 5 reads per million in at least half of the samples?
#' @param min.samp Integer; Number of samples in which a feature needs to be covered by at least one read per million.
#'   Defaults to the size of the smallest group of replicates. See \emph{details}.
#' @details
#' In edgeR, it is recommended to remove features without at least 1 read per million in n of the
#'   samples, where n is the size of the smallest group of replicates (determined from the 'groups' vector).
#'
#' @export
diff_expr_filter_counts <-
  function(counts, samp.info, strict=TRUE, min.samp=NULL)
  {
    #### Important: CPMs must be calculated with all counts including the non-aligned ones so as not to distort the results
    #### because: CPM <- apply(countmatrix, 2, function(x) (x/sum(x))*1000000) # uses sum of ALL counts in that sample
    noint_names <- grep("no_feature|ambiguous|too_low_aQual|_lowaqual|_empty|not_aligned|alignment_not_unique", rownames(counts), value=TRUE)
    noint <- rownames(counts) %in% noint_names
    ## In edgeR, it is recommended to remove features without at least 1 read per million in n of the
    ## samples, where n is the size of the smallest group of replicates (determined from the 'groups' vector)
    ### can be overridden by 'min.samp'
    ### if we want to be really strict we only keep a miRNA if there are > 5 reads per million in at least half of the samples
    # TODO: add edgeR filter function as in McElreavey project
    cpms <- edgeR::cpm(counts)
    if (strict) {
      cat("  Using 'strict' filtering...\n")
      keep <- rowSums(cpms > 5) >= ceiling(ncol(counts)/2) & !noint
    } else {
      if(is.null(min.samp)) {
        min.samp <- min(table(samp.info$Groups))
      }
      keep <- rowSums(cpms > 1) >= min.samp & !noint
    }
    cat("  Removing", length(which(!keep)), "and keeping", length(which(keep)), "rows...\n")
    counts <- counts[keep, ]
    ## Visualize and inspect the count table
    colnames(counts) <- basename(colnames(counts))
    return(counts)
  }



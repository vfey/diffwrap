# Functions to read and filter counts
#
# Author: vidal
###############################################################################

#' Function to read counts as produced by htseq-count
#' @param expr.dat \code{character} or \code{list}. String or vector or list of input file paths, or matrix of count values.
#' Allowed are individual text files with counts for one sample each, with gene IDs in the first and counts in the second column or
#' a single counts matrix file containing read counts for all samples with rows corresponding to genes (genomic features) and columns to samples.
#' Negative values or NAs are not allowed and gene IDs are expected in the first column.
#' If \code{miRSEQ=TRUE} this expects the output from CAP-miRSEQ summary script which is also a counts matrix.
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet.
#' @param miRSEQ \code{logical}. Is the input data the output from CAP-miRSEQ summary script?
#' @return A \code{matrix} of raw counts with features in rows and samples in columns, restricted to the
#'   samples listed in \option{samp.info} and ordered as they are there.
#' @examples
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_read_counts(diffwrap_counts, si)
#' dim(counts)
#' @export
diff_expr_read_counts <-
  function(expr.dat, samp.info, miRSEQ = FALSE)
  {
    # order samp.info data frame by SampleNames column
    samp.info <- samp.info[order(samp.info$SampleNames), ]

    if (is.character(expr.dat)) {
      if (length(unlist(expr.dat)) > 1) {
        ## check sample names
        dw_log("Checking expression file names...\n")
        if (is.null(names(expr.dat))) {
          names(expr.dat) <- sub("\\.[a-z]{1,5}$", "", basename(expr.dat))
        }
        if (length(grep("\\.txt$|\\.count$", names(expr.dat)))) {
          names(expr.dat) <- sub("\\.[a-z]{1,5}$", "", names(expr.dat))
        }
        # order input file names by vector names
        expr.dat <- expr.dat[order(names(expr.dat))]
        if (any(names(expr.dat) %in% samp.info$SampleNames)) {
          snam <- which(names(expr.dat) %in% samp.info$SampleNames)
          if (!identical(names(expr.dat)[snam], as.character(samp.info$SampleNames))) {
            stop("Input vector of count file names must be identical to all or a subset of sample names provided in 'samp.info', i.e., it must be the same names in the same order!")
          } else {
            if (nrow(samp.info) < length(expr.dat)) {
              dw_log("  Importing subset of samples...\n")
            }
            expr.dat <- expr.dat[as.character(samp.info$SampleNames)]
          }
        } else {
          stop("Sample names not found in input files")
        }
        # import individual count files
        dw_log("Reading count files for samples:\n")
        dw_log_obj(samp.info$SampleNames)
        counts <- edgeR::readDGE(expr.dat)$counts
      } else if (miRSEQ) {
        # reading output from CAP-miRSEQ summary script
        dw_log("Reading output from CAP-miRSEQ summary script...\n")
        counts <- read.table(expr.dat, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        # assign the mature miRNA id to each row
        row.names(counts) <- counts$Mature.miRNA
        # extract samples used in this comparison
        samples <- make.names(samp.info$SampleNames)
        counts <- counts[, samples]
      } else {
        # reading counts matrix with rows corresponding to genes and columns to samples
        dw_log("Reading counts matrix...\n")
        counts <- read.delim(expr.dat, row.names = 1, check.names = FALSE)
        # here potentially get columns that are not sample names (future functionality)
        if (any(names(counts) %in% samp.info$SampleNames)) {
          snam <- which(names(counts) %in% samp.info$SampleNames)
          if (!identical(names(counts)[snam], as.character(samp.info$SampleNames))) {
            stop("Input vector of count file names must be identical to all or a subset of sample names provided in 'samp.info', i.e., it must be the same names in the same order!")
          } else {
            if (nrow(samp.info) < ncol(counts)) {
              dw_log("  Importing subset of samples...\n")
            }
            dw_log("Reading count files for samples:\n")
            dw_log_obj(samp.info$SampleNames)
            counts <- counts[, as.character(samp.info$SampleNames)]
          }
        } else {
          stop("Sample names not found in input files")
        }
        # import individual count files
      }
    } else if (is.matrix(expr.dat) || is.data.frame(expr.dat)) {
      # count matrix must have feature IDs (e.g., gene symbols) in a specified column or as row names
      if (any(colnames(expr.dat) %in% samp.info$SampleNames)) {
        snam <- which(colnames(expr.dat) %in% samp.info$SampleNames)
        if (!identical(colnames(expr.dat)[snam], as.character(samp.info$SampleNames))) {
          stop("Column names of input matrix must be identical to all or a subset of sample names provided in 'samp.info', i.e., it must be the same names in the same order!")
        } else {
          if (nrow(samp.info) < ncol(expr.dat)) {
            dw_log("  Importing subset of samples...\n")
          }
          expr.dat <- expr.dat[, as.character(samp.info$SampleNames)]
        }
      } else {
        stop("Sample names not found in input column names")
      }
      dw_log("Geting counts from count matrix...\n")
      counts <- edgeR::getCounts(edgeR::DGEList(expr.dat))
    }
    # in case of untypical gene identifiers change those
    ## for "gene:" in front of the Ensembl Gene ID
    rownames(counts) <- sub("^gene:", "", rownames(counts))
    ## Ensure a consistent return type across all input modes.
    ## 'readDGE()' and 'getCounts()' hand back a matrix, whereas 'read.delim()' and
    ## 'read.table()' hand back a data.frame, so without this coercion the class of the
    ## returned object would depend on how the counts happened to be supplied.
    if (!is.matrix(counts)) {
      not.num <- !vapply(counts, is.numeric, logical(1))
      if (any(not.num)) {
        stop("Non-numeric column(s) in the count data: ",
             paste(sQuote(names(counts)[not.num]), collapse = ", "),
             ". Counts must be numeric.", call. = FALSE)
      }
      counts <- as.matrix(counts)
    }
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
#' @return A \code{matrix} of counts with weakly expressed features and the non-informative
#'   \command{htseq-count} summary rows removed.
#' @examples
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_read_counts(diffwrap_counts, si)
#' filtered <- diff_expr_filter_counts(counts, si, strict = TRUE)
#' c(before = nrow(counts), after = nrow(filtered))
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
      dw_log("  Using 'strict' filtering...\n")
      keep <- rowSums(cpms > 5) >= ceiling(ncol(counts)/2) & !noint
    } else {
      if(is.null(min.samp)) {
        min.samp <- min(table(samp.info$Groups))
      }
      keep <- rowSums(cpms > 1) >= min.samp & !noint
    }
    dw_log("  Removing", length(which(!keep)), "and keeping", length(which(keep)), "rows...\n")
    counts <- counts[keep, ]
    ## Visualize and inspect the count table
    colnames(counts) <- basename(colnames(counts))
    return(counts)
  }



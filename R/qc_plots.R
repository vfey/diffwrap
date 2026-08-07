# All QC plots of the diffr pipeline
#
# Author: vidal
###############################################################################

utils::globalVariables(c("PCx", "PCy", "Condition", "Sample", "Ellipse", "x", "y", "Block", "Groups", "SampleName"))

#' Main wrapper function for QC plots
#' @param counts Counts matrix.
#' @param samp.info \code{data.frame}. samp.info object containing information of the project's sample sheet.
#' @param control \code{character}. Name of the control group
#' @param out.l \code{list}. Object returned by \code{diffExpr()} containing all output components from the analysis.
#' @param grp.nam Legend title.
#' @param PC Integer vector of length two or three specifying the principal components to be plotted.
#' @param sample.plot.names Passed to \code{plotMDS()} (\code{labels}): character vector of sample names or labels. Defaults to colnames(d).
#' @param ellipse Logical indicating whether to draw an ellipse around the sample groups.
#' @param ellipse.mapping.groups Optional additional grouping for ellipse drawing. Overrides sample groups.
#'   Use for selective highlighting of user-defined sample groups.
#' @param ellipse.grp.nam Not implemented.
#' @param label.samples Logical; should the points be annotated with sample name labels?
#' @param geom.point.size Numeric passed to \code{geom_point} giving the point size.
#' @param label.font.size Numeric passed to \code{geom_text_repel} giving the label font size.
#' @param plot.ellipse.legend Logical; should the ellipse legend be plotted. \code{NA}, the default, will plot it if any aesthetics are mapped.
#' @param circle \code{logical}. Draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param varname.size \code{double}. size of the text for variable names
#' @param var.axes \code{logical}. draw arrows for the variables?
#' @param pairs Factor of identifiers specifying paired samples for paired or other block designs, or batch effects.
#' @param pairs.name Legend title for paired or block design variables in MDS ggplot.
#' @param gene.selection Character passed to \code{plotMDS()} specifying the mode to select genes for comparisons.
#' @param n Passed to \code{plotMDS()} (\code{top}): Integer; number of top genes used to calculate pairwise distances.
#' @param type Character; one of "both", "uncorrected", "pseudo-corrected" describing which values should be plotted.
#'   "uncorrected" will plot the input counts matrix while "pseudo-corrected" will plot pseudo counts for
#'   blocked designs (e.g., paired samples or batch factors).
#' @param analysis.name Character used in the plot title and the output file name if the plot is saved to a PDF
#' @param out.dir Character; path where to save PDF. Required; no default is used so that nothing is written to the working directory unintentionally.
#' @return The input list \option{out.l}, with a \code{QCplots} element added (or extended) holding the
#'   generated quality control plots as \code{ggplot} objects, named after the plot type and the value of
#'   \option{type}.
#' @examples
#' \donttest{
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' lcpm <- edgeR::cpm(counts, log = TRUE)
#' out.l <- diff_expr_QC_plots(counts = lcpm, samp.info = si, control = "control",
#'                             out.l = list(), grp.nam = "Group",
#'                             sample.plot.names = colnames(lcpm),
#'                             analysis.name = "demo", out.dir = tempdir())
#' names(out.l$QCplots)
#' }
#' @export
diff_expr_QC_plots <-
  function(counts, samp.info, control, out.l, grp.nam=NULL, PC=c(1,2,3), sample.plot.names=NULL,
           ellipse=TRUE, ellipse.mapping.groups=NULL, ellipse.grp.nam=NULL, label.samples=TRUE,
           geom.point.size=2, label.font.size = 5, plot.ellipse.legend=NA, circle=TRUE, varname.size=0, var.axes=FALSE,
           pairs=NULL, pairs.name=NULL, gene.selection="common", n=500, type=NULL, analysis.name=NULL, out.dir)
  {
    dw_log("  Doing PCA...\n")
    PCA <- diff_expr_PCA(counts=counts, n=n)
    groups <- relevel(samp.info$Groups, ref=control)
    samp.name <- samp.info$SampleNames
    main <- paste(type, analysis.name, sep="_")
    main <- gsub("\\.{1,}", "_", make.names(main))
    if (!length(out.l$QCplots)) out.l$QCplots <- list()
    dw_log("  Plotting...\n")
    pdf_file <- file.path(out.dir, paste0(Sys.Date(), "_", main, "_QC_plots.pdf"))
    dw_log("   Saving plot to", pdf_file, "...\n")
    pdf(pdf_file, width=11, height=11)
    qc.dev <- grDevices::dev.cur()
    main <- paste(type, analysis.name)
    dw_log("   MDS ggplot...\n")
    out.l$QCplots[[paste0(type, "_MDS")]] <- diff_expr_ggplot_mds(counts=counts, samp.name=samp.name, groups=groups, grp.nam=grp.nam, pairs=pairs,
                                                   pairs.name=pairs.name, gene.selection=gene.selection, dim.plot=PC[1:2], main=main)
    dw_log("   PCA ggbiplot...\n")
    out.l$QCplots[[paste0(type, "_PCAbiplot")]] <- diff_expr_PCA_ggbiplot(PCA=PCA, groups=groups, grp.nam=grp.nam, ellipse=ellipse, circle=circle,
                                                           varname.size=varname.size, var.axes=var.axes, main=main)
    dw_log("   PCA ggplot...\n")
    out.l$QCplots[[paste0(type, "_PCAlabelledPlot")]] <- diff_expr_PCA_ggplot(PCA=PCA, samp.name=sample.plot.names, groups=groups, grp.nam=grp.nam, PC=c(1,2),
                                                               main=main, ellipse=ellipse, ellipse.mapping.groups=ellipse.mapping.groups, ellipse.grp.nam=ellipse.grp.nam, label.samples=label.samples,
                                                               geom.point.size=geom.point.size, label.font.size = label.font.size, plot.ellipse.legend=plot.ellipse.legend)
    dw_log("   PCA 3d scatterplot...\n")
    diff_expr_3d_scatterplot(PCA=PCA, samp.name=sample.plot.names, groups=groups, grp.nam=grp.nam, PC=PC[1:3], main=main)
    dw_log("   PCA cluster dendrogram...\n")
    diff_expr_dendro_plot(counts=counts, groups=groups, grp.nam=grp.nam, main=main, col.grps=TRUE)
    dw_dev_off(qc.dev)
    dw_log("  Plotting finished.\n")
    return(out.l)
  }

#' Wrapper around `limma::plotMDS` to generate a MDS plot
#' @param d Passed to \code{plotMDS()} (\code{x}): Any data object that can be coerced to a matrix of log-expression values,
#'   for example an ExpressionSet or an EList. Rows represent genes or genomic features while columns represent samples.
#' @param groups Factor of sample groups for colouring and legend.
#' @param n Passed to \code{plotMDS()} (\code{top}): Integer; number of top genes used to calculate pairwise distances.
#' @param sample.plot.names Passed to \code{plotMDS()} (\code{labels}): character vector of sample names or labels. Defaults to colnames(d).
#' @param analysis.name Character used in the plot title and the output file name if the plot is saved to a PDF
#' @param do.pdf Logical indicating whether a PDF should be produced.
#' @param out.dir Character; path where to save PDF. Required; no default is used so that nothing is written to the working directory unintentionally.
#' @seealso [plotMDS()]
#' @return No return value. Called for its side effect of drawing a multidimensional scaling plot,
#'   optionally into a PDF file below \option{out.dir}.
#' @examples
#' \donttest{
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' d <- edgeR::calcNormFactors(edgeR::DGEList(counts, group = groups))
#' diff_expr_mds_plot(d, groups = groups, do.pdf = FALSE, out.dir = tempdir())
#' }
#' @export
diff_expr_mds_plot <-
  function(d, groups, n=500, sample.plot.names=NULL, analysis.name=NULL, do.pdf=FALSE, out.dir)
  {
    if (!is.null(sample.plot.names)) {
      dw_log("  *** Using custom sample labels: ***\n  ", head(sample.plot.names), "\n")
    }
    mds.dev <- NULL
    if (do.pdf) {
      pdf(file.path(out.dir, paste0(analysis.name, "_MDS_plot.pdf")), width=11, height=11)
      mds.dev <- grDevices::dev.cur()
    } else {
      # only when drawing to the caller's device do we change (and must restore) its par;
      # capturing par before a pdf() would force a device open and disturb the plot output
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)
    }
    par(mar = c(6,6,5,3))
    limma::plotMDS(d, top=n, labels=sample.plot.names, main=paste0("MDS plot for '", analysis.name, "' normalised DGEList"), col = rainbow(length(levels(groups)))[factor(groups)])
    legend("bottomright", legend=levels(groups), pch=15, col=rainbow(length(levels(groups))))
    if (do.pdf) {
      dw_dev_off(mds.dev)
    }
  }

#' Function to generate a MDS plot using `ggplot2`
#' @param counts Counts matrix.
#' @param samp.name Optional sample names to be used in the plot, given as character vector.
#' @param groups Sample groups for plot annotation as character vector or factor.
#' @param grp.nam Legend title.
#' @param pairs Factor of identifiers specifying paired samples for paired or other block designs, or batch effects.
#' @param pairs.name Legend title for paired or block design variables in MDS ggplot.
#' @param gene.selection Character passed to \code{plotMDS()} specifying the mode to select genes for comparisons.
#' @param dim.plot Integer vector of length two passed to \code{plotMDS()} specifying the principal components to be plotted.
#' @param main Plot title.
#' @seealso [plotMDS()]
#' @return A \code{ggplot} object containing the multidimensional scaling plot.
#' @examples
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' lcpm <- edgeR::cpm(counts, log = TRUE)
#' g <- diff_expr_ggplot_mds(lcpm, samp.name = colnames(lcpm), groups = groups)
#' class(g)
#' @export
diff_expr_ggplot_mds <-
  function(counts, samp.name, groups, grp.nam=NULL, pairs=NULL, pairs.name=NULL, gene.selection="common", dim.plot=c(1,2), main=NULL)
  {
    dw_log("    Preparing data...", "\n")
    mds <- limma::plotMDS(counts, gene.selection=gene.selection, dim.plot=dim.plot, plot=FALSE)
    dat <- data.frame(SampleName=samp.name, Groups=groups, x=mds$x, y=mds$y)
    dw_log("    Plotting...", "\n")
    g <- ggplot(dat, aes(x, y, colour=Groups))
    if (!is.null(pairs)) {
      dat <- data.frame(SampleName=samp.name, Groups=groups, Block=pairs, x=mds$x, y=mds$y)
      if (is.factor(dat$Block) && length(levels(dat$Block)) <= 6) {
        dw_log("\n     Block design using a discrete variable. Shaping points by blocking variable...\n")
        g <- ggplot(dat, aes(x, y, colour=Groups, shape=Block))
        if (!is.null(pairs.name)) {
          g <- g + scale_shape_discrete(name=pairs.name)
        }
      } else if (is.numeric(dat$Block)) {
        dw_log("\n     Block design using a continuous variable. Sizing points by blocking variable...\n")
        g <- ggplot(dat, aes(x, y, colour=Groups, size=Block))
        if (!is.null(pairs.name)) {
          g <- g + scale_size_continuous(name=pairs.name)
        }
      }
    }
    g <- g + geom_point()
    g <- g + ggrepel::geom_text_repel(
      data = dat,
      aes(label = SampleName),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      show.legend = FALSE)
    if (!is.null(grp.nam)) {
      g <- g + scale_color_discrete(name=grp.nam)
    }
    g <- g + ggtitle(paste0("Multi-dimensional scaling (", main, ")"))
    print(g)
    return(g)
  }

## PCA plots
#' Function to generate a PCA biplot using `ggbiplot.n`, a version of `ggbiplot` from https://github.com/vqv/ggbiplot.
#' @param PCA List of class \code{prcomp} or just a list having a component \code{x} that contains the rotated variables from \code{prcomp}.
#' @param groups Sample groups for plot annotation as character vector or factor.
#' @param grp.nam Legend title.
#' @param ellipse Logical indicating whether to draw an ellipse around the sample groups.
#' @param circle \code{logical}. Draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param varname.size \code{double}. size of the text for variable names
#' @param var.axes \code{logical}. draw arrows for the variables?
#' @param main Plot title.
#' @param fix.aspect \code{logical}. Should the aspect ratio of the x- and y-axes be kept constant for different plot sizes?
#' @param tweak \code{logical}. Should the plot theme be tweaked? Will apply values set in axes.title.size, legend.text.size and legend.title.size. Defaults to TRUE.
#' @param ... Arguments passed to \code{ggbiplot.n()}.
#' @seealso [ggbiplot.n()]
#' @return A \code{ggplot} object containing the PCA biplot. The plot is returned rather than drawn, so
#'   it has to be printed to appear on a device.
#' @examples
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' pca <- diff_expr_PCA(edgeR::cpm(counts, log = TRUE), n = 100)
#' g <- diff_expr_PCA_ggbiplot(pca, groups = groups)
#' class(g)
#' @export
diff_expr_PCA_ggbiplot <-
  function(PCA, groups, grp.nam=NULL, ellipse=TRUE, circle=TRUE, varname.size=0, var.axes=FALSE, main=NULL, fix.aspect=FALSE, tweak=FALSE, ...)
  {
    main <- paste("Biplot for PCA (", main, ")")
    dw_log("    Plotting...", "\n")
    g <- ggbiplot.n(PCA, var.scale = 1, obs.scale = 1, groups = groups, grp.nam=grp.nam, ellipse = ellipse, circle = circle, varname.size = varname.size, var.axes = var.axes, main=main, fix.aspect=fix.aspect, tweak=tweak, ...)
    print(g)
    return(g)
  }

#' Function to generate an ordinary two-dimensional PCA plot using `ggplot2`
#' @param PCA List of class \code{prcomp} or just a list having a component \code{x} that contains the rotated variables from \code{prcomp}.
#' @param samp.name Optional sample names to be used in the plot, given as character vector.
#' @param groups Sample groups for plot annotation as character vector or factor.
#' @param grp.nam Legend title.
#' @param PC Integer vector of length two specifying the principal components to be plotted.
#' @param main Plot title.
#' @param ellipse Logical indicating whether to draw an ellipse around the sample groups.
#' @param ellipse.mapping.groups Optional additional grouping for ellipse drawing. Overrides sample groups.
#'   Use for selective highlighting of user-defined sample groups.
#' @param ellipse.grp.nam Not implemented.
#' @param label.samples Logical; should the points be annotated with sample name labels?
#' @param geom.point.size Numeric passed to \code{geom_point} giving the point size.
#' @param label.font.size Numeric passed to \code{geom_text_repel} giving the label font size.
#' @param plot.ellipse.legend Logical; should the ellipse legend be plotted. \code{NA}, the default, will plot it if any aesthetics are mapped.
#' @param do.plot Logical; should the plot be printed to the graphics device? Defaults to \code{TRUE}.
#' @seealso [stat_ellipse()]
#' @return A \code{ggplot} object containing the labelled PCA scatterplot.
#' @examples
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' pca <- diff_expr_PCA(edgeR::cpm(counts, log = TRUE), n = 100)
#' g <- diff_expr_PCA_ggplot(pca, samp.name = NULL, groups = groups, do.plot = FALSE)
#' class(g)
#' @export
diff_expr_PCA_ggplot <-
  function(PCA, samp.name=NULL, groups, grp.nam=NULL, PC=c(1,2), main=NULL, ellipse = TRUE, ellipse.mapping.groups=NULL, ellipse.grp.nam=NULL,
           label.samples = TRUE, geom.point.size = 2, label.font.size = 5, plot.ellipse.legend=NA, do.plot = TRUE)
  {

    dw_log("    Preparing data...", "\n")
    if (is.null(samp.name)) {
      samp.n <- rownames(PCA$x)
    } else if (length(samp.name)==1 && is.na(samp.name)) {
      samp.n <- ""
    } else if (identical(rownames(PCA$x), names(samp.name))) {
      samp.n <- samp.name
    } else {
      samp.n <- NULL
      warning("If not a named character vector of the same length and with identical names as 'rownames(PCA$x)' 'samp.name' can only be 'NULL' to use the row names of the principal components matrix or 'NA' to be empty.")
    }

    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2), 1)
    dataGG <- data.frame(PCx=PCA$x[, PC[1]], PCy=PCA$x[, PC[2]], Condition=groups, Sample = { if (!is.null(samp.n)) { samp.n } else { "" }})

    if (is.null(grp.nam)) {
      grp.nam <- "Group"
    }
    if (is.null(ellipse.grp.nam)) {
      eg <- "Ellipse"
    }

    #Additional grouping for ellipses if provided
    if (ellipse && !is.null(ellipse.mapping.groups)) {
      dataGG$Ellipse <- ellipse.mapping.groups
      eg <- ellipse.grp.nam
    } else {
      dataGG$Ellipse <- groups
      eg <- grp.nam
    }
    ## determine minimum group size to generate feedback on potentially missing ellipses
    grp.size <- plyr::daply(dataGG, "Ellipse", nrow)

    dw_log("    Plotting...\n")
    g <- ggplot(data=dataGG, aes(PCx, PCy, color=Condition))
    g <- g + geom_point(size = geom.point.size)
    g <- g + ggtitle(paste0("PCA (", main, ")"))
    g <- g + labs(x=paste0("PC", PC[1], ": ", round(percentVar[PC[1]], 4), "% variance explained"),
                  y=paste0("PC", PC[2], ": ", round(percentVar[PC[2]], 4), "% variance explained"))
    g <- g + scale_colour_brewer(name=grp.nam, type="qual", palette=3, direction=1) # allows a maximum of 12 colours
    g <- g + theme_bw(base_size = 10) # if more than 10 groups yellow will be used which is not easily readable on white bg, so change bg to grey
    if (length(levels(dataGG$Condition))+length(levels(dataGG$Ellipse)) > 10) {
      g <- g + theme(panel.background = element_rect(fill = "#f4f4f4"))
    }

    g <- g + theme(panel.border = element_blank(),
                   axis.line = element_line(color='black'),
                   panel.grid.major = element_line(linewidth = 0.2),
                   panel.grid.minor = element_line(linewidth = 0.2))

    if ( ellipse ) {
      if (min(grp.size<4)) {
        dw_log("    # NOTE: Too few samples in one or more groups. Some or all ellipses may not be plotted.\n")
      }
      dw_log("     Adding ellipse...\n     Using", sQuote(eg), "for ellipse mapping\n")
      if (is.null(ellipse.mapping.groups)) {
        dw_log("      Plotting ellipses around main sample groups...\n")
        g <- g + stat_ellipse(type="t", show.legend=plot.ellipse.legend) #assumes a multivariate t-distribution
      } else {
        dw_log("      Plotting ellipses around second factor groups...\n")
        g <- g + stat_ellipse(mapping = aes(PCx, PCy, linetype=Ellipse), type = "t", inherit.aes = FALSE, show.legend=plot.ellipse.legend)
      }
    }

    if (label.samples){
      if (!is.null(samp.n)) {
        g <- g + ggrepel::geom_text_repel(aes(label=Sample),
                                          data=dataGG,
                                          size = label.font.size,
                                          box.padding = unit(0.35, "lines"),
                                          point.padding = unit(0.3, "lines"),
                                          show.legend = FALSE)
      }
    }
    if (do.plot) {
      print(g)
    }
    return(g)
  }

#' Function to generate a 3D scatterplot
#' @param PCA List of class \code{prcomp} or just a list having a component \code{x} that contains the rotated variables from \code{prcomp}.
#' @param samp.name Optional sample names to be used in the plot, given as character vector.
#' @param groups Sample groups for plot annotation as character vector or factor.
#' @param grp.nam Legend title.
#' @param PC Integer vector of length three specifying the principal components to be plotted.
#' @param main Plot title.
#' @return No return value. Called for its side effect of drawing a three-dimensional PCA scatterplot on
#'   the current graphics device.
#' @examples
#' \donttest{
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' pca <- diff_expr_PCA(edgeR::cpm(counts, log = TRUE), n = 100)
#' diff_expr_3d_scatterplot(pca, groups = groups)
#' }
#' @export
diff_expr_3d_scatterplot <-
  function(PCA, samp.name=NULL, groups, grp.nam=NULL, PC=c(1,2,3), main=NULL)
  {
    if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
      stop("Package ", sQuote("scatterplot3d"),
           " must be installed to draw the 3D PCA scatterplot.", call. = FALSE)
    }
    if (ncol(PCA$x) < max(PC[1:3])) {
      dw_log("    # NOTE: fewer than 3 principal components available; skipping 3D scatterplot.\n")
      return(invisible(NULL))
    }
    dw_log("    Preparing data...", "\n")
    if (is.null(samp.name)) {
      samp.n <- rownames(PCA$x)
    } else if (length(samp.name)==1 && is.na(samp.name)) {
      samp.n <- ""
    } else if (identical(rownames(PCA$x), names(samp.name))) {
      samp.n <- samp.name
    } else {
      samp.n <- NULL
      warning("If not a named character vector of the same length and with identical names as 'rownames(PCA$x)' 'samp.name' can only be 'NULL' to use the row names of the principal components matrix or 'NA' to be empty.")
    }

    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2), 1)
    dat <- data.frame(PCx=PCA$x[, PC[1]], PCy=PCA$x[, PC[2]], PCz=PCA$x[, PC[3]], Condition=groups)
    cond <- levels(dat$Condition)
    col_seq <- rev(rep(c(seq.int(from=2, to=12, by=2), seq.int(from=1, to=11, by=2)), ceiling(length(cond)/12))[seq_len(length(cond))])
    cols <- RColorBrewer::brewer.pal(n=12, name="Paired")[col_seq]
    plcol <- rep(cols[1], nrow(dat))
    names(plcol) <- dat$Condition
    if (is.null(grp.nam)) {
      grp.nam <- "Condition"
    }
    for (i in 2:length(cond)) {
      plcol[dat$Condition == cond[i]] <- cols[i]
    }
    dw_log("    Plotting...\n     scatterplot...\n")
    s3d <- scatterplot3d::scatterplot3d(dat[,1:3],        # x y and z axis
                                        color=plcol, pch=19,        # circle color indicates no. of cylinders
                                        type="h", lty.hplot=2,       # lines to the horizontal plane
                                        scale.y=.75,                 # scale y axis (reduce by 25%)
                                        main=paste0("3-D Scatterplot for PCA (", main, ")"))
    dw_log("     point labels...\n")
    s3d.coords <- s3d$xyz.convert(dat[,1:3])
    text(s3d.coords$x, s3d.coords$y,     # x and y coordinates
         labels=samp.n,       # text to plot
         pos=4, cex=.5)                  # shrink text 50% and place to right of points)
    # add the legend
    dw_log("     legend...\n")
    legend("topleft", inset=.05,      # location and inset
           bty="n", cex=.5,              # suppress legend box, shrink text 50%
           title=grp.nam,
           legend=unique(names(plcol)), fill=unique(plcol))
  }

#' Function to generate dendrogram plots based on hierarchical clustering
#' @param counts Counts matrix.
#' @param groups Sample groups for plot annotation as character vector or factor.
#' @param grp.nam Legend title.
#' @param main Plot title.
#' @param col.grps Logical indicating whether to colour the dendrogram by groups.
#' @return No return value. Called for its side effect of drawing a hierarchical clustering dendrogram on
#'   the current graphics device.
#' @examples
#' \donttest{
#' si <- diff_expr_get_samp_info(diffwrap_samp_info, "SampleName", "Group")
#' counts <- diff_expr_filter_counts(diff_expr_read_counts(diffwrap_counts, si), si)
#' groups <- stats::relevel(si$Groups, ref = "control")
#' diff_expr_dendro_plot(edgeR::cpm(counts, log = TRUE), groups = groups)
#' }
#' @export
diff_expr_dendro_plot <-
    function(counts, groups, grp.nam=NULL, main=NULL, col.grps=FALSE)
  {

    # test if needed packages are installed
    if (col.grps && !requireNamespace("dendextend", quietly = TRUE)) {
      stop(
        paste("Package", sQuote("dendextend"), "must be installed to colour the dendrogram by groups."),
        call. = FALSE
      )
    }

    dw_log("    Hierarchical clustering...", "\n")
    hc <- hclust(dist(t(counts)))
    plot(hc, main=paste("Hierarchical Clustering (", main, ")"))
    dw_log("    Generating coloured dendrogram...", "\n")
    dend <- as.dendrogram(hc)
    main <- paste("Hierarchical Clustering (", main, ")")
    if (col.grps) {
      dendextend::labels_colors(dend) <- as.numeric(groups)[hc$order]
      main <- paste0("Hierarchical Clustering (", main, ")\n[coloured by ", grp.nam, "]")
    }
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar=c(8,6,6,4))
    plot(dend, main=main)
  }



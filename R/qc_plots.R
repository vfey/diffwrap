# All QC plots of the diffr pipeline
# 
# Author: vidal
###############################################################################


#' Main wrapper function for QC plots
#' @export
diff_expr_QC_plots <-
		function(counts, samp.info, control, out.l, grp.nam=NULL, PC=c(1,2,3), sample.plot.names=NULL,
				ellipse=TRUE, ellipse.mapping.groups=NULL, ellipse.grp.nam=NULL, label.samples=TRUE,
				geom.point.size=2, label.font.size = 5, plot.ellipse.legend=NA, circle=TRUE, varname.size=0, var.axes=FALSE,
				pairs=NULL, pairs.name=NULL, gene.selection="common", n=500, type=NULL, analysis.name=NULL, out.dir=".")
{
	cat("  Doing PCA...\n")
	PCA <- diff_expr_PCA(counts=counts, n=n)
	cat("  done\n")
	groups <- relevel(samp.info$Groups, ref=control)
	samp.name <- samp.info$SampleNames
	main <- paste(type, analysis.name, sep="_")
	main <- gsub("\\.{1,}", "_", make.names(main))
	out.l$QCplots <- list()
	cat("  Plotting...\n")
	pdf_file <- file.path(out.dir, paste0(Sys.Date(), "_", main, "_QC_plots.pdf"))
	cat("   Saving plot to", pdf_file, "...\n")
	pdf(pdf_file, width=11, height=11)
	main <- paste(type, analysis.name)
	cat("   MDS ggplot...\n")
	out.l$QCplots[["MDS"]] <- diff_expr_ggplot_mds(counts=counts, samp.name=samp.name, groups=groups, grp.nam=grp.nam, pairs=pairs,
			pairs.name=pairs.name, gene.selection=gene.selection, dim.plot=PC[1:2], main=main)
	cat("   done.\n   PCA ggbiplot...\n")
	out.l$QCplots[["PCAbiplot"]] <- diff_expr_PCA_ggbiplot(PCA=PCA, groups=groups, grp.nam=grp.nam, ellipse=ellipse, circle=circle,
			varname.size=varname.size, var.axes=var.axes, main=main)
	cat("   done.\n   PCA ggplot...\n")
	out.l$QCplots[["PCAlabelledPlot"]] <- diff_expr_PCA_ggplot(PCA=PCA, samp.name=sample.plot.names, groups=groups, grp.nam=grp.nam, PC=c(1,2),
			main=main, ellipse=ellipse, ellipse.mapping.groups=ellipse.mapping.groups, ellipse.grp.nam=ellipse.grp.nam, label.samples=label.samples,
			geom.point.size=geom.point.size, label.font.size = label.font.size, plot.ellipse.legend=plot.ellipse.legend)
	cat("   done.\n   PCA 3d scatterplot...\n")
	diff_expr_3d_scatterplot(PCA=PCA, samp.name=sample.plot.names, groups=groups, grp.nam=grp.nam, PC=PC[1:3], main=main)
	cat("   done.\n   PCA cluster dendrogram...\n")
	diff_expr_dendro_plot(counts=counts, groups=groups, grp.nam=grp.nam, main=main)
	cat("   done.\n")
	dev.off()
	cat("  Plotting finished.\n")
	return(out.l)
}

#' Function to generate a MDS plot using `limma::plotMDS`
#' @export
diff_expr_mds_plot <-
		function(d, groups, n=500, sample.plot.names=NULL, analysis.name=NULL, do.pdf=FALSE, out.dir=".")
{
	if (!is.null(sample.plot.names)) {
		cat("  *** Using custom sample labels: ***\n  ", head(sample.plot.names), "\n")
	}
	if (do.pdf) {
		pdf(file.path(out.dir, paste0(analysis.name, "_MDS_plot.pdf")), width=11, height=11)
	}
	par(mar = c(6,6,5,3))
	plotMDS(d, top=n, labels=sample.plot.names, main=paste0("MDS plot for '", analysis.name, "' normalised DGEList"), col = rainbow(length(levels(groups)))[factor(groups)])
	legend("bottomright", legend=levels(groups), pch=15, col=rainbow(length(levels(groups))))
	if (do.pdf) {
		dev.off()
	}
}

#' Function to generate a MDS plot usig `ggplot2`
#' @export
diff_expr_ggplot_mds <-
		function(counts, samp.name, groups, grp.nam=NULL, pairs=NULL, pairs.name=NULL, gene.selection="common", dim.plot=c(1,2), main=NULL)
{
	cat("    Preparing data...")
	mds <- plotMDS(counts, gene.selection=gene.selection, dim.plot=dim.plot, plot=FALSE)
	dat <- data.frame(SampleName=samp.name, Groups=groups, x=mds$x, y=mds$y)
	cat("done\n    Plotting...")
	g <- ggplot(dat, aes(x, y, colour=Groups))
	if (!is.null(pairs)) {
		dat <- data.frame(SampleName=samp.name, Groups=groups, Block=pairs, x=mds$x, y=mds$y)
		if (is.factor(dat$Block) && length(levels(dat$Block)) <= 6) {
			cat("\n     Block design using a discrete variable. Shaping points by blocking variable...\n")
			g <- ggplot(dat, aes(x, y, colour=Groups, shape=Block))
			if (!is.null(pairs.name)) {
				g <- g + scale_shape_discrete(name=pairs.name)
			}
		} else if (is.numeric(dat$Block)) {
			cat("\n     Block design using a continuous variable. Sizing points by blocking variable...\n")
			g <- ggplot(dat, aes(x, y, colour=Groups, size=Block))
			if (!is.null(pairs.name)) {
				g <- g + scale_size_continuous(name=pairs.name)
			}
		}
	}
	g <- g + geom_point()
	g <- g + geom_text_repel(
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
	cat("    done\n")
	return(g)
}

## PCA plots
#' Function to generate a PCA biplot using `medseqr::ggbiplot`
#' @export
diff_expr_PCA_ggbiplot <-
		function(PCA, groups, grp.nam=NULL, ellipse=TRUE, circle=TRUE, varname.size=0, var.axes=FALSE, main=NULL, fix.aspect=FALSE, tweak=FALSE, ...)
{
	main <- paste("Biplot for PCA (", main, ")")
	cat("    Plotting...")
	g <- medseqr::ggbiplot.n(PCA, var.scale = 1, obs.scale = 1, groups = groups, grp.nam=grp.nam, ellipse = ellipse, circle = circle, varname.size = varname.size, var.axes = var.axes, main=main, fix.aspect=fix.aspect, tweak=tweak, ...)
	print(g)
	cat("done\n")
	return(g)
}

#' Function to generate an ordinary two-dimensional PCA plot using `ggplot2`
#' @export
diff_expr_PCA_ggplot <-
		function(PCA, samp.name=NULL, groups, grp.nam=NULL, PC=c(1,2), main=NULL, ellipse = TRUE, ellipse.mapping.groups=NULL, ellipse.grp.nam=NULL,
				label.samples = TRUE, geom.point.size = 2, label.font.size = 5, plot.ellipse.legend=NA)
{
	
	cat("    Preparing data...")
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
	grp.size <- plyr::daply(dataGG, .(Ellipse), nrow)
	
	cat("done\n    Plotting...\n")
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
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.2))
	
	if ( ellipse ) {
		if (min(grp.size<4)) {
			cat("    # NOTE: Too few samples in one or more groups. Some or all ellipses may not be plotted.\n")
		}
		cat("     Adding ellipse...\n     Using", sQuote(eg), "for ellipse mapping\n")
		if (is.null(ellipse.mapping.groups)) {
			cat("      Plotting ellipses around main sample groups...\n")
			g <- g + stat_ellipse(type="t", show.legend=plot.ellipse.legend) #assumes a multivariate t-distribution
		} else {
			cat("      Plotting ellipses around second factor groups...\n")
			g <- g + stat_ellipse(mapping = aes(PCx, PCy, linetype=Ellipse), type = "t", inherit.aes = F, show.legend=plot.ellipse.legend)
		}
	}
	
	if (label.samples){
		if (!is.null(samp.n)) {
			g <- g + geom_text_repel(aes(label=Sample),
					data=dataGG,
					size = label.font.size,
					box.padding = unit(0.35, "lines"),
					point.padding = unit(0.3, "lines"),
					show.legend = F)
		}
	}
	print(g)
	cat("done\n")
	return(g)
}

#' Function to generate a 3D scatterplot
#' @export
diff_expr_3d_scatterplot <-
		function(PCA, samp.name=NULL, groups, grp.nam=NULL, PC=c(1,2,3), main=NULL)
{
	cat("    Preparing data...")
	if (is.null(samp.name)) {
		samp.n <- rownames(PCA$x)
	} else if (is.na(samp.name)) {
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
	col_seq <- rev(rep(c(seq.int(from=2, to=12, by=2), seq.int(from=1, to=11, by=2)), ceiling(length(cond)/12))[1:length(cond)])
	cols <- brewer.pal(n=12, name="Paired")[col_seq]
	plcol <- rep(cols[1], nrow(dat))
	names(plcol) <- dat$Condition
	if (is.null(grp.nam)) {
		grp.nam <- "Condition"
	}
	for (i in 2:length(cond)) {
		plcol[dat$Condition == cond[i]] <- cols[i]
	}
	cat("done\n    Plotting...\n     scatterplot...\n")
	s3d <- scatterplot3d(dat[,1:3],        # x y and z axis
			color=plcol, pch=19,        # circle color indicates no. of cylinders
			type="h", lty.hplot=2,       # lines to the horizontal plane
			scale.y=.75,                 # scale y axis (reduce by 25%)
			main=paste0("3-D Scatterplot for PCA (", main, ")"))
	cat("     point labels...\n")
	s3d.coords <- s3d$xyz.convert(dat[,1:3])
	text(s3d.coords$x, s3d.coords$y,     # x and y coordinates
			labels=samp.n,       # text to plot
			pos=4, cex=.5)                  # shrink text 50% and place to right of points)
# add the legend
	cat("     legend...\n")
	legend("topleft", inset=.05,      # location and inset
			bty="n", cex=.5,              # suppress legend box, shrink text 50%
			title=grp.nam,
			legend=unique(names(plcol)), fill=unique(plcol))
	cat("    done\n")
}

#' Function to generate dendrogram plots based on hierarchical clustering
#' @export
diff_expr_dendro_plot <-
		function(counts, groups, grp.nam=NULL, main=NULL)
{
	cat("    Hierarchical clustering...")
	hc <- hclust(dist(t(counts)))
	plot(hc, main=paste("Hierarchical Clustering (", main, ")"))
	cat("done\n    Generating coloured dendrogram...")
	dend <- as.dendrogram(hc)
	labels_colors(dend) <- as.numeric(groups)[hc$order]
	par(mar=c(8,6,6,4))
	plot(dend, main=paste0("Hierarchical Clustering (", main, ")\n[coloured by ", grp.nam, "]"))
	cat("done\n")
}



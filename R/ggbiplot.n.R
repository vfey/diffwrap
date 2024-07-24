# Modified version of ggbiplot from https://github.com/vqv/ggbiplot to allow for more
# user control over graphical parameters such as font size or line width.
#
# Author: vidal
###############################################################################

utils::globalVariables(c("xvar", "yvar", "muted", "varname", "angle", "hjust"))

#' Make a biplot of PCA output data using ggplot2.
#'
#' @param pcobj an object returned by prcomp() or princomp()
#' @param choices \code{length-two numeric}. which PCs to plot; default: 1:2
#' @param scale \code{length-one numeric}. covariance biplot (scale = 1) [default], form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param obs.scale \code{length-one numeric}. scale factor to apply to observations; default: 1-scale
#' @param var.scale \code{length-one numeric}. scale factor to apply to variables; default: scale
#' @param pc.biplot \code{logical}. for compatibility with biplot.princomp()
#' @param groups \code{factor}. optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param grp.nam \code{character}. optional name of the grouping variable to be used as legend title
#' @param ellipse \code{logical}. draw a normal data ellipse for each group?
#' @param ellipse.prob \code{double}. size of the ellipse in Normal probability
#' @param labels \code{character}. optional vector of labels for the observations
#' @param labels.size \code{double}. size of the text used for the labels
#' @param alpha \code{double}. alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param circle \code{logical}. Draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param circle.prob \code{double}. size of the circle in Normal probability
#' @param var.axes \code{logical}. draw arrows for the variables?
#' @param varname.size \code{double}. size of the text for variable names
#' @param varname.adjust \code{double}. adjustment factor the placement of the variable names, >= 1 means farther from the arrow
#' @param varname.abbrev \code{logical}. whether or not to abbreviate the variable names
#' @param point.size \code{double}. expansion factor for point size; uses \code{rel()} internally
#' @param axes.title.size \code{double}. expansion factor for axes title sizes; uses \code{rel()} internally
#' @param legend.text.size \code{double}. expansion factor for legend text size; uses \code{rel()} internally
#' @param legend.title.size \code{double}. expansion factor for legend title size; uses \code{rel()} internally
#' @param ellipse.lwd \code{double}. expansion factor for ellipse line width; uses \code{rel()} internally
#' @param main \code{character}. Plot title. If \code{NULL} (default) no title will be added.
#' @param fix.aspect \code{logical}. Should the aspect ratio of the x- and y-axes be kept constant for different plot sizes?
#' @param tweak \code{logical}. Should the plot theme be tweaked? Will apply values set in axes.title.size, legend.text.size and legend.title.size. Defaults to TRUE.
#' @param tidy \code{logical}. If TRUE, \code{theme_minimal} will be applied.
#' @param ... currently not in use
#' @return The final plot object returned by ggplot.
#' @export
ggbiplot.n <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
		obs.scale = 1 - scale, var.scale = scale, groups = NULL, grp.nam=NULL,
		ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3,
		alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69,
		varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE,
		point.size = 1, axes.title.size = 1, legend.text.size = 1, legend.title.size = 1,
		ellipse.lwd=1, main=NULL, fix.aspect=TRUE, tweak=TRUE, tidy=TRUE, ...)
{
	stopifnot(length(choices) == 2)
	if (inherits(pcobj, "prcomp")) {
		nobs.factor <- sqrt(nrow(pcobj$x) - 1)
		d <- pcobj$sdev
		u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
		v <- pcobj$rotation
	}
	else if (inherits(pcobj, "princomp")) {
		nobs.factor <- sqrt(pcobj$n.obs)
		d <- pcobj$sdev
		u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
		v <- pcobj$loadings
	}
	else if (inherits(pcobj, "PCA")) {
		nobs.factor <- sqrt(nrow(pcobj$call$X))
		d <- unlist(sqrt(pcobj$eig)[1])
		u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
		v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord),
								1]), FUN = "/")
	}
	else if (inherits(pcobj, "lda")) {
		nobs.factor <- sqrt(pcobj$N)
		d <- pcobj$svd
		u <- predict(pcobj)$x/nobs.factor
		v <- pcobj$scaling
		d.total <- sum(d^2)
	}
	else {
		stop("Expected a object of class prcomp, princomp, PCA, or lda")
	}
	choices <- pmin(choices, ncol(u))
	df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale,
					FUN = "*"))
	v <- sweep(v, 2, d^var.scale, FUN = "*")
	df.v <- as.data.frame(v[, choices])
	names(df.u) <- c("xvar", "yvar")
	names(df.v) <- names(df.u)
	if (pc.biplot) {
		df.u <- df.u * nobs.factor
	}
	r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
	v.scale <- rowSums(v^2)
	df.v <- r * df.v/sqrt(max(v.scale))
	if (obs.scale == 0) {
		u.axis.labs <- paste("standardized PC", choices, sep = "")
	}
	else {
		u.axis.labs <- paste("PC", choices, sep = "")
	}
	u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)",
					100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
	if (!is.null(labels)) {
		df.u$labels <- labels
	}
	if (!is.null(groups)) {
		df.u$groups <- groups
	}
	if (varname.abbrev) {
		df.v$varname <- abbreviate(rownames(v))
	}
	else {
		df.v$varname <- rownames(v)
	}
	df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
	df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
	g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) +
			ylab(u.axis.labs[2])
	if (fix.aspect) {
		cat("\n     Note: Fixing aspect ratio of axes.\n    ")
		g <- g + coord_equal()
	}
	if (tweak) {
		cat(" Note: Tweaking size of text elements.\n    ")
		g <- g + theme(axis.title.x=element_text(size=rel(axes.title.size)), axis.title.y=element_text(size=rel(axes.title.size)), legend.text=element_text(size=rel(legend.text.size)),
				legend.title=element_text(size=rel(legend.title.size)))
	}
	if (var.axes) {
		if (circle) {
			theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi,
							length = 50))
			circle <- data.frame(xvar = r * cos(theta), yvar = r *
							sin(theta))
			g <- g + geom_path(data = circle, color = muted("white"),
					size = 1/2, alpha = 1/3)
		}
		g <- g + geom_segment(data = df.v, aes(x = 0, y = 0,
						xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2,
								"picas")), color = muted("red"))
	}
	if (!is.null(df.u$groups) && ellipse) {
		theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
		circle <- cbind(cos(theta), sin(theta))
		ell <- plyr::ddply(df.u, "groups", function(x) {
					if (nrow(x) <= 2) {
						return(NULL)
					}
					sigma <- var(cbind(x$xvar, x$yvar))
					mu <- c(mean(x$xvar), mean(x$yvar))
					ed <- sqrt(qchisq(ellipse.prob, df = 2))
					data.frame(sweep(circle %*% chol(sigma) * ed, 2,
									mu, FUN = "+"), groups = x$groups[1])
				})
		names(ell)[1:2] <- c("xvar", "yvar")
		g <- g + geom_path(data = ell, aes(color = groups, group = groups), size=ellipse.lwd)
	}
	if (var.axes) {
		g <- g + geom_text(data = df.v, aes(label = varname,
						x = xvar, y = yvar, angle = angle, hjust = hjust),
				color = "darkred", size = varname.size)
	}
	if (!is.null(df.u$labels)) {
		if (!is.null(df.u$groups)) {
			g <- g + geom_point(aes(color = groups), alpha = alpha, size=point.size)
			if (!is.null(grp.nam)) {
				g <- g + scale_color_discrete(name=grp.nam)
			}
			g <- g + geom_label_repel(
					aes(label = labels, color = groups),
					size = labels.size,
					box.padding = unit(0.35, "lines"),
					point.padding = unit(0.3, "lines"),
					show.legend = F
			)

#			g <- g + geom_text(aes(label = labels, color = groups),
#					size = labels.size)
			if (!is.null(grp.nam)) {
				g <- g + scale_color_discrete(name=grp.nam)
			}
		}
		else {
			g <- g + geom_point(alpha = alpha, size=point.size)
			g <- g + geom_label_repel(
					aes(label = labels),
					size = labels.size,
					box.padding = unit(0.35, "lines"),
					point.padding = unit(0.3, "lines"),
					show.legend = F
			)
#			g <- g + geom_text(aes(label = labels), size = labels.size)
		}
	}
	else {
		if (!is.null(df.u$groups)) {
			g <- g + geom_point(aes(color = groups), alpha = alpha, size=point.size)
			if (!is.null(grp.nam)) {
				g <- g + scale_color_discrete(name=grp.nam)
			}
		}
		else {
			g <- g + geom_point(alpha = alpha, size=point.size)
		}
	}
	if(!is.null(main)) {
		g <- g + ggtitle(main)
	}
	if (tidy) {
		g <- g + theme_minimal()
	}
	return(g)
}


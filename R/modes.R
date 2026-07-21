# Mode resolution: the single place where (pairs, block, do.voom) are interpreted.
#
# Two orthogonal axes resolve into one immutable spec:
#   pairs.role : how the `pairs` column enters the model  ("none" | "fixed" | "block")
#   engine     : test framework                           ("edger" | "limma")
# Valid modes are the legal cells of role x engine, with one constraint:
#   role == "block" forces engine == "limma" (edgeR has no duplicateCorrelation).
#
# Author: vidal
###############################################################################

# --- Registry -----------------------------------------------------------------
# Adding a scenario = add ONE entry here (and, if it needs a genuinely new
# mechanism, implement that in make_design / make_contrasts / diff_expr_fit).
#   design       : which design matrix to build      ("means" ~0+groups | "additive" ~pairs+groups)
#   contrasts    : how contrasts are obtained         ("explicit" makeContrasts | "coef" | "coef+explicit")
#   correlation  : use limma::duplicateCorrelation?   (TRUE/FALSE)
#   force_engine : NULL, or an engine this role pins  (e.g. "limma")
.dw_roles <- list(
  none  = list(design = "means",    contrasts = "explicit",      correlation = FALSE, force_engine = NULL),
  fixed = list(design = "additive", contrasts = "coef+explicit", correlation = FALSE, force_engine = NULL),
  block = list(design = "means",    contrasts = "explicit",      correlation = TRUE,  force_engine = "limma")
  # To add e.g. a paired-GLM-with-correlation role later:
  #   blockglm = list(design = "means", contrasts = "explicit", correlation = TRUE, force_engine = "edger")
  # then teach diff_expr_fit() how to honour correlation under the edgeR engine.
)

# --- Constructor --------------------------------------------------------------
new_dw_mode <- function(role, engine, notes = character()) {
  spec <- .dw_roles[[role]]
  structure(
    list(
      role        = role,
      engine      = engine,
      design      = spec$design,
      contrasts   = spec$contrasts,
      correlation = spec$correlation,
      do.voom     = identical(engine, "limma"),
      notes       = notes
    ),
    class = "diffwrap_mode"
  )
}

# --- Resolver -----------------------------------------------------------------
# `pairs` here is the argument as passed to diffExpr() (a column name or NULL),
# resolved BEFORE it is turned into a factor downstream.
dw_resolve_mode <- function(pairs = NULL, block = FALSE, do.voom = FALSE) {
  notes <- character()

  # --- axis 1: role ---
  role <- if (isTRUE(block)) {
    "block"
  } else if (!is.null(pairs)) {
    "fixed"
  } else {
    "none"
  }

  # --- validation: fail loud rather than take a wrong path silently ---
  if (isTRUE(block) && is.null(pairs)) {
    stop("'block = TRUE' requires a 'pairs' column to use as the correlation block.",
         call. = FALSE)
  }

  # --- axis 2: engine (with the one hard constraint) ---
  forced <- .dw_roles[[role]]$force_engine
  if (!is.null(forced)) {
    engine <- forced
    if (identical(forced, "limma") && !isTRUE(do.voom)) {
      notes <- c(notes, "block design -> voom enforced (edgeR has no duplicateCorrelation equivalent)")
    }
  } else {
    engine <- if (isTRUE(do.voom)) "limma" else "edger"
  }

  new_dw_mode(role = role, engine = engine, notes = notes)
}

# --- Log-friendly rendering ---------------------------------------------------
# format() is the primitive the log writer consumes; print() just displays it.
#' @export
format.diffwrap_mode <- function(x, ...) {
  role_lab   <- c(none = "unpaired", fixed = "paired (fixed effect)", block = "blocked (correlation)")[x$role]
  engine_lab <- c(edger = "edgeR GLM", limma = "limma/voom")[x$engine]
  lines <- c(
    sprintf("Resolved analysis mode: %s | %s", role_lab, engine_lab),
    sprintf("  design      : %s", if (x$design == "means") "~0 + groups (means model)" else "~pairs + groups (intercept)"),
    sprintf("  contrasts   : %s", x$contrasts),
    sprintf("  correlation : %s", if (x$correlation) "duplicateCorrelation" else "none")
  )
  if (length(x$notes)) lines <- c(lines, paste0("  note        : ", x$notes))
  lines
}

#' @export
print.diffwrap_mode <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

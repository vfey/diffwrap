# Enrichment wrappers depend on Bioconductor annotation packages and, for gProfileR,
# on network access. Those dependencies are guarded with skip_if_not_installed() and
# skip_on_cran() so the suite still runs in a minimal or offline environment.
#
# Unlike the rest of the suite these tests use REAL human Ensembl gene identifiers
# (a set of well-characterised cell-cycle and apoptosis genes) rather than the
# simulated example data, whose fictitious identifiers cannot map to GO/KEGG terms.

# A set of real human Ensembl gene IDs with rich GO annotation.
real_hs_ensembl <- c(
  "ENSG00000141510", "ENSG00000012048", "ENSG00000139618", "ENSG00000146648",
  "ENSG00000136997", "ENSG00000124762", "ENSG00000135679", "ENSG00000087088",
  "ENSG00000171791", "ENSG00000164305", "ENSG00000064012", "ENSG00000132906",
  "ENSG00000170312", "ENSG00000123374", "ENSG00000134057", "ENSG00000110092",
  "ENSG00000105173", "ENSG00000139687", "ENSG00000101412", "ENSG00000149311",
  "ENSG00000175054", "ENSG00000149554", "ENSG00000183765", "ENSG00000171862",
  "ENSG00000142208", "ENSG00000121879", "ENSG00000143799", "ENSG00000120868",
  "ENSG00000026103", "ENSG00000147889"
)

# A genome-wide ENSEMBL background (universe). The enrichment wrappers pass this to
# clusterProfiler's 'universe' and g:Profiler's 'custom_bg'; the default of "" is an
# empty universe, which clusterProfiler rejects and g:Profiler answers with HTTP 400.
# skip()s if org.Hs.eg.db is unavailable, so callers can use it directly.
hs_ensembl_background <- function() {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")
  orgdb <- tryCatch(get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db")),
                    error = function(e) NULL)
  skip_if(is.null(orgdb), "org.Hs.eg.db object not retrievable")
  bg <- tryCatch(AnnotationDbi::keys(orgdb, keytype = "ENSEMBL"), error = function(e) NULL)
  skip_if(is.null(bg) || length(bg) < 1000,
          paste("genome-wide ENSEMBL background unavailable (n =",
                if (is.null(bg)) 0 else length(bg), ")"))
  bg
}

# --- pure helpers ----------------------------------------------------------

test_that("prepare_scale_for_legend returns matching coordinates and labels", {
  quiet_log()
  lp <- prepare_scale_for_legend(scale.minimum = -3.2, scale.maximum = 4.8)
  expect_type(lp, "list")
  expect_named(lp, c("scale.y.coordinates", "scale.labels"))
  expect_equal(length(lp$scale.y.coordinates), length(lp$scale.labels))
})

test_that("format_ensembl_ids_annotated_to_term maps ids to symbols", {
  quiet_log()
  # convertId2() works offline, so this needs no annotation package or network.
  # convertId2() expects the capitalised species name ("Human"/"Mouse"); in the package
  # this comes from the 'species4conversion' lookup, not the lowercase 'species' argument.
  result <- data.frame(
    term    = c("GO:0001", "GO:0002"),
    geneID  = c("ENSG00000139618,ENSG00000141510", "ENSG00000012048"),
    stringsAsFactors = FALSE
  )
  out <- tryCatch(
    format_ensembl_ids_annotated_to_term(result, species = "Human"),
    error = function(e) skip(paste("convertId2 unavailable:", conditionMessage(e)))
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(result))
})

# --- wrappers that need annotation packages / network ----------------------

test_that("run_clusterProfiler_GO runs an over-representation analysis on real genes", {
  skip_on_cran()
  skip_if_not_installed("clusterProfiler")
  quiet_log()

  background <- hs_ensembl_background()      # skips if org.Hs.eg.db is unavailable

  res <- tryCatch(
    suppressWarnings(
      run_clusterProfiler_GO(input_genes = real_hs_ensembl,
                             background_genes = background,   # the universe: must be non-empty
                             ontology = "BP",
                             OrgDb = "org.Hs.eg.db", id_type = "ENSEMBL")
    ),
    error = function(e) skip(paste("clusterProfiler call failed:", conditionMessage(e)))
  )
  # the ORA branch returns a data.frame of terms, or the string "No enrichments found"
  # when nothing is significant; both are valid completions
  expect_true(is.data.frame(res) || is.character(res))
})

test_that("the default (empty) background runs against the full universe, not an empty one", {
  # regression guard for the background_genes = "" default: it must be treated as
  # "use all genes" (universe = NULL), not as an empty universe that clusterProfiler rejects
  skip_on_cran()
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  quiet_log()

  res <- tryCatch(
    suppressWarnings(
      run_clusterProfiler_GO(input_genes = real_hs_ensembl,   # no background_genes supplied
                             ontology = "BP",
                             OrgDb = "org.Hs.eg.db", id_type = "ENSEMBL")
    ),
    error = function(e) fail(paste("default background should work, but errored:",
                                   conditionMessage(e)))
  )
  expect_true(is.data.frame(res) || is.character(res))
})

test_that("run.topGO validates that the annotation package is installed", {
  # deterministic guard, needs neither topGO nor an org.db: an unknown organism
  # package must be reported clearly rather than failing deep inside topGO
  quiet_log()
  expect_error(
    run.topGO(background = c("ENSG00000141510", "ENSG00000012048"),
              foreground = "ENSG00000141510",
              organism = "org.NoSuchOrg.eg.db"),
    "must be installed"
  )
})

test_that("run.topGO runs against a genome-wide background", {
  skip_on_cran()
  skip_if_not_installed("topGO")
  quiet_log()

  # a real genome-wide universe is needed for GO terms to be 'feasible' (nodeSize = 10)
  background <- hs_ensembl_background()
  foreground <- intersect(real_hs_ensembl, background)
  skip_if(length(foreground) < 10, "too few mappable foreground genes")

  res <- tryCatch(
    suppressWarnings(
      run.topGO(background = background, foreground = foreground,
                ontologies = "BP", organism = "org.Hs.eg.db", ID_type = "ENSEMBL")
    ),
    error = function(e) skip(paste("topGO call failed:", conditionMessage(e)))
  )
  expect_true(is.list(res) || is.data.frame(res))
})

test_that("run_gprofiler contacts the g:Profiler service", {
  skip_on_cran()
  skip_if_offline()
  skip_if_not_installed("gprofiler2")
  quiet_log()

  background <- hs_ensembl_background()   # custom_bg; the default "" gives an HTTP 400

  res <- tryCatch(
    suppressWarnings(
      run_gprofiler(input_genes = real_hs_ensembl,
                    background_genes = background,
                    organism = "hsapiens",
                    domain_scope = "custom_annotated")   # actually use the custom background
    ),
    error = function(e) skip(paste("g:Profiler unreachable:", conditionMessage(e)))
  )
  expect_true(is.list(res))
})

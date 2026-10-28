# Run the simulation study.
#   Rscript benchmarks/simulation/run.R point 300 10   # study, replicates, cores
#   Rscript benchmarks/simulation/run.R mse 200 10      # bootstrap B = 100 (methods.R)
#   Rscript benchmarks/simulation/run.R point 20 10 B,C /tmp/pilot  # scenarios, output dir
# Uses the installed package (R CMD INSTALL . --preclean). One file per
# scenario in benchmarks/simulation/results/; existing files are skipped.
`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x
args <- commandArgs(TRUE)
study <- args[1] %||% "point"
R <- as.integer(args[2] %||% "300")
cores <- as.integer(args[3] %||% "10")
only <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else NULL

suppressPackageStartupMessages({
  library(flexsae)
  library(parallel)
})
# load the comparison packages in the parent: every replicate is a fresh fork
# (mc.preschedule = FALSE), so otherwise each fit timing includes loading them
for (pkg in c("sae", "nlfh")) if (requireNamespace(pkg, quietly = TRUE)) loadNamespace(pkg)
source("benchmarks/simulation/dgp.R")
source("benchmarks/simulation/methods.R")

scenarios <- if (study == "point") point_scenarios() else mse_scenarios()
if (!is.null(only)) scenarios <- Filter(function(s) s$id %in% only, scenarios)
out_dir <- if (length(args) >= 5) args[5] else file.path("benchmarks/simulation/results", study)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

RNGkind("L'Ecuyer-CMRG")
for (i in seq_along(scenarios)) {
  sc <- scenarios[[i]]
  file <- file.path(out_dir, paste0(sc$id, ".rds"))
  if (file.exists(file)) next
  design <- make_design(sc, seed = 1000 + i)
  set.seed(20261028 + i)
  streams <- vector("list", R)
  s <- .Random.seed
  for (r in seq_len(R)) {
    streams[[r]] <- s
    s <- nextRNGStream(s)
  }
  t0 <- Sys.time()
  reps <- mclapply(seq_len(R), function(r) {
    assign(".Random.seed", streams[[r]], envir = globalenv())
    rep <- draw_replicate(sc, design)
    res <- suppressWarnings(fit_methods(rep$data, model_covariates(sc$shape), mse = study == "mse",
                                        seed = r))
    list(theta = rep$theta, res = res)
  }, mc.cores = cores, mc.preschedule = FALSE)
  saveRDS(list(scenario = sc, design = design, R = R, reps = reps,
               session = list(flexsae = as.character(packageVersion("flexsae")), R = R.version.string,
                              time = Sys.time() - t0)), file)
  cat(sprintf("%s %s: %.1f min\n", format(Sys.time(), "%H:%M"), sc$id,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

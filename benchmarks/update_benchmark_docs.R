# Summarise the simulation results into inst/extdata (read by
# vignettes/simulation.Rmd). Run after benchmarks/simulation/run.R:
#   Rscript benchmarks/update_benchmark_docs.R
source("benchmarks/simulation/summarise.R")
source("benchmarks/simulation/dgp.R")

point <- summarise_dir("benchmarks/simulation/results/point", summarise_point)
mse <- summarise_dir("benchmarks/simulation/results/mse", summarise_mse)
# boost_sae() runs separately on the same replicates (check S12 in DISCREPANCIES.md):
#   Rscript benchmarks/simulation/run.R point 300 10 all benchmarks/simulation/results/boost \
#     fh_linear,boost_replace,boost_correct
if (dir.exists("benchmarks/simulation/results/boost")) {
  boost <- summarise_dir("benchmarks/simulation/results/boost", summarise_point)
  point <- rbind(point, boost[boost$method != "fh_linear", ])
}

describe <- function(scs) {
  do.call(rbind, lapply(scs, function(s) {
    data.frame(scenario = s$id, shape = s$shape, D = s$D, ratio = s$ratio, psi = s$psi,
               v_dist = s$v_dist, unsampled = s$unsampled, lambda = s$lambda, stringsAsFactors = FALSE)
  }))
}
files <- list.files("benchmarks/simulation/results/point", full.names = TRUE)
info <- readRDS(files[1])$session
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  point = point, mse = mse,
  scenarios = describe(point_scenarios()), mse_scenarios = describe(mse_scenarios()),
  settings = list(R_point = max(point$R), R_mse = max(mse$R), B = 100,
                  flexsae = info$flexsae, R_version = info$R, date = Sys.Date())
), "inst/extdata/simulation.rds")
cat("written inst/extdata/simulation.rds\n")

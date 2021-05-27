library(ascotsims)
library(parallel)
library(tidyverse)

print_out <- function(...) cat(sprintf(...), sep = '', file = stdout())

RNGkind("L'Ecuyer-CMRG")
set.seed(6562375)

ncores <- parallel::detectCores(logical = F) - 1

sims <- 10000

# Specify which configurations to simulate
cfg <- list(
  b = list(
    c(qlogis(0.2), 0, 0, 0, 0, log(1.10), 0),
    c(qlogis(0.2), 0, 0, 0, 0, log(1/1.10), 0),
    c(qlogis(0.2), 0, 0, 0, 0, log(1/1.25), 0),
    c(qlogis(0.2), 0, 0, 0, 0, log(1/1.50), 0),
    c(qlogis(0.2), 0, 0, 0, 0, log(1/2.00), 0)),
  n_seq = seq(400, 5000, 200),
  S0 = diag(c(100, rep(1, 6))),
  effective_thres = 0.99,
  futility_thres = 0.95,
  superior_thres = 0.99,
  inferior_thres = 0.01,
  delta = log(1.1)
)
cfgtbl <- enframe(cfg) %>%
  spread(name, value) %>%
  unnest(c(b, delta, effective_thres, futility_thres, inferior_thres, superior_thres)) %>%
  mutate(config = row_number())

# Setup output storage
res_par <- vector('list', nrow(cfgtbl))
res_arm <- vector('list', nrow(cfgtbl))
res_trt <- vector('list', nrow(cfgtbl))
res_trial <- vector('list', nrow(cfgtbl))
res_result <- vector('list', nrow(cfgtbl))
res_trt_intrm <- vector('list', nrow(cfgtbl))

begin <- Sys.time()

print_out("Running fin_eff_domB.R\n")
for(z in cfgtbl$config) {
  print_out("Running config: %3s\n", z)
  res <- mclapply(1:sims, function(j) {
    ascot_trial5(
      b = cfgtbl$b[[z]],
      n_seq = cfgtbl$n_seq[[z]],
      S0 = cfgtbl$S0[[z]],
      delta = cfgtbl$delta[z],
      effective_thres = cfgtbl$effective_thres[z],
      futility_thres = cfgtbl$futility_thres[z],
      inferior_thres = cfgtbl$inferior_thres[z],
      superior_thres = cfgtbl$superior_thres[z],
      min_any = F,
      early_t = 0,
      delay = 0,
      return_all =  TRUE)
  }, mc.cores = ncores)

  res_par[[z]] <- tibble_par_quantities(res, mc.cores = ncores) %>%
    mutate(config = z)
  res_arm[[z]] <- tibble_arm_quantities(res, mc.cores = ncores) %>%
    mutate(config = z)
  res_trt[[z]] <- tibble_trt_quantities(res, mc.cores = ncores) %>%
    mutate(config = z)
  res_trial[[z]] <- tibble_trial_quantities(res) %>%
    mutate(config = z)
  res_result[[z]] <- tibble_result_quantities(res) %>%
    mutate(config = z)
  res_trt_intrm[[z]] <- tibble_trt_quantities(res, final = F, mc.cores = ncores) %>%
    mutate(config = z)
}

end <- Sys.time()

# Group storage
res_par <- bind_rows(res_par, .id = "config") %>%
  group_by(config) %>%
  nest()
res_arm <- bind_rows(res_arm, .id = "config") %>%
  group_by(config) %>%
  nest()
res_trt <- bind_rows(res_trt, .id = "config") %>%
  group_by(config) %>%
  nest()
res_trial <- bind_rows(res_trial, .id = "config") %>%
  group_by(config) %>%
  nest()
res_result <- bind_rows(res_result, .id = "config") %>%
  group_by(config) %>%
  nest()
res_trt_intrm <- bind_rows(res_trt_intrm, .id = "config") %>%
  group_by(config) %>%
  nest()

saveRDS(
  list(
    cores = ncores,
    time = end - begin,
    sims = sims,
    cfg = cfgtbl,
    res_par = res_par,
    res_arm = res_arm,
    res_trt = res_trt,
    res_trial = res_trial,
    res_result = res_result,
    res_trt_intrm = res_trt_intrm),
  "~/out_files/ascot/ascot_sims5/fin_eff_domB.rds")

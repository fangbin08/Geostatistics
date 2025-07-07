library(RTMB)
library(tidyverse)
library(fmesher)

# ----------------------------
# 1. Load model data
# ----------------------------
# Set paths
path_pred <- '/Users/binfang/Documents/Household_variable_mapping/predictions/'
path_output <- '/Users/binfang/Documents/Household_variable_mapping/csv_files_accum/'

dat_read <- read.csv(paste0(path_pred, 'model_floors.csv'), na.strings = "") %>%
  drop_na()

# Subset to selected countries 
iso3_list <- c("COD", "BDI", "RWA", "UGA", "TZA")
dat_read <- dat_read[substr(dat_read$cluster, 1, 3) %in% iso3_list, ]

# Convert cluster to factor
dat_read$cluster <- as.factor(dat_read$cluster)

# Target variable as ordered factor
dat_read$target_var <- factor(dat_read$target_var, levels = c("Flr_Fin", "Flr_Rud", "Flr_Nat"))

# ----------------------------
# 2. Build SPDE mesh using coordinates
# ----------------------------
dat_read_sample <- dat_read %>% group_by(cluster) %>% slice_sample(n=1)
coords <- as.matrix(dat_read_sample[, c("Longitude", "Latitude")])
bnd <- fm_nonconvex_hull(coords, convex = -0.3)

write.csv(coords, paste0(path_pred, '/coords_cf.csv'), row.names = FALSE)

mesh <- fm_mesh_2d(
  loc = coords,
  boundary = bnd,
  max.edge = c(0.2, 2)
)

plot(mesh, main = "")
points(coords, pch = 21, bg = 1, col = "white", cex = 1.5)

spde <- fm_fem(mesh)
projObs <- fm_basis(mesh, loc = coords)

# ----------------------------
# 3. Prepare data for RTMB with spatial random effects
# ----------------------------
num_cols <- c("time", "Latitude", "Longitude", "accessibility", 
              "aridity", "cropland", "density", "elevation", 
              "evi", "footprint", "gdp", "growing", "hdi",
              "irrigation", "lst_range", "night_light",
              "pasture", "pet", "water_dist")
dat_read[num_cols] <- scale(dat_read[num_cols])

climate_dummies <- model.matrix(~ factor(climate) - 1, dat_read)
land_dummies    <- model.matrix(~ factor(land) - 1, dat_read)
urban_dummies   <- model.matrix(~ factor(urban) - 1, dat_read)

xmat <- as.matrix(cbind(dat_read[, num_cols], climate_dummies, land_dummies, urban_dummies))

dat <- list(
  y = as.integer(dat_read$target_var),
  x = xmat,
  projObs = projObs,
  spde = spde
)

# ----------------------------
# 4. Initialize parameters
# ----------------------------
n_class <- length(unique(dat$y))

par <- list(
  beta = rep(0, ncol(dat$x)),
  thbase = 0,
  logthdelta = rep(log(1), n_class - 1),
  omega = rep(0, mesh$n),
  logTau = 0,
  logKappa = 0,
  logsd = 0
)

map_list <- list(
  logthdelta = factor(c(NA, 1)),
  logsd = factor(NA)
)

# ----------------------------
# 5. Define joint negative log-likelihood with spatial random field
# ----------------------------
invlogit <- function(x) 1 / (1 + exp(-x))

jnll <- function(par) {
  "[<-" <- ADoverload("[<-")
  local <- getAll(par, dat)
  
  th <- numeric(n_class - 1)
  th[1] <- thbase
  for (j in 2:(n_class - 1)) {
    th[j] <- th[j - 1] + exp(logthdelta[j - 1])
  }
  
  Q <- exp(logTau) * (exp(logKappa)^4 * spde$c0 + 2 * exp(logKappa)^2 * spde$g1 + spde$g2)
  
  jnll_val <- -dgmrf(omega, Q = Q, log = TRUE)
  
  spatial_effect <- as.vector(projObs %*% omega)
  eta <- as.vector(x %*% beta) + spatial_effect
  
  logits <- invlogit(outer(th, eta, "-"))
  
  probs <- rbind(logits[1, ], 
                 apply(logits[-1, , drop = FALSE], 2, diff),
                 1 - logits[n_class - 1, ])
  
  prob_obs <- probs[cbind(dat$y, seq_along(dat$y))]
  
  jnll_val <- jnll_val - sum(log(prob_obs + 1e-12))
  return(jnll_val)
}

# ----------------------------
# 6. Compile and fit
# ----------------------------
obj <- MakeADFun(
  func = jnll,
  par = par,
  data = dat,
  random = "omega",
  map = map_list,
  control = list(newton.loop = 1),
  checkParameterOrder = FALSE
)

cat("Initial function value:", obj$fn(obj$par), "\n")
cat("Initial gradient norm:", sqrt(sum(obj$gr(obj$par)^2)), "\n")

fit <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 500, iter.max = 500, rel.tol = 1e-6))

sdr <- sdreport(obj)
pl <- as.list(sdr, "Est")

saveRDS(pl, paste0(path_pred, 'pl_spatial.rds'))

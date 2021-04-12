# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("PrevMap") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Cleaned react data
react <- readr::read_csv("data/processed/react_clean.csv")

# NON-SPATIAL MODELS -----------------------------------------------------------

# Covariates to include in the model
cov_names <- c("round", "elog_comprate")

f <- as.formula(paste("prev ~", paste0(cov_names, collapse = "+")))
f_linear <- as.formula(paste("elogit ~", paste0(cov_names, collapse = "+")))

# Linear model with emp. logit prev
fit_linear <- lm(f_linear, data = react)
summary(fit_linear)
pearson_residl <- residuals(fit_linear, type = "pearson")


# Binomial model 
fit_binom <- glm(f, family = "binomial", data = react, weights = tested)
summary(fit_binom)
pearson_residb <- residuals(fit_binom, type = "pearson")


# Compute variogram on residuals for each round
for (i in unique(react$round)) {
  id_round <- which(react$round == i)
  coords <- react[id_round, c("utm_x", "utm_y")]
  
  ggvario(coords, data = pearson_residl[id_round], xlab = "Distance (km)",
          show_nbins = F, nsim = 1000) +
    labs(subtitle = "Calculated from the Pearson residuals of a linear model
                     fitted to the empirical logit of prevalence.")
  ggsave(paste0("figs/empvariol_", gsub(" ", "", i), ".pdf"), w = 8, h = 6)
  
  ggvario(coords, data = pearson_residb[id_round], xlab = "Distance (km)",
          show_nbins = F, nsim = 1000) +
    labs(subtitle = "Calculated from the Pearson residuals of a binomial model.")
  ggsave(paste0("figs/empvariob_", gsub(" ", "", i), ".pdf"), w = 8, h = 6)
  
}

# SPATIAL MODEL ----------------------------------------------------------------

# Fit a geostatistical model for every round

react <- as.data.frame(react)

for (i in unique(react$round)) {
  
  message("Fitting model for ", i)
  
  # Subset  to round i
  temp <- react[react$round == i, ]
  
  # Start with a linear geostatistical model for the emp. logit 
  geo_linear <- linear.model.MLE(formula = update(f_linear, ~ . -round - elog_comprate),
                                 coords =  ~ utm_x + utm_y,
                                 data = temp, start.cov.pars = c(10, 1),
                                 kappa = 0.5, messages = T, method = "nlminb")
  
  print((summary(geo_linear, l = F)))
  
  # Fitting binomial geostatistical model
  cmcmc <- control.mcmc.MCML(n.sim = 60000, burnin = 10000, thin = 5)
  
  par0 <- as.numeric(coef(geo_linear))
  p <- length(par0) - 3
  f_geo <- update(f_linear, positive ~ . - round - elog_comprate)
  
  init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
  geo_binomial <- binomial.logistic.MCML(formula = f_geo, units.m = ~ tested,
                                         coords = ~ utm_x + utm_y,
                                         data = temp,
                                         par0 = par0,
                                         control.mcmc = cmcmc, 
                                         kappa = 0.5,
                                         start.cov.pars = init_cov_pars,
                                         method = "nlminb", messages = T)
  print(summary(geo_binomial, l = F))
  saveRDS(geo_binomial, 
          file = paste0("output/models/geobinom_", gsub(" ", "", i), ".rds"))
}

# Extract estimated parameters of the spatil GP for all rounds
rounds <- unique(react$round)
sigma2 <- rep(NA, length(rounds))
phi <- rep(NA, length(rounds))
tau2 <- rep(NA, length(rounds))

for (i in 1:length(rounds)) {
  mod <- readRDS(file = paste0("output/models/geobinom_", gsub(" ", "", rounds[i]), ".rds"))
  sigma2[i] <- as.numeric(coef(mod)[ncol(mod$D) + 1])
  phi[i] <- as.numeric(coef(mod)[ncol(mod$D) + 2])
  tau2[i] <- as.numeric(coef(mod)[ncol(mod$D) + 3])
}

# SPATIO - TEMPORAL MODEL ------------------------------------------------------

# Generate a numeric variable for time
times <- c("2020-05-01", "2020-06-01", 
           "2020-06-19", "2020-07-07",
           "2020-07-24", "2020-08-11",
           "2020-08-20", "2020-09-08")

times <- as.Date(times)
midpoints <- round(diff(times)[seq(1, length(times), by = 2)] / 2)

times <- as.numeric(times[seq(1, length(times), by = 2)] + midpoints)
times <- times - times[1] + 1

react$time <- rep(times, each = 315)

# Initialise model parameters
par0 <- as.numeric(coef(fit_binom))
p <- length(par0)
psi0 <- 15
par0 <- c(par0, c(mean(sigma2), mean(phi), mean(tau2)), psi0)
init_cov_pars <- c(par0[p + 2],  par0[p + 3] / par0[p + 1], par0[p + 4])

# Fit spatio-temporal binomial geostatistical model with
# separable double-matern space time correlation function
f_geo <- update(f, positive ~ .)
geo_binomial <- binomial.logistic.MCML(formula = f_geo, 
                                       units.m = ~ tested,
                                       coords = ~ utm_x + utm_y,
                                       times = ~ time, 
                                       sst.model = "DM",
                                       data = react,
                                       par0 = par0,
                                       control.mcmc = cmcmc, 
                                       kappa = 0.5,
                                       kappa.t = 0.5,
                                       start.cov.pars = init_cov_pars,
                                       method = "nlminb", messages = T)
print(summary(geo_binomial, l = F))
geo_binomial$formula <- f_geo
saveRDS(geo_binomial, "output/models/geotime_binomial.rds")

# Function to create grid of prediction and 
# project the population raster

elogit <- function(num, den) {
  log((num + 0.5) / (den - num + 0.5))
}

create_grid <- function(resolution, study_area, pop, cutoff = 0, filename = NULL) {
  crs <- st_crs(study_area)$proj4string
  pred_grid <- raster(res = resolution, ext = extent(study_area), crs = crs, val = 1)
  res_pred <- res(pred_grid)[1]
  res_pop <- res(pop)[1] * 110.567
  fct <- floor(res_pred / res_pop)
  if(fct >= 2) pop <- aggregate(pop, fun = sum, fact = fct) 
  if(crs != projection(pop)) pop <- projectRaster(pop, pred_grid,
                                                  method = "ngb")
  
  compare_cond <- compareRaster(pred_grid, pop, res = T, orig = T, 
                                stopiffalse = F, showwarning = T)
  stopifnot(compare_cond)    
  pop_mask <- pop
  pop_mask[pop_mask <= cutoff] <- NA
  # writeRaster(pop, filename)
  pred_raster = mask(pred_grid, pop_mask)
  out <- list(pop = pop, 
              raster = pred_raster,
              coords = rasterToPoints(pred_raster)[, 1:2])
  if(!is.null(filename)) saveRDS(out, file = filename)
  return(out) 
}


plotRaster <- function(x, ...) rasterVis::levelplot(x, margin = F, ...)

# Convert epsg to epsg KM
epsgKM <- function(x) {
  crs <- st_crs(x)
  proj4KM <- gsub(pattern = "+.units=m", replacement = "+units=km", 
                  crs$proj4string)
  return(proj4KM)
}

# Envelope for variogram
variog_envelope <- function (geodata, coords = geodata$coords, data = geodata$data, 
                             obj.variog, nsim = 99, save.sim = FALSE, messages) 
{
  call.fc <- match.call()
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  obj.variog$v <- NULL
  if ((is.matrix(data) | is.data.frame(data))) 
    if (ncol(data) > 1) 
      stop("envelops can be computed for only one data set at once")
  if (!is.null(obj.variog$estimator.type)) 
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    if (abs(obj.variog$lambda) < 1e-04) 
      data <- log(data)
    else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
  }
  xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ xmat + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim, "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data) {
    return(data[sample(1:n.data)])
  }
  simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, 
                       n.data = n.data)
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      temp <- .C("binit", as.integer(obj.variog$n.data), 
                 as.double(as.vector(coords[, 1])), as.double(as.vector(coords[, 
                                                                               2])), as.double(as.vector(sim)), as.integer(nbins), 
                 as.double(as.vector(obj.variog$bins.lim)), as.integer(estimator.type == 
                                                                         "modulus"), as.double(max(obj.variog$u)), as.double(cbin), 
                 vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin), 
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  }
  else {
    variog.vbin <- function(x, ...) {
      variog(geodata = geodata, 
             data = x, uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type, 
             nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist, 
             pairs.min = obj.variog$pairs.min, direction = obj.variog$direction, 
             tolerance = obj.variog$tolerance, messages.screen = FALSE,...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (save.sim == FALSE) 
    simula$data <- NULL
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], 
                  v.upper = limits[2, ])
  if (save.sim) 
    res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}

# Calculate and plot the variogram
ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords))/3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "distance", 
                    show_nbins = T) {
  require(geoR)
  require(ggplot2)
  coords <- as.matrix(coords)
  min_dist <- min(dist(coords))
  if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
  envmc <- variog_envelope(coords = coords, data = data, 
                           obj.variog = empvario, nsim = nsim, messages = F)
  dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                        lowemp = envmc$v.lower, upemp = envmc$v.upper, 
                        nbins = empvario$n)
  p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
    geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = .3) +
    geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
    scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                       breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
    scale_y_continuous(name = "semivariance", 
                       #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1), 
                       limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
    ggtitle("Empirical semivariogram") +
    theme_classic()
  p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
  if(show_nbins) p2 else p1
}

check_dups <- function(data, cols) {
  data$row_id <- 1:nrow(data)
  dups <- data %>%
    group_by_at(cols) %>% 
    filter(n() > 1)  %>% 
    ungroup()
  if(nrow(dups) == 0) {
    dups <- NULL
  } else {
    dups$dup_id <- unclass(factor(apply(dups[, cols], 1, paste0, collapse = "")))
    dups <- dups[order(dups$dup_id), ]
  }
  return(dups)
}


joint_weighted_prev <- function(joint_data) {
  # joint_data = df object with 2 + N columns:
  #   - aggregation ID
  #   - population (pop)
  #   - N columns, each one is a sample from the joint predictive distribution
  #     (on the logit scale)
  
  require(dplyr)
  require(tidyr)
  
  joint_list <- split(joint_data, f = joint_data[, 1])
  
  wpop_prev <- function(x) {
    pop <- t(as.matrix(x$pop))
    prev <- as.matrix(x[, - c(1:2)])
    num <- as.numeric(pop %*% prev)
    den <- sum(pop)
    prev_iu <- num / den
    return(prev_iu)
  }
  joint_wprev <- do.call(rbind, lapply(joint_list, wpop_prev))
  joint_wprev <- as_tibble(cbind(ID = row.names(joint_wprev), 
                                 as.data.frame(joint_wprev)))
  names(joint_wprev) <- c(names(joint_data)[1], 1:(ncol(joint_data) - 2))
  
  joint_wprev <- joint_wprev %>% 
    tidyr::gather(key = "sample_id", value = "prev", -1)
  
  return(joint_wprev)
}

create_tab <- function(x) {
  pars <- as.numeric(x$estimate)
  ci_up <- pars + qnorm(0.975) * sqrt(diag(x$covariance))
  ci_low <- pars - qnorm(0.975) * sqrt(diag(x$covariance))
  ci <- cbind(ci_low, ci_up) 
  estimates <- PrevMap:::coef.PrevMap(x) 
  p <- length(estimates) - 4
  estimates <- round(estimates, 3)
  
  ci[p + 3, ] <- ci[p + 3, ] + ci[p + 1, ]
  ci[(p + 1):length(estimates), ] <- exp(ci[(p + 1):length(estimates), ])
  ci <- round(ci, 3)
  ci <- apply(ci, 1, function(x) paste0("(", x[1], ", ", x[2], ")"))
  tab <- tibble(Parameter = names(estimates), Estimate = estimates, CI = ci)
  return(tab)
}

sptime_pred <- function (object, grid.pred, time.pred, predictors = NULL, 
                         control.mcmc, type = "marginal", 
                         S.sim = NULL, save.sim = TRUE,
                         plot.correlogram = TRUE, messages = TRUE) {
  
  if (length(predictors) > 0 && class(predictors) != "data.frame") 
    stop("'predictors' must be a data frame with columns' names matching those in the data used to fit the model.")
  if (length(predictors) > 0 && any(is.na(predictors))) 
    stop("missing values found in 'predictors'.")
  
  p <- object$p <- ncol(object$D)
  n.pred <- nrow(grid.pred)
  coords <- object$coords
  if (object$p == 1) {
    predictors <- matrix(1, nrow = n.pred)
  } else {
    if (length(dim(predictors)) == 0) 
      stop("covariates at prediction locations should be provided.")
    predictors <- as.matrix(model.matrix(delete.response(terms(formula(object$formula))), data = predictors))
    if (nrow(predictors) != nrow(grid.pred)) 
      stop("the provided values for 'predictors' do not match the number of prediction locations in 'grid.pred'.")
    if (ncol(predictors) != ncol(object$D)) 
      stop("the provided variables in 'predictors' do not match the number of explanatory variables used to fit the model.")
  }
  
  # Extract estimated model coefficients
  beta <- object$estimate[1:p]
  sigma2 <- exp(object$estimate[p + 1])
  phi <- exp(object$estimate[p + 2])
  if (length(object$fixed.rel.nugget) == 0) {
    tau2 <- sigma2 * exp(object$estimate[p + 3])
  } else {
    tau2 <- object$fixed.rel.nugget * sigma2
  }
  psi <- exp(object$estimate[p + 4])
  
  # Calcuate mu for observed locations and times
  mu <- object$D %*% beta
  
  # Generate spatial correlation matrix
  U <- as.matrix(dist(coords))
  R.s <- exp(- U / phi)
  
  # Generate time correlation matrix
  U.t <- as.matrix(dist(object$times))
  R.t <- exp(- U.t / psi)
  
  # Generat full spatio-temporal var covariance matrix
  Sigma <- sigma2 * R.s * R.t
  diag(Sigma) <- diag(Sigma) + tau2
  
  S.sim <- S.sim
  if (is.null(S.sim)) {
    S.sim.res <- Laplace.sampling(mu, Sigma, object$y, 
                                  object$units.m, control.mcmc, 
                                  plot.correlogram = plot.correlogram, 
                                  messages = messages)
    S.sim <- S.sim.res$samples
  }
  if(save.sim) saveRDS(S.sim, paste0(paste("laplace_sampling", Sys.Date(), sep = "_"), ".rds"))
  n.samples <- dim(S.sim)[1]
  
  
  U.pred.coords <- as.matrix(proxy::dist(grid.pred, coords))
  R.spred <- exp(- U.pred.coords / phi)
  
  U.pred.times <- as.matrix(proxy::dist(time.pred, object$times))  
  R.tpred <- exp(- U.pred.times / psi)
    
  C <- sigma2 * R.spred * R.tpred 
  Sigma.inv <- solve(Sigma)
  A <- C %*% Sigma.inv
  mu.pred <- as.numeric(predictors %*% beta)
  
  mu.cond <- mu.pred + A %*% (t(S.sim) - as.numeric(mu))
  if(type == "marginal") {
    if (messages) cat("Calculating", type, "predictions", "\n")
    sd.cond <- sqrt(sigma2 - rowSums(A * C))
    eta.sim <- matrix(rnorm(n.pred * n.samples, mu.cond, sd.cond), n.pred, n.samples)
  } else {
    if (messages) cat("Calculating joint predictions\n")
    U.spred <- as.matrix(dist(grid.pred))
    U.tpred <- as.matrix(dist(time.pred))
    Sigma.pred <- sigma2 * exp(- U.spred / phi) * exp(- U.tpred / psi)
    Sigma.cond <- Sigma.pred - A %*% t(C)
    eta.sim <- t(matrix(rnorm(n.pred * n.samples), ncol = n.pred) %*% chol(Sigma.cond) + t(mu.cond))
    rownames(eta.sim) <- NULL
  }
  out <- list()
  out$D_pred <- mu.pred
  out$samples <- eta.sim
  out$coords <- grid.pred
  out$times <- time.pred
  return(out)
}

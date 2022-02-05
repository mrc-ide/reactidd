
# Save this file as `R/plot_p_spline_prev.R`

#' Plotting function for the P_spline model
#'
#' @export
#' @param X date vector.
#' @param p_splinefit fit of the model to the same set of data using reactidd::stan_p_spline()
#' @param target_dist_between_knots sets the number of days between adjacent knots (default = 5)
#' @param spline_degree sets the degree of the splines (default = 3)
#' @param link_function sets the link function used for the model (logit for REACT-1, log for phe case data)
#' @param n_mean shape parameter of the gamma distribution used for the generation time
#' @param b_mean  rate parameter of the gamma distribution used for the generation time
#' @param min_date_num the minimum date for the period over which average growth rate should be calculated.
#' @param max_date_num the maximum date for the period over which average growth rate should be calculated.
#' @param label informative label to label the estimates produced
#' @return A single row of a dataframe that includes estimates of growth rate, R and doubling/halving times with 95% credible intervals
#'
estimate_p_spline_average_R <- function(X, p_spline_fit, target_dist_between_knots = 5, spline_degree = 3, link_function="logit",
                                        min_date_num, max_date_num, n_mean =2.29, b_mean =0.36 ,label="Estimate"){

  inv_logit <- function(Num){
    1/(1+exp(-Num))
  }
  min_date_num <- as.numeric(min_date_num)
  max_date_num <- as.numeric(max_date_num)
  ff <- rstan::extract(p_spline_fit)
  X <- as.numeric(X)
  X <- seq(min(X),max(X),by=1)

  min_date_numeric <- min(X)
  max_date_numeric <- max(X)
  num_knots <- ceiling((max_date_numeric- min_date_numeric)/target_dist_between_knots)+7
  days_per_knot <- (max_date_numeric - min_date_numeric)/(num_knots -7)
  num_basis <- num_knots + spline_degree - 1
  num_data <- length(X)
  knots <- unname(seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, length.out = num_knots))



  X_new <- seq(min(X)-3*days_per_knot, max(X)+3*days_per_knot, 0.1)
  B_true <- splines::bs(X_new, df=num_basis, degree=spline_degree, intercept = TRUE)
  B_true <- t(predict(B_true, X))
  Y_array <- array(data=NA, dim=c(nrow(ff$a), length(X)))

  #a0<-mean(ff$a0)
  for(i in seq_len(nrow(ff$a))){
    a <- array(NA, num_basis)

    for(j in seq_len(length(a))){
      a[j] <- ff$a[i,j]
    }
    Y_array[i,] <- as.vector(a%*%B_true)
  }

  growth_rate_list <- c()
  delta_x <- max_date_num - min_date_num
  min_index <- match(min_date_num, X)
  max_index <- match(max_date_num, X)

  for(i in seq_len(nrow(Y_array))){
    if(link_function=="logit"){
      delta_y <- log(inv_logit(Y_array[i,min_index]))-log(inv_logit(Y_array[i,max_index]))
    }
    else if(link_function=="log"){
      delta_y <- Y_array[i,min_index]-Y_array[i,max_index]
    }

    growth_rate_list <- c(growth_rate_list, -delta_y/delta_x)

  }

  growth_rate <- growth_rate_list

  r<-quantile(growth_rate, c(0.025,0.5,0.975))
  prob_r_pos <- length(growth_rate[growth_rate>0])*100/ length(growth_rate)
  R <- (1+r/b_mean)**n_mean
  t_half <- log(2)/r

  return(data.frame(label = label,
                    growth_rate = r[[2]],
                    growth_rate_lb = r[[1]],
                    growth_rate_ub = r[[3]],
                    R = R[[2]],
                    R_lb = R[[1]],
                    R_ub = R[[3]],
                    t_half_double = t_half[[2]],
                    t_half_double_lb = t_half[[1]],
                    t_half_double_ub = t_half[[3]],
                    prob_r_pos = prob_r_pos))
}




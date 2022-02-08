#' Simulate shedding data
#'
#' @export
#' @param N_3tests The number of indiviudals with 3 tests to simulate (default = 662).
#' @param N_2tests The number of indiviudals with 2 tests to simulate (default = 212).
#' @param min_T minimum date for second test after first date (default = 6)
#' @param max_T maximum date for second/third test after first date (default = 17)
#' @param sens the true value of P0 the sensitivity of the test (default = 0.79)
#' @param k the exponential decay rate of swab-positivity (default = 0.07126)
#' @param tau time delay before exponential decay (default = 0)
#' @return A data.frame  objectwith indiviudal level data for N_3tests + N_2tests individuals
#'

shedding_simulate_data <- function(N_3tests=662, N_2tests=212,min_T=6, max_T=17,sens=0.79, k=0.07126, tau=0){

  df <- data.frame()
  dif <- (max_T - min_T)/2

  for(i in seq_len(N_3tests)){

    # Random time for test to occur
    t1 <- round(runif(1,min_T,min_T+dif))
    t2 <- round(runif(1,min_T+dif, max_T))
    # True positive at those days?
    p1 <- rbinom(1, 1, min(exp(-k*(t1-tau)), 1) )
    p2gp1 <- min(exp(-k*(t2-tau)), 1)/min(exp(-k*(t1-tau)), 1)
    p2 <- p1 * rbinom(1, 1, p2gp1 )
    # Test positive at those days?
    e1 <- rbinom(1, 1, p1*sens)
    e2 <- rbinom(1, 1, p2*sens)

    row_df <- data.frame(passcode=i,
                         time_shed0 = 0,
                         time_shed1 = t1,
                         time_shed2 = t2,
                         estbinres_shed0 = 1,
                         estbinres_shed1 = e1,
                         estbinres_shed2 = e2)
    df <- rbind(df, row_df)
  }

  for(i in seq_len(N_2tests)){

    # Random time for test to occur
    t1 <- round(runif(1,min_T,max_T))
    # True positive at those days?
    p1 <- rbinom(1, 1, min(exp(-k*(t1-tau)), 1) )
    # Test positive at those days?
    e1 <- rbinom(1, 1, p1*sens)

    row_df <- data.frame(passcode=i+N_3tests,
                         time_shed0 = 0,
                         time_shed1 = t1,
                         time_shed2 = NA,
                         estbinres_shed0 = 1,
                         estbinres_shed1 = e1,
                         estbinres_shed2 = NA)
    df <- rbind(df, row_df)
  }


  return(df)
}

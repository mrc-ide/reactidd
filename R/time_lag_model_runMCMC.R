# Save this file as `R/time_lag_model_runMCMC.R

#' MCMC model fitting - run the model to get whole posterior
#'
#' @export
#' @param init_pars
#' @param N number of iterations
#' @param ensemble scaling factor parameter
#' @param dates
#' @param ENG
#' @return Average log-likelihood value over spline posterior specific values of beta and alpha

time_lag_model_runMCMC <- function(init_pars, N, ensemble, dates, ENG){
  chain <- data.frame(beta = init_pars[1],
                      alpha = init_pars[2])
  current_pars <- as.numeric(chain[1,])
  p_current <- reactidd::time_lag_model_my_prior(current_pars) + reactidd::time_lag_model_likelihood(ensemble,dates,ENG, beta=current_pars[1], alpha = current_pars[2])
  for(i in 1:N){
    print(i)
    current_pars <- as.numeric(chain[i,])
    new_pars <- reactidd::time_lag_model_param_new(current_pars)
    if(reactidd::time_lag_model_my_prior(new_pars)>-10000){
      p_proposal <- reactidd::time_lag_model_my_prior(new_pars) + reactidd::time_lag_model_likelihood(ensemble, dates, ENG, beta = new_pars[1], alpha = new_pars[2])
      p_accept <- p_proposal- p_current
      accept <- runif(1)< exp(p_accept)
    }
    else{
      accept =FALSE
    }

    if(accept){
      chain_add <- data.frame(beta = new_pars[1],
                              alpha = new_pars[2])
      chain<- rbind(chain,chain_add)
      p_current <- p_proposal
    }
    else{
      chain_add <- data.frame(beta = current_pars[1],
                              alpha = current_pars[2])
      chain<- rbind(chain,chain_add)
    }

  }
  return(chain)
}

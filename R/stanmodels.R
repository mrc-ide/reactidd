# Generated by rstantools.  Do not edit by hand.

# names of stan models
stanmodels <- c("b_splines_actual", "b_splines_actual_phe", "b_splines_actual_weighted", "b_splines_actual_weighted_incidence", "b_splines_actual_weighted_two_variants", "breakpoint_complex", "breakpoint_simple", "linear", "linear_phe", "linear_weighted", "shedding_exp", "shedding_exp_tdelay")

# load each stan module
Rcpp::loadModule("stan_fit4b_splines_actual_mod", what = TRUE)
Rcpp::loadModule("stan_fit4b_splines_actual_phe_mod", what = TRUE)
Rcpp::loadModule("stan_fit4b_splines_actual_weighted_mod", what = TRUE)
Rcpp::loadModule("stan_fit4b_splines_actual_weighted_incidence_mod", what = TRUE)
Rcpp::loadModule("stan_fit4b_splines_actual_weighted_two_variants_mod", what = TRUE)
Rcpp::loadModule("stan_fit4breakpoint_complex_mod", what = TRUE)
Rcpp::loadModule("stan_fit4breakpoint_simple_mod", what = TRUE)
Rcpp::loadModule("stan_fit4linear_mod", what = TRUE)
Rcpp::loadModule("stan_fit4linear_phe_mod", what = TRUE)
Rcpp::loadModule("stan_fit4linear_weighted_mod", what = TRUE)
Rcpp::loadModule("stan_fit4shedding_exp_mod", what = TRUE)
Rcpp::loadModule("stan_fit4shedding_exp_tdelay_mod", what = TRUE)

# instantiate each stanmodel object
stanmodels <- sapply(stanmodels, function(model_name) {
  # create C++ code for stan model
  stan_file <- if(dir.exists("stan")) "stan" else file.path("inst", "stan")
  stan_file <- file.path(stan_file, paste0(model_name, ".stan"))
  stanfit <- rstan::stanc_builder(stan_file,
                                  allow_undefined = TRUE,
                                  obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name,
                            model_cppcode = stanfit$cppcode)
  # create stanmodel object
  methods::new(Class = "stanmodel",
               model_name = stanfit$model_name,
               model_code = stanfit$model_code,
               model_cpp = stanfit$model_cpp,
               mk_cppmodule = function(x) get(paste0("model_", model_name)))
})

data {
  int num_data;             // number of data points
  int Y[num_data];
  vector[num_data] X;
}


parameters {
  real alpha;  // intercept
  real beta;  // gradient
  real<lower=0> phi;
}


model {
  //Likelihood
  Y ~ neg_binomial_2_log(X*beta+alpha, phi);
}

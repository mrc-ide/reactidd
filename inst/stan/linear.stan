data {
  int num_data;             // number of data points
  int Y[num_data];
  int N[num_data];
  vector[num_data] X;
}


parameters {
  real alpha;  // intercept
  real beta;  // gradient
}


model {
  //Likelihood
  Y ~ binomial(N, exp(X*beta+alpha));
}

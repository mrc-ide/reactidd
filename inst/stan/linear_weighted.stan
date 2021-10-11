data {
  int num_data;             // number of data points
  row_vector[num_data] Y;
  row_vector[num_data] N;
  vector[num_data] X;
}


parameters {
  real alpha;  // intercept
  real beta;  // gradient
}

transformed parameters {
  vector[num_data] Y_hat;
  Y_hat = exp((X*beta)+alpha);
}

model {
  //Likelihood
  target+= Y*log(Y_hat) + (N-Y)*log(1-Y_hat);
}

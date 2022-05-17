functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
      //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
    b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
          (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
          (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
        w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}


data {
  int num_data;			// number of dates	
  int num_knots;		// number of knots
  vector[num_knots] knots;	// the sequence of knots
  int spline_degree;		// the degree of spline (is equal to order - 1)
  row_vector[num_data] Y;		// Data for group 1
  row_vector[num_data] N;
  row_vector[num_data] YL;		// Data for group 2
  row_vector[num_data] NL;
  real X[num_data];		// Dates
}

transformed data {
  int num_basis = num_knots + spline_degree - 1;					// total number of B-splines
  matrix[num_basis, num_data] B;							// matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots;					// set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
}


parameters {
  row_vector[num_basis] a_raw1;
  row_vector[num_basis] a_raw2;
  real<lower=0> tau;
  real<lower=0> eta;
}


transformed parameters {
  row_vector[num_basis] a1;
  row_vector[num_basis] a2;
  vector[num_data] Y_hat1;
  vector[num_data] Y_hat2;
  a1[1] = a_raw1[1];
  a1[2] = a_raw1[2];
  a2[1] = a_raw2[1];
  a2[2] = a_raw2[2];

  for (i in 3:num_basis){
    a1[i] = 2*a1[i-1] - a1[i-2] + a_raw1[i]*tau;
    a2[i] = 2*a2[i-1] - a2[i-2] + a_raw2[i]*tau;
  }
  Y_hat1 = inv_logit(to_vector(a1*B));
  Y_hat2 = inv_logit(to_vector(a2*B));


}

model {
  vector[num_data] prop;
  for (i in 1:num_data){
    prop[i] = Y_hat1[i]/(Y_hat1[i]+Y_hat2[i]);
  }

  // Priors
  a_raw1[3:num_basis] ~ normal(0,1);
  a_raw2[3:num_basis] ~ normal(0,1);

  tau ~ inv_gamma(0.0001, 0.0001);
  eta ~ inv_gamma(0.0001, 0.0001);

  // Mixed Effects
  a_raw1[3:num_basis] ~ normal(a_raw2[3:num_basis], eta);
  

  //Likelihood
  target+= Y*log(Y_hat1+Y_hat2) + (N-Y)*log(1-Y_hat1-Y_hat2);
  target+= YL*log(prop) + (NL-YL)*log(1-prop);
}

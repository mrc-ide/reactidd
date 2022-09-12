// Author: Oliver Eales
// Github: Eales96
data {
  int num_data;           	// number of data points
  int num_all;           	// number of days
  row_vector[num_data] Y;	// Weighed positives
  row_vector[num_data] N;	// Weighted totals
  vector[num_all] Vax1;	// Proportion vaccinated 1 dose
  vector[num_all] Vax2;	// Proportion vaccinated 2 dose
  vector[num_all] Delta;	// Proportion infected Delta
  real X[num_data];		// Dates (as numerics)
  real max_X;
  real X_all[num_all];
  int indexes[num_data];
  real restriction_dates[5];
  real gen;
}

transformed data{
  real tau1;
  real tau2;
  real tau3;
  real tau4;
  real tau5;

  tau1 = restriction_dates[1];
  tau2 = restriction_dates[2];
  tau3 = restriction_dates[3];
  tau4 = restriction_dates[4];
  tau5 = restriction_dates[5];
}


parameters {
  real<lower=0> prev0;
  vector<lower=-0.5, upper=0.5>[6] beta;			// growth factors
  real<lower=0, upper=14> delay;
  real<lower=0, upper=1> vef1;
  real<lower=0, upper=vef1> vef2;
  real<lower=1., upper=10> delta_trans;
}


transformed parameters {
  vector[num_all] beta_int;		// intrinsic beta at time t
  vector[num_all] beta_t;		// actual beta at time t
  vector[num_all] prev_t;		// prev_t
  vector[num_data] Y_hat;		//
  
  for(i in 1:num_all){
    if(X_all[i] < tau1+delay){
      beta_int[i] = exp(beta[1]);
    }
    else if(X_all[i] < tau2+delay){
      beta_int[i] = exp(beta[2]);
    }
    else if(X_all[i] < tau3+delay){
      beta_int[i] = exp(beta[3]);
    }
    else if(X_all[i] < tau4+delay){
      beta_int[i] = exp(beta[4]);
    }
    else if(X_all[i] < tau5+delay){
      beta_int[i] = exp(beta[5]);
    }
    else{
      beta_int[i] = exp(beta[6]);
    }
  }


  beta_t = beta_int .*( (1-Delta+Delta*delta_trans) .*(1-Vax1-Vax2+Vax1*vef1+Vax2*vef2));
  
  prev_t[1] = prev0;

  for(i in 2:num_all){
    prev_t[i] = prev_t[i-1]*pow(beta_t[i-1], 1/gen);
  }

  Y_hat = prev_t[indexes];
  
}

model {
  //priors

  //Likelihood
  target+= Y*log(Y_hat) + (N-Y)*log(1-Y_hat);
}

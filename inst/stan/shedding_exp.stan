data {
  int N_11;            			// number of participants with two results
  vector[N_11] time1_11;		// time of first repeat test

  int N_10;            			// number of participants with two results
  vector[N_10] time1_10;		// time of first repeat test

  int N_111;            		// number of participants with three results
  vector[N_111] time1_111;		// time of first repeat test
  vector[N_111] time2_111;		// time of second repeat test

  int N_101;            		// number of participants with three results
  vector[N_101] time1_101;		// time of first repeat test
  vector[N_101] time2_101;		// time of second repeat test

  int N_110;            		// number of participants with three results
  vector[N_110] time1_110;		// time of first repeat test
  vector[N_110] time2_110;		// time of second repeat test

  int N_100;            		// number of participants with three results
  vector[N_100] time1_100;		// time of first repeat test
  vector[N_100] time2_100;		// time of second repeat test


}

parameters {
  real<lower=0, upper=1> sens;

  real<lower=0> k;

}



model {
  target+=log(sens*exp(-time1_11*k));

  target+=log( (1-sens)*exp(-time1_10*k) + (1-exp(-time1_10*k)) );

  target+=log( sens*sens*exp(-time2_111*k) );

  target+=log( sens*(1-sens)*exp(-time2_110*k) + sens*(exp(-time1_110*k) - exp(-time2_110*k)) );

  target+=log( (1-sens)*sens*exp(-time2_101*k) );

  target+=log( (1-sens)*(1-sens)*exp(-time2_100*k) + (1-sens)*(exp(-time1_100*k) - exp(-time2_100*k)) +(1-exp(-time1_100*k)) );

}

generated quantities{
  real log_lik[N_11+N_10+N_111+N_110+N_101+N_100];
  for(i in 1:N_11){
    log_lik[i] = log(sens*exp(-time1_11[i]*k));
  }
  for(i in (1+N_11):(N_11+N_10)){
    log_lik[i] = log( (1-sens)*exp(-time1_10[i-N_11]*k) + (1-exp(-time1_10[i-N_11]*k)) );
  }
  for(i in (1+N_11+N_10):(N_11+N_10+N_111)){
    log_lik[i] =log( sens*sens*exp(-time2_111[i-N_11-N_10]*k) );
  }
  for(i in (1+N_11+N_10+N_111):(N_11+N_10+N_111+N_110)){
    log_lik[i] =log( sens*(1-sens)*exp(-time2_110[i-N_11-N_10-N_111]*k) + sens*(exp(-time1_110[i-N_11-N_10-N_111]*k) - exp(-time2_110[i-N_11-N_10-N_111]*k)) );
  }
  for(i in (1+N_11+N_10+N_111+N_110):(N_11+N_10+N_111+N_110+N_101)){
    log_lik[i] = log( (1-sens)*sens*exp(-time2_101[i-N_11-N_10-N_111-N_110]*k) );
  }
  for(i in (1+N_11+N_10+N_111+N_110+N_101):(N_11+N_10+N_111+N_110+N_101+N_100)){
    log_lik[i] = log( (1-sens)*(1-sens)*exp(-time2_100[i-N_11-N_10-N_111-N_110-N_101]*k) + (1-sens)*(exp(-time1_100[i-N_11-N_10-N_111-N_110-N_101]*k) - exp(-time2_100[i-N_11-N_10-N_111-N_110-N_101]*k)) +(1-exp(-time1_100[i-N_11-N_10-N_111-N_110-N_101]*k)) );
  }
}

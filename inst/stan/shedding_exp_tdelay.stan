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

  real<lower=0,upper=10> t_d;
}

transformed parameters {
  vector[N_11] prop1_11;		
  vector[N_10] prop1_10;		
  vector[N_111] prop1_111;	
  vector[N_111] prop2_111;		
  vector[N_101] prop1_101;		
  vector[N_101] prop2_101;		
  vector[N_110] prop1_110;		
  vector[N_110] prop2_110;		
  vector[N_100] prop1_100;		
  vector[N_100] prop2_100;		
  
  for(i in 1:N_11) {
    prop1_11[i] = fmin(1, exp(-k*(time1_11[i]-t_d)) );
  }

  for(i in 1:N_10){
    prop1_10[i] = fmin(1, exp(-k*(time1_10[i]-t_d)) );
  }

  for(i in 1:N_111){
    prop1_111[i] = fmin(1, exp(-k*(time1_111[i]-t_d)) );
    prop2_111[i] = fmin(1, exp(-k*(time2_111[i]-t_d)) );
  }
  
  for(i in 1:N_110){
    prop1_110[i] = fmin(1, exp(-k*(time1_110[i]-t_d)) );
    prop2_110[i] = fmin(1, exp(-k*(time2_110[i]-t_d)) );
  }

  for(i in 1:N_101){
    prop1_101[i] = fmin(1, exp(-k*(time1_101[i]-t_d)) );
    prop2_101[i] = fmin(1, exp(-k*(time2_101[i]-t_d)) );
  }

  for(i in 1:N_100){
    prop1_100[i] = fmin(1, exp(-k*(time1_100[i]-t_d)) );
    prop2_100[i] = fmin(1, exp(-k*(time1_100[i]-t_d)) );
  }
}



model {
  target+=log(sens*prop1_11);

  target+=log( (1-sens)*prop1_10 + (1-prop1_10) );

  target+=log( sens*sens*prop2_111 );

  target+=log( sens*(1-sens)*prop2_110 + sens*(prop1_110 - prop2_110) );

  target+=log( (1-sens)*sens*prop2_101 );

  target+=log( (1-sens)*(1-sens)*prop2_100 + (1-sens)*(prop1_100 - prop2_100) +(1-prop1_100) );

}

generated quantities{
  real log_lik[N_11+N_10+N_111+N_110+N_101+N_100];
  for(i in 1:N_11){
    log_lik[i] = log(sens*prop1_11[i]);
  }
  for(i in (1+N_11):(N_11+N_10)){
    log_lik[i] = log( (1-sens)*prop1_10[i-N_11] + (1-prop1_10[i-N_11]) );
  }
  for(i in (1+N_11+N_10):(N_11+N_10+N_111)){
    log_lik[i] = log( sens*sens*prop2_111[i-N_11-N_10] );
  }
  for(i in (1+N_11+N_10+N_111):(N_11+N_10+N_111+N_110)){
    log_lik[i] =log( sens*(1-sens)*prop2_110[i-N_11-N_10-N_111] + sens*(prop1_110[i-N_11-N_10-N_111] - prop2_110[i-N_11-N_10-N_111]) );
  }
  for(i in (1+N_11+N_10+N_111+N_110):(N_11+N_10+N_111+N_110+N_101)){
    log_lik[i] = log( (1-sens)*sens*prop2_101[i-N_11-N_10-N_111-N_110] );
  }
  for(i in (1+N_11+N_10+N_111+N_110+N_101):(N_11+N_10+N_111+N_110+N_101+N_100)){
    log_lik[i] = log( (1-sens)*(1-sens)*prop2_100[i-N_11-N_10-N_111-N_110-N_101] + (1-sens)*(prop1_100[i-N_11-N_10-N_111-N_110-N_101] - prop2_100[i-N_11-N_10-N_111-N_110-N_101]) +(1-prop1_100[i-N_11-N_10-N_111-N_110-N_101]) );
  }
}

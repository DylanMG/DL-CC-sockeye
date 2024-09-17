data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops
  vector[N] logRS; //the response variable
  int<lower=1> pop[N]; //vector indexing populaitons
  vector[J] SmaxPR; //Smax estimates for each pop
  vector[J] SmaxPR_SD; //Smax error for each pop
  matrix[N,J] X; //predictor matrix of spawners
  matrix[N,J] int_X; //matrix to turn on/off the intercept
  //data for out of sample predictions
  int<lower=1> N_OOS; //N obs for OOS to predict
  vector[N_OOS] y_OOS; //log(recruits per spawner)
  matrix[N_OOS,J] X_OOS; //spawners in time T
  matrix[N_OOS,J] int_X_OOS; //spawners in time T
  int<lower=1> pop_OOS; //numeric population factor level for OOS estimate
}

transformed data{
  vector[J] logbeta_pr;
  vector[J] logbeta_pr_sig;

  for(i in 1:J){
    logbeta_pr_sig[i] = sqrt(log(1+((1/SmaxPR_SD[i])*(1/SmaxPR_SD[i]))/((1/SmaxPR[i])*(1/SmaxPR[i])))); //this converts sigma on the untransformed scale to a log scale
    logbeta_pr[i] = log(1/SmaxPR[i])-0.5*logbeta_pr_sig[i]*logbeta_pr_sig[i]; //convert smax prior to per capita slope - transform to log scale with bias correction
  }
}

parameters {
  vector<lower=0>[J] log_a_pop; //populations' alpha estimates
  vector[J] log_b_pop; //populations' beta estimates
  vector<lower=0>[J] sigma_pop; //error
}

transformed parameters {
  vector[J] b_pop;
  b_pop = exp(log_b_pop);
}

model {
  //priors
  log_a_pop ~ normal(1.5,2);
  log_b_pop ~ normal(logbeta_pr, logbeta_pr_sig); 
  sigma_pop ~ normal(1,1);
  
  //likelihood
  logRS ~ normal(int_X*log_a_pop - X*b_pop, sigma_pop[pop]);
}

generated quantities{
  vector[N_OOS]  log_lik; 
  for (i in 1:N_OOS){log_lik[i] = normal_lpdf(y_OOS[i] | int_X_OOS[i]*log_a_pop - X_OOS[i]*b_pop , sigma_pop[pop_OOS]);
  }
  vector[J] Smax_pop; //calc Smax directly
  for (i in 1:J){Smax_pop[i] = 1/b_pop[i]; 
  }
  array[N] real logRS_rep; //posterior predictive check
  logRS_rep = normal_rng(int_X*log_a_pop - X*b_pop, sigma_pop[pop]);
}
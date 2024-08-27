data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops
  int<lower=1> pop[N]; //vector indexing populations
  matrix[N, J] X; //predictor matrix of spawners
  vector[N] logRS; //the response variable
  vector[J] SmaxPR; //Smax estimates for each pop
  vector[J] SmaxPR_SD; //Smax error for each pop
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
  real log_a; //global alpha
  real<lower=0> sd_a_pop; //variance of a means around pops
  vector<lower=0>[J] log_a_pop; //pop specific a's bound by 0

  vector[J] log_b_pop; //slopes for each pop 
  
  real<lower=0> sigma; //global variance term
  vector<lower=0>[J] sigma_pop; //pop specific sigmas bound by 0
  real<lower=0> sd_sigma_pop; //variance of a means around pops
}

transformed parameters {
  vector[J] b_pop;
  b_pop = exp(log_b_pop);
}

model {
  //priors
  log_a ~ normal(1.5,2);
  sd_a_pop ~ gamma(2,3);  
  log_a_pop ~ normal(log_a, sd_a_pop);
  
  log_b_pop ~ normal(logbeta_pr, logbeta_pr_sig); 
  
  //variance priors
  sigma ~ normal(1,1); //variance priors
  sd_sigma_pop ~ normal(0,1);
  sigma_pop ~ normal(sigma, sd_sigma_pop);
  
  //likelihood model
  logRS ~ normal(log_a_pop[pop] - X*b_pop, sigma_pop[pop]);
}

//helper to get log_lik from extract_log_lik help file
generated quantities{
  vector[N]  log_lik; 
  for (i in 1:N){log_lik[i] = normal_lpdf(logRS[i] | log_a_pop[pop[i]] - X[i]*b_pop, sigma_pop[pop[i]]);
  }
    vector[J] Smax_pop; 
  for (i in 1:J){Smax_pop[i] = 1/b_pop[i]; //calc Smax directly
  }
}

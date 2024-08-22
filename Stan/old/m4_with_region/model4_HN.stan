data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops
  int<lower=1> K; //n regions
  int<lower=1> pop[N]; //vector indexing populations
  int<lower=1> region[N]; //vector indexing regions
  int<lower=1> reg_pop[J]; //region-pop index
  matrix[N, J] X; //predictor matrix of spawners
  vector[N] logRS; //the response variable
  vector[J] SmaxPR; //Smax estimates for each pop
  vector[J] SmaxPR_SD; //Smax error for each pop
}

parameters {
  real log_a; //global alpha
  vector<lower=0>[J] log_a_pop; //population alpha deviations
  vector[K] log_a_region; //regional alpha deviations
  vector<lower=0>[J] b_pop; //slope sfor each pop 
  real<lower=0> sigma; //main variance term
  real<lower=0> sd_a_pop; //variance of means for around pops
  real<lower=0> sd_a_region; //variance of means for region
}

model {
  //priors
  log_a ~ normal(0,1);
  log_a_region ~ normal(0, sd_a_region);
  log_a_pop ~ normal(log_a + log_a_region[reg_pop], sd_a_pop);
  b_pop ~ normal(pow(SmaxPR, -1), pow(SmaxPR_SD, -1));
  
  //variance priors
  sigma ~ gamma(2,3); 
  sd_a_pop ~ normal(0,1);
  sd_a_region ~ normal(0,1);

  //likelihood model
  logRS ~ normal(log_a_pop[pop] - X*b_pop, sigma);
}

//helper to get log_lik from extract_log_lik help file
generated quantities{
  vector[N]  log_lik; 
  for (i in 1:N){log_lik[i] = normal_lpdf(logRS[i] | log_a_pop[pop[i]] - X[i]*b_pop, sigma);
  }
}

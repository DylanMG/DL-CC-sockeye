data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops 
  int<lower=1> K; //n regions
  matrix[N, J] X; //predictor matrix of spawners 
  vector[N] logRS; //outcome variable
  int<lower=1> pop[N]; //vector indexing populaitons
  int<lower=1> region[N]; //vector indexing regions
  int<lower=1> reg_pop[J]; //region-pop index
}

parameters {
  real log_a; //global alpha
  vector<lower=0>[J] log_a_pop; //population alpha deviations
  vector[K] log_a_region; //regional alpha deviations
  vector[J] log_b_pop; //logged slopes
  real<lower=0> sigma; //error
  real<lower=0> sd_a_pop; //variance of means for pops
  real<lower=0> sd_a_region; //variance of means for region
}

transformed parameters {
  vector[J] b_pop; 
  b_pop = exp(log_b_pop); //prevents b (density dependence) from being negative
}

model {
  //priors
  log_a ~ gamma(3,2);
  log_a_region ~ normal(0, sd_a_region);
  log_a_pop ~ normal(log_a + log_a_region[reg_pop], sd_a_pop);
  log_b_pop ~ normal(-12, 3);
  
  //variance priors
  sigma ~ gamma(2,3);
  sd_a_pop ~ gamma(2,3);
  sd_a_region ~ gamma(2,3);

  logRS ~ normal(log_a_pop[pop] - X*b_pop, sigma); 
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i] = normal_lpdf(logRS[i] | log_a_pop[pop[i]] - X[i]*b_pop, sigma);
  }
}

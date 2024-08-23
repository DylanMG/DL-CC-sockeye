data {
  int<lower=1> N; //length of data
  int<lower=1> J; //number of pops in the model 
  vector[N] logRS; //outcome variable
  int<lower=1> pop[N]; //vector indexing populaitons
  matrix[N, J] X; //predictor matrix of spawners
}

parameters {
  real log_a; //global alpha
  vector<lower=0>[J] log_a_pop; //pop specific as bound by 0
  real<lower=0> sd_a_pop; //variance of a means around pops
  vector[J] log_b_pop; //slopes for each pop 
  
  real<lower=0> sigma; //global variance term
  vector<lower=0>[J] sigma_pop; //pop specific sigmas bound by 0
  real<lower=0> sd_sigma_pop; //variance of a means around pops
}

transformed parameters {
  vector[J] b_pop; 
  b_pop = exp(log_b_pop); //prevents b (density dependence) from being negative
}

model {
  log_a ~ normal(1.5,2);//priors
  sd_a_pop ~ gamma(2,3);  
  log_a_pop ~ normal(log_a, sd_a_pop);
  log_b_pop ~ normal(-12, 3);
  
  sigma ~ normal(1,1); //variance priors
  sd_sigma_pop ~ normal(0,1);
  sigma_pop ~ normal(sigma, sd_sigma_pop);

  logRS ~ normal(log_a_pop[pop] - X*b_pop, sigma_pop[pop]); //likliehood
}                                                       

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i] = normal_lpdf(logRS[i]|log_a_pop[pop[i]] - X[i]*b_pop, sigma_pop[pop[i]]);
  }
  vector[J] Smax_pop; //calc Smax directly
  for (i in 1:J){Smax_pop[i] = 1/b_pop[i]; 
  }
  array[N] real logRS_rep; //posterior predictive check
  logRS_rep = normal_rng(log_a_pop[pop] - X*b_pop, sigma_pop[pop]);
}

data {
  int<lower=1> N; //length of data
  int<lower=1> J; //number of pops in the model 
  vector[N] logRS; //outcome variable
  int<lower=1> pop[N]; //vector indexing populaitons
  matrix[N, J] X; //predictor matrix of spawners
}

parameters {
  real log_a; //global alpha
  vector<lower=0>[J] log_a_pop; //pop specific a's bound by 0
  vector[J] log_a_dev; //pop specific deviation from global
  vector[J] log_b_pop; //vector of fixed slope priors
  
  real<lower=0> sigma; //error
  real sd_a_pop; //variance of a means around pops
}

transformed parameters {
  vector[J] b_pop; 
  b_pop = exp(log_b_pop); //prevents b (density dependence) from being negative
}

model {
  log_a ~ gamma(3,2); 
  log_a_dev ~ normal(0,1); //pop level deviaitons from the global alpha
  log_a_pop ~ normal(log_a + log_a_dev, sd_a_pop);
  log_b_pop ~ normal(-12, 3);
  
  sigma ~ gamma(2,3);  //global sigma
  sd_a_pop ~ gamma(2,3);

  logRS ~ normal(log_a_pop[pop] - X*b_pop, sigma);
}                                                       

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i] = normal_lpdf(logRS[i]|log_a_pop[pop[i]] - X[i]*b_pop, sigma);
  }
}

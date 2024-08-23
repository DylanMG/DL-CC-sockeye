data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops
  vector[N] logRS; //the response variable
  int<lower=1> pop[N]; //vector indexing populaitons
  matrix[N,J] X; //predictor matrix of spawners
  matrix[N,J] int_X; //matrix to turn on/off the intercept so alphas stay independent
}

parameters {
  vector<lower=0>[J] log_a_pop; //intercepts
  vector[J] log_b_pop; //population slopes
  vector<lower=0>[J] sigma_pop; //error
}

transformed parameters {
  vector[J] b_pop; 
  b_pop = exp(log_b_pop); //prevents b (density dependence) from being negative
}

model {
  //priors
  log_a_pop ~ normal(1.5,2);
  log_b_pop ~ normal(-12, 3);
  sigma_pop ~ normal(1,1);

  //likelihood
  logRS ~ normal(int_X*log_a_pop - X*b_pop, sigma_pop[pop]);
}

//helper to get log_lik from extract_log_lik help file
generated quantities{
  vector[N]  log_lik; 
  for (i in 1:N){log_lik[i] = normal_lpdf(logRS[i] | int_X[i]*log_a_pop - X[i]*b_pop , sigma_pop);
  }
  vector[J] Smax_pop; //calc Smax directly
  for (i in 1:J){Smax_pop[i] = 1/b_pop[i]; 
  }
  array[N] real logRS_rep; //posterior predictive check
  logRS_rep = normal_rng(int_X*log_a_pop - X*b_pop, sigma_pop[pop]);
}

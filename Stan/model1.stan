data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops
  vector[N] logRS; //the response variable
  matrix[N,J] X; //predictor matrix of spawners
  matrix[N,J] int_X; //matrix to turn on/off the intercept
}

parameters {
  vector<lower=0>[J] log_a_pop; //intercepts
  vector[J] log_b_pop; //population slopes
  real<lower=0> sigma; //error
}

transformed parameters {
  vector[J] b_pop; 
  b_pop = exp(log_b_pop); //prevents b (density dependence) from being negative
}

model {
  //priors
  log_a_pop ~ gamma(3,2);
  log_b_pop ~ normal(-12, 3);
  sigma ~ gamma(2,3);

  //likelihood
  logRS ~ normal(int_X*log_a_pop - X*b_pop, sigma);
}

//helper to get log_lik from extract_log_lik help file
generated quantities{
  vector[N]  log_lik; 
  for (i in 1:N){log_lik[i] = normal_lpdf(logRS[i] | int_X[i]*log_a_pop - X[i]*b_pop , sigma);
  }
}

data {
  int<lower=1> N; //length of data
  int<lower=1> J; //n pops
  vector[N] logRS; //the response variable
  vector[J] SmaxPR; //Smax estimates for each pop
  vector[J] SmaxPR_SD; //Smax error for each pop
  matrix[N,J] X; //predictor matrix of spawners
  matrix[N,J] int_X; //matrix to turn on/off the intercept
}

parameters {
  vector<lower=0>[J] log_a_pop; //populations' alpha estimates
  vector<lower=0>[J] b_pop; //populations' beta estimates
  real<lower=0> sigma; //error
}


model {
  //priors
  log_a_pop ~ gamma(3,2);
  b_pop ~ normal(pow(SmaxPR, -1), pow(SmaxPR_SD, -1));

  sigma ~ gamma(2,3); 
  //likelihood model
  logRS ~ normal(int_X*log_a_pop - X*b_pop, sigma);
}

generated quantities{
  vector[N]  log_lik; 
  for (i in 1:N){log_lik[i] = normal_lpdf(logRS[i] | int_X[i]*log_a_pop - X[i]*b_pop, sigma);
  }
}


data {
  int n_schools;
  array[n_schools] real y;
  vector<lower = 0>[n_schools] sigma;
}


parameters {
  real mu;
  real<lower = 0> tau;
  vector[n_schools] theta;
}

model {
  mu ~ normal(5, 3);
  tau ~ normal(0, 10);

  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}

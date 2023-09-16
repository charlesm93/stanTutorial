
data {
  int n_schools;
  array[n_schools] real y;
  vector<lower = 0>[n_schools] sigma;
}

parameters {
  real mu;
  real<lower = 0> tau;
  vector[n_schools] eta;
  // vector[n_schools] theta;
}

transformed parameters {
  vector[n_schools] theta = mu + tau * eta;
}

model {
  mu ~ normal(5, 3);
  tau ~ normal(0, 10);

  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}


data {
  int n_schools;
  array[n_schools] real y;
  vector<lower = 0>[n_schools] sigma;
}

parameters {
  real mu;
  real<lower = 0> tau;
  // vector[n_schools] theta;
}

model {
  mu ~ normal(5, 3);
  tau ~ normal(0, 10);

  y ~ normal(mu, sqrt(square(tau) + square(sigma)));  // p(y | mu, tau)
}

generated quantities {
  real theta[n_schools];
  for (i in 1:n_schools) {
    real conjugate_variance =  1 / (1 / square(sigma[i]) + 1 / square(tau));
    real conjugate_mean =
      (y[i] / square(sigma[i]) + mu / square(tau)) * conjugate_variance;

    theta[i] = normal_rng(conjugate_mean, sqrt(conjugate_variance));
  }  // p(theta | y, mu, tau)
}

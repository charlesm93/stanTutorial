
data {
  int<lower = 1> nEvent;
  int<lower = 1> nObs;
  array[nObs] int<lower = 1> iObs;

  // Event schedule
  array[nEvent] int<lower = 1> cmt;
  array[nEvent] int evid;
  array[nEvent] int addl;
  array[nEvent] int ss;
  array[nEvent] real amt;
  array[nEvent] real time;
  array[nEvent] real rate;
  array[nEvent] real ii;

  vector<lower = 0>[nObs] cObs;
}

transformed data {
  int nTheta = 5;
  int nCmt = 3;
}


parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> VC;
  real<lower = 0> VP;
  real<lower = 0> ka;
  real<lower = 0> sigma;
}

transformed parameters {
  array[nTheta] real theta = {CL, Q, VC, VP, ka};
  matrix<lower = 0>[nCmt, nEvent]
    mass = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);

  row_vector<lower = 0>[nEvent] concentration = mass[2, ] ./ VC;
}

model {
  // priors
  CL ~ lognormal(log(10), 0.25); 
  Q ~ lognormal(log(15), 0.5);
  VC ~ lognormal(log(35), 0.25);
  VP ~ lognormal(log(105), 0.5);
  ka ~ lognormal(log(2.5), 1);
  sigma ~ normal(0, 1);

  // likelihood
  cObs ~ lognormal(log(concentration[iObs]), sigma);
}

generated quantities {
  array[nObs] real concentrationObsPred 
    = exp(normal_rng(log(concentration[iObs]), sigma));
    
  vector[nObs] log_lik;
  for (i in 1:nObs) log_lik[i]
    = lognormal_lpdf(cObs[i] | log(concentration[iObs[i]]), sigma);
}

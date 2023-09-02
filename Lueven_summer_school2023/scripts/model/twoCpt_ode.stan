
functions {
  vector system(real time, vector y, array[] real theta,
                array[] real x_r, array[] int x_i) {
    real CL = theta[1];
    real Q = theta[2];
    real VC = theta[3];
    real VP = theta[4];
    real ka = theta[5];

    vector[3] dydt;
    dydt[1] = - ka * y[1];
    dydt[2] = ka * y[1] - (CL + Q) / VC * y[2] + Q / VP * y[3];
    dydt[3] = Q / VC * y[2] - Q / VP * y[3];

    return dydt;
  }
}

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
  // vector[nObs] logCObs = log(cObs);
  int nTheta = 5;
  int nCmt = 3;
  
  // NOTE: optional control parameters for pmx_solve_rk45.
  real rel_tol = 1e-2;
  real abs_tol = 1e-2;
  real rel_tol_IS = 1e-10;
  real abs_tol_IS = 1e-10;
  int max_num_steps = 1000;
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
  matrix<lower = 0>[nCmt, nEvent] mass;

  // mass = pmx_solve_rk45(system, nCmt,
  //                       time, amt, rate, ii, evid, cmt, addl, ss, theta);

  // Torsten function with control parameters for ODEs.
  mass = pmx_solve_rk45(system, nCmt,
                        time, amt, rate, ii, evid, cmt, addl, ss, theta,
                        rel_tol, abs_tol, max_num_steps);

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
  
  // PSIS diagnostic for ODE tolerance
  matrix<lower = 0>[nCmt, nEvent] massIS
    = pmx_solve_rk45(system, nCmt,
                     time, amt, rate, ii, evid, cmt, addl, ss, theta,
                     rel_tol_IS, abs_tol_IS, max_num_steps);
  
  row_vector<lower = 0>[nEvent] concentrationIS = massIS[2, ] ./ VC;
  
  vector[nObs] log_lik;
  for (i in 1:nObs) log_lik[i]
    = lognormal_lpdf(cObs[i] | log(concentration[iObs[i]]), sigma);
  
  vector[nObs] log_lik_IS;
  for (i in 1:nObs) log_lik_IS[i]
    = lognormal_lpdf(cObs[i] | log(concentrationIS[iObs[i]]), sigma);
  
  real log_ratios = sum(log_lik_IS - log_lik);
}

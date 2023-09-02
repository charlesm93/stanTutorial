
data {
  int<lower = 1> nEvent;
  int<lower = 1> nObs;
  array[nObs] int<lower = 1> iObs;
  int<lower = 1> nSubjects;
  array[nSubjects] int<lower = 1> start;
  array[nSubjects] int<lower = 1> end;

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
  int nIIV = 5;
}

parameters {
  // Population parameters
  real<lower = 0> CL_pop;
  real<lower = 0> Q_pop;
  real<lower = 0> VC_pop;
  real<lower = 0> VP_pop;
  real<lower = 0> ka_pop;

  // Inter-individual variability
  vector<lower = 0>[nIIV] omega;
  matrix[nSubjects, nTheta] eta;
  real<lower = 0> sigma;
}

transformed parameters {
  array[nSubjects, nTheta] real theta;
  vector<lower = 0>[nTheta] 
    theta_pop = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, ka_pop});
  row_vector<lower = 0>[nEvent] concentration;
  row_vector<lower = 0>[nObs] concentrationObs;
  matrix<lower = 0>[nCmt, nEvent] mass;

  for (j in 1:nSubjects) {
    theta[j, ] = to_array_1d(exp(omega .* eta'[, j]) .* theta_pop);

    mass[, start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]],
                                               amt[start[j]:end[j]],
                                               rate[start[j]:end[j]],
                                               ii[start[j]:end[j]],
                                               evid[start[j]:end[j]],
                                               cmt[start[j]:end[j]],
                                               addl[start[j]:end[j]],
                                               ss[start[j]:end[j]],
                                               theta[j, ]);

    concentration[start[j]:end[j]] = 
                      mass[2, start[j]:end[j]] / theta[j, 3];
  }

  concentrationObs = concentration[iObs];
}

model {
  // priors
  CL_pop ~ lognormal(log(10), 0.25); 
  Q_pop ~ lognormal(log(15), 0.5);
  VC_pop ~ lognormal(log(35), 0.25);
  VP_pop ~ lognormal(log(105), 0.5);
  ka_pop ~ lognormal(log(1), 0.25);
  sigma ~ normal(0, 1);
  omega ~ normal(0, 0.2);

  // hierarchical prior
  for (j in 1:nSubjects) eta[j, ] ~ normal(0, 1);

  // likelihood
  cObs ~ lognormal(log(concentrationObs), sigma);
}


generated quantities {
  // predictions for existing patients
  array[nObs] real concentrationObsPred
    = lognormal_rng(log(concentrationObs), sigma);

  // predictions for new patients
  array[nObs] real cObsNewPred;
  matrix<lower = 0>[nCmt, nEvent] massNew;
  array[nSubjects, nTheta] real thetaNew;
  row_vector<lower = 0>[nEvent] concentrationNew;
  row_vector<lower = 0>[nObs] concentrationObsNew;

  for (j in 1:nSubjects) {
    thetaNew[j, ] = lognormal_rng(log(theta_pop), omega);

    massNew[, start[j]:end[j]]
      = pmx_solve_twocpt(time[start[j]:end[j]],
                         amt[start[j]:end[j]],
                         rate[start[j]:end[j]],
                         ii[start[j]:end[j]],
                         evid[start[j]:end[j]],
                         cmt[start[j]:end[j]],
                         addl[start[j]:end[j]],
                         ss[start[j]:end[j]],
                         thetaNew[j, ]);

      concentrationNew[start[j]:end[j]]
        = massNew[2, start[j]:end[j]] / thetaNew[j, 3];

      concentrationObsNew = concentrationNew[iObs];
  }

  cObsNewPred = lognormal_rng(log(concentrationObsNew), sigma);
}

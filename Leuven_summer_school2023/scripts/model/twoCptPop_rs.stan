
functions {
  real partial_sum(array[] int subject, int start_subject, int end_subject,
                   array[] int start, array[] int end,
                   array[] real time, array[] real amt, array[] real rate,
                   array[] real ii, array[] int evid, array[] int cmt,
                   array[] int addl, array[] int ss, int nCmt,
                   array[] int iObs, real sigma, vector cObs,
                   vector omega, matrix eta, vector theta_pop) {
    // NOTE: I used the fact there is the same number of events
    // and observations for each patient to simplify the bookkeeping.
    
    // index off-set due to the fact we may not start with the first patient
    int index_offset = start[start_subject] - 1;
    
    int nEvents_sub = end[end_subject] - start[start_subject] + 1;
    int nSubjects_sub = size(subject);
    int nObs_sub = nEvents_sub - nSubjects_sub;
    int nTheta = 5;
    
    array[nSubjects_sub, nTheta] real theta;
    for (j in 1:nSubjects_sub) {
      theta[j, ] = to_array_1d(
        exp(omega .* eta'[, j + start_subject - 1]) .* theta_pop);
    }
    
    
  // for (j in 1:nSubjects) {
  //   theta[j, ] = to_array_1d(theta_pop); 
  //   // to_array_1d(exp(omega .* eta'[, j]) .* theta_pop);
  // }
    
    row_vector[nEvents_sub] concentration;
    row_vector[nObs_sub] concentrationObs;
    matrix[nCmt, nEvents_sub] mass;
    
    for (j in 1:nSubjects_sub) {
      int start_offset = start[subject[j]] - index_offset;
      int end_offset = end[subject[j]] - index_offset;
      
      mass[, start[j]:end[j]]
        = pmx_solve_twocpt(time[start_offset:end_offset],
                           amt[start_offset:end_offset],
                           rate[start_offset:end_offset],
                           ii[start_offset:end_offset],
                           evid[start_offset:end_offset],
                           cmt[start_offset:end_offset],
                           addl[start_offset:end_offset],
                           ss[start_offset:end_offset],
                           theta[j, ]);
      
      concentration[start[j]:end[j]]
        = mass[2, start[j]:end[j]] / theta[j, 3];
    }

    int nObs_patient = 19;
    int startObs = (start_subject - 1) * nObs_patient + 1;
    int endObs = startObs + nObs_patient * nSubjects_sub - 1;

    concentrationObs
      = concentration[iObs[1:(nObs_patient * nSubjects_sub)]];

    vector[nObs_patient * nSubjects_sub] cObs_sub = cObs[startObs:endObs];

    return (lognormal_lpdf(cObs_sub | log(concentrationObs), sigma));
  }
}

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
  
  // index for partial_sum function
  array[nSubjects] int subject;
  for (i in 1:nSubjects) subject[i] = i;
  
  int grain_size = 1;
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
  // array[nSubjects, nTheta] real theta;
  vector<lower = 0>[nTheta] 
    theta_pop = to_vector({CL_pop, Q_pop, VC_pop, VP_pop, ka_pop});
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
  target += reduce_sum(partial_sum, subject, grain_size,
                       start, end,
                       time, amt, rate, ii, evid, cmt, addl, ss, nCmt,
                       iObs, sigma, cObs, omega, eta, theta_pop);
}


generated quantities {
  // predictions for existing patients
  // array[nObs] real concentrationObsPred
  //   = lognormal_rng(log(concentrationObs), sigma);

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

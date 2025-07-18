data {
  int<lower=1> N; // total number of observations
  array[N] int Y; // response variable
  int<lower=1> m; // number of basis
  // data for spline 
  // spline basis function matrices for control group
  matrix[N, m] S;
}



parameters {
  real Intercept; // temporary intercept for centered predictors
  // parameters for spline s(age_std, by = trt, bs = "cs")0
  // scalar for the spline coefficient
  real<lower = 0> alpha;
  vector[m] gamma; // unscaled spline coefficient in the control group
  // unscaled spline coefficients in the treatment group
}

transformed parameters {
  // actual spline coefficients for control group
  vector[m] s_0;
  

  // compute actual spline coefficients
  s_0 = alpha * gamma;
  
}

model {
  // define sum_gamma, sum_c, yhat, and ZS in the model block
  real sum_gamma;

  sum_gamma = sum(gamma);

  

  alpha ~ exponential(1);
  Intercept ~ normal(0,3);
  
  // random walk prior for spline coefficient
  
  for (j in 2:m) {
    gamma[j] ~ normal(gamma[j-1],1);
  }
  sum_gamma ~ normal(0, 0.01*m);
  
  
  
  target += bernoulli_logit_lpmf(Y | Intercept + S*s_0);

}

generated quantities {
  
}

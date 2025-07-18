// fit a monotonic spline model to each arm seperately

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
  real alpha;
  vector[m-1] gamma; // unscaled spline coefficient on the log scale
}

transformed parameters {
  // actual spline coefficients for control group
  vector[m] b;
  vector[m] coef; 
  vector[m] s_0;
  
  b = append_row(0, gamma);
  coef = softmax(b);
  

  // compute actual spline coefficients
  s_0 = alpha * coef;
  
}

model {
  // define yhat, and ZS in the model block

  
  
  profile("priors"){
    target += normal_lupdf(Intercept |0, 3);
    target += normal_lupdf(alpha |0, 3);
    
    target += logistic_lupdf(gamma | 0, 0.5);
  }
  // AR(1) prior for spline coefficient
  profile("likelihood"){
    target += bernoulli_logit_lpmf(Y | Intercept + S*s_0);
  }
  
}

generated quantities {
  
}
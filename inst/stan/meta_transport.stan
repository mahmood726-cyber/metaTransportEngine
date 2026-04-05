data {
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1> T;
  vector[S] y;
  vector<lower=0>[S] v;
  matrix[S, P] Z;
  matrix[T, P] ZT;
}

parameters {
  real beta0;
  vector[P] beta;
  real<lower=0> tau;
  vector[S] theta_raw;
}

transformed parameters {
  vector[S] mu;
  vector[S] theta;

  mu = beta0 + Z * beta;
  theta = mu + tau * theta_raw;
}

model {
  beta0 ~ normal(0, 10);
  beta ~ normal(0, 2.5);
  tau ~ normal(0, 0.5);
  theta_raw ~ normal(0, 1);

  y ~ normal(theta, sqrt(v));
}

generated quantities {
  vector[T] theta_T;
  vector[T] theta_T_pred;

  for (t in 1:T) {
    theta_T[t] = beta0 + ZT[t] * beta;
    theta_T_pred[t] = normal_rng(theta_T[t], tau);
  }
}

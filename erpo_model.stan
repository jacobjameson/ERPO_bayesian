
data {
  int<lower=1> N; int<lower=1> n_states; int<lower=1> n_years;
  int<lower=1> KP; int<lower=1> Kx;

  int<lower=0> y[N];
  vector[N] l1; vector[N] l2;
  matrix[N,KP]  P;  matrix[N,KP]  P1;  matrix[N,KP]  P2;
  vector[N] offset; vector[N] offset1; vector[N] offset2;
  matrix[N,Kx] X;
  matrix[N,n_years] T; matrix[N,n_years] T1; matrix[N,n_years] T2;
  int<lower=1,upper=n_states> S[N];
  real<lower=0> prior_sd;
}
parameters {
  real alpha;
  real<lower=0,upper=1>  delta1;
  real<lower=-1,upper=1> delta2;
  vector[KP] beta;
  vector[Kx] gamma_X;
  vector[n_years] zeta;
  real<lower=0> invphi;
  vector[n_states] state_raw;
}
transformed parameters {
  vector[n_states] state_eff = state_raw - mean(state_raw);
  real phi = 1.0 / invphi;
  vector[N] log_mu =
        alpha
      + delta1 .* (l1 - offset1 - P1*beta)
      + delta2 .* (l2 - offset2 - P2*beta)
      + P*beta + X*gamma_X + T*zeta + offset;

  for (i in 1:N) log_mu[i] += state_eff[S[i]];
}
model {
  alpha  ~ normal(0, sqrt(10));
  delta1 ~ normal(0.5, 1);
  delta2 ~ normal(0,   1);

  beta     ~ normal(0, prior_sd);
  gamma_X  ~ normal(0, 0.1);
  zeta     ~ normal(0, 0.2);
  invphi   ~ normal(0, 0.1) T[0,];
  state_raw~ normal(0, 0.5);

  y ~ neg_binomial_2_log(log_mu, phi);
}
generated quantities {
  vector[6] IRR;
  IRR[1] = exp(beta[1]);                     //  year 0
  IRR[2] = exp(beta[1] + 0.2*beta[2]);       //  year 1
  IRR[3] = exp(beta[1] + 0.4*beta[2]);       //  year 2
  IRR[4] = exp(beta[1] + 0.6*beta[2]);       //  year 3
  IRR[5] = exp(beta[1] + 0.8*beta[2]);       //  year 4
  IRR[6] = exp(beta[1] +       beta[2]);     //  >=5 yr
}


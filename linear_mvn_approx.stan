functions {
  // numerically robust covariance estimate; ym is expected to be the
  // column wise mean
  matrix robust_cov(matrix y, row_vector ym) {
    matrix[cols(y),cols(y)] cov;

    //cov = yc' * yc /(rows(y) - 1);
    cov = crossprod(y - rep_matrix(ym, rows(y))) / (rows(y) - 1.0);

    return(0.5 * (cov + cov'));
  }

  real mvn_approx_lpdf(vector y_prime_bar, vector x_prime,
                       int J_prime,
                       vector delta,
                       vector mu_a, real beta,
                       vector sigma_a, real sigma_y,
                       matrix xi  // N(0,1) as K x J_tilde
                       ) {
    int J_tilde = cols(xi);
    int K       = rows(delta);
    int T_prime = rows(y_prime_bar);
    vector[T_prime] M_tilde;
    matrix[T_prime,T_prime] Sigma_tilde;
    vector[K] mu_a_tilde;
    matrix[K,J_tilde] a_tilde;
    matrix[T_prime,J_tilde] y_tilde;

    mu_a_tilde = mu_a + delta;

    a_tilde = rep_matrix(mu_a_tilde, J_tilde) + diag_pre_multiply(sigma_a, xi);

    for (k in 1:T_prime) {
      y_tilde[k] = a_tilde[1] + a_tilde[2]*x_prime[k] + beta*x_prime[k]^2;
      M_tilde[k] = mean(y_tilde[k]);
    }
    // simulated cov matrix + residual error term
    Sigma_tilde = robust_cov(y_tilde', to_row_vector(M_tilde)) +
      diag_matrix(rep_vector( square(sigma_y), T_prime));

    return multi_normal_lpdf(y_prime_bar| M_tilde, Sigma_tilde/J_prime);
  }

  vector colMeans(matrix y) {
    vector[cols(y)] m;
    for (i in 1:cols(y))
      m[i] = mean(col(y, i));
    return(m);
  }

}
data {
  int<lower=-1,upper=1> fit_all;
  int<lower= 0,upper=1> fit_local;
  int J;                         // #people in y
  int T;                         // #time points in y
  int K;                         // #parameters in delta
  int K_phi;                     // #parameters in phi
  matrix[J,T] y;
  int J_prime;                   // #people in second dataset
  int T_prime;                   // #time points in second dataset
  matrix[J_prime,T_prime] y_prime;
  vector[T] x;
  vector[T_prime] x_prime;
  vector[K_phi] mu_phi_p;        // prior mean
  cov_matrix[K_phi] Sigma_phi_p; // prior variance
  vector[K] mu_delta_p;          // prior mean
  cov_matrix[K] Sigma_delta_p;   // prior variance
  int J_tilde;                   // number of fake patients used for approximate log-lik
  int C;
  int<lower=1, upper=C> CHAIN_ID;
  matrix[K, 2*J_tilde] xi[C];         // N(0,1) random numbers used for integration
}
transformed data {
  vector[T_prime] ones         = rep_vector(1, T_prime);
  matrix[T_prime,K] X_prime    = append_col(ones, x_prime);
  cov_matrix[T_prime] identity = diag_matrix(ones);
  int K_phi_mean = K_phi - 3;
  cholesky_factor_cov[K_phi_mean] L_Sigma_phi_p   = cholesky_decompose(Sigma_phi_p[1:K_phi_mean,1:K_phi_mean]);
  cholesky_factor_cov[K]          L_Sigma_delta_p = cholesky_decompose(Sigma_delta_p);
  int J_tilde_chain = (CHAIN_ID % 2 == 0) ? 2*J_tilde : J_tilde;
  matrix[K, J_tilde_chain] xi_tilde = xi[CHAIN_ID, 1:K, 1:J_tilde_chain];
  vector[T_prime] y_prime_bar = colMeans(y_prime);

  print("Running chain ", CHAIN_ID, " with J_tilde_chain = ", J_tilde_chain, "...");
}
parameters {
  matrix[K,J] eta_a;
  real<lower=0> sigma_y;
  vector[K] mu_a;
  vector<lower=0>[K] sigma_a;    // Now assuming indep
  real beta;                     // Shared parameter
  vector[K] delta;
  matrix[K,J_prime] eta_prime_a;
}
transformed parameters {
  matrix[K,J] a;
  matrix[K,J_prime] a_prime;
  vector[K_phi] phi;
  real y_log;
  real y_prime_log;
  
  a = rep_matrix(mu_a, J) + diag_pre_multiply(sigma_a, eta_a);
  a_prime = rep_matrix(mu_a + delta, J_prime) + diag_pre_multiply(sigma_a, eta_prime_a);
  
  phi[1] = mu_a[1];
  phi[2] = mu_a[2];
  phi[3] = beta;
  phi[4] = log(sigma_a[1]);
  phi[5] = log(sigma_a[2]);
  phi[6] = log(sigma_y);

  {
    matrix[J,T] y_pred;
    for (j in 1:J)
      for (k in 1:T)
        y_pred[j,k] = a[1,j] + a[2,j]*x[k] + beta*x[k]^2;
  
    y_log = normal_lpdf(to_vector(y)| to_vector(y_pred), sigma_y);
  }

  if (!fit_local) {
    if (fit_all == 1) {
      // include raw data from prime data set
      matrix[J_prime,T_prime] y_prime_pred;
      for (j in 1:J_prime)
        for (k in 1:T_prime)
          y_prime_pred[j,k] = a_prime[1,j] + a_prime[2,j]*x_prime[k] + beta*x_prime[k]^2;
      
      y_prime_log = normal_lpdf(to_vector(y_prime)| to_vector(y_prime_pred), sigma_y);
    } else if (fit_all == -1) {
      // integrated approach
      vector[T_prime] y_prime_pred_bar;
      matrix[T_prime,T_prime] Sigma_y_prime_bar;
      Sigma_y_prime_bar = (X_prime*diag_matrix(sigma_a .* sigma_a)*X_prime' + sigma_y^2*identity)/J_prime;
      for (k in 1:T_prime)
        y_prime_pred_bar[k] = mu_a[1] + delta[1] + (mu_a[2] + delta[2])*x_prime[k] + beta*x_prime[k]^2;
      
       y_prime_log = multi_normal_lpdf(y_prime_bar| y_prime_pred_bar, Sigma_y_prime_bar);
    } else {
      // approximate approach
      y_prime_log = mvn_approx_lpdf(y_prime_bar| x_prime, J_prime, delta, mu_a, beta, sigma_a, sigma_y, xi_tilde);
    }
  } else {
    y_prime_log = 0;
  }
}
model {
  target += y_log;
  target += y_prime_log;

  target += normal_lpdf(to_vector(eta_a)|       0, 1);
  target += normal_lpdf(to_vector(eta_prime_a)| 0, 1);

  // assign prior to all means
  target += multi_normal_cholesky_lpdf(phi[1:K_phi_mean]| mu_phi_p[1:K_phi_mean], L_Sigma_phi_p  );
  target += multi_normal_cholesky_lpdf(delta|             mu_delta_p,             L_Sigma_delta_p);

  // assign priors to variance components on the natural scale, always
  // assume these uncorrelated
  target += normal_lpdf(sigma_a[1]| 0, sqrt(Sigma_phi_p[K_phi_mean+1,K_phi_mean+1]));
  target += normal_lpdf(sigma_a[2]| 0, sqrt(Sigma_phi_p[K_phi_mean+2,K_phi_mean+2]));
  target += normal_lpdf(sigma_y   | 0, sqrt(Sigma_phi_p[K_phi_mean+3,K_phi_mean+3]));
}

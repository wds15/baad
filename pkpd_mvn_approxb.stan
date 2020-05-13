functions {

  // transform VA into transformed VA
  real valogit(real va)      { return( logit(va) );      }
  // inverse transformed to original
  real inv_valogit(real tva) { return(inv_logit(tva));  }

  // takes sampling parameters and produces model prediction for all
  // patients
  matrix evaluate_model(vector x, int[] DRUG,
			real Lalpha_0, real Lalpha_s, real lkappa,
                        real lEmax, real delta,
                        vector eta_Lalpha_0, vector eta_lkappa,
                        real sigma_Lalpha_0, real sigma_lkappa) {
    int J = num_elements(eta_Lalpha_0);
    int T = num_elements(x);
    vector[J] Lalpha_0_j = Lalpha_0 + sigma_Lalpha_0 * eta_Lalpha_0;
    vector[J] lkappa_j   = lkappa   + sigma_lkappa * eta_lkappa;
    matrix[J,T] mu;
    vector[3] Emax_s;
    
    Emax_s[1] = 0.;
    Emax_s[2] = exp(lEmax);
    Emax_s[3] = exp(lEmax + delta);
    
    // a + (b-a) * exp(-kout*x) = a * (1 - exp(-kout*x)) + b * exp(-kout*x)
    for(j in 1:J) {
      real kout;
      real Lalpha_s_j;
      kout = exp(-lkappa_j[j]);
      Lalpha_s_j = Lalpha_s + Emax_s[DRUG[j]+1];
      for(k in 1:T)
        mu[j,k] = inv_valogit(Lalpha_s_j + (Lalpha_0_j[j]-Lalpha_s_j) * exp(-kout * x[k]));
    }

    return(mu);
  }

  void pretty_print(matrix x) {
    if (rows(x) == 0) {
      print("empty matrix");
      return;
    }
    for (m in 1:rows(x)) {
      row_vector[cols(x)] rv;
      for (n in 1:cols(x))
	rv[n] = round(1000*x[m,n])/1000.;
      print("row ", m, " = ", rv);
    }
  }

  vector colMeans(matrix y) {
    vector[cols(y)] m;
    for (i in 1:cols(y))
      m[i] = mean(col(y, i));
    return(m);
  }

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
                      int DRUG_prime,
                      vector delta,
                      real Lalpha_0, real Lalpha_s, real lkappa, real lEmax,
                      real sigma_Lalpha_0, real sigma_lkappa, real sigma_y,
                      matrix xi) {
    int J_tilde = cols(xi);
    int T_prime = num_elements(x_prime);
    vector[T_prime] y_prime_mean;
    matrix[T_prime,T_prime] y_prime_var;
    int DRUG[J_tilde] = rep_array(DRUG_prime, J_tilde);
    matrix[J_tilde,T_prime] y_prime_sim;

    y_prime_sim = evaluate_model(x_prime, DRUG,
                                 Lalpha_0, Lalpha_s, lkappa, lEmax, delta[1],
                                 to_vector(xi[1]), to_vector(xi[2]),
                                 sigma_Lalpha_0, sigma_lkappa);

    y_prime_mean = colMeans(y_prime_sim);
    // simulated cov matrix + residual error term
    y_prime_var = robust_cov(y_prime_sim, to_row_vector(y_prime_mean)) +
      diag_matrix(rep_vector(square(sigma_y), num_elements(x_prime)));

    return multi_normal_lpdf(y_prime_bar| y_prime_mean, y_prime_var/J_prime);
  }

}
data {
  int<lower=1> K;                         // #parameters in delta
  int<lower=K> K_phi;                     // #parameters in phi
  int<lower=0> J;
  int<lower=0> T;
  int<lower=0> J_prime;
  int<lower=0> T_prime;
  matrix[J,T] y;
  matrix[J_prime,T_prime] y_prime;
  vector[T] x;
  vector[T_prime] x_prime;
  int<lower=0> J_tilde;
  int<lower=0,upper=2> DRUG[J+J_prime];      // drug administered (0 for Placebo, 1 for directly fit and 2 for fitting via summaries)
  int<lower=0,upper=1> fit_all;
  int<lower=0,upper=1> fit_local;

  vector[K_phi] mu_phi_p;          // prior mean
  vector[K]     mu_delta_p;        // prior delta mean
  cov_matrix[K_phi] Sigma_phi_p;   // prior variance
  cov_matrix[K]     Sigma_delta_p; // prior delta variance
  int C;
  int<lower=1,upper=C> CHAIN_ID;   // the CHAIN_ID is passed into this block by rstan
  matrix[2, 2*J_tilde] xi[C];      // N(0,1) random numbers used for integration
}
transformed data {
  int K_phi_mean = K_phi - 3;
  cholesky_factor_cov[K_phi_mean] L_Sigma_phi_p   = cholesky_decompose(Sigma_phi_p[1:K_phi_mean,1:K_phi_mean]);
  cholesky_factor_cov[K]          L_Sigma_delta_p = cholesky_decompose(Sigma_delta_p);
  int J_tilde_chain = (CHAIN_ID % 2 == 0) ? 2*J_tilde : J_tilde;
  matrix[2, J_tilde_chain] xi_tilde = xi[CHAIN_ID, 1:2, 1:J_tilde_chain];
  vector[T_prime] y_prime_bar = colMeans(y_prime);
  int<lower=0,upper=2> DRUG_local[J] = segment(DRUG, 1, J);
  int<lower=0,upper=2> DRUG_prime[J_prime] = segment(DRUG, J+1, J_prime);

  print("Running chain ", CHAIN_ID, " with J_tilde_chain = ", J_tilde_chain, "...");
}
parameters {
  // population parameters
  real Lalpha_0;              // logit-VA latent at baseline
  real Lalpha_s;              // alpha_s = kin_s/kout
  real lkappa;                // kout = 1/exp(ltau)
  real lEmax;                 // emax for reference drug
  vector[K] delta;            // change in Emax for second drug

  real<lower=0> sigma_Lalpha_0;
  real<lower=0> sigma_lkappa;

  vector[J] eta_Lalpha_0;
  vector[J] eta_lkappa;
  vector[J_prime] eta_Lalpha_0_prime;
  vector[J_prime] eta_lkappa_prime;

  real<lower=0> sigma_y;            // prediction noise term
}
transformed parameters {
  real y_log;
  real y_prime_log;
  vector[K_phi] phi;

  phi[1] = Lalpha_0;
  phi[2] = Lalpha_s;
  phi[3] = lkappa;
  phi[4] = lEmax;
  phi[5] = log(sigma_Lalpha_0);
  phi[6] = log(sigma_lkappa);
  phi[7] = log(sigma_y);

  {
    // we do not need to store y_pred (and y_prime_pred) for each
    // iteration such that we put these into this local block
    matrix[J,T] y_pred;

    y_pred       = evaluate_model(x, DRUG_local,
                                  Lalpha_0, Lalpha_s, lkappa, lEmax, delta[1],
                                  eta_Lalpha_0, eta_lkappa,
                                  sigma_Lalpha_0, sigma_lkappa);

    y_log = normal_lpdf(to_vector(y)| to_vector(y_pred), sigma_y);

    // now get log-lik contribution for y_prime which depends on the
    // scenario
    if(!fit_local) {
      // we are including data from the prime data set
      if(fit_all == 1) {
        // by individual data...
        matrix[J_prime,T_prime] y_prime_pred;
        y_prime_pred = evaluate_model(x_prime, DRUG_prime, 
                                      Lalpha_0, Lalpha_s, lkappa, lEmax, delta[1],
                                      eta_Lalpha_0_prime, eta_lkappa_prime,
                                      sigma_Lalpha_0, sigma_lkappa);
    
        y_prime_log = normal_lpdf(to_vector(y_prime)| to_vector(y_prime_pred), sigma_y);
      } else {
        // .. or by mean data only
        y_prime_log = mvn_approx_lpdf(y_prime_bar| x_prime, J_prime, DRUG_prime[1], delta,
                                      Lalpha_0, Lalpha_s, lkappa, lEmax, 
                                      sigma_Lalpha_0, sigma_lkappa, sigma_y,
                                      xi_tilde);
      }
    } else {
      // or not at all (prime data ignored)
      y_prime_log = 0;
    }
  }

}
model {
  target += y_log;
  target += y_prime_log;

  target += normal_lpdf(eta_Lalpha_0      | 0, 1);
  target += normal_lpdf(eta_lkappa        | 0, 1);
  target += normal_lpdf(eta_Lalpha_0_prime| 0, 1);
  target += normal_lpdf(eta_lkappa_prime  | 0, 1);

  // assign prior to all means
  target += multi_normal_cholesky_lpdf(phi[1:K_phi_mean]| mu_phi_p[1:K_phi_mean], L_Sigma_phi_p  );
  target += multi_normal_cholesky_lpdf(delta            | mu_delta_p,             L_Sigma_delta_p);
  
  // assign priors to variance components on the natural scale, always
  // assume these are uncorrelated a priori
  target += normal_lpdf(sigma_Lalpha_0| 0, sqrt(Sigma_phi_p[K_phi_mean+1,K_phi_mean+1]));
  target += normal_lpdf(sigma_lkappa  | 0, sqrt(Sigma_phi_p[K_phi_mean+2,K_phi_mean+2]));
  target += normal_lpdf(sigma_y       | 0, sqrt(Sigma_phi_p[K_phi_mean+3,K_phi_mean+3]));
}
generated quantities {
}

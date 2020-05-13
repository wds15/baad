functions {
  // utility functions
  int rle_elem_count(int[] set) {
    int U;
    U = num_elements(set) == 0 ? 0 : 1;
    for(i in 2:num_elements(set)) {
      if(set[i-1] != set[i])
        U = U + 1;
    }
    return(U);
  }
  
  int[] rle_int(int[] set) {
    int res[rle_elem_count(set)];
    int c;
    res[1] = 1;
    c = 1;
    for(i in 2:num_elements(set)) {
      if(set[i-1] == set[i]) {
        res[c] = res[c] + 1;
      } else {
        c = c + 1;
        res[c] = 1;
      }
    }
    return(res);
  }

  // count number times elem appears in test set
  int count_elem(int[] test, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
      if(test[i] == elem)
	count = count + 1;
    return(count);
  }
  
  // count number times elems appears in test set
  int[] count_elems(int[] test, int[] elems) {
    int counts[num_elements(elems)];
    for(i in 1:num_elements(elems))
      counts[i] = count_elem(test, elems[i]);
    return(counts);
  }

  // find elements in test which are equal to elem
  int[] which_elem(int[] test, int elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
      if(test[i] == elem) {
	res[ci] = i;
	ci = ci + 1;
      }
    return(res);
  }

  // find indices of elements in x which corresond to the first match
  // in y (see match in R). If no match is found, set to 0 (boolean
  // false)
  int[] match_int(int[] x, int[] y) {
    int res[num_elements(x)];
    int Nx;
    int Ny;
    Nx = num_elements(x);
    Ny = num_elements(y);
    // do a brute-force seach for the (first) indices of elements in x
    // have in y
    for(i in 1:Nx) {
      int j;
      res[i] = 0;
      j = 0;
      while(j != Ny) {
	j = j + 1;
	if(y[j] == x[i]) {
	  res[i] = j;
	  j = Ny;
	}
      }
    }
    return(res);
  }

  // compare a vector of ints to a scalar and return the boolean
  // vector which is true if element in v is equal to s
  int[] int_vs_eq(int[] v, int s) {
    int res[num_elements(v)];
    for(i in 1:num_elements(v))
      res[i] = v[i] == s;
    return(res);
  }

  // given a repeated-length-element structure, return unpacked
  // representation, i.e. element S[i] is repeated N[i] times
  int[] unpack_rle_int(int[] N, int[] S) {
    int res[sum(N)];
    int k;
    k = 1;
    for(i in 1:num_elements(N)) {
      for(j in 1:N[i])
	res[k + j - 1] = S[i];
      k = k + N[i];
    }
    return(res);
  }
  
  int[] unpack_rle_id(int[] N) {
    int range[num_elements(N)];
    for(i in 1:num_elements(N)) range[i] = i;
    return(unpack_rle_int(N, range));
  }

  int[] prepare_next_counts(int[] id, int[] evid, real[] time) {
    int j;
    int count;
    int N;
    int Nnext[size(id)];
    N = size(id);
    j = -1;
    for(i in 1:N) {
      int k;
      k = N-i+1;
      if(id[k] != j) {
	// new patient, reset
	count = 0;
	j = id[k];
      }
      // new dose, reset
      if(evid[k] == 1)
	count = 0;
      // we just had a dose and the time-points are the same => make
      // this observation then count=0
      if(count == 1 && time[k] == time[k+1])
	count = 0;
      Nnext[k] = count;
      count = count + 1;
    }
    return(Nnext);
  }
  
  matrix prepare_pk(int[] id, int[] evid, int[] drug, real[] amt, real[] time) {
    vector[2] k_e;
    vector[2] kDa;
    real state_lc;
    real last_state_lc;
    real state_time;
    real lconcM[size(amt)];
    real Ke;
    int j;
    int N;
    matrix[size(id),4] lconcT;
    N = size(id);
    k_e[1] = log(2.) / 9.;      // Ranibizumab
    k_e[2] = log(2.) / 9.;      // Aflibercept
    // from below ref (kilo Daltons or g/mol)
    // Austin J Clin Ophthalmol. 2014;1(3): 1016.
    kDa[1] = 48.;               // Ranibizumab
    kDa[2] = 115.;              // Aflibercept
    state_lc = -25;
    // convert dose from mg/mL to nM
    for(i in 1:N) {
      if(evid[i] == 1 && amt[i] > 0)
        lconcM[i] = log((amt[i]/4.) * 1e3 / kDa[drug[i]]);
      else
        lconcM[i] = -25;
    }
    j = 0;
    for(i in 1:N) {
      if(id[i] != j) {
	// new subject
	if(time[i] != 0)
	  reject("Patient id = ", id[i], " does not start with a time==0 record which must be given, can be any type of record.");
	j = id[i];
	state_lc = -25;
        last_state_lc = -25;
	state_time = time[i];
	Ke = k_e[drug[j]];
      }
      lconcT[i,1] = time[i];
      lconcT[i,2] = state_lc - Ke * (time[i] - state_time);
      lconcT[i,3] = Ke;
      lconcT[i,4] = (last_state_lc + lconcT[i,2])/2;

      // in case we have a dosing event, inject it and update the
      // reference time
      if(evid[i] == 1 && amt[i] > 0) {
	state_lc = log_sum_exp(state_lc - Ke * (time[i] - state_time), lconcM[i] );
        lconcT[i,2] = state_lc;
	state_time = time[i];
      }
      
      last_state_lc = lconcT[i,2];
    }

    return(lconcT);
  }
  

  real[] drug_disease_stim_kinL_Et_ode(real t,
                                       real[] y,
                                       real[] theta,
                                       real[] x_r,
                                       int[] x_i) {
    real dydt[1];         // return array
    real y_pred;          // state: y[1]
    
    real emax_0;          // theta[1]
    real lec_50;          // theta[2]
    real r180;            // theta[3]
    real beta;            // theta[4]
    real k_in;            // theta[5]
    real k_out;           // theta[6]

    real start_t;         // x_r[1]
    real lconc_0;         // x_r[2]
    real K;               // x_r[3]
    real hill;            // x_r[4]
    
    real lconcentration;
    real emax;
    real stim;

    y_pred = y[1];
    
    emax_0 = theta[1];
    lec_50 = theta[2];
    r180 = theta[3];
    beta = theta[4];
    k_in = theta[5];
    k_out = theta[6];

    start_t = x_r[1];
    lconc_0 = x_r[2];
    K = x_r[3];
    hill = x_r[4];
    //hill = 1;

    // decay is relative to system start time
    lconcentration = lconc_0 - K * (t - start_t);
    emax = emax_0 * (r180 + (1 - r180) * exp(-beta * t / 30.));
    stim = emax * inv_logit(hill * (lconcentration - lec_50) );

    // standard stim of kin turn-over
    //dydt[1] = k_in * (1 + stim) - k_out * y_pred;

    // modified turn-over
    dydt[1] = k_in - k_out * (y_pred - stim);
    
    return dydt;
  }

  vector drug_disease_stim_kinL_Et(real[] time,
                                   real g0,          
				   data matrix lconc,
				   int[] Nnext,
                                   real[] params,
                                   data real hill) {
    vector[size(time)] gv;
    real state_g[1];
    int x_i[0];
    int i = 2;
    
    // we assume that time[1] == 0 !!! This is checked in the
    // transformed data block
    state_g[1] = g0;

    // for t=0 the function value is equal to the initial condition
    gv[1] = g0;

    // now go stepwise through data set
    while(i <= size(time)) {
      if(time[i-1] == time[i]) {
	gv[i] = state_g[1];
        i = i + 1;
      } else {
	// run integrate_ode until we find the next dose or the end of the record
	real res[Nnext[i] + 1,1];
	res = integrate_ode_rk45(drug_disease_stim_kinL_Et_ode,
                                  state_g, time[i-1],
                                  segment(time, i, Nnext[i] + 1),
                                  params,
                                  to_array_1d(append_col(lconc[i-1], hill)),
                                  x_i);

	for(l in 0:Nnext[i]) {
	  gv[i+l] = res[l+1,1];
        }
	state_g[1] = res[Nnext[i] + 1,1];
	i = i + Nnext[i];
      }
    }
    return gv;
  }

  // use step-wise analytical solution to ODE as approximation,
  // i.e. the time-changing parts are step-wise approximated as
  // constant
  vector drug_disease_stim_kinL_Et_approx(real[] time,
                                          real g0,
                                          data vector lconc,
                                          real[] theta,
                                          data real hill) {
    vector[size(time)] gv;
    real emax_0;          // theta[1]
    real lec_50;          // theta[2]
    real r180;            // theta[3]
    real beta;            // theta[4]
    real k_in;            // theta[5]
    real k_out;           // theta[6]
    real emax_cur;
   
    emax_0 = theta[1];
    lec_50 = theta[2];
    r180 = theta[3];
    beta = theta[4];
    k_in = theta[5];
    k_out = theta[6];

    // we assume that time[1] == 0 !!! This is checked in the
    // transformed data block

    // for t=0 the function value is equal to the initial condition
    gv[1] = g0;
    emax_cur = emax_0;
    
    // now go stepwise through data set
    for(i in 2:size(time)) {
      real emax_next;
      emax_next = emax_0 * (r180 + (1 - r180) * exp(-beta * time[i] / 30.));
      if(time[i-1] == time[i]) {
	gv[i] = gv[i-1];
      } else {
        // calculate the time-changing constant between current and
        // the next time-point where we want the solution; we then
        // take the average
        real stim;
        real asymp;
        stim = (emax_cur + emax_next)/2. * inv_logit(hill * (lconc[i] - lec_50) );
        // NOTE: NON-STANDARD TURN-OVER, INSTEAD of kin * (1+stim) we
        // USE (kin/kout+stim)
        //asymp = k_in * (1+stim)/k_out;
        asymp = (k_in/k_out+stim);
        gv[i] = asymp + (gv[i-1] - asymp) * exp(-k_out * (time[i] - time[i-1]));
      }
      emax_cur = emax_next;
    }
    return gv;
  }

  // transform BCVA into transformed BCVA
  real    valogit(real va)   { return( logit(va/100.) ); }
  vector Vvalogit(vector va) { return( logit(va/100.) ); }
  // inverse transformed to original
  real    inv_valogit(real   tva) { return(100.*inv_logit(tva)); }
  vector Vinv_valogit(vector tva) { return(100.*inv_logit(tva)); }

  // evaluate model for a given drug (given via EC50+Emax) and pbo
  vector evaluate_model(int use_ode,
                        int[] id, int[] M, int[] pbo, int[] grp, vector age, real[] time,
                        data matrix lconcT,
                        int[] Nnext,
                        data real hill,
                        real Lalpha_0,
                        real lalpha_s,
                        real lkout,
                        real lEC50,
                        real LEmax_0,
                        real lEmax_s,
                        real lk_Emax,
                        vector beta,
                        vector eta_Lalpha_0,
                        vector eta_lkout,
                        vector eta_lEmax_s,
                        int[]  drug,
                        vector delta) {
    int N;
    int J;
    int y_index;
    vector[size(id)] ypred;
    real k_Emax;
    real base_low;
    real base_delta;
    real valogit_ref;

    N = size(id);
    J = size(M);

    k_Emax = exp(lk_Emax);

    // truncation between 25 and 70 at baseline such that we cut from
    // 20 to 75 the mean. this corresponds to 25 to 69 on the
    // transformed scale
    base_low   = valogit(20.);
    base_delta = valogit(75.) - base_low;
    valogit_ref = valogit(55.);
    
    y_index = 1;
    for (j in 1:J) {
      // individual level parameters are set in this loop
      real g0;
      real Lg0s;
      real g0_j;
      real Lg0s_j;
      real kin_j;
      real kout_j;
      real alpha_s_j;
      real Emax_0_j;
      real Emax_s_j;
      real gamma_s_j;
      real lEmax_s_drug;
      real lEC50_drug;

      lEmax_s_drug = drug[y_index] == 1 ? lEmax_s : lEmax_s + delta[1];
      lEC50_drug   = drug[y_index] == 1 ? lEC50   : lEC50   + delta[2];

      Lg0s = Lalpha_0;
      // allow study specific intercepts in relation to the grand mean
      // and at baseline
      if(grp[y_index] > 1)
        Lg0s = Lg0s + beta[1 + grp[y_index] - 1];
      Lg0s_j = Lg0s + eta_Lalpha_0[j];
      g0   = base_low + base_delta * inv_logit(Lg0s);
      g0_j = base_low + base_delta * inv_logit(Lg0s_j);

      alpha_s_j = g0 - exp(lalpha_s);
      kout_j = exp(lkout + eta_lkout[j]);
      kin_j  = alpha_s_j * kout_j;

      Emax_s_j = exp(lEmax_s_drug + eta_lEmax_s[j] + beta[1] * (g0_j-valogit_ref));

      gamma_s_j = inv_logit(-LEmax_0);
      Emax_0_j  = Emax_s_j / gamma_s_j;

      if(pbo[y_index]) { // analytic solution for placebo (no drug) 
        for(i in 1:M[j]) {
          real va;
          va = alpha_s_j + (g0_j - alpha_s_j) * exp(-kout_j * time[y_index + i - 1]);
          ypred[y_index + i - 1] = inv_valogit(va);
        }
      } else {
        vector[M[j]] g;
        real params[6];
        params[1] = Emax_0_j;
        params[2] = lEC50_drug;
        params[3] = gamma_s_j;
        params[4] = k_Emax;
        params[5] = kin_j;
        params[6] = kout_j;
        if(use_ode) {
          g = drug_disease_stim_kinL_Et(segment(time, y_index, M[j]),
                                         g0_j,          
                                         block(lconcT, y_index, 1, M[j], 3),
                                         segment(Nnext, y_index, M[j]), 
                                         params,
                                         hill);

        } else {
          g = drug_disease_stim_kinL_Et_approx(segment(time, y_index, M[j]),
                                                g0_j,          
                                                sub_col(lconcT, y_index, 4, M[j]),
                                                params,
                                                hill);
        }
        for(i in 1:M[j])
          ypred[y_index + i - 1] = inv_valogit(g[i]);
      }
      y_index = y_index + M[j];
    }
    
    return(ypred);
  }

  // numerically robust covariance estimate with given column means
  matrix robust_cov(matrix y, vector y_col_means) {
    matrix[cols(y),cols(y)] cov;
    //matrix[rows(y),cols(y)] yc;
    
    // center column wise and take crossproduct
    //yc = y - rep_matrix(to_row_vector(y_col_means), rows(y));
    //cov = yc' * yc /(rows(y) - 1);
    cov = crossprod(y - rep_matrix(to_row_vector(y_col_means), rows(y))) /(rows(y) - 1.);
    return(0.5 * (cov + cov'));
  }

  vector colMeans(matrix y) {
    vector[cols(y)] m;
    for (i in 1:cols(y))
      m[i] = mean(col(y, i));
    return(m);
  }

  real mvn_approx_log(matrix y_prime_bar,
                      int[] J_prime,
                      int[] grp_prime, int[] sarm_prime,
                      int[] M_prime,
                      data matrix lconc_prime,
                      int[] Nnext_prime, int[] mdv_prime,
                      int[] drug_prime, vector age_mean_prime, vector age_sd_prime,
                      matrix[] xi,
                      data real hill,
                      real Lalpha_0,
                      real lalpha_s,
                      real lkout,
                      real lEC50,
                      real LEmax_0,
                      real lEmax_s,
                      real lk_Emax,
                      vector beta,
                      real sigma_Lalpha_0,
                      real sigma_lkout,
                      real sigma_lEmax_s,
                      vector delta,
                      real sigma_y) {
    int S;
    int A;
    int O;
    int J_tilde;
    //real time_prime[rows(lconc_prime)];
    int s_index;
    real log_weight;

    //time_prime = to_array_1d(col(lconc_prime, 1));

    J_tilde = rows(xi[1]);
    S = size(J_prime);
    s_index = 1;
    //G = max(grp_prime);
    A = max(sarm_prime);
    O = cols(y_prime_bar);

    log_weight = 0;
    
    // now simulate S trials and get for each trial a log_weight
    for (i in 1:S) {
      matrix[J_tilde, O] y_prime_sim;
      vector[O] y_prime_mean;
      matrix[O,O] y_prime_var;
      int M[1];
      int pbo[M_prime[i]];
      int grp[M_prime[i]];
      int id[M_prime[i]];
      int Nnext_s[M_prime[i]];
      int obs_ind[count_elem(segment(mdv_prime, s_index, M_prime[i]), 0)];
      int drug_prime_s[M_prime[i]];

      obs_ind = which_elem(segment(mdv_prime, s_index, M_prime[i]), 0);
      
      pbo = rep_array(0, M_prime[i]);
      grp = rep_array(grp_prime[s_index], M_prime[i]);
      id  = rep_array(1, M_prime[i]);
      drug_prime_s = rep_array(drug_prime[s_index], M_prime[i]);
      M[1] = M_prime[i];
      Nnext_s = segment(Nnext_prime, s_index, M_prime[i]);

      // Speed could be improved if we could use accumulators for
      // mean, variance and covariance as we could then avoid
      // allocation of the huge y_prime_sim array.
      
      // loop over each patient, one by one
      for (j in 1:J_tilde) {
        vector[1] eta_Lalpha_0_sim;
        vector[1] eta_lkout_sim;
        vector[1] eta_lEmax_s_sim;
        vector[M_prime[i]] y_prime_sim_raw;

        eta_Lalpha_0_sim[1] = sigma_Lalpha_0 * xi[i,j,1]; // take random number for respective arm i, row j and column 1
        eta_lkout_sim[1] = sigma_lkout * xi[i,j,2];       // take random number for respective arm i, row j and column 2
        eta_lEmax_s_sim[1] = sigma_lEmax_s * xi[i,j,3];
 
        y_prime_sim_raw = evaluate_model(0,
                                          id, M, pbo, grp,
                                          rep_vector( age_mean_prime[i], M_prime[i]),
                                          to_array_1d(sub_col(lconc_prime, s_index, 1, M_prime[i])),
                                          //segment(to_array_1d(col(lconc_prime, 1)), s_index, M_prime[i]),
                                          block(lconc_prime, s_index, 1, M_prime[i], 4),
                                          Nnext_s,
                                          hill,
                                          Lalpha_0,
                                          lalpha_s,
                                          lkout,
                                          lEC50,
                                          LEmax_0,
                                          lEmax_s,
                                          lk_Emax,
                                          beta,
                                          eta_Lalpha_0_sim,
                                          eta_lkout_sim,
                                          eta_lEmax_s_sim,
                                          drug_prime_s,
                                          delta);

        y_prime_sim[j] = to_row_vector(y_prime_sim_raw[obs_ind]);
      }

      y_prime_mean = colMeans(   y_prime_sim );

      // get simulated cov matrix and add in residual variance
      y_prime_var  = robust_cov( y_prime_sim, y_prime_mean ) + diag_matrix(rep_vector(square(sigma_y), O));

      log_weight = log_weight + multi_normal_log(to_vector(y_prime_bar[i]), y_prime_mean, y_prime_var/J_prime[i]);
      
      s_index = s_index + M_prime[i];
    }
    
    return log_weight;
  }

  int countValidObs(vector dv, int[] mdv, real[] time, int[] cmt) {
    int count;
    for(i in 1:size(mdv)) {
      if(mdv[i] == 0 && cmt[i] == 2)
	count = count + 1;
    }
    return count;
  }
}
data {
  int<lower=1> N;
  int<lower=0> B;
  int<lower=1,upper=N> id[N];
  
  real<lower=0> time[N];
  vector<lower=0,upper=100>[N] lidv;
  int<lower=0> cmt[N];
  int<lower=0, upper=1> mdv[N];
  int<lower=0, upper=2> evid[N];
  real<lower=0> amt[N];

  vector<lower=0, upper=150>[N] age;
  int<lower=0> grp[N];

  int<lower=0, upper=1> pbo[N];

  real<lower=0> hill;

  int<lower=0, upper=1> use_ode;
  int<lower=0, upper=1> calc_ipred;

  int<lower=1> K_phi;              // #parameters in phi
  vector[K_phi] mu_phi;            // prior mean
  vector<lower=0>[K_phi] sigma_phi;   // prior sigma

  // HEP2
  int<lower=0> S_prime;            // # of to-be included mean-study arms
  int<lower=0> T_prime;            // # of time points in external data set
  int<lower=0> K;                  // #parameters in delta
  vector[K]     mu_delta;          // prior delta mean
  vector<lower=0>[K] sigma_delta;  // prior delta sigma

  matrix[S_prime,T_prime] y_prime_bar; // mean external data
  int<lower=0> J_prime[S_prime];       // #patients per study arm external
  int<lower=0> J_tilde;                // #fake patients to simulate
  int<lower=0> N_prime;

  real<lower=0> time_prime[N_prime];
  int<lower=0, upper=1>  mdv_prime[N_prime];
  int<lower=0, upper=2> evid_prime[N_prime];
  real<lower=0> amt_prime[N_prime];
  vector<lower=0, upper=150>[S_prime] age_mean_prime;
  vector<lower=0, upper=150>[S_prime] age_sd_prime;
  int<lower=1> drug_prime[N_prime];

  int<lower=max(grp)>  grp_prime[N_prime];
  int<lower=1> sarm_prime[N_prime];

  int<lower=1> C;
  int<lower=1, upper=C> chain;
  matrix[2*J_tilde,3] xi[C,S_prime];         // N(0,1) random numbers used for integration

  // transforms to the raw scale
  vector[K_phi] scale_phi;
  vector[K_phi] location_phi;
}
transformed data {
  int J;
  int M[max(id)];
  real weight_pbo;
  real weight_trt;
  vector<lower=0>[N] odv;
  vector[count_elem(mdv[which_elem(pbo, 0)], 0)] odv_trt;
  vector[count_elem(mdv[which_elem(pbo, 1)], 0)] odv_pbo;
  int O;
  int Ot;
  int Op;
  int ind_obs[count_elem(mdv, 0)];
  int ind_obs_trt[count_elem(mdv[which_elem(pbo, 0)], 0)];
  int ind_obs_pbo[count_elem(mdv[which_elem(pbo, 1)], 0)];
  int drug[N];
  vector[N] lconcM;            // log-doses divided by V in nM/mL for all patients, V=4
  int Nnext[N];                // counts for each entry the number of data lines until the next dose event
  matrix[N,4] lconcT;          // log-concentration and time at all times
  int G;                       // number of different studies (MARINA is the reference and does not count)
  int P;
  int M_prime[S_prime];
  int Nnext_prime[N_prime];    // counts for each entry the number of data lines until the next dose event
  matrix[N_prime,4] lconc_prime;// log-concentration and time at all times
  int J_tilde_chain;
  matrix[chain % 2 == 0 ? J_tilde : 2*J_tilde, 3] xi_chain[S_prime];
  
  // total number of observations
  O = count_elem(mdv, 0);
  ind_obs = which_elem(mdv, 0);
  odv = rep_vector(0, N);

  {
    int i_trt;
    int i_pbo;
    
    // grab elements from dv which will enter the likelihood
    i_trt = 1;
    i_pbo = 1;

    for (i in 1:O) {
      odv[ind_obs[i]] = lidv[ind_obs[i]];
      if (pbo[ind_obs[i]] == 0) {
        ind_obs_trt[i_trt] = ind_obs[i];
        i_trt = i_trt + 1;
      } else {
        ind_obs_pbo[i_pbo] = ind_obs[i];
        i_pbo = i_pbo + 1;
      }
    }
  }

  Ot = num_elements(ind_obs_trt);
  Op = num_elements(ind_obs_pbo);

  // take data on the original scale
  odv_trt = odv[ind_obs_trt];
  odv_pbo = odv[ind_obs_pbo];

  // number of patients
  J = max(id);
  print("no of patients J = ", J);
  
  // build up auxilary ragged array structure indices
  // contains the number of entries in the data vector for a given
  // patient
  M = rle_int(id);
  print("avg no of records per patient = ", mean(to_vector(M)));

  // count the number of infusions and IV bolus events for reporting
  print("avg no/dose of IV boluses per patient: ", 1.0*count_elem(evid, 1)/J, " / ", sum(amt)/count_elem(evid, 1) );

  // the internal data set is only Lucentis
  drug = rep_array(1, N);
  
  // calculate weights for treatment data and for pbo data to adjust
  // for the hugely different sample sizes
  {
    real Nt;
    real Np;
    Nt = rle_elem_count(id[ind_obs_trt]); // # of treated patients
    Np = rle_elem_count(id[ind_obs_pbo]); // # of control patients
    weight_trt = sqrt((Nt+Np)/Nt * 0.5);
    weight_pbo = sqrt((Nt+Np)/Np * 0.5);

    weight_pbo = weight_pbo / weight_trt;
    weight_trt = 1.;

    print("# of treated patients: ", Nt);
    print("# of control patients: ", Np);
    print("weighting treatment log-lik with ", square(weight_trt));
    print("weighting control   log-lik with ", square(weight_pbo));
  }

  G = max(grp);

  P = calc_ipred ? O : 0;

  lconcT = prepare_pk(id, evid, rep_array(1, N), amt, time);

  // prepare a vector which counts for each patient always how many
  // observations the next dose is away, doing so by counting from the
  // back; if at the same time point first an observation and then a
  // dosing is scheduled, the counting is such that the dose is at the
  // observation (to stop the integrator there)
  Nnext = prepare_next_counts(id, evid, time);

  if(K_phi != 12+B)
    reject("K_phi must be ", 12+B, " !");

  // prepare data structures for external data
  if(S_prime != 0) {
    M_prime = rle_int(sarm_prime);
    Nnext_prime = prepare_next_counts(sarm_prime, evid_prime, time_prime);
    lconc_prime = prepare_pk(sarm_prime, evid_prime, drug_prime, amt_prime, time_prime);
  
    print("Running chain ", chain, "...");

    J_tilde_chain = chain % 2 == 0 ? J_tilde : 2*J_tilde;

    print("Using ", J_tilde_chain, " fake patients.");
  
    xi_chain = xi[chain, 1:S_prime, 1:J_tilde_chain];
  }
}
parameters {
  vector[K_phi] phi_raw;
  vector[J] xi_Lalpha_0;      // log-BCVA latent at baseline
  vector[J] xi_lkout;
  vector[J] xi_lEmax_s;
  vector[K] delta;
}
transformed parameters {
  vector[K_phi] phi;
  real ldva_m;
  real ldva_r;
  real Lalpha_0;              // log-BCVA latent at baseline
  real lalpha_s;              // alpha_s = kin_s/kout = g0 - exp(lalpha_s)  per subject
  real lkout;                 // log(kout)
  real lEC50;                 // log EC50
  real LEmax_0;               // emax as delta BCVA from in ss
  real lEmax_s;               // logit for Emax_0 * frac = Emax_s
  real lk_Emax; 
  vector[B] beta;             // regression coefficiencts
  real sigma_eta_Lalpha_0;  
  real sigma_eta_lkout;
  real sigma_eta_lEmax_s;
  real sigma_y_pbo;            // noise term, placebo
  real sigma_y_trt;            // noise term, treatment

  // shift and scale such that raw parameters are 0 centered and have
  // scale 1
  phi = location_phi + scale_phi .* phi_raw;

  Lalpha_0 = phi[1];
  ldva_m = phi[2];
  lkout = mu_phi[3];
  lEC50 = phi[4];
  LEmax_0 = phi[5];
  ldva_r = phi[6];
  lalpha_s = ldva_m - ldva_r;
  lEmax_s  = ldva_m + ldva_r;
  lk_Emax = phi[7];
  beta = phi[8:(8+B-1)];
  sigma_eta_Lalpha_0 = exp(phi[8+B]);
  sigma_eta_lkout = exp(phi[8+B+1]);
  sigma_eta_lEmax_s = exp(phi[8+B+2]);
  sigma_y_trt = exp(phi[8+B+3]);
  sigma_y_pbo = exp(phi[8+B+4]);
}
model {
  vector[N] ypred;
  vector[Ot] ypred_obs_trt;
  vector[Op] ypred_obs_pbo;
  vector[J] eta_Lalpha_0;
  vector[J] eta_lkout;
  vector[J] eta_lEmax_s;

  // log-normals on the means
  phi[1:(8+B-1)]   ~ normal(mu_phi[1:(8+B-1)], sigma_phi[1:(8+B-1)] );
  // half-normals on the variance components
  exp(phi[(8+B):K_phi]) ~ normal(0, sigma_phi[(8+B):K_phi] );
  // add Jacobian correction since we apply the prior on the natural
  // scale, but sample on the log-scale for variance components
  target += phi[(8+B):K_phi];
  
  delta ~ normal(mu_delta, sigma_delta);

  xi_Lalpha_0 ~ normal(0., 1.);
  xi_lkout    ~ normal(0., 1.);
  xi_lEmax_s  ~ normal(0., 1.);

  eta_Lalpha_0 = xi_Lalpha_0 * sigma_eta_Lalpha_0;
  eta_lkout = xi_lkout * sigma_eta_lkout;
  eta_lEmax_s = xi_lEmax_s * sigma_eta_lEmax_s;

  ypred = evaluate_model(use_ode, id, M, pbo, grp, age, time,
                          lconcT, Nnext,
                          hill,
                          Lalpha_0,
                          lalpha_s,
                          lkout,
                          lEC50,
                          LEmax_0,
                          lEmax_s,
                          lk_Emax,
                          beta,
                          eta_Lalpha_0,
                          eta_lkout,
                          eta_lEmax_s,
                          drug,
                          delta
                          );

  ypred_obs_trt = ypred[ind_obs_trt];
  ypred_obs_pbo = ypred[ind_obs_pbo];

  // vectorized likelihood
  odv_trt ~ normal(ypred_obs_trt, sigma_y_trt);
  odv_pbo ~ normal(ypred_obs_pbo, sigma_y_pbo);

  // likelihood for summary data
  if(S_prime != 0) {
    target += mvn_approx_log(y_prime_bar, J_prime, grp_prime, sarm_prime,
                             M_prime, lconc_prime, Nnext_prime,
                             mdv_prime, drug_prime, age_mean_prime, age_sd_prime,
                             xi_chain,
                             hill,
                             Lalpha_0,
                             lalpha_s,
                             lkout,
                             lEC50,
                             LEmax_0,
                             lEmax_s,
                             lk_Emax,
                             beta,
                             sigma_eta_Lalpha_0,
                             sigma_eta_lkout,
                             sigma_eta_lEmax_s,
                             delta,
                             sigma_y_trt
                             );
  }
  
}
generated quantities {
  vector[P] ipred;
  real gss;
  real gss_t;
  real kappa;
  real lambda;
  real Emax_eff;
    
  gss    = inv_valogit(valogit(55) - exp(lalpha_s));
  gss_t  = inv_valogit(valogit(55) - exp(lalpha_s) + exp(lEmax_s));
  kappa  = exp(-lkout)/365;
  lambda = inv_logit(-LEmax_0);
  Emax_eff = gss_t/gss - 1.;
  
  if(calc_ipred) {
    vector[N] ypred;
    vector[J] eta_Lalpha_0;
    vector[J] eta_lkout;
    vector[J] eta_lEmax_s;
    
    eta_Lalpha_0 = xi_Lalpha_0 * sigma_eta_Lalpha_0;
    eta_lkout = xi_lkout * sigma_eta_lkout;
    eta_lEmax_s = xi_lEmax_s * sigma_eta_lEmax_s;
    
    ypred = evaluate_model(use_ode, id, M, pbo, grp, age, time,
                            lconcT, Nnext,
                            hill,
                            Lalpha_0,
                            lalpha_s,
                            lkout,
                            lEC50,
                            LEmax_0,
                            lEmax_s,
                            lk_Emax,
                            beta,
                            eta_Lalpha_0,
                            eta_lkout,
                            eta_lEmax_s,
                            drug,
                            delta);
    
    ipred = ypred[ind_obs];
  }
}

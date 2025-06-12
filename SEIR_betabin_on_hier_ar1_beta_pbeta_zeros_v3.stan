data {
  int<lower=1> N_C;
  int<lower=0> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT-1,N_C] real<lower=0> ii;
  array[N_C] int<lower=0> first;

  }
  
parameters {

  matrix[TT, N_C] u_t_logit_eta;
  matrix[TT, N_C] w_t_logit_eta;
  matrix[TT, N_C] v_t_logit_eta;
  matrix[TT, N_C] raw_log_beta_mat;


  //real<lower=-1,upper=1> phi;
  vector[TT] z;
  real<lower=0> sigma; 
  real<lower=0> sig_beta;
  real p_raw;  // unconstrained
  //real phi_p_raw;  // concentration parameter for beta-binomial approximation
  real<lower=0> phi_p; // concentration parameter for beta-binomial approximation
  array[N_C] real i0_raw;  
  array[N_C] real eta_raw; // unconstrained
  array[N_C] real gamma_raw; // unconstrained
  real v_raw; // unconstrained
  real rho_ir_raw; // Spatial range parameter
  real rho_ei_raw; // Spatial range parameter
}

transformed parameters {
  matrix[TT, N_C] s_t;
  matrix[TT, N_C] se_t;
  matrix[TT, N_C] ei_t;
  matrix[TT, N_C] i_t;
  matrix[TT, N_C] e_t = rep_matrix(0,TT,N_C);
  matrix[TT, N_C] ir_t;
  matrix[TT, N_C] u_t;
  matrix[TT, N_C] v_t;
  matrix[TT, N_C] w_t;
  vector[TT] log_beta;
  matrix[TT, N_C] log_beta_mat;
  matrix[TT, N_C] beta_mat;
  real<lower=0,upper=0> rho_se;

  //real<lower=0> phi_p = exp(phi_p_raw);
  real u_t_mean;
  real v = inv_logit(v_raw);
  real phi = 2 * v - 1;
  real p = inv_logit(p_raw);
  real rho_ei = inv_logit(rho_ei_raw);
  real rho_ir = inv_logit(rho_ir_raw);
  array[N_C] real<lower=0,upper=1> i0 = inv_logit(i0_raw);
  array[N_C] real<lower=0,upper=1> gamma = inv_logit(gamma_raw);
  array[N_C] real<lower=0,upper=1> eta = inv_logit(eta_raw);
  rho_se = 0;


  se_t=rep_matrix(0,TT,N_C);
  ir_t=rep_matrix(0,TT,N_C);
  ei_t=rep_matrix(0,TT,N_C);

  for (ct in 1:N_C){
  i_t[1:(first[ct]-1),ct] = rep_vector(0,first[ct]-1);
  i_t[first[ct],ct]=i0[ct];

  s_t[1:(first[ct]-1),ct] = rep_vector(1,first[ct]-1);
  s_t[first[ct],ct] = 1-i0[ct];

  ei_t[first[ct],ct] = i0[ct];
  }
 

  log_beta[1] = sigma / sqrt(1 - phi^2) * z[1];
  for (n in 2:TT) {
    log_beta[n] = phi * log_beta[n - 1] + sigma * z[n];
  }
  for (t in 1:TT) {
    for (ct in 1:N_C) {
      log_beta_mat[t, ct] = log_beta[t] + sig_beta * raw_log_beta_mat[t, ct];
      beta_mat[t, ct] = exp(log_beta_mat[t, ct]);

    }
  }


  for (ct in 1:N_C) {
    for (n in 2:TT) {
      if (n > first[ct]) {
            // Begin conditional block
            //print("n: ", n, "ct: ", ct, "i_t[n-1,ct]: ", i_t[n-1,ct], "e_t[n-1,ct]: ", e_t[n-1,ct], "s_t[n-1,ct]: ", s_t[n-1,ct]);
            u_t_mean = exponential_cdf(beta_mat[n - 1,ct] * i_t[n - 1, ct] | 1);
            u_t[n-1, ct] = inv_logit((u_t_logit_eta[n-1, ct]) * (sqrt((1 - rho_se) / (pop_size[ct] * s_t[n - 1, ct] * u_t_mean * (1 - u_t_mean)) + rho_se / (u_t_mean * (1 - u_t_mean)))) + logit(u_t_mean));

            if (ei_t[n-1, ct] == 0) {
              v_t[n-1, ct] = inv_logit((v_t_logit_eta[n-1, ct]) * (sqrt((1 - rho_ei) / (pop_size[ct] * e_t[n - 1, ct] * eta[ct] * (1 - eta[ct])) + rho_ei / (eta[ct] * (1 - eta[ct])))) + logit(eta[ct]));
              ei_t[n-1, ct] = v_t[n-1, ct] * e_t[n - 1, ct];
            }

            w_t[n-1, ct] = inv_logit((w_t_logit_eta[n-1, ct]) * (sqrt((1 - rho_ir) / (pop_size[ct] * i_t[n - 1, ct] * gamma[ct] * (1 - gamma[ct])) + rho_ir / (gamma[ct] * (1 - gamma[ct])))) + logit(gamma[ct]));
            se_t[n-1, ct] = u_t[n-1, ct] * s_t[n - 1, ct];
            
            ir_t[n-1, ct] = w_t[n-1, ct] * i_t[n - 1, ct];
            s_t[n, ct] = s_t[n - 1, ct] - se_t[n-1, ct];
            if (n > first[ct]+1) {
              e_t[n, ct] = e_t[n - 1, ct] + se_t[n-1, ct] - ei_t[n-1, ct];
            } else {
              e_t[n, ct] = e_t[n - 1, ct] + se_t[n-1, ct];
            }
            i_t[n, ct] = i_t[n - 1, ct] + ei_t[n-1, ct] - ir_t[n-1, ct];
          } 
      }
    }
  }

model {
  // Priors
  rho_ir_raw ~ logistic(0,1);
  rho_ei_raw ~ logistic(0,1); 
  v_raw ~ normal(0,1); //logistic(0,1);

  sigma ~ normal(0, 3); //gamma(1,.01);//
  sig_beta ~ normal(0, 3); //gamma(1,.01);//normal(0, 5);
  phi_p ~ cauchy(0, 10); //gamma(1,.01);//normal(0, 5); implied uniform over [0, inf]
  //phi_p_raw ~ student_t(3, 0, 2.5); // normal(0,3); // concentration parameter for beta-binomial approximation
  i0_raw ~  logistic(0,1); 
  p_raw ~ logistic(0,1);

  gamma_raw  ~ logistic(0,1); 
  eta_raw ~ logistic(0,1);

  to_vector(z) ~ std_normal();
  to_vector(u_t_logit_eta) ~ std_normal();
  to_vector(v_t_logit_eta) ~ std_normal();
  to_vector(w_t_logit_eta) ~ std_normal();
  to_vector(raw_log_beta_mat) ~ std_normal();

  for (ct in 1:N_C) {
    for (i in 1:(TT-1)) {
      if (i >= first[ct]) {
        {
          real n_ii = pop_size[ct] * ei_t[i,ct];
          real mu_ii = p * n_ii;
          real var_ii = n_ii * p * (1 - p) * (n_ii + phi_p) / (1 + phi_p); 
          ii[i,ct] ~ normal(mu_ii, sqrt(var_ii));
        }
      }
    }
  }
}


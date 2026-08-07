data {
  int<lower=1> N_C;
  int<lower=0> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT-1,N_C] int<lower=0> ii;
  array[N_C] int<lower=0> first;
  array[N_C] int<lower=0> last;
  int<lower=1> N_wks_per_period;
}
transformed data {
  array[N_C] int dim_beta_by_c;
  for (c in 1:N_C)
    dim_beta_by_c[c] = last[c] - first[c] + 1;
  int n_beta = sum(dim_beta_by_c);
  int glob_beta_start = min(first);
  int n_glob_beta = TT - glob_beta_start + 1;
  array[N_C] int idx_starts;
  idx_starts[1] = 1;
  for (c in 2:N_C)
    idx_starts[c] = idx_starts[c-1] + dim_beta_by_c[c-1];
  array[N_C] int idx_ends;
  for (c in 1:N_C)
    idx_ends[c] = idx_starts[c] + dim_beta_by_c[c] - 1;
  matrix[TT-1, N_C] pp_real = to_matrix(ii) * diag_matrix(pop_size);
  array[TT] int p_idx;
  for (t in 0:(TT - 1)) {
    p_idx[t+1] = t %/% N_wks_per_period + 1;
  }
  int n_p_pers = max(p_idx);
}
parameters {
  matrix[TT, N_C] u_t_logit_eta;
  matrix[TT, N_C] w_t_logit_eta;
  matrix[TT, N_C] v_t_logit_eta;
  vector[n_beta] raw_log_beta;

  real mu_log_beta;
  vector[n_glob_beta] z;
  real<lower=0> sigma; 
  vector<lower=0>[N_C] sig_beta;
  vector[n_p_pers] eta_p_time;  // unconstrained
  vector[N_C] log_e0_count;
  array[N_C] real eta_raw; // unconstrained
  array[N_C] real gamma_raw; // unconstrained
  real v_raw; // unconstrained
  // real rho_ir_raw; // Spatial range parameter
  // real rho_ei_raw; // Spatial range parameter
}

transformed parameters {
  matrix[TT, N_C] s_t;
  matrix[TT, N_C] se_t;
  matrix[TT, N_C] ei_t;
  matrix[TT, N_C] i_t = rep_matrix(0,TT,N_C);
  matrix[TT, N_C] e_t;
  matrix[TT, N_C] ir_t;
  matrix[TT, N_C] u_t;
  matrix[TT, N_C] v_t;
  matrix[TT, N_C] w_t;
  vector[TT] log_beta;
  matrix[TT, N_C] log_beta_mat;
  matrix[TT, N_C] beta_mat;
  //real<lower=0,upper=0> rho_ei;
  //real<lower=0,upper=0> rho_ir;

  real u_t_mean;
  real v = inv_logit(v_raw);
  real phi = 2 * v - 1;
  matrix[TT, N_C] p;
  real rho_ei = 0; //inv_logit(rho_ei_raw);
  real rho_ir = 0; //inv_logit(rho_ir_raw);
  array[N_C] real<lower=0,upper=1> gamma = inv_logit(gamma_raw);
  array[N_C] real<lower=0,upper=1> eta = inv_logit(eta_raw);
  vector<lower=0, upper=1>[N_C] e0 = exp(log_e0_count) ./ pop_size;
  u_t=rep_matrix(0,TT,N_C);
  v_t=rep_matrix(0,TT,N_C);
  w_t=rep_matrix(0,TT,N_C);
  se_t=rep_matrix(0,TT,N_C);
  ir_t=rep_matrix(0,TT,N_C);
  ei_t=rep_matrix(0,TT,N_C);

 //  i_t[n, ct] = i_t[n - 1, ct] + ei_t[n-1, ct] - ir_t[n-1, ct];


  for (ct in 1:N_C){
    e_t[1:(first[ct]-2),ct] = rep_vector(0,first[ct]-2);
    e_t[first[ct]-1,ct]=e0[ct];

    s_t[1:(first[ct]-2),ct] = rep_vector(1,first[ct]-2);
    s_t[first[ct]-1,ct] = 1-e0[ct];
  }
 
  for (t in 1:TT) {
    for (c in 1:N_C) {
      p[t,c] = inv_logit(eta_p_time[p_idx[t]]);
    }
  }
  

  log_beta[glob_beta_start] = mu_log_beta
                + sigma / sqrt(1 - square(phi)) * z[1];
  {  
    int dummy = 2;
  for (n in (glob_beta_start+1):TT) {
    log_beta[n] = mu_log_beta
                  + phi * (log_beta[n-1] - mu_log_beta)
                  + sigma * z[dummy];
    dummy += 1;
  }
  }
  for (ct in 1:N_C) {
    log_beta_mat[first[ct]:last[ct], ct] = log_beta[first[ct]:last[ct]] 
                                           + sig_beta[ct] * raw_log_beta[idx_starts[ct]:idx_ends[ct]];
    beta_mat[first[ct]:last[ct], ct] = exp(log_beta_mat[first[ct]:last[ct], ct]);
  }


  for (ct in 1:N_C) {
    for (n in 2:TT) {
      if (n >= first[ct] && n <= (last[ct]+1)) {
            // Begin conditional block
            //print("n: ", n, "ct: ", ct, "i_t[n-1,ct]: ", i_t[n-1,ct], "e_t[n-1,ct]: ", e_t[n-1,ct], "s_t[n-1,ct]: ", s_t[n-1,ct]);
            v_t[n, ct] = inv_logit((v_t_logit_eta[n, ct]) * (sqrt((1 - rho_ei) / (pop_size[ct] * e_t[n - 1, ct] * eta[ct] * (1 - eta[ct])) + rho_ei / (eta[ct] * (1 - eta[ct])))) + logit(eta[ct]));
            ei_t[n, ct] = v_t[n, ct] * e_t[n - 1, ct];


            if (n > first[ct]) {
              u_t_mean = exponential_cdf(beta_mat[n - 1,ct] * i_t[n - 1, ct] | 1);
              u_t[n, ct] = inv_logit((u_t_logit_eta[n, ct]) * sqrt(inv(pop_size[ct] * s_t[n - 1, ct] * u_t_mean * (1 - u_t_mean))) + logit(u_t_mean));
                
              w_t[n, ct] = inv_logit((w_t_logit_eta[n, ct]) * (sqrt((1 - rho_ir) / (pop_size[ct] * i_t[n - 1, ct] * gamma[ct] * (1 - gamma[ct])) + rho_ir / (gamma[ct] * (1 - gamma[ct])))) + logit(gamma[ct]));
              se_t[n, ct] = u_t[n, ct] * s_t[n - 1, ct];
              ir_t[n, ct] = w_t[n, ct] * i_t[n - 1, ct];
            }

            
            s_t[n, ct] = s_t[n - 1, ct] - se_t[n, ct];
            e_t[n, ct] = e_t[n - 1, ct] + se_t[n, ct] - ei_t[n, ct];
            i_t[n, ct] = i_t[n - 1, ct] + ei_t[n, ct] - ir_t[n, ct];
          } 
      }
    }
  }

model {
  // Priors
  // rho_ir_raw ~ logistic(0, 1);
  // rho_ei_raw ~ logistic(0, 1);
  v_raw ~ normal(logit(0.9), 0.3);//normal(0,1); //logistic(0,1);
  mu_log_beta ~ normal(0,2);
  sigma ~ normal(0,1); // normal(log(.1), .5); //normal(0,.1);//
  sig_beta ~ normal(0,2); // normal(log(0.03), 0.5);
  log_e0_count ~ normal(log(20), 1.5);
  gamma_raw ~ normal(logit(0.75), 0.5);
  eta_p_time ~ logistic(0, 1);
  
  eta_raw ~ normal(logit(0.75), 1); 

  z ~ std_normal();
  to_vector(u_t_logit_eta) ~ std_normal();
  to_vector(v_t_logit_eta) ~ std_normal();
  to_vector(w_t_logit_eta) ~ std_normal();
  raw_log_beta ~ std_normal();

  for (ct in 1:N_C) {
    for (i in 1:(TT-1)) {
      if (i >= first[ct] && i <= last[ct]) 
        ii[i,ct] ~ poisson(p[i,ct] * ei_t[i,ct] * pop_size[ct]);
    }
  }
}
generated quantities {
  real average_p = sum(to_array_1d(ii)) / sum(ei_t * diag_matrix(pop_size)); 
}

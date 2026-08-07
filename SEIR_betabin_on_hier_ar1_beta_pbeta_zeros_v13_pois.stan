data {
  int<lower=1> N_C;
  int<lower=3> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT - 1, N_C] int<lower=0> ii;
  int<lower=2, upper=TT - 1> model_start;
  int<lower=1> N_wks_per_period;
}

transformed data {
  int n_beta_per_county = TT - model_start;
  int n_beta = N_C * n_beta_per_county;
  int n_glob_beta = TT - model_start + 1;
  array[N_C] int idx_starts;
  array[N_C] int idx_ends;
  array[TT - 1] int p_idx;

  for (c in 1:N_C) {
    idx_starts[c] = (c - 1) * n_beta_per_county + 1;
    idx_ends[c] = c * n_beta_per_county;
  }
  for (t in 1:(TT - 1)) {
    p_idx[t] = (t - 1) %/% N_wks_per_period + 1;
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
  vector[n_p_pers] eta_p_time;

  real mu_log_e0;
  real<lower=0> sd_log_e0;
  vector[N_C] z_e0;

  array[N_C] real eta_raw;
  array[N_C] real gamma_raw;
  real v_raw;
}

transformed parameters {
  matrix[TT, N_C] s_t = rep_matrix(1, TT, N_C);
  matrix[TT, N_C] e_t = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] i_t = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] se_t = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] ei_t = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] ir_t = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] u_t = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] v_t = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] w_t = rep_matrix(0, TT, N_C);
  vector[TT] log_beta = rep_vector(0, TT);
  matrix[TT, N_C] log_beta_mat = rep_matrix(0, TT, N_C);
  matrix[TT, N_C] beta_mat = rep_matrix(0, TT, N_C);
  matrix[TT - 1, N_C] p;

  real v = inv_logit(v_raw);
  real phi = 2 * v - 1;
  real rho_ei = 0;
  real rho_ir = 0;
  array[N_C] real<lower=0, upper=1> gamma = inv_logit(gamma_raw);
  array[N_C] real<lower=0, upper=1> eta = inv_logit(eta_raw);

  vector[N_C] log_e0_count = mu_log_e0 + sd_log_e0 * z_e0;
  vector<lower=0, upper=1>[N_C] e0 = exp(log_e0_count) ./ pop_size;

  for (ct in 1:N_C) {
    e_t[model_start - 1, ct] = e0[ct];
    s_t[model_start - 1, ct] = 1 - e0[ct];
  }

  for (t in 1:(TT - 1)) {
    for (ct in 1:N_C) {
      p[t, ct] = inv_logit(eta_p_time[p_idx[t]]);
    }
  }

  // Centered, stationary AR(1) for the statewide transmission trajectory.
  log_beta[model_start] = mu_log_beta
    + sigma / sqrt(1 - square(phi)) * z[1];
  for (n in (model_start + 1):TT) {
    log_beta[n] = mu_log_beta
      + phi * (log_beta[n - 1] - mu_log_beta)
      + sigma * z[n - model_start + 1];
  }

  for (ct in 1:N_C) {
    log_beta_mat[model_start:(TT - 1), ct] =
      log_beta[model_start:(TT - 1)]
      + sig_beta[ct] * raw_log_beta[idx_starts[ct]:idx_ends[ct]];
    beta_mat[model_start:(TT - 1), ct] =
      exp(log_beta_mat[model_start:(TT - 1), ct]);
  }

  for (ct in 1:N_C) {
    for (n in model_start:TT) {
      real u_t_mean;

      v_t[n, ct] = inv_logit(
        v_t_logit_eta[n, ct]
        * sqrt((1 - rho_ei)
               / (pop_size[ct] * e_t[n - 1, ct]
                  * eta[ct] * (1 - eta[ct]))
               + rho_ei / (eta[ct] * (1 - eta[ct])))
        + logit(eta[ct]));
      ei_t[n, ct] = v_t[n, ct] * e_t[n - 1, ct];

      if (n > model_start) {
        u_t_mean = exponential_cdf(beta_mat[n - 1, ct]
                                   * i_t[n - 1, ct] | 1);
        u_t[n, ct] = inv_logit(
          u_t_logit_eta[n, ct]
          * sqrt(inv(pop_size[ct] * s_t[n - 1, ct]
                     * u_t_mean * (1 - u_t_mean)))
          + logit(u_t_mean));

        w_t[n, ct] = inv_logit(
          w_t_logit_eta[n, ct]
          * sqrt((1 - rho_ir)
                 / (pop_size[ct] * i_t[n - 1, ct]
                    * gamma[ct] * (1 - gamma[ct]))
                 + rho_ir / (gamma[ct] * (1 - gamma[ct])))
          + logit(gamma[ct]));
        se_t[n, ct] = u_t[n, ct] * s_t[n - 1, ct];
        ir_t[n, ct] = w_t[n, ct] * i_t[n - 1, ct];
      }

      s_t[n, ct] = s_t[n - 1, ct] - se_t[n, ct];
      e_t[n, ct] = e_t[n - 1, ct] + se_t[n, ct] - ei_t[n, ct];
      i_t[n, ct] = i_t[n - 1, ct] + ei_t[n, ct] - ir_t[n, ct];
    }
  }
}

model {
  v_raw ~ normal(logit(0.9), 0.3);
  mu_log_beta ~ normal(0, 2);
  sigma ~ normal(0, 1);
  sig_beta ~ normal(0, 2);

  // Partial pooling is on the initial exposed count scale.
  mu_log_e0 ~ normal(log(20), 2);
  sd_log_e0 ~ normal(0, 1);
  z_e0 ~ std_normal();

  gamma_raw ~ normal(logit(0.75), 0.5);
  eta_p_time ~ logistic(0, 1);
  eta_raw ~ normal(logit(0.75), 1);

  z ~ std_normal();
  to_vector(u_t_logit_eta) ~ std_normal();
  to_vector(v_t_logit_eta) ~ std_normal();
  to_vector(w_t_logit_eta) ~ std_normal();
  raw_log_beta ~ std_normal();

  for (ct in 1:N_C) {
    for (t in model_start:(TT - 1)) {
      ii[t, ct] ~ poisson(p[t, ct] * ei_t[t, ct] * pop_size[ct]);
    }
  }
}

generated quantities {
  real expected_detectable = 0;
  real observed_cases = 0;
  real average_p;

  for (ct in 1:N_C) {
    for (t in model_start:(TT - 1)) {
      observed_cases += ii[t, ct];
      expected_detectable += ei_t[t, ct] * pop_size[ct];
    }
  }
  average_p = observed_cases / expected_detectable;
}

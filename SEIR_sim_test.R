library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(MASS)
library(VGAM)


N_C <- 11
TT <- 30



pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
#pop=read_csv("./utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
pop_size = pop$Population_2020[counties]





ii=readRDS("div_dat_ex_1008.rds")#[[2]]

first=apply(matrix(as.numeric(as.matrix(ii)!=0),ncol=N_C),2,function(x) which(x==1)[1])


last<- unlist(1:N_C %>% purrr::map(function(col) {
  w <- which(ii[,col] == 0)
  w=w[which(w>first[col])][1]-1
  w=ifelse(is.na(w),TT-1,w)
}))

tt <- cmdstan_model("SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v7.stan")
dat <- 
  list(
    ii = ii,#as.matrix(d1)[,counties],
    TT = TT,
    N_C = N_C,
    pop_size = pop_size,
    #D=dist[counties,counties]/10,
    first=first,
    last=last,
    min_first=min(first)
  )

init=\() {list(u_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  v_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  w_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  raw_log_beta_mat = matrix(runif(TT*N_C, -1, 1), TT, N_C),#matrix(rnorm(TT*N_C, 0, 1), TT, N_C),
                                  p_raw = runif(1, -.25,.25),
                                  #kappa = runif(1, 0,1),
                                  phi_p = runif(1, 100, 200),
                                  v_raw=runif(1,1.5,2.5),
                                  z=rnorm(TT,0,.25),
                                  mu_log_beta = rnorm(1, 0, .05),
                                  sigma = runif(1, .1, .25),
                                  sig_beta = runif(1, .05, .1),
                                  i0_raw = runif(N_C,-3,-1),#rbeta(N_C, 0.01*50, 0.99*50),
                                  #rho_si = runif(1, 0.0001, 0.005),
                                  rho_ei_raw = runif(1, -0.1, 0.1),
                                  rho_ir_raw = runif(1, -0.1, 0.1),
                                  gamma_raw = runif(N_C, 0,.5),
                                  eta_raw = runif(N_C,.25,.75))}
#saveRDS(dat,file="dat2.rds")
#dat=readRDS("dat2.rds")
#seed(123)
init=list(init(),init(),init(),init(),init(),init(),init(),init(),init(),init())
#inits <- replicate(15, init_fn(), simplify = FALSE)

fit = tt$sample(data = dat, chains = 10,
                 adapt_delta = 0.99,
                 max_treedepth = 16,
                 init = init,
                 iter_warmup = 1500,
                 iter_sampling = 1000, parallel_chains = 10)#,


#fit$summary("p")

diagn=fit$diagnostic_summary();diagn

init_vals <- fit$init()

init_tbl <- tibble(
  chain     = seq_along(init_vals),
  sigma     = map_dbl(init_vals, "sigma"),
  sig_beta  = map_dbl(init_vals, "sig_beta"),
  phi_p     = map_dbl(init_vals, "phi_p"),
  p_raw     = map_dbl(init_vals, "p_raw"),
  v_raw     = map_dbl(init_vals, "v_raw"),
  mu_log_beta = map_dbl(init_vals, "mu_log_beta"),
  beta_min = unlist(map(init_vals, function(x) min(x[["raw_log_beta_mat"]]))),
  beta_max = unlist(map(init_vals, function(x) max(x[["raw_log_beta_mat"]])))
  # etc… pull out any other scalars you care about
)
init_tbl %>% mutate(div=diagn$num_divergent)


np_fit <- nuts_params(fit)
mcmc_pairs(fit$draws(c("sigma","sig_beta","phi_p","p_raw","eta","gamma")), np = np_fit, 
           pars = c("sigma","sig_beta","phi_p","p_raw","eta[1]","gamma[1]"),
            off_diag_args = list(size = 0.75),condition=pairs_condition(chains=list(1:5,6:10)))
 


library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(MASS)
library(VGAM)

iter=as.numeric(commandArgs(T))
#2,500
#1,1000

chn=2#2#as.numeric(commandArgs(T))
fiter=sample(1:500,1)

#set.seed(123)
#beta <- 1.05
#TT <- 100

#counties=c(1,2,3,8)#c(1,2,3,4,5,6,7,8,10,11,12) # dist/10

fit=as_cmdstan_fit("./tst_folder_sig/SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v10-202511131059-1-3f135e.csv")

#fit2=as_cmdstan_fit(paste0("./tst_folder_sig/SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v6-202508251617-",chn,"-b39fc5.csv"))

#fit$output_files()=paste0(getwd(),"/",new_csv)

#fit <- read_cmdstan_csv(new_csv)


N_C <- 11
N_C_use=11
TT <- 30


rho_se = 0
rho_ei = fit$draws("rho_ei", format = "draws_matrix")[fiter]
rho_ir = fit$draws("rho_ir", format = "draws_matrix")[fiter]
phi = fit$draws("phi", format = "draws_matrix")[fiter]
sig_beta = fit$draws("sig_beta", format = "draws_matrix")[fiter]
mu_log_beta = fit$draws("mu_log_beta", format = "draws_matrix")[fiter]
#first=sample(3:8,N_C,replace=T);first

#pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
pop=read_csv("./utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
pop_size = pop$Population_2020[counties][1:N_C_use]
#pop=1e4*sample(1:N_C)#pop[counties,]

#counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
#pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
#pop=pop$Population_2020[counties]#1e4*sample(5:N_C)#pop[counties,]
#if(iter %in% 1:100) {p=.15} else if (iter %in% 101:200) {p=.5} else {p=.85}
if(iter %in% 1:200) {p=.15} else if (iter %in% 201:400) {p=.5} else {p=.25}


#fit$draws("p", format = "draws_matrix")[fiter]
sigma = fit$draws("sigma", format = "draws_matrix")[fiter]

#log_beta <- numeric(TT-1)
#log_beta[1] <- rnorm(1,0,sigma/sqrt(1-phi^2));log_beta[1]

#for (t in 2:(TT-1)) {
#  log_beta[t] <- phi * log_beta[t - 1] + rnorm(1,0,sigma)
#}

#log_beta_mat=matrix(0,TT-1,N_C)
#for (t in 1:(TT-1)) {
#  log_beta_mat[t,] <- log_beta[t] + sig_beta * rnorm(N_C,0,1)
#}

beta=matrix(fit$draws("beta_mat", format = "draws_matrix")[fiter,], nrow = TT, ncol = N_C)[,1:N_C_use]

#exp(log_beta) %>% plot
#beta=exp(log_beta_mat)

#saveRDS(fit$summary("beta_mat")$mean %>% matrix(TT,N_C),file="beta.rds")
#beta=readRDS('beta.rds')

# create an empty list
gamma <- as.numeric(fit$draws("gamma", format = "draws_matrix")[fiter,])[1:N_C_use]#readRDS('gamma.rds')#runif(N_C, min = 0.59, max = 0.87)

eta <- as.numeric(fit$draws("eta", format = "draws_matrix")[fiter,])[1:N_C_use]#readRDS('eta.rds')#runif(N_C, min = 0.29, max = 0.88)

#E_0 <- round(.009*pop_size)#ample(10:20,N_C,TRUE)#round(c(11902,6600,1120,737))#
i0= as.numeric(fit$draws("i0", format = "draws_matrix")[fiter,])[1:N_C_use]#readRDS('i0.rds')#runif(N_C,.0003,.01)
I_0 <- round(i0*pop_size)#

first=readRDS("first.rds")[1:N_C_use]#sample(3:8,N_C,replace=T);first
#phi_p = as.numeric(fit$draws("phi_p", format = "draws_matrix")[fiter])#readRDS('phi_p.rds')#runif(N_C, min = 0.05, max = 0.5)

S_0 <- pop_size - I_0
R <- S <- I <- E <- matrix(0,TT, N_C_use)
ii <- SE <- IR <- EI <- matrix(0,TT-1, N_C_use)
p_detect <- p
for (i in 1:N_C_use){
S[1:(first[i]-1),i]= 0     
S[first[i],i] <- S_0[i]

I[1:(first[i]-1),i] = 0
I[first[i],i] <- I_0[i]

EI[first[i],i]=I_0[i]
}


#imp_rate=50


for (ct in  1:N_C_use) {
    for (t in 2:(TT)){
    if (t > first[ct]) {
     SE[t-1,ct]=rbetabinom(1,S[t-1,ct],prob=1-exp(- sum(beta[t-1,ct] * I[t-1,ct]/ pop_size[ct])),rho=rho_se) 
     if (EI[t-1,ct] == 0) {
          EI[t-1,ct]=rbetabinom(1,E[t-1,ct],eta[ct],rho=rho_ei) 
     } 
     IR[t-1,ct]=rbetabinom(1,I[t-1,ct],gamma[ct],rho=rho_ir)
     E[t,ct]=E[t-1,ct]+SE[t-1,ct]-(EI[t-1,ct]*as.numeric(t> (first[ct]+1))) #+ rpois(1,pop_size[ct]/1000)
     I[t,ct]=I[t-1,ct]+EI[t-1,ct]-IR[t-1,ct]
     S[t,ct]=S[t-1,ct]-SE[t-1,ct]
     R[t,ct]=R[t-1,ct]+IR[t-1,ct]
  }
  }
}



ii=matrix(rbinom(N_C_use*(TT-1),EI,p_detect),nrow=TT-1);ii

#ii=ii[,-11]
ii[is.na(ii)] <- 0
N_C=ncol(ii)
#ii_save
#ii=ii[,-c(8,11)] # remove counties with no data
last<- unlist(1:N_C %>% purrr::map(function(col) {
  w <- which(ii[,col] == 0)
  w=w[which(w>first[col])][1]-1
  w=ifelse(is.na(w),TT-1,w)
}))

# run_len <- 20      # <-- choose how many consecutive time points define the tail
# case_cut <- 75    # <-- "less than 100 cases"

# last <- purrr::map_int(1:N_C_use, function(col) {

#   x <- ii[, col]

#   # ---- existing rule: stop at first 0 after first[col] ----
#   w0 <- which(x == 0)
#   w0 <- w0[w0 > first[col]][1] - 1
#   w0 <- ifelse(is.na(w0), TT - 1, w0)

#   # ---- NEW: only allow "low-run truncation" after the peak ----
#   peak_idx <- which.max(x)   # first occurrence of the max

#   # start searching AFTER the peak (and after first[col], if that matters)
#   start_idx <- max(first[col], peak_idx) + 1
#   end_idx   <- TT - 1

#   if (start_idx > end_idx) return(w0)

#   low <- x[start_idx:end_idx] < case_cut

#   r <- rle(low)
#   ends   <- cumsum(r$lengths)
#   starts <- ends - r$lengths + 1

#   hit <- which(r$values & r$lengths >= run_len)[1]

#   wlow <- if (is.na(hit)) {
#     TT - 1
#   } else {
#     # truncate right before the low-run begins
#     (start_idx + starts[hit] - 1) - 1
#   }

#   min(w0, wlow)
# })

#ii  %>% as_tibble() %>% mutate(x=1:(TT-1)) %>% gather(Key, Value,-x) %>%
#ggplot(aes(x=x, y=Value)) + geom_line() + facet_wrap(~Key,scales="free_y")


#apply(ii, 2, function(col) {
# w <- which(col != 0)
#  if (length(w)) tail(w, 1) else NA_integer_
#})

#ii_0=rbinom(N_C,I_0,p_detect)

#for (i in 1:N_C){
#ii[first[i],i]=ii_0[i]
#}
#saveRDS(ii,file="ii.rds")
#ii=readRDS("ii.rds")
tt <- cmdstan_model("SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v10.stan")
dat <- 
  list(
    ii = ii,#as.matrix(d1)[,counties],
    TT = TT,
    N_C = dim(ii)[2],
    pop_size = pop_size,
    #D=dist[counties,counties]/10,
    first=first,
    last=last,
    min_first=min(first)
  )

init = \() {list(u_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  v_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  w_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  #raw_log_beta_mat = matrix(runif(TT*N_C, -1, 1), TT, N_C),#matrix(rnorm(TT*N_C, 0, 1), TT, N_C),
                                  raw_log_beta = rnorm(sum(dat$last - dat$first + 1)),
                                  #log_phi_p = rnorm(1,4,1),
                                  p_raw = runif(1, -.25,.25),
                                  #kappa = runif(1, 0,1),
                                  v_raw=runif(1,1.5,2.5),
                                  z=rnorm(TT-min(first) + 1,0,.25),
                                  mu_log_beta = rnorm(1, 0, .05),
                                  sigma = runif(1, .1, .25),
                                  sig_beta = runif(1, .1,.25),
                                  i0_raw = runif(N_C,-9,-7),#rbeta(N_C, 0.01*50, 0.99*50),
                                  #rho_si = runif(1, 0.0001, 0.005),
                                  rho_ei_raw = runif(1, -0.1, 0.1),
                                  rho_ir_raw = runif(1, -0.1, 0.1),
                                  gamma_raw = runif(N_C, 0,.5),
                                  eta_raw = runif(N_C,.25,.75))}
#saveRDS(dat,file="dat2.rds")
#dat=readRDS("dat2.rds")
#seed(123)
init=list(init(),init(),init(),init())

fit = tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.99,
                 max_treedepth = 18,
                 init = init,
                 iter_warmup = 500,
                 iter_sampling = 500, parallel_chains = 4,
                 step_size=.0009)#,


diagn=fit$diagnostic_summary()





# sv=fit$draws("log_beta") %>% posterior::as_draws_rvars() 
# #sv=fit$draws("log_beta") %>% posterior::as_draws_rvars() 

# sv$beta_mat[,c(1,2,3,4)]

# np_fit <- nuts_params(fit)
# png("tst_np.png")
# mcmc_pairs(fit$draws(c("sigma","sig_beta","p","mu_log_beta","i0","rho_ei_raw","rho_ir_raw","phi")), np = np_fit, 
#            pars = c("sigma","sig_beta","p","mu_log_beta","i0[1]","phi"),
#             off_diag_args = list(size = 0.75))
# dev.off()


#  z_t_d <- fit$draws("ei_t", format = "draws_array") |> posterior::as_draws_rvars()
#  z_t_d <- z_t_d$ei_t
#  qpt025 <- quantile(z_t_d,0.025)
#  qpt975 <- quantile(z_t_d,0.975)
#  pop4=pop_size
#  obs=rbind(data.frame(value="mean",sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
#  data.frame(value="lwr",sweep(qpt025[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
#  data.frame(value="upr",sweep(qpt975[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")))

#  obs=obs %>% gather(County,Cases,-value) %>% pivot_wider(names_from=value,values_from=Cases) %>% unnest() %>% mutate(date=rep(1:(TT-1),N_C))







#qplot(exp(rlogis(1e6,0,1)))

# Performance for simulation
# Average absolute error p, coverage, CI width  
# Average absolute error of betas by county, coverage, CI width
# Average absolute error EI by county, coverage, CI width
# Average absolute error of I by county, coverage, CI width

# Change p for each simulation using the same chn iter combo

#p

p1=abs(p-fit$summary("p")$mean)
p2=between(p, fit$summary("p")$q5, fit$summary("p")$q95)
p3=fit$summary("p")$q95 - fit$summary("p")$q5

#Beta Performance
beta1=abs(matrix(fit$summary("beta_mat")$mean,nrow=TT) - beta)
beta2=matrix(between(as.vector(beta),fit$summary("beta_mat")$q5, fit$summary("beta_mat")$q95),nrow=TT)
beta3=matrix(fit$summary("beta_mat")$q95 - fit$summary("beta_mat")$q5,nrow=TT)

# EI Performance
EI1=unlist(1:N_C %>% purrr::map(function(x) mean(abs(matrix(fit$summary("ei_t")$mean,nrow=TT)[(first[x]):(last[x]),x] - (EI/pop_size[x])[(first[x]):(last[x]),x]))))
EI2=unlist(1:N_C %>% purrr::map(function(x) mean(between((EI/pop_size[x])[(first[x]):(last[x]),x],
    as.vector(matrix(fit$summary("ei_t")$q5,nrow=TT)[(first[x]):(last[x]),x]),
    as.vector(matrix(fit$summary("ei_t")$q95,nrow=TT)[(first[x]):(last[x]),x])))))
EI3=unlist(1:N_C %>% purrr::map(function(x) {mean(
    matrix(fit$summary("ei_t")$q95,nrow=TT)[(first[x]):(last[x]),x] -
    matrix(fit$summary("ei_t")$q5,nrow=TT)[(first[x]):(last[x]),x])}))


# I Performance
I1=unlist(1:N_C %>% purrr::map(function(x) mean(abs(matrix(fit$summary("i_t")$mean,nrow=TT)[(first[x]):(last[x]+1),x] - (I/pop_size[x])[(first[x]):(last[x]+1),x]))))
I2=unlist(1:N_C %>% purrr::map(function(x) mean(between((I/pop_size[x])[(first[x]):(last[x]+1),x],
    as.vector(matrix(fit$summary("i_t")$q5,nrow=TT)[(first[x]):(last[x]+1),x]),
    as.vector(matrix(fit$summary("i_t")$q95,nrow=TT)[(first[x]):(last[x]+1),x])))))
I3=unlist(1:N_C %>% purrr::map(function(x) {mean(
    matrix(fit$summary("i_t")$q95,nrow=TT)[(first[x]):(last[x]+1),x] -
    matrix(fit$summary("i_t")$q5,nrow=TT)[(first[x]):(last[x]+1),x])}))

saveRDS(list(dat,diagn,init,p1,p2,p3,beta1,beta2,beta3,EI1,EI2,EI3,I1,I2,I3,fit$summary()),file=paste0("./chn3/SEIR_sim",chn,"_",iter,".rds"))
# Save R_hat

#sum_fit <- fit$summary()


###
# draws_df$divergent__

# draws_array <- fit$draws(inc_warmup = FALSE, format = "draws_array")

# # mcmc_pairs will colour points with divergent__ == 1 in red
# mcmc_pairs(
#   draws_array,
#   pars = c("mu_log_beta", "phi", "sigma")
# )

# np_fit <- nuts_params(fit)
# png("tst23r.png")
# mcmc_pairs(fit$draws(c("mu_log_beta", "phi", "sigma","sig_beta","p","beta_mat[1,1]")), np = np_fit, pars = c("mu_log_beta", "phi", "sigma","sig_beta","p","beta_mat[1,1]"),
#             off_diag_args = list(size = 0.75))
# dev.off()





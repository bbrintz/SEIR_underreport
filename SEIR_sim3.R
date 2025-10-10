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
fiter=500#500#sample(1:1500,1)

#set.seed(123)
#beta <- 1.05
#TT <- 100

#counties=c(1,2,3,8)#c(1,2,3,4,5,6,7,8,10,11,12) # dist/10

fit=as_cmdstan_fit(paste0("./tst_folder_sig/SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v6-202509052035-",chn,"-1d5eb8.csv"))


#fit2=as_cmdstan_fit(paste0("./tst_folder_sig/SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v6-202508251617-",chn,"-b39fc5.csv"))


fiter=fiter
#fit$output_files()=paste0(getwd(),"/",new_csv)

#fit <- read_cmdstan_csv(new_csv)


N_C <- 11
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
pop_size = pop$Population_2020[counties]
#pop=1e4*sample(1:N_C)#pop[counties,]

#counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
#pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
#pop=pop$Population_2020[counties]#1e4*sample(5:N_C)#pop[counties,]
if(iter %in% 1:100) {p=.15} else if (iter %in% 101:200) {p=.5} else {p=.85}


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

beta=matrix(fit$draws("beta_mat", format = "draws_matrix")[fiter,], nrow = TT, ncol = N_C)

#exp(log_beta) %>% plot
#beta=exp(log_beta_mat)

#saveRDS(fit$summary("beta_mat")$mean %>% matrix(TT,N_C),file="beta.rds")
#beta=readRDS('beta.rds')

# create an empty list
gamma <- as.numeric(fit$draws("gamma", format = "draws_matrix")[fiter,])#readRDS('gamma.rds')#runif(N_C, min = 0.59, max = 0.87)

eta <- as.numeric(fit$draws("eta", format = "draws_matrix")[fiter,])#readRDS('eta.rds')#runif(N_C, min = 0.29, max = 0.88)

#E_0 <- round(.009*pop_size)#ample(10:20,N_C,TRUE)#round(c(11902,6600,1120,737))#
i0= as.numeric(fit$draws("i0", format = "draws_matrix")[fiter,])#readRDS('i0.rds')#runif(N_C,.0003,.01)
I_0 <- round(i0*pop_size)#

first=readRDS("first.rds")

phi_p = as.numeric(fit$draws("phi_p", format = "draws_matrix")[fiter])#readRDS('phi_p.rds')#runif(N_C, min = 0.05, max = 0.5)

S_0 <- pop_size - I_0
R <- S <- I <- E <- matrix(0,TT, N_C)
ii <- SE <- IR <- EI <- matrix(0,TT-1, N_C)
p_detect <- p
for (i in 1:N_C){
S[1:(first[i]-1),i]= 0     
S[first[i],i] <- S_0[i]

I[1:(first[i]-1),i] = 0
I[first[i],i] <- I_0[i]

EI[first[i],i]=I_0[i]
}


#imp_rate=50


for (ct in  1:N_C) {
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



ii=matrix(rbetabinom(N_C*(TT-1),EI,p_detect,rho=1/(1+phi_p)),nrow=TT-1);ii


ii[is.na(ii)] <- 0

last<- unlist(1:N_C %>% purrr::map(function(col) {
  w <- which(ii[,col] == 0)
  w=w[which(w>first[col])][1]-1
  w=ifelse(is.na(w),TT-1,w)
}))

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
tt <- cmdstan_model("SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v6.stan")
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
                                  
                                  p_raw = runif(1, -.25,.25),
                                  #kappa = runif(1, 0,1),
                                  phi_p = runif(1, 10, 1000),
                                  v_raw=runif(1,.1,.25),
                                  z=rnorm(TT,0,.25),
                                  mu_log_beta = rnorm(1, 0, .25),
                                  sigma = runif(1, .05, .5),
                                  sig_beta = runif(1, .25, .75),
                                  i0_raw = runif(N_C,-3,-1),#rbeta(N_C, 0.01*50, 0.99*50),
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
                 max_treedepth = 16,
                 init = init,
                 iter_warmup = 1500,
                 iter_sampling = 1000, parallel_chains = 4)#,


diagn=fit$diagnostic_summary()


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
beta1=apply(abs(matrix(fit$summary("beta_mat")$mean,nrow=TT) - beta),2,mean)
beta2=apply(matrix(between(as.vector(beta),fit$summary("beta_mat")$q5, fit$summary("beta_mat")$q95),nrow=TT),2,mean)
beta3=apply(matrix(fit$summary("beta_mat")$q95 - fit$summary("beta_mat")$q5,nrow=TT),2,mean)

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

saveRDS(list(dat,diagn,init,p1,p2,p3,beta1,beta2,beta3,EI1,EI2,EI3,I1,I2,I3,fit$summary()),file=paste0("./chn2fiter500_v2/SEIR_sim",chn,"_",iter,".rds"))
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
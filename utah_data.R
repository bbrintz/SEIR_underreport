library(tidyverse) 
library(cmdstanr)
library(bayesplot)
library(posterior)
library(here)



#set working directory to current file path
case = read_csv("./data/time_series_covid19_confirmed_US.csv") %>% filter(`Province_State`=="Utah")
test = read_csv("./data/time_series_covid19_US_testing_by_state.csv") %>% filter(state=="UT")
vax = read_csv("./data/COVID-19_Vaccinations_in_the_United_States_County.csv") %>% filter(Recip_State=="UT")
test=test %>% mutate(date=mdy(date)) %>% arrange(date) %>% dplyr::select(date,tests_combined_total)
vax=vax %>% mutate(date=mdy(Date)) %>% arrange(date) %>% dplyr::select(date,vax=Series_Complete_Yes) %>% 
  group_by(date) %>% summarize(vax=sum(vax))

case=case %>% dplyr::select(-c(1:5,7:11)) %>% gather("date","cases",-Admin2) %>% mutate(date=mdy(date)) %>% group_by(Admin2,date) %>% summarize(cases=sum(as.numeric(cases)))
case=case %>% ungroup() %>% mutate(date=(date = 7 * (as.numeric(date - min(date)) %/% 7) + min(date))) %>% 
  group_by(date,Admin2) %>% summarize(cases=sum(cases)) %>%#,tests=sum(tests_combined_total-lag(tests_combined_total),na.rm=T))#,
                               #vax=sum(vax-lag(vax),na.rm=T))
                               #filter(date>=ymd("2020-04-15")) %>% 
                               filter(date<ymd("2021-06-02"))
                               #filter(date>=ymd("2020-10-06")) %>% filter(date<ymd("2021-02-02"))


case=case %>% group_by(Admin2) %>% mutate(new_cases=cases-lag(cases)) %>% arrange(Admin2,date) %>% filter(!is.na(new_cases))# %>% filter(Admin2 %in% readRDS("final_counties.rds")) %>%# filter(any(cases>1000)) %>%
#saveRDS(unique(dat_final$Admin2),"final_counties.rds")


dat_final=case %>% filter(Admin2 %in% readRDS("final_counties.rds"))#group_by(Admin2) %>% filter(any(cases>100))
pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))

# get the distance between each county (Admin2) in dat_final
dist=matrix(0,nrow=length(unique(dat_final$Admin2)),ncol=length(unique(dat_final$Admin2)))
for(i in 1:12){
  for(j in 1:12){
    dist[i,j]=geosphere::distHaversine(as.matrix(pop %>% dplyr::select(lon=Longitude,lat=Latitude))[c(i,j),])/1000
  }
} 

d1=dat_final %>% filter(date<ymd("2020/8/19")) %>% 
ungroup() %>% rename(County="Admin2") %>% left_join(pop,by="County") %>% arrange(desc(Population_2020),date) %>%
dplyr::select(-Population_2020,-Latitude,-Longitude,-cases) %>% pivot_wider(names_from=County,values_from=new_cases) %>%
dplyr::select(-date)
d1[6,9]=1

tt <- cmdstan_model("SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v10.stan")#cmdstan_model("SEIR_betabin_vary_beta_nospat.stan")#"stoch_beta_spatial_SI_utah_betabin.stan")


#d1=d1 %>% dplyr::select(`Salt Lake`,Utah,Davis,`Weber-Morgan`)

    TT = nrow(d1)+1

#1,2,3,8 for 4 counties # dist/10
counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
#counties=c(1,2,3,8)
N_C = length(counties)#ncol(d1)

ii=as.matrix(d1)[,counties]

first=apply(matrix(as.numeric(as.matrix(d1)[,counties]!=0),ncol=N_C),2,function(x) which(x==1)[1])
last<- unlist(1:N_C %>% purrr::map(function(col) {
  w <- which(ii[,col] == 0)
  w=w[which(w>first[col])][1]-1
  w=ifelse(is.na(w),TT-1,w)
}))
dat <- 
  list(
    t0 = min(first), 
    ii = as.matrix(d1)[,counties],
    TT = TT,
    N_C = N_C,
    pop_size = pop$Population_2020[counties],
    D=dist[counties,counties]/10,
    first=first,
    last=last,
    min_first=min(first)
  )

dat$ii %>% as_tibble %>% gather(County,Cases) %>% mutate(date=rep(1:(TT-1),N_C),first=rep(first,each=TT-1)) %>% 
  ggplot(aes(x=date,y=Cases)) + geom_line() + geom_vline(aes(xintercept=first),color="red") + facet_wrap(~County,scales="free_y") + theme_bw() + ylab("New Cases") + xlab("Date")

? cmdstanr::cmdstan_model

set.seed(123)

# No e0 
# Calculate I0, but get rid of V_t at first time point 
# First time point, I0 is EI so use I0 for detection thing
fit = tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.90,
                 max_treedepth = 14,
                 init = \() {list(u_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  v_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  w_t_logit_eta = matrix(rnorm(TT*N_C, 0,1), TT, N_C),
                                  #raw_log_beta_mat = matrix(runif(TT*N_C, -1, 1), TT, N_C),#matrix(rnorm(TT*N_C, 0, 1), TT, N_C),
                                  raw_log_beta = rnorm(sum(dat$last - dat$first + 1)),
                                  #log_phi_p = rnorm(1,4,1),
                                  p_raw = runif(1, -.25,.25),
                                  #kappa = runif(1, 0,1),
                                  phi_p = runif(1, 100, 200),
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
                                  eta_raw = runif(N_C,.25,.75))},#rbeta(N_C, 0.5 * 4, 0.5 * 4))},
                 iter_warmup = 500,#1500,
                 iter_sampling = 500, parallel_chains = 4,
                 output_dir = paste0(getwd(),"/tst_folder_sig/"),
                 step_size=.0009
                 )



#fit=readRDS("fit_utah_10.rds")

#new_csv=unlist(1:10 %>% purrr::map(~paste0("SEIR_betabin_hier_ar1_beta_zeros_v2-202502041329-",.,"-5031f7.csv")))
#new_csv=unlist(1:10 %>% purrr::map(~paste0("SEIR_betabin_hier_ar1_beta_zeros_v2-202502121643-",.,"-2f66c5.csv")))
#new_csv1=unlist(1:10 %>% purrr::map(~paste0("./ignore/SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v6-202508250811-",.,"-a5dbaf.csv")))
#new_csv2=unlist(1:10 %>% purrr::map(~paste0("./beta_bin_results/SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v2-202504300918-",.,"-56a7f2.csv")))

new_csv1=unlist(1:10 %>% purrr::map(~paste0("./tst_folder_sig/SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v8-202510131149-",.,"-6a1619.csv")))




fit$diagnostic_summary()
fit$summary("p")

np_fit <- nuts_params(fit)

np_fit %>% 
  filter(Parameter == "stepsize__") %>%
  mutate(stage = case_when(Iteration <= iter_warmup ~ "warmup", TRUE ~ "sampling")) %>%
  group_by(Chain, stage) %>%
  summarize(median_stepsize = median(Value), .groups="drop")

mcmc_pairs(fit$draws(c("p","phi","sigma","sig_beta","rho_ei","i0","mu_log_beta")), np = np_fit, pars = c("p","phi","sigma","sig_beta","mu_log_beta","i0[1]"),
                  off_diag_args = list(size = 0.75))

#fit$output_files()=paste0(getwd(),"/",new_csv)

#fit <- read_cmdstan_csv(new_csv)
fit1=as_cmdstan_fit(new_csv1)
fit=fit1
sv=fit$draws("beta_mat") %>% posterior::as_draws_rvars() 
#sv=fit$draws("log_beta") %>% posterior::as_draws_rvars() 

sv$beta_mat[,c(1,2,3,4)]


fit2=as_cmdstan_fit(new_csv2)



ppost=rbeta(15000,fit1$draws("p",format="draws_matrix")[,1]*fit1$draws("phi_p",format="draws_matrix")[,1],
(1-fit1$draws("p",format="draws_matrix")[,1])*fit1$draws("phi_p",format="draws_matrix")) 
qplot(ppost)
qplot(fit1$draws("p",format="draws_matrix")[,1])

fit2$summary(c("lp__"))

log_lik1 <- fit1$draws("lp__", format = "matrix")
loo1 <- loo(fit1)
log_lik2 <- fit2$draws("lp__", format = "matrix")
loo2 <- loo(log_lik2)

waic1 <- waic(log_lik1)

loo_compare(loo1, loo2)

fit_sum=fit$summary() %>% as_tibble 
fit_sum %>% filter(str_detect(variable,"^i_t")) %>% mutate(var=substr(variable,5,10 )) %>% 
mutate(var=gsub("]","",var)) %>% 
separate(var,into=c("Time","County")) %>% 
mutate(Time=as.numeric(Time),County=as.numeric(County)) %>% 
ggplot(aes(x=Time,y=mean)) + geom_point() + facet_wrap(~County)

draws=fit$draws()
df_draws=as_draws_df(draws)
write.csv(df_draws %>% as_tibble,file="utah_draws_020325.csv")
saveRDS(df_draws %>% as_tibble,file="utah_draws_020325.rds")


fit$diagnostic_summary()
fit$sampler_diagnostics()
fit$time()


fit$init() %>% purrr::map(~.$eta[1])

write.csv(fit$summary(c("p","phi","sigma","sig_beta","rho_ei","rho_ir")) %>% as_tibble %>%
select(-mad,-ess_bulk) %>% mutate_if(is.numeric,~round(.,3)),file="betabin_ests.csv")

fit$summary(c("p","phi","sigma","sig_beta","phi_p","rho_ei","rho_ir"))
fit$summary(c("phi_p"))



#betas=fit$draws("beta_mat", format = "draws_array") |> posterior::as_draws_rvars()

fit$draws("log_beta") %>% as_tibble %>% gather(Key,Value) %>%
  extract(Key, into = c("variable", "index"), regex = "(.*)\\[(\\d+)\\]", convert = TRUE) %>%
  group_by(index)%>% summarize(Val=mean(Value),Low=quantile(Value,0.025),High=quantile(Value,0.975)) %>% ungroup() %>%
  mutate_all(~as.numeric(.)) %>%
  ggplot(aes(x=index,y=Val,ymin=Low,ymax=High)) + geom_line() + geom_ribbon(alpha=.25)

betas=fit$draws("beta_mat") %>% as_tibble %>% gather(Key,Value) %>%
  extract(Key, into = c("n1", "n2"), regex = ".*\\[(\\d+),(\\d+)\\].*") %>%
  group_by(n1,n2) %>% summarize(Val=mean(Value),Low=quantile(Value,0.025),High=quantile(Value,0.975)) %>% ungroup() %>%
  mutate_all(~as.numeric(.))
###
log_betas=fit$draws("log_beta") %>% as_tibble %>% gather(Key,Value) %>%
  extract(Key, into = c("variable", "index"), regex = "(.*)\\[(\\d+)\\]", convert = TRUE) %>%
  group_by(index)%>% summarize(Val=mean(Value),Low=quantile(Value,0.025),High=quantile(Value,0.975)) %>%
  filter(index>9,index<29)

betas=fit$draws("beta_mat") %>% as_tibble %>% gather(Key,Value) %>%
  extract(Key, into = c("n1", "n2"), regex = ".*\\[(\\d+),(\\d+)\\].*") %>%
  group_by(n1,n2) %>% summarize(Val=mean(Value),Low=quantile(Value,0.025),High=quantile(Value,0.975)) %>% ungroup() %>%
  mutate_all(~as.numeric(.)) %>%  filter(n1>9,n1<29)

log_betas %>% ggplot(aes(x=index,y=exp(Val),ymin=exp(Low),ymax=exp(High))) + geom_line() + geom_ribbon(alpha=.25) +
  #geom_hline(yintercept=0,linetype="dashed") + geom_vline(aes(xintercept=29)) +
  geom_point(data=betas,aes(x=n1,y=Val)) +
  ylab("log beta") + xlab("County") + theme_bw() + ylim(0,5)
  
log_betas_time=rbind(log_betas %>% mutate(time=1),log_betas %>% mutate(time=2))

betas_time=rbind(betas %>% rename(index=n1) %>% mutate(Val=exp(log_betas$Val[match(index,log_betas$index)])) %>% mutate(time=1),
betas %>% rename(index=n1) %>% mutate(time=2))

saveRDS(list(log_betas_time,betas_time),file="betas_time.rds")


p <- ggplot() +
  geom_line(data = log_betas_time, aes(x = index, y = exp(Val)), size = 1) +
  geom_ribbon(data = log_betas, aes(x = index, ymin = exp(Low), ymax = exp(High)),
              alpha = 0.25, fill = "blue") +
  # Use the same transformation for y if desired (here I use exp(Val) so they align with the line)
  geom_point(data = betas_time, aes(x = index, y = Val), color = "red", size = 2) +
  labs(x = "County", y = "transformed log beta") +
  theme_bw() +
  ylim(0, 5) +
  transition_time(time)

###

fit$draws("log_beta") %>% as_tibble %>% gather(Key,Value) %>%
  extract(Key, into = c("variable", "index"), regex = "(.*)\\[(\\d+)\\]", convert = TRUE) %>%
  group_by(index)%>% summarize(Val=mean(Value)) %>% pull(Val) %>% saveRDS("log_beta_nobetabinom.rds")


quartz()
betas %>% filter(n1 > (first[n2]+1),n1<30) %>% ggplot(aes(x=n1,y=Val,ymin=Low,ymax=High))  + geom_line() + geom_ribbon(alpha=.25) + facet_wrap(~n2,scales="free_y")
ggsave("utah_beta_020325.png")
fit$draws("p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1500,10),group=key,color=key)) + geom_line()

fit$draws("phi_p") %>% as_tibble %>% gather() %>% ggplot(aes(x=value,group=key,fill=key)) + geom_density(alpha=.25) + theme_bw() +
scale_fill_viridis_d()

library(ggplot2)
library(gganimate)
library(dplyr)
library(tidyr)
library(viridis)

plt <- rbind(
  data.frame(variable = "Prior", p = plogis(rnorm(1e6, 0, 1)), Time = 1),
  fit$draws("p") %>% 
    as_tibble() %>% 
    pivot_longer(everything(), names_to = "variable", values_to = "p") %>% 
    mutate(Time = 2)
) %>%
  ggplot(aes(x = p, group = variable, fill = variable)) +
  geom_density(alpha = 0.25) +
  # Optionally add the theoretical density overlay here if desired
  # stat_function(
  #   fun = function(p) dnorm(qlogis(p))/(p*(1 - p)),
  #   color = "#2a1414", size = 1, xlim = c(0, 1)
  # ) +
  xlab("Clinical Detection Rate") + 
  ylab("Density") +
  ggtitle("Posterior Distribution of Clinical Detection Rate") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  transition_states(Time, transition_length = 2, state_length = 1) +
  ease_aes('linear') +
  view_follow(fixed_x = FALSE, fixed_y = FALSE)    

# Render the animation, setting parameters like number of frames and fps
anim <- animate(plt, nframes = 100, fps = 10, renderer = gifski_renderer())
# Save the animation to a file
anim_save("clinical_detection_rate.gif", animation = anim)

fit$draws("p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,10),group=key,color=key)) + geom_line()


ggsave("utah_p_dens_020325.png")
fit$draws("sigma") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1500,10),group=key,color=key)) + geom_line()
ggsave("utah_sigma_020325.png")
fit$draws("sig_beta") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1500,10),group=key,color=key)) + geom_line()
ggsave("utah_sig_beta_020325.png")
fit$draws("phi") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1500,10),group=key,color=key)) + geom_line()
ggsave("utah_phi_020325.png")

fit$draws("phi_p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,10),group=key,color=key)) + geom_line()


fit$draws("eta") %>% as_tibble %>% gather() %>% extract(key, into = c("n1"), regex = ".*\\[(\\d+)].*") %>% 
group_by(n1) %>% summarize(value=mean(value)) %>% ungroup() %>% mutate(n1=as.numeric(n1)) %>%
ggplot(aes(y=value,x=n1)) + geom_point()
ggsave("utah_eta_020325.png")

fit$draws("gamma") %>% as_tibble %>% gather() %>% extract(key, into = c("n1"), regex = ".*\\[(\\d+)].*") %>% 
group_by(n1) %>% summarize(value=mean(value)) %>% ungroup() %>% mutate(n1=as.numeric(n1)) %>%
ggplot(aes(y=value,x=n1)) + geom_point()
ggsave(utah_gamma_020325.png)

z_t_d <- fit1$draws("ei_t", format = "draws_array") |> posterior::as_draws_rvars()
z_t_d <- z_t_d$ei_t
qpt025 <- quantile(z_t_d,0.025)
qpt975 <- quantile(z_t_d,0.975)

pop4=dat$pop_size


obs=rbind(data.frame(value="mean",sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")),
data.frame(value="lwr",sweep(qpt025[1,-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")),
data.frame(value="upr",sweep(qpt975[1,-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")))

obs=obs %>% gather(County,Cases,-value) %>% pivot_wider(names_from=value,values_from=Cases) %>% unnest() %>% mutate(date=rep(1:(TT-1),N_C)) #%>%
#filter(County %in% unique(truth$County)[1:4])#



#truth=SI %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
#mutate(County= gsub("V", "X", County)) #%>% filter(County %in% unique(truth$County)[1:4])

fit$summary("p")

truth=dat$ii %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) 
truth$County


unique(truth$County)
quartz()

dts=d1 %>% select(Date=date) %>% mutate(date=1:29)
obs %>% mutate(County=rep(unique(truth$County),each=29),date=date+1) %>% filter(date<30) %>%
left_join(truth,by=c("County","date")) %>% 
ggplot(aes(x=date,y=mean,ymin=lwr,ymax=upr)) + geom_point(color="blue") + 
geom_ribbon(alpha=.75,fill="lightblue") +
facet_wrap(~County,scales="free_y") + 
geom_point(data=dat$ii %>% as_tibble %>% gather(County,Cases) %>% mutate(date=rep(1:(TT-1),N_C),lwr=1,upr=1),
aes(x=date,y=Cases),color="black") + geom_line() + theme_bw() + ylab("New Cases") + xlab("Date") + 
scale_x_continuous(breaks=seq(1, 29, by=1),labels=format(dts$Date,"%m/%d/%y")) + 
theme(axis.text.x = element_text(angle = 90, hjust = 0,size=8)) + 
theme(legend.position="none")

ggsave("utah_102125_bb.png",width=12,height=8)

write.csv(data.frame(mu, k, a, b), file="rbeta_params.csv")



mu=.6
k=5
a=mu*k
b=k-a

rbeta(1e6,a,b) %>% hist(breaks=100,freq=FALSE)


plogis(rlogis(1e6,0,1)) %>% hist(breaks=100,freq=FALSE)
exp(rnorm(1e6,0,1)) %>% hist(breaks=100,freq=FALSE)



qplot(rgamma(1e6,shape=1,rate=.01))


# assume you’ve already done:
# tt <- cmdstan_model("SEIR_betabin_on_hier_ar1_beta_pbeta_zeros_v4.stan")
# dat  <- (your data list)

# 1) do a prior‐only run
prior_fit <- tt$sample(
  data          = dat,
  chains        = 500,
  iter_warmup   = 0,
  iter_sampling = 1,
  refresh       = 0,
  fixed_param   = TRUE
)



# 2) extract the draws for beta_mat[1,1] and beta_mat[2,2]
draws_df <- as_draws_df(prior_fit$draws())

hist(draws_df$`beta_mat[1,2]`)

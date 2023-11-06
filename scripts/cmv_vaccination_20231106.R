##########################################################################################
# R script for inference of CMV transmission model parameters using serological data     #
# (PIENTER2) and  data from a prospective cohort study (CROCUS; PhD thesis M Korndewal). #
# Code is based on earlier work with Chris van Dorp and Sophia de Jong.                  #
# Code is copyrighted by Michiel van Boven and licensed under BSD (3-clause).            #      
##########################################################################################

library(rstan)
library(loo)
library(readr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# import CMV data
cmvdata <- read_csv("cmvdata.csv")

# filter non-western subjects
cmvdata = cmvdata[cmvdata$nl == '1',]

# filter ages <6 months
cmvdata = cmvdata[cmvdata$lftinmnd2 > 6,]

# select relevant columns
Ages = cmvdata$LFTB
Titers = cmvdata$boxcoxvalues
N = length(Ages)

# censoring and spike
RightCensor = 3.41
Spike = -2.91
Censor = rep(0,N)
for (i in 1:N){
  if (Titers[i] > RightCensor){
    Censor[i] = 1
  } 
  else if (Titers[i] <= Spike){
    Censor[i] = 2
  }
} 

# make gender numeric
Gender= rep(0,N)
Gender[cmvdata$gesl2 == 'male'] = 1 # female = 0

# import contact matrix
ContactData = read_csv("contact_intensities_aggregated.csv")
MM = ContactData$mMM # length 256 = 16*16
FM = ContactData$mFM
MF = ContactData$mMF
FF = ContactData$mFF
Contact_MM = matrix(MM,nrow=16,ncol=16)
Contact_FM = matrix(FM,nrow=16,ncol=16)
Contact_MF = matrix(MF,nrow=16,ncol=16)
Contact_FF = matrix(FF,nrow=16,ncol=16)

# preliminaries
DeltaA = 5
A = 16
BirthContribution = c(0, 0, 0, 0.0101, 0.0893, 0.2680, 0.3788, 0.2173, 0.0351, 0.0014, 0, 0, 0, 0, 0, 0) #www.cbs.nl
num_knots = 5 
spline_degree = 3
knots = seq(from = 0, to = DeltaA*A, by = DeltaA*A/(num_knots-1))
ts = 1:(DeltaA*A)
mode = 0 # 0=normal sampling, 1=WBIC sampling

# data dictionary and initial values for Stan
data_values = list(
  'N' = N,
  'A' = A,
  'DeltaA' = DeltaA,
  'Titers' = Titers,
  'Ages' = Ages,
  'Censor' = Censor,
  'BirthContribution' = BirthContribution,
  'RightCensor' = RightCensor,
  'Contact_MM' = Contact_MM,
  'Contact_FM' = Contact_FM,
  'Contact_MF' = Contact_MF,
  'Contact_FF' = Contact_FF,
  'MuS' = -1.68,
  'MuL' = 1.03,
  'MuB' = 2.6,
  'SigmaS' = 1 / 6.27,
  'SigmaL' = 1 / 1.52,
  'SigmaB' = 1 / 2.03,
  'Penalty' = 1e4,
  'Gender' = Gender,
  'ts' = ts,
  'numbertestedinfants' = 31484,
  'numbercCMVinfants' = 154,
  'spline_degree' = spline_degree,
  'num_knots' = num_knots,
  'knots' = knots,
  'reducinf' = 1.0,
  'mode' = mode
) 

# initial values for parameters 
parameter_inits = function(){
  return(list(beta1 = 0.002,
              beta2 = 0.02,  
              lambda_f = c(0.015,0.019,0.021,0.016,0.011,0.012,0.013,0.013,0.013,0.012,0.011,0.009,0.0081,0.00083,0.0071,0.0058),
              lambda_m = c(0.014,0.019,0.017,0.013,0.0091,0.01,0.012,0.014,0.013,0.011,0.0094,0.0083,0.009,0.0076,0.0066,0.0065),
              z = 0.4,
              probLtoB = 0.5,
              S0 = 0.80,
              qcCMV = 0.05,
              nu = 0.4,
              a_raw = matrix(c(0.002,0.002,0.01,0.01,0.02,0.03,0.03,0.002,0.002,0.01,0.01,0.02,0.03,0.03), nrow = 2, ncol = num_knots + spline_degree - 1)
        ) 
  )
}

# estimate parameters with Stan
# main analysis with gamma(2,50) priors for spline weights
fit = stan(file = 'cmv_22082019.stan', 
  data = data_values, 
  init = parameter_inits, 
  iter = 2000, 
  warmup = 1000, 
  thin = 20, 
  chains = 20,
  control = list(adapt_delta = 0.97, max_treedepth = 15)
)

# plot traces
traceplot(
  fit,
  pars = c(
    "beta1",
    "beta2",
    "z",
    "S0",
    "nu",
    "qcCMV",
    "probLtoB",
    "rho_f[10]",
    "rho_f[40]",
    "rho_f[70]"
  )
)
traceplot(fit, pars=c("rho_f"))
traceplot(fit, pars=c("rho_m"))
traceplot(fit, pars=c("a_raw"))
traceplot(fit, pars=c("lambda_f"))
traceplot(fit, pars=c("lambda_m"))


# pairplots
pairs(
  fit,
  pars = c(
    "beta1",
    "beta2",
    "z",
    "S0",
    "nu",
    "probLtoB",
    "qcCMV"
  )
)

# some diagnostics
stan_diag(fit)

# some parameter outputs 
print(
  fit,
  pars = c(
    "beta1",
    "beta2",
    "z",
    "S0",
    "nu",
    "probLtoB",
    "rho_m",
    "rho_f",
    "a_raw"
  ),
  digits = 5
)

# extract parameters
params = rstan::extract(fit)

# select sample with highest log posterior likelihood (MAP)
hist(params$lp__, nclass = 30)
max(params$lp__)
posmax <- as.numeric(which(params$lp__ == max(params$lp__)))
MAP <- as.data.frame(fit)[posmax, c("beta1",
                                    "beta2",
                                    "z",
                                    "S0",
                                    "nu",
                                    "probLtoB",
                                    "qcCMV")] 
MAP

# save fit
save(fit, file = "cmv_22082019.rda")

# checks
rm(fit)
load("xxx.rda")
load(file = "high reactivation_10052019.rda")
print(fit) 

# calculate WAIC
LL = extract_log_lik(fit, parameter_name = 'log_lik')
waic(LL) # now obsolete
loo(LL) # better
write.csv(LL, file = "logliks.csv", row.names = FALSE)

# calculate WBIC; note set mode=1 in separate run
print(fit, pars="log_like", digits=5) 

# write selected output
output = as.data.frame(fit, pars=c("beta1", "beta2", 
                                   "probLtoB", 
                                   "z", "S0", "nu", "qcCMV",
                                   "rho_m", "rho_f",
                                   "lambda_m", "lambda_f", 
                                   "S_m", "S_f", "L_m", "L_f", "B_m", "B_f", "a_raw")) 

write.csv(output, file = "high reactivation (10052019).csv", row.names = FALSE)

# plot of reactivation rates
col2rgb('grey80', alpha=TRUE) 
greytrans <- rgb(204,204,204, 127, maxColorValue=255)

rho_med_f <- array(NA, DeltaA*A)
rho_ub_f <- array(NA, DeltaA*A)
rho_lb_f <- array(NA, DeltaA*A)
rho_med_m <- array(NA, DeltaA*A)
rho_ub_m <- array(NA, DeltaA*A)
rho_lb_m <- array(NA, DeltaA*A)
for (j in 1:80) {
  rho_med_f[j] <- median(params$rho_f[,j])
  rho_lb_f[j] <- quantile(params$rho_f[,j],probs = 0.025)
  rho_ub_f[j] <- quantile(params$rho_f[,j],probs = 0.975)
  rho_med_m[j] <- median(params$rho_m[,j])
  rho_lb_m[j] <- quantile(params$rho_m[,j],probs = 0.025)
  rho_ub_m[j] <- quantile(params$rho_m[,j],probs = 0.975)
} 
df_rho <- data.frame(list(age = ts, median_f = rho_med_f, lower_f = rho_lb_f, upper_f = rho_ub_f, 
                          median_m = rho_med_m, lower_m = rho_lb_m, upper_m = rho_ub_m))

plot_rho <- ggplot(data = df_rho, aes(x = age-0.5)) + 
  geom_line(aes(y = median_f, color = "female")) +
  geom_line(aes(y = median_m, color = "male")) +
  scale_colour_manual("", values = c("female"="purple", "male"="brown")) +
  geom_ribbon(aes(ymin = lower_f, ymax = upper_f), alpha = 0.3, fill = "purple") +
  geom_ribbon(aes(ymin = lower_m, ymax = upper_m), alpha = 0.3, fill = "brown") +
  labs(x = "Age (yr)" , y = expression(paste("Reactivation rate (", yr^-1, ")", sep = ""))) +  
  scale_x_continuous(limits = c(0, 80), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(0, 0.075), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  theme(legend.background = element_rect(fill="transparent"), legend.position=c(0.13,0.95)) 
 
plot_rho 

# plot female prevalences
S_med_f <- array(NA, DeltaA*A)
S_ub_f <- array(NA, DeltaA*A)
S_lb_f <- array(NA, DeltaA*A)
L_med_f <- array(NA, DeltaA*A)
L_ub_f <- array(NA, DeltaA*A)
L_lb_f <- array(NA, DeltaA*A)
B_med_f <- array(NA, DeltaA*A)
B_ub_f <- array(NA, DeltaA*A)
B_lb_f <- array(NA, DeltaA*A)

for (j in 1:80) {
  S_med_f[j] <- median(params$S_f[,j])
  S_lb_f[j] <- quantile(params$S_f[,j], probs = 0.025)
  S_ub_f[j] <- quantile(params$S_f[,j], probs = 0.975)
  L_med_f[j] <- median(params$L_f[,j])
  L_lb_f[j] <- quantile(params$L_f[,j], probs = 0.025)
  L_ub_f[j] <- quantile(params$L_f[,j], probs = 0.975)
  B_med_f[j] <- median(params$B_f[,j])
  B_lb_f[j] <- quantile(params$B_f[,j], probs = 0.025)
  B_ub_f[j] <- quantile(params$B_f[,j], probs = 0.975)
} 

df_prev_f <- data.frame(list(age = ts, 
                             median_S_f = S_med_f, lower_S_f = S_lb_f, upper_S_f = S_ub_f, 
                             median_L_f = L_med_f, lower_L_f = L_lb_f, upper_L_f = L_ub_f, 
                             median_B_f = B_med_f, lower_B_f = B_lb_f, upper_B_f = B_ub_f
)
) 

plot_prev_f <- ggplot(data = df_prev_f, aes(x= age-0.5)) + 
  geom_line(aes(y = median_S_f, color = "S")) +
  geom_line(aes(y = median_L_f, color = "L")) +
  geom_line(aes(y = median_B_f, color = "B")) +
  scale_colour_manual("", values = c("S"="blue", "L"="purple", "B"="red")) +
  geom_ribbon(aes(ymin = lower_S_f, ymax = upper_S_f), alpha = 0.3, fill = "blue") +
  geom_ribbon(aes(ymin = lower_L_f, ymax = upper_L_f), alpha = 0.3, fill = "purple") + 
  geom_ribbon(aes(ymin = lower_B_f, ymax = upper_B_f), alpha = 0.3, fill = "red") +
  labs(x = "Age (yr)" , y = "Prevalence") +  
  scale_x_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  theme(legend.background = element_rect(fill="transparent"), legend.position=c(0.9,0.90)) 

plot_prev_f

# plot male prevalences
S_med_m <- array(NA, DeltaA*A)
S_ub_m <- array(NA, DeltaA*A)
S_lb_m <- array(NA, DeltaA*A)
L_med_m <- array(NA, DeltaA*A)
L_ub_m <- array(NA, DeltaA*A)
L_lb_m <- array(NA, DeltaA*A)
B_med_m <- array(NA, DeltaA*A)
B_ub_m <- array(NA, DeltaA*A)
B_lb_m <- array(NA, DeltaA*A)

for (j in 1:(DeltaA*A)) {
  S_med_m[j] <- median(params$S_m[,j])
  S_lb_m[j] <- quantile(params$S_m[,j], probs = 0.025)
  S_ub_m[j] <- quantile(params$S_m[,j], probs = 0.975)
  L_med_m[j] <- median(params$L_m[,j])
  L_lb_m[j] <- quantile(params$L_m[,j], probs = 0.025)
  L_ub_m[j] <- quantile(params$L_m[,j], probs = 0.975)
  B_med_m[j] <- median(params$B_m[,j])
  B_lb_m[j] <- quantile(params$B_m[,j], probs = 0.025)
  B_ub_m[j] <- quantile(params$B_m[,j], probs = 0.975)
} 

df_prev_m <- data.frame(list(age = ts, 
                             median_S_m = S_med_m, lower_S_m = S_lb_m, upper_S_m = S_ub_m, 
                             median_L_m = L_med_m, lower_L_m = L_lb_m, upper_L_m = L_ub_m, 
                             median_B_m = B_med_m, lower_B_m = B_lb_m, upper_B_m = B_ub_m
) 
)

plot_prev_m <- ggplot(data = df_prev_m, aes(x= age-0.5)) + 
  geom_line(aes(y = median_S_m, color = "S")) +
  geom_line(aes(y = median_L_m, color = "L")) +
  geom_line(aes(y = median_B_m, color = "B")) +
  scale_colour_manual("", values = c("S"="blue", "L"="purple", "B"="red")) +
  geom_ribbon(aes(ymin = lower_S_m, ymax = upper_S_m), alpha = 0.3, fill = "blue") +
  geom_ribbon(aes(ymin = lower_L_m, ymax = upper_L_m), alpha = 0.3, fill = "purple") + 
  geom_ribbon(aes(ymin = lower_B_m, ymax = upper_B_m), alpha = 0.3, fill = "red") +
  labs(x = "Age (yr)" , y = "Prevalence") +  
  scale_x_continuous(limits = c(0, 80), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  theme(legend.background = element_rect(fill="transparent"), legend.position=c(0.9,0.92)) 

plot_prev_m

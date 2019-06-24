library(INLA)
library(ggplot2)
library(knitr)
library(kableExtra)
library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- get(load("./assignment1.RData"))
# show head of the data
kable(head(data)) %>% kable_styling()

# add numeric column super_region_code to data for easier calculation
data <- data %>% mutate(super_region_code = as.numeric(as.factor(data$Super_region_name_1)))

# show super_region_name - super_region_code
kable(data.frame(super_region_name = unique(data$Super_region_name_1), 
                 super_region_code = as.numeric(as.factor(unique(data$Super_region_name_1))))) %>% kable_styling()

# scatter plot - global simple linear model
ggplot(data, aes(x=logSAT, y=logPM25)) + 
  geom_point(aes(color = Super_region_name_1)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed", color ="black")+ labs(color='The Seven Super Regions') 

# scatter plot - fit a linear line for each region
scatter_plot <- ggplot(data, aes(x=logSAT, y=logPM25)) + 
  geom_point(aes(color = Super_region_name_1)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed", color ="black") +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(color = Super_region_name_1))


# tiny R square for slr for each region
output = data.frame(super_region = rep(NA,7), R2 = rep(NA,7))
dat = data.frame(logpm25 = data$logPM25, logSAT = data$logSAT,
                 super_region = data$super_region_code, super_name = data$Super_region_name_1)
for (i in 1:7){
  dat_local = dat %>% filter(super_region == i)
  output$R2[i] = summary(lm(logpm25 ~ logSAT, data = dat_local))$r.squared
  output$super_region[i] = unique(dat_local$super_name)
}
print(output, digits=2)


# Plot: prior with default priors
set.seed(seed = 10000011)

sigma_ <- 1 / sqrt(rgamma(1, 1, rate = 100))
beta0j <- rnorm(185, 0, sigma_) # 185 countries 
beta0k <- rnorm(7, 0, sigma_) # 7 super regions
beta1j <- rnorm(185, 0, sigma_) 
beta1k <- rnorm(7, 0, sigma_)
beta0 <- rnorm(1, 0, 1000)
beta1 <- rnorm(1, 0, 100)

Nsim <- length(data$country_code_1)

data1 <- data.frame(log_pm25=data$logPM25, 
                    sim=beta0 + beta0j[data$country_code_1] + beta0k[data$super_region_code] +
                      (beta1 + beta1j[data$country_code_1] + beta1k[data$super_region_code]) * data$logSAT +
                      rnorm(Nsim, mean = 0, sd = sigma_))


xysim_labs <- labs(
  x = expression(paste("Observed ", log(PM[2.5]))),
  y = "Simulated data"
)

ggplot(data=data1, aes(x = log_pm25, y = sim)) + 
  geom_point(alpha = 0.1, color = "red") + 
  xysim_labs


# Plot: prior predictive with weakly informative priors
set.seed(seed = 10000011)
sigma <- abs(rnorm(1, 0, 1))
beta0j <- rnorm(185, 0, sigma) # 185 countries 
beta0k <- rnorm(7, 0, sigma) # 7 super regions
beta1j <- rnorm(185, 0, sigma) 
beta1k <- rnorm(7, 0, sigma)
beta0 <- rnorm(1, 0, 1)
beta1 <- rnorm(1, 0, 1)

data2 <- data.frame(log_pm25=data$logPM25, 
                    sim = beta0 + beta0j[data$country_code_1] + beta0k[data$super_region_code] +
                      (beta1 + beta1j[data$country_code_1] + beta1k[data$super_region_code]) * data$logSAT +
                      rnorm(Nsim, mean = 0, sd = sigma))

ggplot(data2, aes(x = log_pm25, y = sim)) +
  geom_point(alpha = 0.1) + 
  xysim_labs




# full model
g = "./world.adj"
full_Super_region_inter = data$super_region_code
full_Super_region_slope = data$super_region_code
full_country_code_inter = data$country_code_1
full_country_code_slope = data$country_code_1

fullmod = logPM25 ~ logSAT +
  f(full_country_code_inter, model="bym2", graph=g,
    hyper = list(prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_inter, model="iid",hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_country_code_slope, logSAT, model="bym2", graph=g,
    hyper = list(prec = list( prior="logtnormal", param=c(0, 1))))

fullresult = inla(fullmod, data = data, family = "gaussian", 
                  control.compute= list(dic=TRUE, cpo = TRUE,config = TRUE), 
                  control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))



# iid model
iid_Super_region_inter = data$super_region_code
iid_Super_region_slope = data$super_region_code

iidmod = logPM25 ~ logSAT + 
  f(iid_Super_region_inter, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(iid_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1))))

iidresult = inla(iidmod, data = data, family = "gaussian", 
                 control.compute= list(dic=TRUE, cpo = TRUE,config = TRUE),
                 control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))



# simple linear regression
linear_mod = logPM25 ~ logSAT
linearresult = inla(linear_mod, data = data, family = "gaussian", 
                    control.compute= list(dic=TRUE, cpo = TRUE,config = TRUE),
                    control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))




# pit for the models
hist(fullresult$cpo$pit, main = "PIT for full Model")
hist(iidresult$cpo$pit, main = "PIT for iid model")
hist(linearresult$cpo$pit, main = "PIT for SLR model")

# cpo for the models
cpofull = -0.5*sum(log(fullresult$cpo$cpo))
cpoiid = -0.5*sum(log(iidresult$cpo$cpo))
cposlr = -0.5*sum(log(linearresult$cpo$cpo))

# super region name - super region code - number of observations
kable(data %>% group_by(Super_region_name_1)%>% count() %>% 
  left_join(data.frame(Super_region_name_1 = unique(data$Super_region_name_1), 
             super_region_code = as.numeric(as.factor(unique(data$Super_region_name_1)))),by = "Super_region_name_1") )





# prepare data for cv1
set.seed(2);sample(c(1,3,4,5,7), 2) # produce 1,4, Randomly put 500 from 1 and 100 from 4 to test set. 
set.seed(1000);takeout1 = sample(seq(1:518),500)
test_data_wout_1 = (data %>% filter(super_region_code == 1))[takeout1,]
set.seed(1000);takeout4 = sample(seq(1:273),100)
test_data_wout_4 = (data %>% filter(super_region_code == 4))[takeout4,]

test_data_wout_14 = rbind(test_data_wout_1, test_data_wout_4)

train_data_wout_14 = rbind((data %>% filter(super_region_code != 1, super_region_code != 4)), 
                           (data %>% filter(super_region_code == 4))[-takeout4,],
                           (data %>% filter(super_region_code == 1))[-takeout1,])

# ___________________________________________________cv1____________________________________
# full cv1
full_Super_region_inter = train_data_wout_14$super_region_code
full_Super_region_slope = train_data_wout_14$super_region_code
full_country_code_inter = train_data_wout_14$country_code_1
full_country_code_slope = train_data_wout_14$country_code_1
g = "./world.adj"
fullmod_CV1 = logPM25 ~ logSAT +
  f(full_country_code_inter, model="bym2", graph=g, hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_inter, model="iid",hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_country_code_slope, logSAT, model="bym2", graph=g, hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1))))

fullresult_CV1 = inla(fullmod_CV1, data = train_data_wout_14, family = "gaussian", 
                      control.compute=list(config = TRUE),
                      control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))

# 1000 sample from posterior
set.seed(1000)
samples11 = inla.posterior.sample(n=1000, result = fullresult_CV1) 


# FIND THE INTERCEPT index
intercept_index11 = grep("Intercept", rownames(samples11[[1]]$latent))
intercept_super_index11 = grep("full_Super_region_inter", rownames(samples11[[1]]$latent))
intercept_country_index11 = grep("full_country_code_inter", rownames(samples11[[1]]$latent))[1:185]


## Same for slopes

slope_index11 = grep("logSAT", rownames(samples11[[1]]$latent))
slope_super_index11 = grep("full_Super_region_slope", rownames(samples11[[1]]$latent))
slope_country_index11 = grep("full_country_code_slope", rownames(samples11[[1]]$latent))[1:185]


error11 = 0
for (i in 1:600){
  country = test_data_wout_14[i,5]
  logsat = test_data_wout_14[i,3]
  logpm = test_data_wout_14[i,2]
  super = test_data_wout_14[i,7]
  estimates = rep(NA,1000)
  for (j in 1:1000){
    x = samples11[[j]]
    intercept = (x$latent[intercept_index11] + x$latent[intercept_super_index11][super]
                 + x$latent[intercept_country_index11][country])
    slope = (x$latent[slope_index11] + x$latent[slope_super_index11][super] + x$latent[slope_country_index11][country])
    estimate = intercept + slope*logsat
    estimates[j] =  estimate
  }
  error11 = error11 + (mean(estimates) - logpm)^2
}
error11/600 




# iid CV1
iid_Super_region_inter = train_data_wout_14$super_region_code
iid_Super_region_slope = train_data_wout_14$super_region_code

iidmod_cv1 = logPM25 ~ logSAT  + 
  f(iid_Super_region_inter, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(iid_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1))))

iidresult_cv1 = inla(iidmod_cv1, data = train_data_wout_14, family = "gaussian", 
                     control.compute=list(config = TRUE),
                     control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))
samples21 = inla.posterior.sample(n=1000, result = iidresult_cv1)

# FIND THE INTERCEPT index 
intercept_index21 = grep("Intercept", rownames(samples21[[1]]$latent))
intercept_random_index21 = grep("iid_Super_region_inter", rownames(samples21[[1]]$latent))

## Same for slopes
slope_index21 = grep("logSAT", rownames(samples21[[1]]$latent))
slope_random_index21 = grep("iid_Super_region_slope", rownames(samples21[[1]]$latent))


error21 = 0
for (i in 1:600){
  country = test_data_wout_14[i,5]
  logsat = test_data_wout_14[i,3]
  logpm = test_data_wout_14[i,2]
  super = test_data_wout_14[i,7]
  estimates = rep(NA, 1000)
  for (j in 1:1000){
    x = samples21[[j]]
    inter = x$latent[intercept_index21] + x$latent[intercept_random_index21][super]
    slope = x$latent[slope_index21] + x$latent[slope_random_index21][super]
    estimate = inter + slope*logsat
    estimates[j] = estimate
  }
  error21 = error21 + (mean(estimates) - logpm)^2
}
error21/600




# slr CV1
linear_mod_CV1 = logPM25 ~ logSAT
linearresult_CV1 = inla(linear_mod_CV1, data = train_data_wout_14, family = "gaussian",
                        control.compute=list(config = TRUE),
                        control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))
set.seed(1000); samples31 = inla.posterior.sample(n=1000, result = linearresult_CV1)


# FIND THE INTERCEPT index
intercept_index31 = grep("Intercept", rownames(samples31[[1]]$latent))

## Same for slopes
slope_index31 = grep("logSAT", rownames(samples31[[1]]$latent))

error31 = 0
for (i in 1:600){
  country = test_data_wout_14[i,5]
  logsat = test_data_wout_14[i,3]
  logpm = test_data_wout_14[i,2]
  super = test_data_wout_14[i,7]
  estimates = rep(NA, 1000)
  for (j in 1:1000){
    x = samples31[[j]]
    inter = x$latent[intercept_index31]
    slope = x$latent[slope_index31]
    estimate = inter + slope*logsat
    estimates[j] = estimate
  }
  error31 = error31 + (mean(estimates) - logpm)^2
}
error31/600 
















# prepare data for cv2
set.seed(5);sample(c(1,3,4,5,7), 2) # produce 3,4 

set.seed(10001)
takeout3 = sample(seq(1:308),300) # Randomly put 270 observations from SR4 and 300 observations from SR3 to test set.  
test_data_wout_3 = (data %>% filter(super_region_code == 3))[takeout3,]
takeout4 = sample(seq(1:273),270)
test_data_wout_4 = (data %>% filter(super_region_code == 4))[takeout4,]

test_data_wout_34 = rbind(test_data_wout_3, test_data_wout_4)
View(test_data_wout_34)

train_data_wout_34 = rbind((data %>% filter(super_region_code != 3, super_region_code != 4)), 
                           (data %>% filter(super_region_code == 4))[-takeout4,],
                           (data %>% filter(super_region_code == 3))[-takeout3,])
View(train_data_wout_34)


#  ------------------------------------cv2-------------------------------------------
# cv2 for full model
full_Super_region_inter = train_data_wout_34$super_region_code
full_Super_region_slope = train_data_wout_34$super_region_code
full_country_code_inter = train_data_wout_34$country_code_1
full_country_code_slope = train_data_wout_34$country_code_1
g = "./world.adj"
fullmod_CV2 = logPM25 ~ logSAT +
  f(full_country_code_inter, model="bym2", graph=g, hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_inter, model="iid",hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_country_code_slope, logSAT, model="bym2", graph=g, hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1))))

fullresult_CV2 = inla(fullmod_CV2, data = train_data_wout_34, family = "gaussian", 
                      control.compute=list(config = TRUE),
                      control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))

set.seed(1000);samples12 = inla.posterior.sample(n=1000, result = fullresult_CV2) # 1000 sample from posterior 

intercept_index12 = grep("Intercept", rownames(samples12[[1]]$latent))
intercept_super_index12 = grep("full_Super_region_inter", rownames(samples12[[1]]$latent))
intercept_country_index12 = grep("full_country_code_inter", rownames(samples12[[1]]$latent))[1:185]

slope_index12 = grep("logSAT", rownames(samples12[[1]]$latent))
slope_super_index12 = grep("full_Super_region_slope", rownames(samples12[[1]]$latent))
slope_country_index12 = grep("full_country_code_slope", rownames(samples12[[1]]$latent))[1:185]


error12 = 0
for (i in 1:570){
  country = test_data_wout_34[i,5]
  logsat = test_data_wout_34[i,3]
  logpm = test_data_wout_34[i,2]
  super = test_data_wout_34[i,7]
  estimates = rep(NA, 1000)
  for (j in 1:1000){
    x = samples12[[j]]
    intercept = (x$latent[intercept_index12] + x$latent[intercept_super_index12][super]
                 + x$latent[intercept_country_index12][country])
    slope = (x$latent[slope_index12] + x$latent[slope_super_index12][super] + x$latent[slope_country_index12][country])
    estimate = intercept + slope*logsat
    estimates[j] = estimate
  }
  error12 = error12 + (mean(estimate) - logpm)^2
}
error12/570 









# iid cv2

iid_Super_region_inter = train_data_wout_34$super_region_code
iid_Super_region_slope = train_data_wout_34$super_region_code

iidmod_cv2 = logPM25 ~ logSAT  + 
  f(iid_Super_region_inter, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(iid_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1))))

iidresult_cv2 = inla(iidmod_cv2, data = train_data_wout_34, family = "gaussian", 
                     control.compute=list(config = TRUE),
                     control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))
set.seed(1000);samples22 = inla.posterior.sample(n=1000, result = iidresult_cv2)


# FIND THE INTERCEPT!!! 
intercept_index22 = grep("Intercept", rownames(samples22[[1]]$latent))
intercept_random_index22 = grep("iid_Super_region_inter", rownames(samples22[[1]]$latent))

## Same for slopes
slope_index22 = grep("logSAT", rownames(samples22[[1]]$latent))
slope_random_index22 = grep("iid_Super_region_slope", rownames(samples22[[1]]$latent))


error22 = 0
for (i in 1:570){
  country = test_data_wout_34[i,5]
  logsat = test_data_wout_34[i,3]
  logpm = test_data_wout_34[i,2]
  super = test_data_wout_34[i,7]
  estimates = rep(NA, 1000)
  for (j in 1:1000){
    x = samples22[[j]]
    inter = x$latent[intercept_index22] + x$latent[intercept_random_index22][super]
    slope = x$latent[slope_index22] + x$latent[slope_random_index22][super]
    estimate = inter + slope*logsat
    estimates[j] = estimate
  }
  error22 = error22 + (mean(estimates) - logpm)^2
}
error22/570 





# slr cv2
linear_mod_CV2 = logPM25 ~ logSAT
linearresult_CV2 = inla(linear_mod_CV2, data = train_data_wout_34, family = "gaussian",
                        control.compute=list(config = TRUE),
                        control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))
set.seed(1000);samples32 = inla.posterior.sample(n=1000, result = linearresult_CV2)


# FIND THE INTERCEPT index
intercept_index32 = grep("Intercept", rownames(samples32[[1]]$latent))

## Same for slopes
slope_index32 = grep("logSAT", rownames(samples32[[1]]$latent))


error32 = 0
for (i in 1:570){
  country = test_data_wout_34[i,5]
  logsat = test_data_wout_34[i,3]
  logpm = test_data_wout_34[i,2]
  super = test_data_wout_34[i,7]
  estimates = rep(NA, 1000)
  for (j in 1:1000){
    x = samples32[[j]]
    inter = x$latent[intercept_index32]
    slope = x$latent[slope_index32]
    estimate = inter + slope*logsat
    estimates[j] = estimate
  }
  error32 = error32 + (mean(estimates) - logpm)^2
}
error32/570





# prepare data for cv3
set.seed(21);sample(c(1,3,4,5,7), 2) # produce 3,5

set.seed(100011); takeout3 = sample(seq(1:308),300)
test_data_wout_3 = (data %>% filter(super_region_code == 3))[takeout3,]
set.seed(100011); takeout5 = sample(seq(1:464),300)
test_data_wout_5 = (data %>% filter(super_region_code == 5))[takeout5,]

test_data_wout_35 = rbind(test_data_wout_3, test_data_wout_5)


train_data_wout_35 = rbind((data %>% filter(super_region_code != 3, super_region_code != 5)), 
                           (data %>% filter(super_region_code == 5))[-takeout5,],
                           (data %>% filter(super_region_code == 3))[-takeout3,])
View(train_data_wout_35)
View(test_data_wout_35)
#-------------------------------------cv3-------------------------------------
# full cv3
full_Super_region_inter = train_data_wout_35$super_region_code
full_Super_region_slope = train_data_wout_35$super_region_code
full_country_code_inter = train_data_wout_35$country_code_1
full_country_code_slope = train_data_wout_35$country_code_1

fullmod_CV3 = logPM25 ~ logSAT +
  f(full_country_code_inter, model="bym2", graph=g, hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_inter, model="iid",hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(full_country_code_slope, logSAT, model="bym2", graph=g, hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1))))

fullresult_CV3 = inla(fullmod_CV3, data = train_data_wout_35, family = "gaussian", 
                      control.compute=list(config = TRUE),
                      control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))
samples13 = inla.posterior.sample(n=1000, result = fullresult_CV3) # 1000 sample from posterior 

# FIND THE INTERCEPT index
intercept_index13 = grep("Intercept", rownames(samples13[[1]]$latent))
intercept_super_index13 = grep("full_Super_region_inter", rownames(samples13[[1]]$latent))
intercept_country_index13 = grep("full_country_code_inter", rownames(samples13[[1]]$latent))[1:185]


## Same for slopes

slope_index13 = grep("logSAT", rownames(samples13[[1]]$latent))
slope_super_index13 = grep("full_Super_region_slope", rownames(samples13[[1]]$latent))
slope_country_index13 = grep("full_country_code_slope", rownames(samples13[[1]]$latent))[1:185]


error13 = 0
for (i in 1:600){
  country = test_data_wout_35[i,5]
  logsat = test_data_wout_35[i,3]
  logpm = test_data_wout_35[i,2]
  super = test_data_wout_35[i,7]
  estimates = rep(NA,1000)
  for (j in 1:1000){
    x = samples13[[j]]
    intercept = (x$latent[intercept_index13] + x$latent[intercept_super_index13][super]
                 + x$latent[intercept_country_index13][country])
    slope = (x$latent[slope_index13] + x$latent[slope_super_index13][super] + x$latent[slope_country_index13][country])
    estimate = intercept + slope*logsat
    estimates[j] =  estimate
  }
  error13 = error13 + (mean(estimates) - logpm)^2
}
error13/600


# iid cv3

iid_Super_region_inter = train_data_wout_35$super_region_code
iid_Super_region_slope = train_data_wout_35$super_region_code

iidmod_cv3 = logPM25 ~ logSAT  + 
  f(iid_Super_region_inter, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1)))) + 
  f(iid_Super_region_slope, logSAT, model="iid", hyper = list(
    prec = list( prior="logtnormal", param=c(0, 1))))

iidresult_cv3 = inla(iidmod_cv3, data = train_data_wout_35, family = "gaussian", 
                     control.compute=list(config = TRUE),
                     control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))
set.seed(1000);samples23 = inla.posterior.sample(n=1000, result = iidresult_cv3)


# FIND THE INTERCEPT index 
intercept_index23 = grep("Intercept", rownames(samples23[[1]]$latent))
intercept_random_index23 = grep("iid_Super_region_inter", rownames(samples23[[1]]$latent))

## Same for slopes
slope_index23 = grep("logSAT", rownames(samples23[[1]]$latent))
slope_random_index23 = grep("iid_Super_region_slope", rownames(samples23[[1]]$latent))

mse_iid_cv1 = 0
error23 = 0
for (i in 1:600){
  country = test_data_wout_35[i,5]
  logsat = test_data_wout_35[i,3]
  logpm = test_data_wout_35[i,2]
  super = test_data_wout_35[i,7]
  estimates = rep(NA, 1000)
  for (j in 1:1000){
    x = samples23[[j]]
    inter = x$latent[intercept_index23] + x$latent[intercept_random_index23][super]
    slope = x$latent[slope_index23] + x$latent[slope_random_index23][super]
    estimate = inter + slope*logsat
    estimates[j] = estimate
  }
  error23 = error23 + (mean(estimates) - logpm)^2
}
error23/600 


# slr cv3

linear_mod_CV3 = logPM25 ~ logSAT
linearresult_CV3 = inla(linear_mod_CV3, data = train_data_wout_35, family = "gaussian", 
                        control.compute=list(config = TRUE),
                        control.fixed = list(mean = 0, prec=1, mean.intercept=0, prec.intercept=1))


# FIND THE INTERCEPT index
intercept_index33 = grep("Intercept", rownames(samples33[[1]]$latent))

## Same for slopes
slope_index33 = grep("logSAT", rownames(samples33[[1]]$latent))

error33 = 0
for (i in 1:600){
  country = test_data_wout_35[i,5]
  logsat = test_data_wout_35[i,3]
  logpm = test_data_wout_35[i,2]
  super = test_data_wout_35[i,7]
  estimates = rep(NA, 1000)
  for (j in 1:1000){
    x = samples33[[j]]
    inter = x$latent[intercept_index33]
    slope = x$latent[slope_index33]
    estimate = inter + slope*logsat
    estimates[j] = estimate
  }
  error33 = error33 + (mean(estimates) - logpm)^2
}
error33/600 


errorfull = (error11 + error12 +error13)/600/3

erroriid = (error21 + error22 + error23)/570/3

errorslr = (error31 + error32 + error33)/1800


data.frame(model=c("full","iid","SLR"), MSE=c(errorfull, erroriid, errorslr)) %>% kable() %>% kable_styling()








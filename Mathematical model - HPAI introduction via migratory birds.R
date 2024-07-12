###Highly Pathogenic Avian Influenza transmission at the waterfowl-backyard poultry farm interface###
start_time <- Sys.time()
n <- 1000 #Number of iterations (Monte Carlo)
set.seed(145689)

#Packages needed
library(beepr)
library(deSolve)
library(EnvStats)
library(ggplot2)
library(gridExtra)
library(epiR)
library(RColorBrewer)
library(hrbrthemes)
library(scales)
library(IRdisplay)
library(mc2d)
library(TruncatedDistributions)
library(ggpubr)
library(rpart)				        # Popular decision tree algorithm
library(rattle)					# Fancy tree plot
library(rpart.plot)				# Enhanced tree plots
library(RColorBrewer)				# Color selection for fancy tree plot
library(caret)

initial.time=Sys.time() #Stores the time the simulation begins (this is then used to determine how long it takes for the model tu run)

#=============Initial Values======================
#Populations
Nmbar<-750 #Resident mallards
Np<-192    #Poultry farms
Nw<-118    #Migratory mute swans
Nm<-2250   #Migratory mallards

prev.res.mal<- 0 #Prevalence of infection among resident mallards at the beginning of the simulation

AI_MultiSpecies.init <-c(
  E_init=0,     #Environment
  S_w_init=0,   #Susceptible mute swan
  I_w_init=0,   #Infectious mute swan
  D_w_init=0,   #Dead mute swan
  
  S_m_init=0,   #Susceptible migratory mallards
  I_m_init=0,   #Infectious migratory mallards
  R_m_init=0,   #Recovered migratory mallards
  D_m_init=0,   #Dead migratory mallards
  
  S_mbar_init=Nmbar*(1-prev.res.mal), #Susceptible resident mallards
  I_mbar_init=Nmbar*(prev.res.mal),   #Infectious resident mallards
  R_mbar_init=0,                      #Recovered resident mallards
  D_mbar_init=0,                      #Dead resident mallards
  
  S_p_init=Np,                        #Infection-free poultry farms
  I_p_init=0,                         #HPAI-affected poultry farms
  D_p_init=0                          #Culled poultry farms
 ) 

#=============Summary=========================
#This model assess the transmission dynamics of highly pathogenic avian influenza (HPAI) among wild birds
#and backyard poultry farms. Its main objective is to understand how timing of migration and the duration of the stopover period of mallards
#influences the probability of HPAI infection in backyard poultry farms. The model consider 4 important animal components: 
#1) migratory mute swans, 2) migratory mallards, 3) resident mallards, and 4) backyard poultry farms. 
#Additionally, it considers an environmental compartment where HPAI survives for a short period.
#The model accounts for intras-pecific direct transmission of HPAI in wild birds (frequency-dependent) and indirect transmission
#between HPAI in the environment and birds.
#Model unit of time = day
#=============Parameters======================


epsilon_w <-rpert(n, min=10^3.23, mode=10^3.23, max=10^3.23)*2031 #HPAI shedding rate by mute swans (virions per g of feces)

epsilon_m <- rpert(n, min=10^2, mode=10^2.5, max=10^3.3)*48  #HPAI shedding rate by mallards swans (virions per g of feces)

r<- 1.23 #Rate of HPAI removal from the environment

rho_sw <- 1  #Proportion of migratory mute swans arriving as susceptible (100%)

day_arrival_muteswans <- as.integer(rpert(n, min=0, max=60, mode=30)) #Arrival day of swan fall migration (Sept 1st to Oct 31). The model starts on Sept 1st (day 0)

migdays_w <- 1 #Length of the migration event for mute swans (over how many days mute swans migrate in/out)

b_w <- Nw #Total number of swans that migrate into the area over "migdays_w"

transmission_rate<- 0.50 #Transmission of HPAI between mallards (per bird per day)

density_times_lower_than_mallards = 2.5 #How much lower is population density in mute swans compared to mallards

beta_ww <- transmission_rate/density_times_lower_than_mallards #HPAI transmission rate between mute swans (per bird per day)

mu_w <- rpert(n, min=0.0007, max=0.0009, mode=0.0008) #Natural mortality rate for mute swans (per bird per day)

alpha <- 10^-10.0 #Indirect transmission rate between an infectious dose of HPAI in the environment and wild birds or poultry (per virion per bird  per day)

eta_w <- 10^0.95 #HPAI infectious dose for mute swans (virion)

rho_iw <- 1-rho_sw #Proportion of migratory mute swans that arrived HPAI-infected to Croatia (considered to be zero)

d_w <- 1/rpert(n, min = 3.0, max = 5.5, mode = 4.5)  #Death rate from HPAI infection in mute swans (per bird per day)

day_departure_muteswans <- as.integer(rpert(n, min=153, max=211, mode=182)) #Departure day of swan spring migration (February 1st to March 31st) #min=153, max=211, mode=182

rho_sm <-  0.99 #Proportion of migratory mallards that arrived susceptible to Croatia

day_arrival_mallards <- as.integer(rpert(n, min=30, max=62, mode=46)) #Arrival day of mallard fall migration (Oct 1 to Nov 2) Need to allow for a stopover of <= 28 days 

migdays_m <- 1 #assumed to be the same coming in and leaving.

b_m <- Nm #b_m represents the total number of mallards that immigrated into the area over migdays_m

beta_mm <- transmission_rate #HPAI transmission rate between mallards (per bird per day)

eta_m <- runif(n,min=10^3,max=10^4) #HPAI infectious dose for mallards (virion)

mu_mal <- rpert(n, min=0.00086, max=0.0039, mode=0.0018) #Natural mortality rate for mallards (per bird per day)

rho_im <- 1-rho_sm #Proportion of migratory mallards that arrive infected to Croatia (a single individual introduces the infection)

gamma_mal <- 1/rpert(n, min = 4.2, max = 7.2, mode = 5.4) #Recovery rate from HPAI infection in mallards (per bird per day)

dead_m <- rpert(n, min = 0, max = 1, mode = 0) #Probability of death from HPAI infection in mallards

d_m <- dead_m/rpert(n, min = 4.6, max = 7.6, mode = 5.7)  #Death rate from HPAI infection in mallards (per bird per day)

duration_mallard_stopover <- as.integer(runif(n, min=7, max=28)) #Length of time mallards spend at the staging site during the fall (between 7-28 days)

day_departure_mallards <- day_arrival_mallards + duration_mallard_stopover #Departure day of mallard fall migration (arrival + stay of 7 to 28 days)

#calibrated value: 
beta_mp<-  0.00000052125 #HPAI transmission rate from  mallards to backyard poultry farms (per farm per day)

d_p<- 1/rpert(n, min = 2, max = 5, mode = 3.4) #Rate of backyard poultry  farm culling (per farm per day)

#=============Model equations======================

PiecewiseFlux<-function(t, t1, t2, b) {ifelse(t1<=floor(t %% 365) && floor(t %% 365) <=t2,b,0)} #This function is used to establish the breeding season in resident mallards

AI_MultiSpecies.dyn<- function(t, var, par) {
  
  # Parameters 
  epsilon_w<- par[1];
  epsilon_m<- par[2];
  day_arrival_muteswans <- par[3];
  mu_w <- par[4];
  d_w <- par[5];
  day_arrival_mallards <- par[6];
  eta_m<- par[7];
  mu_mal <- par[8];
  gamma_mal <- par[9];
  d_m <- par[10];
  d_p <- par[11];


  #State values
  E = var[1];
  S_w = var[2];
  I_w = var[3];
  D_w = var[4];
  S_m = var[5];
  I_m = var[6];
  R_m = var[7];
  D_m = var[8];
  S_mbar = var[9];
  I_mbar = var[10];
  R_mbar = var[11];
  D_mbar = var[12];
  S_p = var[13];
  I_p = var[14];
  D_p = var[15];

  
# Ordinary differential equations

# Mute Swans
  dE =  epsilon_w*I_w + epsilon_m*(I_m+I_mbar) - E*r 
  
  dS_w = rho_sw * PiecewiseFlux(t, day_arrival_muteswans, day_arrival_muteswans+migdays_w - 1, b_w/migdays_w) - 
    ifelse(S_w + I_w == 0, 0, (beta_ww*I_w*S_w)/(S_w+I_w)) - mu_w*S_w - alpha*S_w*(E/eta_w)
  
  dI_w = rho_iw * PiecewiseFlux(t, day_arrival_muteswans, day_arrival_muteswans+migdays_w - 1, b_w/migdays_w)  + 
    ifelse(S_w + I_w == 0, 0, (beta_ww*I_w*S_w)/(S_w+I_w)) + alpha*S_w*(E/eta_w)  - (mu_w + d_w)*I_w  
  
  dD_w =  d_w*I_w
  
# Migratory mallards
  dS_m = rho_sm * PiecewiseFlux(t, day_arrival_mallards, day_arrival_mallards+migdays_m -1 , b_m/migdays_m ) - 
    ifelse(S_m + I_m + R_m == 0, 0, (beta_mm * (I_m + I_mbar)*S_m)/(S_m+I_m+R_m+S_mbar+I_mbar+R_mbar)) - alpha*S_m*(E/eta_m) - mu_mal*S_m
  
  dI_m = rho_im * PiecewiseFlux(t, day_arrival_mallards, day_arrival_mallards+migdays_m -1 , b_m/migdays_m ) + 
    ifelse(S_m + I_m + R_m == 0, 0, (beta_mm * (I_m + I_mbar)*S_m)/(S_m+I_m+R_m+S_mbar+I_mbar+R_mbar)) + alpha*S_m*(E/eta_m) - (mu_mal + gamma_mal + d_m)*I_m  
  
  dR_m = (1 - rho_sm - rho_im)*PiecewiseFlux(t, day_arrival_mallards, day_arrival_mallards+migdays_m -1, b_m/migdays_m) + gamma_mal*I_m - mu_mal*R_m 
  
  dD_m =  d_m*I_m
  
# Resident mallards
  dS_mbar = (if(t >= 242 & t <= 364) {
    2.985 * mu_mal * (S_mbar+I_mbar+R_mbar)
  } 
  else 0)  - (beta_mm*(I_m + I_mbar)*S_mbar)/(S_m+I_m+R_m+S_mbar+I_mbar+R_mbar) - alpha*S_mbar*(E/eta_m) - mu_mal*S_mbar
  
  dI_mbar = (beta_mm*(I_m + I_mbar)*S_mbar)/(S_m+I_m+R_m+S_mbar+I_mbar+R_mbar) + alpha*S_mbar*(E/eta_m) - (gamma_mal + mu_mal + d_m)*I_mbar

  dR_mbar =  gamma_mal*I_mbar - mu_mal*R_mbar 
  
  dD_mbar = d_m*I_mbar
  
# Poultry
  dS_p = - beta_mp*(I_m+I_mbar)*S_p 
  
  dI_p =  beta_mp*(I_m+I_mbar)*S_p - (d_p)*I_p
  
  dD_p = d_p*I_p


  
  return(list(c(dE, dS_w, dI_w, dD_w, dS_m, dI_m, dR_m, dD_m, dS_mbar, dI_mbar, dR_mbar, dD_mbar, dS_p, dI_p, dD_p 
                )))
}

#=====================================================================================================
# Run the model

AI_MultiSpecies.sol=NULL
E=NULL
S_w=NULL
I_w=NULL
D_w=NULL
S_m=NULL
I_m=NULL
R_m=NULL
D_m=NULL
S_mbar=NULL
I_mbar=NULL
R_mbar=NULL
D_mbar=NULL
S_p=NULL
I_p=NULL
D_p=NULL


#=====================================================================================================
# Monte Carlo simulation
#=====================================================================================================

for(i in  1:n)
{
  AI_MultiSpecies.par<- c(
    epsilon_w[i],
    epsilon_m[i],
    day_arrival_muteswans[i],
    mu_w[i],
    d_w[i],
    day_arrival_mallards[i],
    eta_m[i],
    mu_mal[i],
    gamma_mal[i],
    d_m[i],
    d_p[i]
  )

  AI_MultiSpecies.t <- seq(0,365, 1) # Change the 1*365 to add multiple years to the model
  
  eventdat <- data.frame(var = c(rep (c("S_m_init", "I_m_init", "R_m_init","S_w_init", "I_w_init"),1)),
                         time = c( rep(day_departure_mallards[i], 3), rep(day_departure_muteswans[i],2)),
                         value = c(rep(0,5)),
                         method = c(rep("rep", 5))
  ) 
  eventdat
  
  AI_MultiSpecies.sol [[i]] <- lsode(AI_MultiSpecies.init, AI_MultiSpecies.t, AI_MultiSpecies.dyn, 
                                     AI_MultiSpecies.par, events = list(data = eventdat))
}

#This tell us how long it took for the model to run
difftime(Sys.time(), initial.time, units='mins') 

colnames(AI_MultiSpecies.sol[[1]])

removed<- which((sapply(1:length(AI_MultiSpecies.sol), function(i) dim(AI_MultiSpecies.sol[[i]])[1])<366) == T) #Identify iterations where the model run into an error. 

ModelOutput<-AI_MultiSpecies.sol
for(i in 1:length(ModelOutput)){
  if(nrow(ModelOutput[[i]])<365) 
    ModelOutput[[i]]<-NULL
    
} ###Iterations with errors are removed

sapply(1:length(ModelOutput), function(i) dim(ModelOutput[[i]])[1]) 
sum(sapply(1:length(ModelOutput), function(i) dim(ModelOutput[[i]])[1])<366)

##########################################
# Evaluate each compartment for diagnosis
##########################################

colnames(ModelOutput[[1]])

times=0:365
plot<- as.data.frame(ModelOutput) 

#=====================================================================================================
# Extract Monte Carlo results
#===========================================================================
inf_w<- NULL
inf_m<- NULL
inf_mbar<- NULL
env<- NULL
deaths_w <- NULL
deaths_m<-NULL
deaths_mbar<- NULL
deaths_p<- NULL
residprevI_mbar <- NULL
overlap <- NULL
wdays <- NULL
mdays <- NULL
mbardays <- NULL
pdays <- NULL
firstdetect <- NULL
infected_farm <- NULL

for (i in 1:length(ModelOutput)) {
  inf_w[i] <- ModelOutput[[i]][365,4]
  inf_m[i] <- ModelOutput[[i]][365,7]
  inf_mbar[i] <- ModelOutput[[i]][365,11]
  env[i]<- ModelOutput[[i]][365, 2]
  deaths_w[i] <- ModelOutput[[i]][365,5]  #Cumulative annual deaths in swans
  deaths_m[i] <- ModelOutput[[i]][365,9] #Cumulative annual deaths in mig mallards
  deaths_mbar[i] <- ModelOutput[[i]][365,13] #Cumulative annual deaths in resident mallards
  deaths_p[i] <- ModelOutput[[i]][365,15]+ModelOutput[[i]][365,16] #* (4717/192) #Cumulative annual deaths in poultry #ModelOutput[[i]][365,15] 
                                                               #is multiplied by the average number of chicken in those farms to get the number of chicken deaths 
                                                               #4717/192 is for JASTREBARSKO only.
  infected_farm[i] <-ModelOutput[[i]][365,15]+ModelOutput[[i]][365,16]
  overlap[i] <- (if (day_arrival_mallards[i] >= day_arrival_muteswans[i]) {
    duration_mallard_stopover[i]
  } 
  else if(day_departure_mallards[i] > day_arrival_muteswans[i]) {
    day_departure_mallards[i]-day_arrival_muteswans[i]
  }
  else if(day_departure_mallards[i] <= day_arrival_muteswans[i]) {
    0
  }
  )
  
 }

output <- as.data.frame(cbind(overlap, deaths_w, deaths_m, deaths_mbar, infected_farm,residprevI_mbar))


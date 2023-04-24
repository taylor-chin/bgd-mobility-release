# metapopulation SEIR model
# adapted from https://science.sciencemag.org/content/368/6490/489.full
# code: https://github.com/SenPei-CU/COVID-19
# http://mathesaurus.sourceforge.net/octave-r.html


# requires pop data to run from utils.R file
spatial_seir <- function(x, M, pop, ts, pop0) {
  
  num_loc = length(pop)
  dt = 1 #timestep = 1 day
  tmstep = 1 #integrate forward for one day
  num_states = 7
  
  Sidx = seq(from =1, to = num_states*num_loc, by = num_states)
  Eidx = seq(from =2, to = num_states*num_loc, by = num_states)
  Isidx = seq(from =3, to = num_states*num_loc, by = num_states)
  Iaidx = seq(from =4, to = num_states*num_loc, by = num_states)
  obsidx = seq(from =5, to = num_states*num_loc, by = num_states)
  totobsidx = seq(from =6, to = num_states*num_loc, by = num_states)  # added
  leaveiidx = seq(from =7, to = num_states*num_loc, by = num_states)  # added
  
  betaidx = num_states*num_loc + 1
  muidx = num_states*num_loc + 2
  thetaidx = num_states*num_loc + 3
  Zidx = num_states*num_loc + 4
  alphaidx = num_states*num_loc + 5
  Didx = num_states*num_loc + 6
  
  S = matrix(0, num_loc, tmstep + 1)
  E = matrix(0, num_loc, tmstep + 1)
  Is = matrix(0, num_loc, tmstep + 1)
  Ia = matrix(0, num_loc, tmstep + 1)
  Incidence = matrix(0, num_loc, tmstep + 1)
  obs = matrix(0, num_loc, tmstep + 1)
  Totobs_count = matrix(0, num_loc, tmstep + 1) # added
  totobs = matrix(0, num_loc, tmstep + 1) # added
  Leavei_count = matrix(0, num_loc, tmstep + 1) # added
  leavei = matrix(0, num_loc, tmstep + 1) # added 
  
  #intialize
  S[,1] = x[Sidx]
  E[,1] = x[Eidx]
  Is[,1] = x[Isidx]
  Ia[,1] = x[Iaidx]
  beta = x[betaidx]   # transmission rate due to symptomatic individuals
  mu = x[muidx]       # transmission rate due to asymp. individuals is multiplied by a factor mu
  theta = x[thetaidx] # mult. factor, which is >1 to reflect underreporting of human movement
  Z = x[Zidx]         # average latency period
  alpha = x[alphaidx] # fraction of symptomatic infections
  D = x[Didx]         # average duration of infection
  
  # We assume that individuals in the Is group do not move between cities,
  # though these individuals can move between cities during the latency period.
  
  # Mij is the daily number of people traveling from city j to city i
  
  # start integration
  tcnt = 0
  for (t in seq(from = (ts+dt), to = (ts+tmstep), by = dt)) {
    
    
    tcnt=tcnt+1
    dt1=dt
    
    #first step
    ESenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (S[,tcnt]/(pop-Is[,tcnt])))
    ESleft = pmin(dt1 * (rep(1, num_loc)*theta) * (S[,tcnt]/(pop-Is[,tcnt])) * colSums(M[[ts]]), dt1*S[,tcnt])
    
    EEenter = dt1* (rep(1,num_loc)*theta) * (M[[ts]] %*% (E[,tcnt]/(pop-Is[,tcnt])))
    EEleft = pmin(dt1 * (rep(1, num_loc)*theta) * (E[,tcnt]/(pop-Is[,tcnt])) * colSums(M[[ts]]), dt1*E[,tcnt])
    
    EIaenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (Ia[,tcnt]/(pop-Is[,tcnt])))
    EIaleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Ia[,tcnt]/(pop-Is[,tcnt])) * colSums(M[[ts]]), dt1*Ia[,tcnt])  
    
    Eexps = dt1*(rep(1,num_loc)*beta) * S[,tcnt] * Is[,tcnt]/pop
    Eexpa = dt1*(rep(1,num_loc)*mu) * (rep(1,num_loc)*beta) * S[,tcnt] * Ia[,tcnt]/pop
    Einfs = dt1*(rep(1, num_loc)*alpha) * E[,tcnt]/(rep(1,num_loc)*Z)
    Einfa = dt1*(rep(1, num_loc)*(1-alpha)) * E[,tcnt]/(rep(1,num_loc)*Z)
    Erecs = dt1*Is[,tcnt]/(rep(1,num_loc)*D)
    Ereca = dt1*Ia[,tcnt]/(rep(1,num_loc)*D)
    
    ESenter = pmax(ESenter,0); ESleft = pmax(ESleft,0)
    EEenter = pmax(EEenter,0); EEleft = pmax(EEleft,0);
    EIaenter = pmax(EIaenter,0); EIaleft= pmax(EIaleft,0);
    Eexps = pmax(Eexps,0); Eexpa = pmax(Eexpa,0);
    Einfs = pmax(Einfs,0); Einfa = pmax(Einfa,0);
    Erecs = pmax(Erecs,0); Ereca = pmax(Ereca,0);
    
    #stochastic version
    ESenter = rpois(num_loc,ESenter); ESleft = rpois(num_loc,ESleft);
    EEenter = rpois(num_loc, EEenter); EEleft = rpois(num_loc, EEleft);
    EIaenter = rpois(num_loc,EIaenter);EIaleft=rpois(num_loc,EIaleft);
    Eexps=rpois(num_loc,Eexps);
    Eexpa=rpois(num_loc,Eexpa);
    Einfs=rpois(num_loc,Einfs);
    Einfa=rpois(num_loc,Einfa);
    Erecs=rpois(num_loc,Erecs);
    Ereca=rpois(num_loc,Ereca);
    
    sk1=-Eexps-Eexpa+ESenter-ESleft;
    ek1=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk1=Einfs-Erecs;
    iak1=Einfa-Ereca+EIaenter-EIaleft;
    ik1i=Einfs;
    totobs1i = Einfs+Einfa;  # added
    leavei1i = EIaleft; # added
    
    #second step
    Ts1 = S[,tcnt]+ sk1/2;
    Te1 = E[,tcnt] + ek1/2;
    Tis1 = Is[,tcnt] + isk1/2;
    Tia1 = Ia[,tcnt] + iak1/2;
    
    
    ESenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (Ts1/(pop-Tis1)))
    ESleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Ts1/(pop-Tis1)) * colSums(M[[ts]]), dt1*Ts1)
    
    EEenter = dt1* (rep(1,num_loc)*theta) * (M[[ts]] %*% (Te1/(pop-Tis1)))
    EEleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Te1/(pop-Tis1)) * colSums(M[[ts]]), dt1*Te1)
    
    EIaenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (Tia1/(pop-Tis1)))
    EIaleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Tia1/(pop-Tis1)) * colSums(M[[ts]]), dt1*Tia1)  
    
    
    Eexps = dt1*(rep(1,num_loc)*beta) * Ts1 * Tis1/pop;
    Eexpa = dt1*(rep(1,num_loc)*mu) * (rep(1,num_loc)*beta) * Ts1 * Tia1/pop;
    Einfs = dt1*(rep(1,num_loc)*alpha) * Te1/(rep(1,num_loc)*Z);
    Einfa = dt1*(rep(1,num_loc)*(1-alpha)) * Te1/(rep(1,num_loc)*Z);
    Erecs = dt1*Tis1/(rep(1,num_loc)*D);
    Ereca = dt1*Tia1/(rep(1,num_loc)*D);
    
    ESenter = pmax(ESenter,0); ESleft = pmax(ESleft,0)
    EEenter = pmax(EEenter,0); EEleft = pmax(EEleft,0);
    EIaenter = pmax(EIaenter,0); EIaleft= pmax(EIaleft,0);
    Eexps = pmax(Eexps,0); Eexpa = pmax(Eexpa,0);
    Einfs = pmax(Einfs,0); Einfa = pmax(Einfa,0);
    Erecs = pmax(Erecs,0); Ereca = pmax(Ereca,0);
    
    #stochastic version
    ESenter = rpois(num_loc,ESenter); ESleft = rpois(num_loc,ESleft);
    EEenter = rpois(num_loc, EEenter); EEleft = rpois(num_loc, EEleft);
    EIaenter = rpois(num_loc,EIaenter);EIaleft=rpois(num_loc,EIaleft);
    Eexps=rpois(num_loc,Eexps);
    Eexpa=rpois(num_loc,Eexpa);
    Einfs=rpois(num_loc,Einfs);
    Einfa=rpois(num_loc,Einfa);
    Erecs=rpois(num_loc,Erecs);
    Ereca=rpois(num_loc,Ereca);
    
    
    sk2=-Eexps-Eexpa+ESenter-ESleft;
    ek2=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk2=Einfs-Erecs;
    iak2=Einfa-Ereca+EIaenter-EIaleft;
    ik2i=Einfs;
    totobs2i = Einfs+Einfa;  # added
    leavei2i = EIaleft; # added
    
    
    # third step
    Ts2 = S[,tcnt]+ sk2/2;
    Te2 = E[,tcnt] + ek2/2;
    Tis2 = Is[,tcnt] + isk2/2;
    Tia2 = Ia[,tcnt] + iak2/2;
    
    
    ESenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (Ts2/(pop-Tis2)))
    ESleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Ts2/(pop-Tis2)) * colSums(M[[ts]]), dt1*Ts2)
    
    EEenter = dt1* (rep(1,num_loc)*theta) * (M[[ts]] %*% (Te2/(pop-Tis2)))
    EEleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Te2/(pop-Tis2)) * colSums(M[[ts]]), dt1*Te2)
    
    EIaenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (Tia2/(pop-Tis2)))
    EIaleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Tia2/(pop-Tis2)) * colSums(M[[ts]]), dt1*Tia2)  
    
    
    Eexps = dt1*(rep(1,num_loc)*beta) * Ts2 * Tis2/pop;
    Eexpa = dt1*(rep(1,num_loc)*mu) * (rep(1,num_loc)*beta) * Ts2 * Tia2/pop;
    Einfs = dt1*(rep(1,num_loc)*alpha) * Te2/(rep(1,num_loc)*Z);
    Einfa = dt1*(rep(1,num_loc)*(1-alpha)) * Te2/(rep(1,num_loc)*Z);
    Erecs = dt1*Tis2/(rep(1,num_loc)*D);
    Ereca = dt1*Tia2/(rep(1,num_loc)*D);
    
    ESenter = pmax(ESenter,0); ESleft = pmax(ESleft,0)
    EEenter = pmax(EEenter,0); EEleft = pmax(EEleft,0);
    EIaenter = pmax(EIaenter,0); EIaleft= pmax(EIaleft,0);
    Eexps = pmax(Eexps,0); Eexpa = pmax(Eexpa,0);
    Einfs = pmax(Einfs,0); Einfa = pmax(Einfa,0);
    Erecs = pmax(Erecs,0); Ereca = pmax(Ereca,0);
    
    #stochastic version
    ESenter = rpois(num_loc,ESenter); ESleft = rpois(num_loc,ESleft);
    EEenter = rpois(num_loc, EEenter); EEleft = rpois(num_loc, EEleft);
    EIaenter = rpois(num_loc,EIaenter);EIaleft=rpois(num_loc,EIaleft);
    Eexps=rpois(num_loc,Eexps);
    Eexpa=rpois(num_loc,Eexpa);
    Einfs=rpois(num_loc,Einfs);
    Einfa=rpois(num_loc,Einfa);
    Erecs=rpois(num_loc,Erecs);
    Ereca=rpois(num_loc,Ereca);
    
    sk3=-Eexps-Eexpa+ESenter-ESleft;
    ek3=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk3=Einfs-Erecs;
    iak3=Einfa-Ereca+EIaenter-EIaleft;
    ik3i=Einfs;
    totobs3i = Einfs+Einfa; # added
    leavei3i = EIaleft; # added
    
    # fourth step
    Ts3 = S[,tcnt]+sk3;
    Te3 = E[,tcnt]+ek3;
    Tis3 = Is[,tcnt]+isk3;
    Tia3 = Ia[,tcnt]+iak3;
    
    ESenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (Ts3/(pop-Tis3)))
    ESleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Ts3/(pop-Tis3)) * colSums(M[[ts]]), dt1*Ts3)
    
    EEenter = dt1* (rep(1,num_loc)*theta) * (M[[ts]] %*% (Te3/(pop-Tis3)))
    EEleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Te3/(pop-Tis3)) * colSums(M[[ts]]), dt1*Te3)
    
    EIaenter = dt1 * (rep(1,num_loc)*theta) * (M[[ts]] %*% (Tia3/(pop-Tis3)))
    EIaleft = pmin(dt1 * (rep(1, num_loc)*theta) * (Tia3/(pop-Tis3)) * colSums(M[[ts]]), dt1*Tia3)  
    
    
    Eexps = dt1*(rep(1,num_loc)*beta) * Ts3 * Tis3/pop;
    Eexpa = dt1*(rep(1,num_loc)*mu) * (rep(1,num_loc)*beta) * Ts3 * Tia3/pop;
    Einfs = dt1*(rep(1,num_loc)*alpha) * Te3/(rep(1,num_loc)*Z);
    Einfa = dt1*(rep(1,num_loc)*(1-alpha)) * Te3/(rep(1,num_loc)*Z);
    Erecs = dt1*Tis3/(rep(1,num_loc)*D);
    Ereca = dt1*Tia3/(rep(1,num_loc)*D);
    
    ESenter = pmax(ESenter,0); ESleft = pmax(ESleft,0)
    EEenter = pmax(EEenter,0); EEleft = pmax(EEleft,0);
    EIaenter = pmax(EIaenter,0); EIaleft= pmax(EIaleft,0);
    Eexps = pmax(Eexps,0); Eexpa = pmax(Eexpa,0);
    Einfs = pmax(Einfs,0); Einfa = pmax(Einfa,0);
    Erecs = pmax(Erecs,0); Ereca = pmax(Ereca,0);
    
    #stochastic version
    ESenter = rpois(num_loc,ESenter); ESleft = rpois(num_loc,ESleft);
    EEenter = rpois(num_loc, EEenter); EEleft = rpois(num_loc, EEleft);
    EIaenter = rpois(num_loc,EIaenter);EIaleft=rpois(num_loc,EIaleft);
    Eexps=rpois(num_loc,Eexps);
    Eexpa=rpois(num_loc,Eexpa);
    Einfs=rpois(num_loc,Einfs);
    Einfa=rpois(num_loc,Einfa);
    Erecs=rpois(num_loc,Erecs);
    Ereca=rpois(num_loc,Ereca);
    
    
    sk4=-Eexps-Eexpa+ESenter-ESleft;
    ek4=Eexps+Eexpa-Einfs-Einfa+EEenter-EEleft;
    isk4=Einfs-Erecs;
    iak4=Einfa-Ereca+EIaenter-EIaleft;
    ik4i=Einfs;
    totobs4i = Einfs+Einfa; # added
    leavei4i = EIaleft; # added
    
    S[,tcnt+1] = S[,tcnt]+round(sk1/6+sk2/3+sk3/3+sk4/6);
    E[,tcnt+1] = E[,tcnt]+round(ek1/6+ek2/3+ek3/3+ek4/6);
    Is[,tcnt+1] = Is[,tcnt] + round(isk1/6+isk2/3+isk3/3+isk4/6);
    Ia[,tcnt+1] = Ia[,tcnt]+round(iak1/6+iak2/3+iak3/3+iak4/6);
    Incidence[,tcnt+1]=round(ik1i/6+ik2i/3+ik3i/3+ik4i/6);
    obs=Incidence[,tcnt+1];
    Totobs_count[,tcnt+1]=round(totobs1i/6+totobs2i/3+totobs3i/3+totobs4i/6); # added
    totobs=Totobs_count[,tcnt+1]; # added
    Leavei_count[,tcnt+1]=round(leavei1i/6+leavei2i/3+leavei3i/3+leavei4i/6); # added
    leavei=Leavei_count[,tcnt+1]; # added
    
    
  }
  
  # update x
  x[Sidx] = S[,tcnt+1];
  x[Eidx] = E[,tcnt+1];
  x[Isidx] = Is[,tcnt+1];
  x[Iaidx] = Ia[,tcnt+1];
  x[obsidx]=obs;
  x[totobsidx]=totobs;  # added
  x[leaveiidx] = leavei;# added
  #update pop
  #print(paste0("First pop output:", pop[1]))
  pop = pop - colSums(M[[ts]])*theta + rowSums(M[[ts]])*theta;
  #print(paste0("Pop output after travel:", pop[1]))
  #print(paste0("Column sum/ppl leaving:",colSums(M[[ts]])[1]))
  #print(paste0("Row sum/ppl entering:",rowSums(M[[ts]])[1]))
  minfrac=0.6;
  pop[which(pop<minfrac*pop0)] = pop0[which(pop<minfrac*pop0)]*minfrac
  #print(paste0("Final pop output:", pop[1]))
  
  return(list(x=x, pop =pop))
}


run_metapop = function(m1.param, # weekday mobility matrix beta distrib parameters
                       m2.param, # weekend mobility matrix beta distrib parameters
                       mobility_source, # name used in df
                       seed_city, # city where outbreak is initialized
                       pop_data, # upazila level population data
                       r0){  # calculate beta based on input r0
  
  # 
  pop_orig_data = pop_data %>% 
    left_join(., district_upa %>% select(District, upazila_code), by = "upazila_code") %>% 
    group_by(District) %>%
    summarise(pop = sum(pop))
  
  #
  num_loc = length(unique(pop_orig_data$District)) # number of locations

  # indices for each compartment + fixed parameters
  num_states = 7
  
  Sidx = seq(from =1, to = num_states*num_loc, by = num_states)
  Eidx = seq(from =2, to = num_states*num_loc, by = num_states)
  Isidx = seq(from =3, to = num_states*num_loc, by = num_states)
  Iaidx = seq(from =4, to = num_states*num_loc, by = num_states)
  obsidx = seq(from =5, to = num_states*num_loc, by = num_states)
  totobsidx = seq(from =6, to = num_states*num_loc, by = num_states)
  leaveiidx = seq(from =7, to = num_states*num_loc, by = num_states)
  betaidx = num_states*num_loc + 1
  muidx = num_states*num_loc + 2
  thetaidx = num_states*num_loc + 3
  Zidx = num_states*num_loc + 4
  alphaidx = num_states*num_loc + 5
  Didx = num_states*num_loc + 6
  
  n_runs = 1000
  x.save.list <- list()
  pop.save.list <- list()
  
  for (run in 1:n_runs){
    
    # draw from beta distribution using saved alpha, beta parameters
    if(mobility_source == "Meta"){ # district level
      draw.m1 = matrix(NA, nrow = 64, ncol = 64) # 64 districts
      draw.m2 = matrix(NA, nrow = 64, ncol = 64) 
      for(i in 1:64){
        for(j in 1:64){
          draw.m1[i,j] = rbeta(1, m1.param$alpha[i,j], m1.param$beta[i,j])
          draw.m2[i,j] = rbeta(1, m2.param$alpha[i,j], m2.param$beta[i,j])
        }
      }
      
      colnames(draw.m1) = colnames(m1.param$alpha)
      rownames(draw.m1) = rownames(m1.param$alpha)
      colnames(draw.m2) = colnames(m2.param$alpha)
      rownames(draw.m2) = rownames(m2.param$alpha)
      
      # normalize, symmetrize, differentiate between 0 and NA
      # this function is from mat_utils.R
      sym.m1 = get_norm_sym_m(draw.m1, pop_input = pop_data, mobility_source = mobility_source)
      sym.m2 = get_norm_sym_m(draw.m2, pop_input = pop_data, mobility_source = mobility_source)
      
      district.m1 = sym.m1
      district.m2 = sym.m2
      
    } else if(mobility_source != "Meta"){
      draw.m1 = matrix(NA, nrow = 544, ncol = 544) # 544 upazilas
      draw.m2 = matrix(NA, nrow = 544, ncol = 544) # 544 upazilas
      for(i in 1:544){
        for(j in 1:544){
          draw.m1[i,j] = rbeta(1, m1.param$alpha[i,j], m1.param$beta[i,j])
          draw.m2[i,j] = rbeta(1, m2.param$alpha[i,j], m2.param$beta[i,j])
        }
      }
      
      colnames(draw.m1) = colnames(m1.param$alpha)
      rownames(draw.m1) = rownames(m1.param$alpha)
      colnames(draw.m2) = colnames(m2.param$alpha)
      rownames(draw.m2) = rownames(m2.param$alpha)
      
      # normalize, symmetrize, differentiate between 0 and NA
      # this function is from mat_utils.R
      sym.m1 = get_norm_sym_m(draw.m1, pop_input = pop_data, mobility_source = mobility_source)
      sym.m2 = get_norm_sym_m(draw.m2, pop_input = pop_data, mobility_source = mobility_source)
      
      # aggregate to districts
      # this function is from mat_utils.R
      district.m1 = get_district_m(sym.m1)
      district.m2 = get_district_m(sym.m2)
    }
    
    district.m1.forpop = district.m1
    
    # zero out diagonal
    diag(district.m1) <- 0
    diag(district.m2) <- 0
    
    m_week <- c(replicate(5, district.m1, simplify = FALSE), replicate(2, district.m2, simplify = FALSE))
    M <- rep(m_week, times = 75)
    num_tsteps = length(M)
    
    # district:
    pop0 = unname(colSums(district.m1.forpop))
    pop = unname(colSums(district.m1.forpop))

    # this only changes vector for op1
    pop0[pop0==0] <- 673503.8 # Jhalokati's population - only missing in op1 data
    pop[pop==0] <- 673503.8 # Jhalokati's population - only missing in op1 data

    #pop0 = pop_orig_data %>% distinct(District, .keep_all = TRUE) %>% arrange(match(District, colnames(district.m1))) %>% pull(pop)
    #pop = pop_orig_data %>% distinct(District, .keep_all = TRUE) %>% arrange(match(District, colnames(district.m1))) %>% pull(pop)
    
    x = rep(NA, length =num_states*num_loc + 6) #vector of initial conditions
    ts = 1 # the current time step for starting the integration
    
    mu = 0.50      # multiplicative factor for the transmission rate of asymptomatic individuals
    alpha = 0.65   # fraction of symptomatic infections
    D = 5  
    
    seedid = match(seed_city, colnames(district.m1))
    
    x[Eidx] = rep(0, num_loc)
    x[((seedid - 1)*num_states)+2] = 500 # seed 500 infected in latent compartment, in metapopulations defined by seedid
    x[Sidx] = pop0 -  x[Eidx]   # everyone except init. infected start as susceptible in the beginning 
    x[Isidx] = rep(0, num_loc)
    x[Iaidx] = rep(0,num_loc)
    x[obsidx] = rep(0,num_loc)
    x[totobsidx] = rep(0,num_loc)
    x[leaveiidx] = rep(0,num_loc)
    x[betaidx] = r0/(alpha*D + (1 - alpha)*mu*D)    # transmission rate due to symptomatic individuals
    x[muidx] = 0.50      #  multiplicative factor for the transmission rate of asymptomatic individuals
    x[thetaidx] = 1     # mult. factor, which is >1 to reflect underreporting of human movement
    x[Zidx] = 3        # average latency period
    x[alphaidx] = 0.65   # fraction of symptomatic infections
    x[Didx] = 5          # average duration of infection
    
    x.save <- matrix(NA, nrow = length(x), ncol = num_tsteps+1)
    pop.save <- matrix(NA, nrow = length(pop), ncol = num_tsteps+1)
    x.save[,1] <- x
    pop.save[,1] <- pop0
    
    for (ts in 1:num_tsteps) {
      sim.tmp <- spatial_seir(x, M, pop, ts, pop0)
      # update values
      x  = x.save[,ts+1] = sim.tmp$x 
      pop = pop.save[,ts+1] = sim.tmp$pop 
      
    }
    
    # save each run 
    
    x.save.list[[run]] = x.save
    pop.save.list[[run]] = pop.save
    
  } # end of runs for loop
  
  # process results ----
  geo_codes = colnames(district.m1)
  
  pop.update = do.call(rbind, pop.save.list) %>%
    as.data.frame() %>%
    mutate(sim = rep(1:n_runs, each = nrow(pop.save.list[[1]]))) %>%
    mutate(geo_code = rep(geo_codes, times = n_runs)) %>%   
    pivot_longer(!c(geo_code, sim), names_to = "time", values_to = "pop") %>%
    mutate(time = as.numeric(sub('.', '', time))) # nrow = 493406
  
  #x.save.list 
  # No cols = time steps
  # Rows = 5 state outputs from spatial_seir function * # locations + 6 parameter values
  
  results = do.call(rbind, x.save.list) %>%
    as.data.frame() %>% 
    mutate(sim = rep(1:n_runs, each = nrow(x.save.list[[1]]))) %>%
    mutate(state = rep(c(rep(c("S", "E", "Is", "Ia", "obs", "totobs", "leavei"), times = num_loc), 
                         "beta", "mu", "theta", "Z", "alpha", "D"), times = n_runs)) %>%
    pivot_longer(!c(sim, state), names_to = "time", values_to = "value") %>%
    mutate(time = as.numeric(sub('.', '', time))) %>% 
    filter(!state %in% c("beta", "mu", "theta", "Z", "alpha", "D")) %>% 
    mutate(geo_code = rep(rep(geo_codes, each = (num_tsteps+1)*num_states), times = n_runs)) %>% 
    left_join(., pop.update, by = c("geo_code", "time", "sim")) %>%
    # parameters:
    mutate(source = mobility_source) %>%
    mutate(r0 = r0) %>%
    mutate(seed_city = seed_city)
    
  # average df - average across sims
  results.avg = results %>%
    group_by(state, time, geo_code, r0, seed_city, source) %>%
    summarise(mean = mean(value, na.rm=TRUE),
              lower = quantile(value, probs = 0.025, na.rm=TRUE),
              upper = quantile(value, probs = 0.975, na.rm=TRUE),
              lower0.1 = quantile(value, probs = 0.1, na.rm=TRUE),
              upper0.9 = quantile(value, probs = 0.9, na.rm=TRUE),
              mean.pop = mean(pop,na.rm= TRUE), 
              lower.pop = quantile(pop, probs = 0.025, na.rm=TRUE),
              upper.pop = quantile(pop, probs = 0.975, na.rm=TRUE)) %>%
    ungroup() 
  
  ## final size and # districts with at least one case metric for barplot
  pop_bang = sum(pop_orig_data$pop)
  
  sum_bydist = results  %>%
    filter(state == "obs") %>%
    group_by(time, source, r0, seed_city, sim) %>%
    # aggregate across geo_code districts
    summarise(value = sum(value))
  
  final_size = sum_bydist %>%
    group_by(source, r0, seed_city, sim) %>%
    arrange(time) %>%
    # cumulative sum
    mutate(final_value = cumsum(value)) %>% 
    filter(time == 526) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(final_size_prop = mean(final_value, na.rm= TRUE)/pop_bang, 
              final_size_prop_lower = quantile(final_value, probs = 0.025, na.rm=TRUE)/pop_bang,
              final_size_prop_upper = quantile(final_value, probs = 0.975, na.rm=TRUE)/pop_bang) %>%
    ungroup()
  
  num_dist_plt2 = results %>%
    filter(state == "obs") %>%
    filter(time <= 30 & value >= 1) %>% 
    group_by(source, r0, seed_city, sim) %>% #
    summarise(num_dist = n_distinct(geo_code)) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(num_dist = mean(num_dist, na.rm= TRUE), 
              num_dist_lower = quantile(num_dist, probs = 0.025, na.rm=TRUE),
              num_dist_upper = quantile(num_dist, probs = 0.975, na.rm=TRUE))

  
  return(list(results, results.avg, sum_bydist, final_size, num_dist_plt2))
}


run_metapop_upa = function(m1.param, # weekday mobility matrix beta distrib parameters
                       m2.param, # weekend mobility matrix beta distrib parameters
                       mobility_source, # name used in df
                       seed_city, # city where outbreak is initialized
                       pop_data, # upazila level population data
                       r0){  # calculate beta based on input r0
  
  # 
  pop_orig_data = pop_data
  
  #
  num_loc = length(unique(pop_orig_data$upazila_code)) # number of locations
  
  # indices for each compartment + fixed parameters
  num_states = 7
  
  Sidx = seq(from =1, to = num_states*num_loc, by = num_states)
  Eidx = seq(from =2, to = num_states*num_loc, by = num_states)
  Isidx = seq(from =3, to = num_states*num_loc, by = num_states)
  Iaidx = seq(from =4, to = num_states*num_loc, by = num_states)
  obsidx = seq(from =5, to = num_states*num_loc, by = num_states)
  moveidx = seq(from =6, to = num_states*num_loc, by = num_states)
  leaveiidx = seq(from =7, to = num_states*num_loc, by = num_states)
  betaidx = num_states*num_loc + 1
  muidx = num_states*num_loc + 2
  thetaidx = num_states*num_loc + 3
  Zidx = num_states*num_loc + 4
  alphaidx = num_states*num_loc + 5
  Didx = num_states*num_loc + 6
  
  n_runs = 100
  
  x.save.list <- list()
  pop.save.list <- list()
  
  for (run in 1:n_runs){
    
    # draw from beta distribution using saved alpha, beta parameters
    draw.m1 = matrix(NA, nrow = 544, ncol = 544) # 544 upazilas
    draw.m2 = matrix(NA, nrow = 544, ncol = 544) # 544 upazilas
    for(i in 1:544){
      for(j in 1:544){
        draw.m1[i,j] = rbeta(1, m1.param$alpha[i,j], m1.param$beta[i,j])
        draw.m2[i,j] = rbeta(1, m2.param$alpha[i,j], m2.param$beta[i,j])
      }
    }
    
    colnames(draw.m1) = colnames(m1.param$alpha)
    rownames(draw.m1) = rownames(m1.param$alpha)
    colnames(draw.m2) = colnames(m2.param$alpha)
    rownames(draw.m2) = rownames(m2.param$alpha)
    
    # normalize, symmetrize, differentiate between 0 and NA
    # this function is from mat_utils.R
    sym.m1 = get_norm_sym_m(draw.m1, pop_input = pop_data, mobility_source = mobility_source)
    sym.m2 = get_norm_sym_m(draw.m2, pop_input = pop_data, mobility_source = mobility_source)
    
    sym.m1[is.na(sym.m1)] <- 0
    sym.m2[is.na(sym.m2)] <- 0
    
    sym.m1.forpop = sym.m1
    
    # zero out diagonal
    diag(sym.m1) <- 0
    diag(sym.m2) <- 0
    
    m_week <- c(replicate(5, sym.m1, simplify = FALSE), replicate(2, sym.m2, simplify = FALSE))
    M <- rep(m_week, times = 75)
    num_tsteps = length(M)
    
    upapop_tofill = pop_orig_data %>% distinct(upazila_code, .keep_all = TRUE) %>% arrange(match(upazila_code, colnames(sym.m1))) %>% pull(pop)
    
    matpop = unname(colSums(sym.m1.forpop))
    pop0 = ifelse(matpop == 0, upapop_tofill, matpop)
    pop = ifelse(matpop == 0, upapop_tofill, matpop)
    
    #pop0 = pop_orig_data %>% distinct(upazila_code, .keep_all = TRUE) %>% arrange(match(upazila_code, colnames(sym.m1))) %>% pull(pop)
    #pop = pop_orig_data %>% distinct(upazila_code, .keep_all = TRUE) %>% arrange(match(upazila_code, colnames(sym.m1))) %>% pull(pop)
    
    x = rep(NA, length =num_states*num_loc + 6) #vector of initial conditions
    ts = 1 # the current time step for starting the integration
    
    mu = 0.50      #multiplicative factor for the transmission rate of asymptomatic individuals
    alpha = 0.65   # fraction of symptomatic infections
    D = 5  
    
    upa_in_seed_district = district_upa %>% filter(District == seed_city) %>% pull(upazila_code)
    top_upa = pop_orig_data %>% filter(upazila_code %in% upa_in_seed_district) %>% arrange(desc(pop)) %>% slice(1:5) %>% pull(upazila_code)
    seedid = match(top_upa, colnames(sym.m1))
    
    x[Eidx] = rep(0, num_loc)
    x[((seedid - 1)*num_states)+2] = 100 # seed 100 infected in latent compartment, in metapopulations defined by seedid
    x[Sidx] = pop0 -  x[Eidx]   # everyone except init. infected start as susceptible in the beginning 
    x[Isidx] = rep(0, num_loc)
    x[Iaidx] = rep(0,num_loc)
    x[obsidx] = rep(0,num_loc)
    x[moveidx] = rep(0,num_loc)
    x[leaveiidx] = rep(0,num_loc)
    x[betaidx] = r0/(alpha*D + (1 - alpha)*mu*D)  # transmission rate due to symptomatic individuals
    x[muidx] = 0.50      #  multiplicative factor for the transmission rate of asymptomatic individuals
    x[thetaidx] = 1     # mult. factor, which is >1 to reflect underreporting of human movement
    x[Zidx] = 3        # average latency period
    x[alphaidx] = 0.65   # fraction of symptomatic infections
    x[Didx] = 5          # average duration of infection
    
    x.save <- matrix(NA, nrow = length(x), ncol = num_tsteps+1)
    pop.save <- matrix(NA, nrow = length(pop), ncol = num_tsteps+1)
    x.save[,1] <- x
    pop.save[,1] <- pop0
    
    for (ts in 1:num_tsteps) {
      sim.tmp <- spatial_seir(x, M, pop, ts, pop0)
      # update values
      x  = x.save[,ts+1] = sim.tmp$x 
      pop = pop.save[,ts+1] = sim.tmp$pop 
      
    }
    
    # save each run 
    
    x.save.list[[run]] = x.save
    pop.save.list[[run]] = pop.save
    
  } # end of runs for loop
  
  # process results ----
  geo_codes = colnames(sym.m1)
  
  pop.update = do.call(rbind, pop.save.list) %>%
    as.data.frame() %>%
    mutate(sim = rep(1:n_runs, each = nrow(pop.save.list[[1]]))) %>%
    mutate(geo_code = rep(geo_codes, times = n_runs)) %>%   
    pivot_longer(!c(geo_code, sim), names_to = "time", values_to = "pop") %>%
    mutate(time = as.numeric(sub('.', '', time))) # nrow = 493406
  
  #x.save.list 
  # No cols = time steps
  # Rows = 5 state outputs from spatial_seir.R function * # locations + 6 parameter values
  
  results = do.call(rbind, x.save.list) %>%
    as.data.frame() %>% 
    mutate(sim = rep(1:n_runs, each = nrow(x.save.list[[1]]))) %>%
    mutate(state = rep(c(rep(c("S", "E", "Is", "Ia", "obs", "move", "leavei"), times = num_loc), 
                         "beta", "mu", "theta", "Z", "alpha", "D"), times = n_runs)) %>%
    pivot_longer(!c(sim, state), names_to = "time", values_to = "value") %>%
    mutate(time = as.numeric(sub('.', '', time))) %>% 
    filter(!state %in% c("beta", "mu", "theta", "Z", "alpha", "D")) %>% 
    mutate(geo_code = rep(rep(geo_codes, each = (num_tsteps+1)*num_states), times = n_runs)) %>% 
    left_join(., pop.update, by = c("geo_code", "time", "sim")) %>%
    # parameters:
    mutate(source = mobility_source) %>%
    mutate(r0 = r0) %>%
    mutate(seed_city = seed_city)
  
  # average df - average across sims
  results.avg = results %>%
    group_by(state, time, geo_code, r0, seed_city, source) %>%
    summarise(mean = mean(value, na.rm=TRUE),
              lower = quantile(value, probs = 0.025, na.rm=TRUE),
              upper = quantile(value, probs = 0.975, na.rm=TRUE),
              lower0.1 = quantile(value, probs = 0.1, na.rm=TRUE),
              upper0.9 = quantile(value, probs = 0.9, na.rm=TRUE),
              mean.pop = mean(pop,na.rm= TRUE)) %>% 
    ungroup()
  
  ## final size and # districts with at least one case metric for barplot
  pop_bang = sum(pop_orig_data$pop)
  
  sum_bydist = results  %>%
    filter(state == "obs") %>%
    group_by(time, source, r0, seed_city, sim) %>%
    # aggregate across geo_code districts
    summarise(value = sum(value))
  
  final_size = sum_bydist %>%
    group_by(source, r0, seed_city, sim) %>%
    arrange(time) %>%
    # cumulative sum
    mutate(final_value = cumsum(value)) %>% 
    filter(time == 526) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(final_size_prop = mean(final_value, na.rm= TRUE)/pop_bang, 
              final_size_prop_lower = quantile(final_value, probs = 0.025, na.rm=TRUE)/pop_bang,
              final_size_prop_upper = quantile(final_value, probs = 0.975, na.rm=TRUE)/pop_bang) %>%
    ungroup()
  
  num_dist_plt = results %>%
    filter(state == "obs") %>%
    filter(time <= 30 & value >= 1) %>% 
    group_by(source, r0, seed_city, sim) %>% #
    summarise(num_dist = n_distinct(geo_code)) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(num_dist = mean(num_dist, na.rm= TRUE), 
              num_dist_lower = quantile(num_dist, probs = 0.025, na.rm=TRUE),
              num_dist_upper = quantile(num_dist, probs = 0.975, na.rm=TRUE))
  
  return(list(results.avg, sum_bydist, final_size, num_dist_plt))
}

# need no uncertainty versions for gravity model
# difference is does not use beta distribution parameters as inputs; differences in matrix processing
run_metapop_nomatuncert = function(m1, # weekday mobility matrix 
                                    m2, # weekend mobility matrix 
                                    mobility_source, # name used in df
                                    seed_city, # city where outbreak is initialized
                                    pop_data, # upazila-level population data
                                    r0){  # calculate beta based on input r0
  # 
  pop_orig_data = pop_data %>% 
    left_join(., district_upa %>% select(District, upazila_code), by = "upazila_code") %>% 
    group_by(District) %>%
    summarise(pop = sum(pop))
  
  #
  num_loc = length(unique(pop_orig_data$District)) # number of locations
  
  # indices for each compartment + fixed parameters
  num_states = 7
  
  Sidx = seq(from =1, to = num_states*num_loc, by = num_states)
  Eidx = seq(from =2, to = num_states*num_loc, by = num_states)
  Isidx = seq(from =3, to = num_states*num_loc, by = num_states)
  Iaidx = seq(from =4, to = num_states*num_loc, by = num_states)
  obsidx = seq(from =5, to = num_states*num_loc, by = num_states)
  totobsidx = seq(from =6, to = num_states*num_loc, by = num_states)
  leaveiidx = seq(from =7, to = num_states*num_loc, by = num_states)
  betaidx = num_states*num_loc + 1
  muidx = num_states*num_loc + 2
  thetaidx = num_states*num_loc + 3
  Zidx = num_states*num_loc + 4
  alphaidx = num_states*num_loc + 5
  Didx = num_states*num_loc + 6
  
  n_runs = 1000
  x.save.list <- list()
  pop.save.list <- list()
  
  for (run in 1:n_runs){
    
    district.m1 <- m1
    district.m2 <- m2
    
    # zero out diagonal
    diag(district.m1) <- 0
    diag(district.m2) <- 0
    
    m_week <- c(replicate(5, district.m1, simplify = FALSE), replicate(2, district.m2, simplify = FALSE))
    M <- rep(m_week, times = 75)
    length(M)
    num_tsteps = length(M)
    
    # this only changes vector for op1
    if(mobility_source == "Gravity model"){
      pop0 = pop_orig_data %>% distinct(District, .keep_all = TRUE) %>% arrange(match(District, colnames(district.m1))) %>% pull(pop)
      pop = pop_orig_data %>% distinct(District, .keep_all = TRUE) %>% arrange(match(District, colnames(district.m1))) %>% pull(pop)
    } else if(mobility_source != "Gravity model") {
      # district:
      pop0 = unname(colSums(m1))
      pop = unname(colSums(m1))
      
      pop0[pop0==0] <- 673503.8 # Jhalokati's population - only missing in op1 data
      pop[pop==0] <- 673503.8 # Jhalokati's population - only missing in op1 data
    }
    
    x = rep(NA, length =num_states*num_loc + 6) #vector of initial conditions
    ts = 1 # the current time step for starting the integration
    
    mu = 0.50      # multiplicative factor for the transmission rate of asymptomatic individuals
    alpha = 0.65   # fraction of symptomatic infections
    D = 5  
    
    seedid = match(seed_city, colnames(district.m1))
    
    x[Eidx] = rep(0, num_loc)
    x[((seedid - 1)*num_states)+2] = 500 # seed in latent compartment, in metapopulations defined by seedid
    x[Sidx] = pop0 -  x[Eidx]   # everyone except init. infected start as susceptible in the beginning 
    x[Isidx] = rep(0, num_loc)
    x[Iaidx] = rep(0,num_loc)
    x[obsidx] = rep(0,num_loc)
    x[totobsidx] = rep(0,num_loc)
    x[leaveiidx] = rep(0,num_loc)
    x[betaidx] = r0/(alpha*D + (1 - alpha)*mu*D)  # transmission rate due to symptomatic individuals
    x[muidx] = 0.50      # multiplicative factor for the transmission rate of asymptomatic individuals
    x[thetaidx] = 1     # mult. factor, which is >1 to reflect underreporting of human movement
    x[Zidx] = 3        # average latency period
    x[alphaidx] = 0.65   # fraction of symptomatic infections
    x[Didx] = 5          # average duration of infection
    
    x.save <- matrix(NA, nrow = length(x), ncol = num_tsteps+1)
    pop.save <- matrix(NA, nrow = length(pop), ncol = num_tsteps+1)
    x.save[,1] <- x
    pop.save[,1] <- pop0
    
    for (ts in 1:num_tsteps) {
      sim.tmp <- spatial_seir(x, M, pop, ts, pop0)
      # update values
      x  = x.save[,ts+1] = sim.tmp$x 
      pop = pop.save[,ts+1] = sim.tmp$pop 
      
    }
    
    # save each run 
    
    x.save.list[[run]] = x.save
    pop.save.list[[run]] = pop.save
    
  } # end of runs for loop
  
  # process results ----
  geo_codes = colnames(district.m1)
  
  pop.update = do.call(rbind, pop.save.list) %>%
    as.data.frame() %>%
    mutate(sim = rep(1:n_runs, each = nrow(pop.save.list[[1]]))) %>%
    mutate(geo_code = rep(geo_codes, times = n_runs)) %>%   
    pivot_longer(!c(geo_code, sim), names_to = "time", values_to = "pop") %>%
    mutate(time = as.numeric(sub('.', '', time))) # nrow = 493406
  
  #x.save.list 
  # No cols = time steps
  # Rows = 5 state outputs from spatial_seir function * # locations + 6 parameter values
  
  results = do.call(rbind, x.save.list) %>%
    as.data.frame() %>% 
    mutate(sim = rep(1:n_runs, each = nrow(x.save.list[[1]]))) %>%
    mutate(state = rep(c(rep(c("S", "E", "Is", "Ia", "obs", "totobs", "leavei"), times = num_loc), 
                         "beta", "mu", "theta", "Z", "alpha", "D"), times = n_runs)) %>%
    pivot_longer(!c(sim, state), names_to = "time", values_to = "value") %>%
    mutate(time = as.numeric(sub('.', '', time))) %>% 
    filter(!state %in% c("beta", "mu", "theta", "Z", "alpha", "D")) %>% 
    mutate(geo_code = rep(rep(geo_codes, each = (num_tsteps+1)*num_states), times = n_runs)) %>% 
    left_join(., pop.update, by = c("geo_code", "time", "sim")) %>%
    # parameters:
    mutate(source = mobility_source) %>%
    mutate(r0 = r0) %>%
    mutate(seed_city = seed_city)
  
  # average df - average across sims
  results.avg = results %>%
    group_by(state, time, geo_code) %>%
    summarise(mean = mean(value, na.rm=TRUE),
              lower = quantile(value, probs = 0.025, na.rm=TRUE),
              upper = quantile(value, probs = 0.975, na.rm=TRUE),
              lower0.1 = quantile(value, probs = 0.1, na.rm=TRUE),
              upper0.9 = quantile(value, probs = 0.9, na.rm=TRUE),
              mean.pop = mean(pop,na.rm= TRUE), 
              lower.pop = quantile(pop, probs = 0.025, na.rm=TRUE),
              upper.pop = quantile(pop, probs = 0.975, na.rm=TRUE)) %>%
    ungroup() %>% # merge with district for mapping %>%
    # parameters
    mutate(source = mobility_source) %>%
    mutate(r0 = r0) %>%
    mutate(seed_city = seed_city)
  
  
  ## final size and # districts with at least one case metric for barplot
  pop_bang = sum(pop_orig_data$pop)
  
  sum_bydist = results  %>%
    filter(state == "obs") %>%
    group_by(time, source, r0, seed_city, sim) %>%
    # aggregate across geo_code districts
    summarise(value = sum(value))
  
  final_size = sum_bydist %>%
    group_by(source, r0, seed_city, sim) %>%
    arrange(time) %>%
    # cumulative sum
    mutate(final_value = cumsum(value)) %>% 
    filter(time == 526) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(final_size_prop = mean(final_value, na.rm= TRUE)/pop_bang, 
              final_size_prop_lower = quantile(final_value, probs = 0.025, na.rm=TRUE)/pop_bang,
              final_size_prop_upper = quantile(final_value, probs = 0.975, na.rm=TRUE)/pop_bang) %>%
    ungroup()
  
  num_dist_plt2 = results %>%
    filter(state == "obs") %>%
    filter(time <= 30 & value >= 1) %>% 
    group_by(source, r0, seed_city, sim) %>% #
    summarise(num_dist = n_distinct(geo_code)) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(num_dist = mean(num_dist, na.rm= TRUE), 
              num_dist_lower = quantile(num_dist, probs = 0.025, na.rm=TRUE),
              num_dist_upper = quantile(num_dist, probs = 0.975, na.rm=TRUE))
  
  return(list(results, results.avg, sum_bydist, final_size, num_dist_plt2))
}

run_metapop_upa_nomatuncert = function(m1, # weekday mobility matrix beta distrib parameters
                                       m2, # weekend mobility matrix beta distrib parameters
                                       mobility_source, # name used in df
                                       seed_city, # city with initialized infecteds
                                       pop_data, # upazila level population data
                                       r0){  # calculate beta based on input r0
  
  # 
  pop_orig_data = pop_data
  
  #
  num_loc = length(unique(pop_orig_data$upazila_code)) # number of locations
  
  # indices for each compartment + fixed parameters
  num_states = 7
  
  Sidx = seq(from =1, to = num_states*num_loc, by = num_states)
  Eidx = seq(from =2, to = num_states*num_loc, by = num_states)
  Isidx = seq(from =3, to = num_states*num_loc, by = num_states)
  Iaidx = seq(from =4, to = num_states*num_loc, by = num_states)
  obsidx = seq(from =5, to = num_states*num_loc, by = num_states)
  moveidx = seq(from =6, to = num_states*num_loc, by = num_states)
  leaveiidx = seq(from =7, to = num_states*num_loc, by = num_states)
  betaidx = num_states*num_loc + 1
  muidx = num_states*num_loc + 2
  thetaidx = num_states*num_loc + 3
  Zidx = num_states*num_loc + 4
  alphaidx = num_states*num_loc + 5
  Didx = num_states*num_loc + 6
  
  n_runs = 100
  x.save.list <- list()
  pop.save.list <- list()
  
  for (run in 1:n_runs){
    
    sym.m1 <- m1
    sym.m2 <- m2
    
    sym.m1[is.na(sym.m1)] <- 0
    sym.m2[is.na(sym.m2)] <- 0
    
    sym.m1.forpop = sym.m1
    
    # zero out diagonal
    diag(sym.m1) <- 0
    diag(sym.m2) <- 0
    
    m_week <- c(replicate(5, sym.m1, simplify = FALSE), replicate(2, sym.m2, simplify = FALSE))
    M <- rep(m_week, times = 75)
    num_tsteps = length(M)
    
    upapop_tofill = pop_orig_data %>% distinct(upazila_code, .keep_all = TRUE) %>% arrange(match(upazila_code, colnames(sym.m1))) %>% pull(pop)
    
    matpop = unname(colSums(sym.m1.forpop))
    
    if(mobility_source == "Gravity model"){
      pop0 = pop_orig_data %>% distinct(upazila_code, .keep_all = TRUE) %>% arrange(match(upazila_code, colnames(sym.m1))) %>% pull(pop)
      pop = pop_orig_data %>% distinct(upazila_code, .keep_all = TRUE) %>% arrange(match(upazila_code, colnames(sym.m1))) %>% pull(pop)
    } else if(mobility_source != "Gravity model") {
      pop0 = ifelse(matpop == 0, upapop_tofill, matpop)
      pop = ifelse(matpop == 0, upapop_tofill, matpop)
    }
    
    x = rep(NA, length =num_states*num_loc + 6) #vector of initial conditions
    ts = 1 # the current time step for starting the integration
    
    mu = 0.50      # multiplicative factor for the transmission rate of asymptomatic individuals
    alpha = 0.65   # fraction of symptomatic infections
    D = 5  
    
    upa_in_seed_district = district_upa %>% filter(District == seed_city) %>% pull(upazila_code)
    top_upa = pop_orig_data %>% filter(upazila_code %in% upa_in_seed_district) %>% arrange(desc(pop)) %>% slice(1:5) %>% pull(upazila_code)
    seedid = match(top_upa, colnames(sym.m1))
    
    x[Eidx] = rep(0, num_loc)
    x[((seedid - 1)*num_states)+2] = 100 # seed infected in latent compartment, in metapopulations defined by seedid
    x[Sidx] = pop0 -  x[Eidx]   # everyone except init. infected start as susceptible in the beginning 
    x[Isidx] = rep(0, num_loc)
    x[Iaidx] = rep(0,num_loc)
    x[obsidx] = rep(0,num_loc)
    x[moveidx] = rep(0,num_loc)
    x[leaveiidx] = rep(0,num_loc)
    x[betaidx] = r0/(alpha*D + (1 - alpha)*mu*D)  # transmission rate due to symptomatic individuals
    x[muidx] = 0.50      #  multiplicative factor for the transmission rate of asymptomatic individuals
    x[thetaidx] = 1     # mult. factor, which is >1 to reflect underreporting of human movement
    x[Zidx] = 3        # average latency period
    x[alphaidx] = 0.65   # fraction of symptomatic infections
    x[Didx] = 5          # average duration of infection
    
    x.save <- matrix(NA, nrow = length(x), ncol = num_tsteps+1)
    pop.save <- matrix(NA, nrow = length(pop), ncol = num_tsteps+1)
    x.save[,1] <- x
    pop.save[,1] <- pop0
    
    for (ts in 1:num_tsteps) {
      sim.tmp <- spatial_seir(x, M, pop, ts, pop0)
      # update values
      x  = x.save[,ts+1] = sim.tmp$x 
      pop = pop.save[,ts+1] = sim.tmp$pop 
      
    }
    
    # save each run 
    
    x.save.list[[run]] = x.save
    pop.save.list[[run]] = pop.save
    
  } # end of runs for loop
  
  # process results ----
  geo_codes = colnames(sym.m1)
  
  pop.update = do.call(rbind, pop.save.list) %>%
    as.data.frame() %>%
    mutate(sim = rep(1:n_runs, each = nrow(pop.save.list[[1]]))) %>%
    mutate(geo_code = rep(geo_codes, times = n_runs)) %>%   
    pivot_longer(!c(geo_code, sim), names_to = "time", values_to = "pop") %>%
    mutate(time = as.numeric(sub('.', '', time))) # nrow = 493406
  
  #x.save.list 
  # No cols = time steps
  # Rows = 5 state outputs from spatial_seir.R function * # locations + 6 parameter values
  
  results = do.call(rbind, x.save.list) %>%
    as.data.frame() %>% 
    mutate(sim = rep(1:n_runs, each = nrow(x.save.list[[1]]))) %>%
    mutate(state = rep(c(rep(c("S", "E", "Is", "Ia", "obs", "move", "leavei"), times = num_loc), 
                         "beta", "mu", "theta", "Z", "alpha", "D"), times = n_runs)) %>%
    pivot_longer(!c(sim, state), names_to = "time", values_to = "value") %>%
    mutate(time = as.numeric(sub('.', '', time))) %>% 
    filter(!state %in% c("beta", "mu", "theta", "Z", "alpha", "D")) %>% 
    mutate(geo_code = rep(rep(geo_codes, each = (num_tsteps+1)*num_states), times = n_runs)) %>% 
    left_join(., pop.update, by = c("geo_code", "time", "sim")) %>%
    # parameters:
    mutate(source = mobility_source) %>%
    mutate(r0 = r0) %>%
    mutate(seed_city = seed_city)
  
  # average df - average across sims
  results.avg = results %>%
    group_by(state, time, geo_code, r0, seed_city, source) %>%
    summarise(mean = mean(value, na.rm=TRUE),
              lower = quantile(value, probs = 0.025, na.rm=TRUE),
              upper = quantile(value, probs = 0.975, na.rm=TRUE), #,
              lower0.1 = quantile(value, probs = 0.1, na.rm=TRUE),
              upper0.9 = quantile(value, probs = 0.9, na.rm=TRUE),
              mean.pop = mean(pop,na.rm= TRUE)) %>% 
    left_join(., district_upa, by =c("geo_code" = "upazila_code")) %>%
    ungroup() # merge with district for mapping %>%
  
  ## final size and # districts with at least one case metric for barplot
  pop_bang = sum(pop_orig_data$pop)
  
  sum_bydist = results  %>%
    filter(state == "obs") %>%
    group_by(time, source, r0, seed_city, sim) %>%
    # aggregate across geo_code districts
    summarise(value = sum(value))
  
  final_size = sum_bydist %>%
    group_by(source, r0, seed_city, sim) %>%
    arrange(time) %>%
    # cumulative sum
    mutate(final_value = cumsum(value)) %>% 
    filter(time == 526) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(final_size_prop = mean(final_value, na.rm= TRUE)/pop_bang, 
              final_size_prop_lower = quantile(final_value, probs = 0.025, na.rm=TRUE)/pop_bang,
              final_size_prop_upper = quantile(final_value, probs = 0.975, na.rm=TRUE)/pop_bang) %>%
    ungroup()
  
  num_dist_plt = results %>%
    filter(state == "obs") %>%
    filter(time <= 60 & value >= 1) %>% 
    group_by(source, r0, seed_city, sim) %>% #
    summarise(num_dist = n_distinct(geo_code)) %>% 
    group_by(source, r0, seed_city) %>%
    summarise(num_dist = mean(num_dist, na.rm= TRUE), 
              num_dist_lower = quantile(num_dist, probs = 0.025, na.rm=TRUE),
              num_dist_upper = quantile(num_dist, probs = 0.975, na.rm=TRUE))
  
  return(list(results.avg, sum_bydist, final_size, num_dist_plt))
}

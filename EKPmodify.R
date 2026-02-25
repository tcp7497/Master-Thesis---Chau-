##########################################################################
# Set of functions defining epidemiological quantities
##########################################################################

# function that assign the length of the infectious period
# Random infectious period (Gamma), in days
# mean_ip and sd_ip are interpretable
infection.period.length <- function(mean_ip = 6, sd_ip = 2){
  shape <- (mean_ip / sd_ip)^2
  scale <- (sd_ip^2) / mean_ip
  rgamma(1, shape = shape, scale = scale)
}

# infectiousness measure - intrinsic generation time distribution
# Gamma infectiousness profile (PDF)
infectiousness.measure <- function(t, mean = 4.5, sd = 2.0){
  if (t < 0) return(0)
  
  shape <- (mean / sd)^2
  scale <- (sd^2) / mean
  
  # density (units 1/day), integrates to 1 over [0, inf)
  return(dgamma(t, shape = shape, scale = scale))
}


# ------------------------------------------------------------------------
# Population Matrix
# ------------------------------------------------------------------------
pop.matrix <- function(n, age.distr = NULL){
  status.matrix <- data.frame(
    infected          = rep(0, n), # 1-infecte, 0-suceptible, -1 recovered
    time.of.infection = rep(NA, n),
    infector          = rep(NA, n), # who infected this person
    Recovery          = rep(Inf, n), # when individual recovers from infection
    ct.rates          = rep(0,n), # ct rate of the individual
    infected_id       = 1:n,
    z_inf             = rep(1, n)   # individual infectivity multiplier
  )
  return(status.matrix)
}

make_network_er <- function(n, p = 0.15){
  A <- matrix(rbinom(n*n, 1, p), n, n)
  diag(A) <- 0
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A
}

# ------------------------------------------------------------------------
# Main Simulation
# ------------------------------------------------------------------------
sim.epi <- function(n, ct.rates, nSeeds, inf.param, pathogen = NULL, t.stop,
                    k = 0.2,  # heterogeneity: Var(z)=1/k
                    nu_mean = 4.5, nu_sd = 2.0,
                    mean_ip = 6, sd_ip = 2,
                    A = NULL){
  
  # create a matrix in which epidemiological states and other infection
  # details are reported
  status.matrix <- pop.matrix(n)
  
  if (is.null(A)) {
    A <- make_network_er(n, p = 0.15)
  }
  neigh_list <- lapply(1:n, function(i) which(A[i, ] == 1))
  
  # Individual-level heterogeneity in infectiousness:
  # z_inf ~ Gamma(shape=k, scale=1/k) => mean=1, var=1/k
  status.matrix$z_inf <- rgamma(n, shape = k, scale = 1/k)
  
  # Contact rates - here it is assumed these are the same (homogeneous population)
  status.matrix$ct.rates <- ct.rates
  
  # Pre-allocation of variables that will keep track of outbreak characteristics
  max_events <- 2 * n + 200 
  time.events <- matrix(NA, nrow = max_events, ncol = 3)
  colnames(time.events) <- c("time", "event", "who")
  event_counter <- 0
  
  infectives <- rep(0, n) 
  current.time <- 1
  
  # Vector storing the next proposed contact time for every individual
  contact.times <- rep(Inf, n)
  
  # Initialization - adding initial cases
  first.cases <- sample(n, nSeeds)
  
  for (first in first.cases){
    status.matrix$infected[first] <- 1 
    status.matrix$time.of.infection[first] <- 0
    status.matrix$Recovery[first] <- current.time + infection.period.length(mean_ip = 6, sd_ip = 2)
    
    event_counter <- event_counter + 1
    time.events[event_counter, ] <- c(current.time, 1, first)
    
    infectives[first] <- 1
    # Generate next contact
    contact.times[first] <- current.time + rexp(1, status.matrix$ct.rates[first])
  }
  
  num_infectives <- length(first.cases)
  
  while(num_infectives > 0 & current.time < t.stop){
    
    # Identify next event times
    # events occur in continuous time - selecting and resolving the ones that occur sooner 
    
    
    # Only look at active infectives to speed up min()
    active_ids <- which(infectives == 1)
    if(length(active_ids) == 0) break
    
    min_ctc_time <- min(contact.times[active_ids])
    min_rec_time <- min(status.matrix$Recovery[active_ids])
    
    if (min_ctc_time < min_rec_time) {
      
      # --- Next Event: Contact ---
      current.time <- min_ctc_time
      
      # Identify infector (handle ties)
      candidates <- active_ids[contact.times[active_ids] == min_ctc_time]
      infector <- if(length(candidates) > 1) sample(candidates, 1) else candidates
      
      # Schedule next contact for this infector
      contact.times[infector] <- current.time + rexp(1, status.matrix$ct.rates[infector])
      
      # Structured mixing
      neighbors <- neigh_list[[infector]]
      if (length(neighbors) == 0) {
        next
      }
      infectee <- sample(neighbors, 1)
      
      # Check if the contacted person is susceptible
      if (status.matrix$infected[infectee] == 0) {
        
        # Calculate acceptance rate using the function directly
        t_since_inf <- current.time - status.matrix$time.of.infection[infector]
        
        acc.rate <- infectiousness.measure(t_since_inf, mean = nu_mean, sd = nu_sd) *
          inf.param *
          status.matrix$z_inf[infector]
        
        # safety: probability must be <= 1
        acc.rate <- min(acc.rate, 1)
        
        # acc.rate should always be less or equal to 1 (it's a probability)
        # sometimes setting could make this value exceed 1 -> check that it does not occur in your system
        
        if (acc.rate > 0 && runif(1) < acc.rate) { #a contact is successful -> infection
          status.matrix$infected[infectee] <- 1
          status.matrix$time.of.infection[infectee] <- current.time
          status.matrix$infector[infectee] <- infector
          status.matrix$Recovery[infectee] <- current.time + infection.period.length(mean_ip = 6, sd_ip = 2)
          
          event_counter <- event_counter + 1
          time.events[event_counter, ] <- c(current.time, 1, infectee)
          
          infectives[infectee] <- 1
          num_infectives <- num_infectives + 1
          
          # Schedule first contact for new infectee
          contact.times[infectee] <- current.time + rexp(1, status.matrix$ct.rates[infectee])
        }
      }
      
    } else {
      # --- Next Event: Recovery ---
      current.time <- min_rec_time
      
      candidates <- active_ids[status.matrix$Recovery[active_ids] == min_rec_time]
      # update the status matrix for the individual who recovers (not infectious anymore-> so no generating contacts)
      for (recovered in candidates) {
        status.matrix$infected[recovered] <- -1
        infectives[recovered] <- 0
        num_infectives <- num_infectives - 1
        contact.times[recovered] <- Inf
        event_counter <- event_counter + 1
        time.events[event_counter, ] <- c(current.time, -1, recovered)
      }
    }
  }
  
  # Trim Matrix
  if(event_counter > 0) time.events <- time.events[1:event_counter, , drop=FALSE]
  
  # ------------------------------------------------------------------------
  # Summary Calculation
  # ------------------------------------------------------------------------
  
  # Create daily counts
  all_days <- floor(time.events[,1] + 0.001)
  max_day <- floor(max(time.events[,1])) + 1
  
  inf_events <- time.events[time.events[,2] == 1, 1]
  rec_events <- time.events[time.events[,2] == -1, 1]
  
  day_counts_inf <- table(factor(floor(inf_events + 0.001), levels = 1:max_day))
  day_counts_rec <- table(factor(floor(rec_events + 0.001), levels = 1:max_day))
  
  daily_inc <- as.numeric(day_counts_inf)
  daily_prev <- (cumsum(daily_inc) + nSeeds) - cumsum(as.numeric(day_counts_rec))
  
  sim.info <- data.frame(
    day = 1:max_day,
    dailyInc = daily_inc,
    dailyPrev = daily_prev
  )
  
  # Weekly Aggregation
  sim.info$week <- ceiling(sim.info$day / 7)
  sim.info$weeklyInc <- ave(sim.info$dailyInc, sim.info$week, FUN=sum)
  
  # Get weekly prevalence (taking the prevalence at the end of the week)
  last_val <- function(x) tail(x, 1)
  weekly_prev_vec <- tapply(sim.info$dailyPrev, sim.info$week, last_val)
  sim.info$weeklyPrev <- weekly_prev_vec[sim.info$week]
  
  return(list(time.events=time.events, status.matrix=status.matrix, sim.info=sim.info))
}
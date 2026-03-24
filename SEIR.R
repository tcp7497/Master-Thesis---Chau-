##########################################################################
# Set of functions defining epidemiological quantities
##########################################################################

# incubation period (Stage E of SEIR model) - Weibull Distribution
incubation.period <- function(mean_e = 5.2, sd_e = 2.8) {
  cv_e <- sd_e / mean_e
  shape_e <- cv_e^-1.086 
  scale_e <- mean_e / gamma(1 + 1/shape_e)
  return(rweibull(1, shape = shape_e, scale = scale_e))
}

# infectiousness measure v(t) (Stage I of SEIR model) - Gamma distribution
infectiousness.measure <- function(t, mean_i = 10.0, sd_i = 3.8){
  if (t < 0) return(0)
  shape_i <- (mean_i / sd_i)^2
  scale_i <- (sd_i^2) / mean_i
  
  # density (units 1/day), integrates to 1 over [0, inf)
  return(dgamma(t, shape = shape_i, scale = scale_i))
}


# ------------------------------------------------------------------------
# Population Matrix
# ------------------------------------------------------------------------
pop.matrix <- function(n, group_assign = NULL){
  if (is.null(group_assign)) 
    group_assign <- rep(1:n, length.out = n) #default: 1 group
  
  
  status.matrix <- data.frame(
    infected          = rep(0, n),   # 1:infectee, 0.5:exposed, 0:suceptible, -1: recovered
    time.of.infection = rep(NA, n),
    onset_time        = rep(NA,n),   # time to start the infection
    infector          = rep(NA, n),  # who infected this person
    Recovery          = rep(Inf, n), # when individual recovers from infection
    ct.rates          = rep(0,n),    # ct rate of the individual
    infected_id       = 1:n,
    z_inf             = rep(1, n),   # individual infectivity multiplier
    group             = group_assign
  )
  return(status.matrix)
}

choose_infectee_structured <- function(n, infector, groups, M){
  gi <- groups[infector]              # group of infector
  probs_group <- M[gi, ]              # row = infector group
  
  # normalize (just in case)
  probs_group <- probs_group / sum(probs_group)
  
  # sample target group
  g_target <- sample(seq_along(probs_group), size = 1, prob = probs_group)
  
  candidates <- which(groups == g_target & seq_len(n) != infector)
  
  # fallback if group has only infector (rare edge case)
  if (length(candidates) == 0) candidates <- setdiff(seq_len(n), infector)
  if (length(candidates) > 1) sample(candidates, size = 1) else candidates
}


# ------------------------------------------------------------------------
# Main Simulation
# ------------------------------------------------------------------------
sim.epi_seir <- function(n, ct.rates, nSeeds, inf.param, t.stop,
                    k = 0.2,  # heterogeneity: Var(z) = 1/k
                    mean_e = 5.2, sd_e = 2.8,
                    mean_i = 10.0, sd_i = 3.8,
                    group_assign = NULL, M = NULL, A = NULL){
  
  # create a matrix in which epidemiological states and other infection
  # details are reported
  status.matrix <- pop.matrix(n, group_assign = group_assign)
  
  # individual-level heterogeneity in infectiousness:
  # z_inf ~ Gamma(shape=k, scale=1/k) => mean=1, var=1/k
  status.matrix$z_inf <- rgamma(n, shape = k, scale = 1/k)
  
  # contact rates - here it is assumed these are the same (homogeneous population)
  status.matrix$ct.rates <- ct.rates
  
  # time when symptoms/infectiousness starts
  status.matrix$onset_time <- NA 
  
  # pre-allocation of variables that will keep track of outbreak characteristics
  max_events <- 5 * n + 200
  time.events <- matrix(NA, nrow = max_events, ncol = 3)
  colnames(time.events) <- c("time", "event", "who")
  event_counter <- 0
  
  # binary vector for those currently in stage I
  infectives <- rep(0, n) 
  current.time <- 0
  
  # vector storing the next proposed contact time for every individual
  contact.times <- rep(Inf, n)
  
  # initial Seeds (S -> E) - adding initial cases
  first.cases <- sample(n, nSeeds)
  
  for (first in first.cases){
    status.matrix$infected[first] <- 0.5 
    status.matrix$time.of.infection[first] <- 0
    
    # Calculate when the person will actually become infectious
    incubation_p <- incubation.period(mean_e, sd_e)
    status.matrix$onset_time[first] <- current.time + incubation_p
    status.matrix$Recovery[first] <- status.matrix$onset_time[first] + 18 # symptomatic period
    
    event_counter <- event_counter + 1
    time.events[event_counter, ] <- c(current.time, 0.5, first)
  }
  
  while(current.time < t.stop){
    
    # identify next event times
    exposed_ids <- which(status.matrix$infected == 0.5)
    min_onset <- if(length(exposed_ids) > 0) 
      min(status.matrix$onset_time[exposed_ids]) else Inf
    
    # events occur in continuous time - selecting and resolving the ones that occur sooner 
    # only look at active infectives to speed up min()
    active_ids <- which(infectives == 1)
    if(length(active_ids) == 0) break
    min_ctc_time <- if(length(active_ids) > 0) 
      min(contact.times[active_ids]) else Inf # The earliest time infected made contact
    min_rec_time <- if(length(active_ids) > 0) 
      min(status.matrix$Recovery[active_ids]) else Inf #The earliest time infected has recovered from the illness
      
    # next event: Contact 
    current.time <- min(min_onset, min_ctc_time, min_rec_time)
    if (current.time > t.stop || current.time == Inf) break
      
      
    if (current.time == min_onset) {
      # stage E -> I (onset of infectiousness)
      candidates <- exposed_ids[which.min(status.matrix$onset_time[exposed_ids])]
      status.matrix$infected[candidates] <- 1
      infectives[candidates] <- 1
      contact.times[candidates] <- current.time + rexp(1, status.matrix$ct.rates[candidates])
      
      event_counter <- event_counter + 1
      time.events[event_counter, ] <- c(current.time, 1, candidates)
        
    } else if (current.time == min_ctc_time) {
      # Contact event
      candidates_inf <- active_ids[which.min(contact.times[active_ids])]
        
    # Select candidates based on scenario  
    if (!is.null(group_assign) && !is.null(M)) {
      # SCENARIO 1: Group-based structured mixing
      infectee <- choose_infectee_structured(
        n, candidates_inf, group_assign, M)
    } else if (!is.null(A)) {
        
      # SCENARIO 2: Network-based mixing
      neighbors <- which(A[candidates_inf, ] == 1)
      infectee <- if(length(neighbors) > 0) sample(neighbors, 1) else next
        
    } else {
      # SCENARIO 3: Free mixing (baseline)
      infectee <- sample(setdiff(1:n, candidates_inf), 1)
    }
      
      
    # if target is susceptible, try to infect
    if (status.matrix$infected[infectee] == 0) {
        
      # calculate acceptance rate using the function directly
      t_since_inf <- current.time - status.matrix$onset_time[candidates_inf]
      acc.rate <- infectiousness.measure(t_since_inf, mean = mean_i, sd = sd_i) *
        inf.param * status.matrix$z_inf[candidates_inf]
        
      # acc.rate should always be less or equal to 1 (it's a probability)
      if (runif(1) < min(acc.rate, 1)) { #a contact is successful -> infection
        status.matrix$infected[infectee] <- 0.5
        status.matrix$time.of.infection[infectee] <- current.time
        status.matrix$onset_time[infectee] <- current.time + incubation.period(mean_e, sd_e)
        status.matrix$Recovery[infectee] <- status.matrix$onset_time[infectee] + 18
          
        event_counter <- event_counter + 1
        time.events[event_counter, ] <- c(current.time, 0.5, infectee)
        }
      }   
          
        # schedule first contact for new infectee
        contact.times[candidates_inf] <- current.time + rexp(1, status.matrix$ct.rates[candidates_inf])
      
    } else {
      # Stage I -> R (Recovery)
      candidates_rec <- active_ids[which.min(status.matrix$Recovery[active_ids])]
     
      # update the status matrix for the individual who recovers (not infectious anymore-> so no generating contacts)
      for (recovered in candidates_rec) {
        status.matrix$infected[recovered] <- -1
        infectives[recovered] <- 0
        contact.times[recovered] <- Inf
        event_counter <- event_counter + 1
        time.events[event_counter, ] <- c(current.time, -1, recovered)
      }
    }
    return(list(time.events = time.events[1:event_counter, ], status.matrix = status.matrix))
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
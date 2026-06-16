# =========================================================================
# Set of functions defining epidemiological quantities
# =========================================================================
# Incubation period (Stage E of SEIR model) - Weibull Distribution
incubation.period <- function(mean_e = 5.2, sd_e = 2.8) {
  cv_e <- sd_e / mean_e
  shape_e <- cv_e^-1.086 
  scale_e <- mean_e / gamma(1 + 1/shape_e)
  return(rweibull(1, shape = shape_e, scale = scale_e))
}

# Infectiousness measure v(t) (Stage I of SEIR model) - Gamma distribution
infectiousness.measure <- function(t, mean_i = 10.0, sd_i = 3.8){
  if (t < 0) return(0)
  shape_i <- (mean_i / sd_i)^2
  scale_i <- (sd_i^2) / mean_i
  
  # Density (units 1/day), integrates to 1 over [0, inf)
  return(dgamma(t, shape = shape_i, scale = scale_i))
}

# Multipliers based on behavior
mod_staff <- c("morning"=1.5, "lunch"=1.0, "afternoon"=1.5, "dinner"=1.0, "evening"=0.8, "night"=0.6)
mod_vol   <- c("morning"=1.2, "lunch"=2.0, "afternoon"=1.0, "dinner"=1.8, "evening"=1.5, "night"=0.5)

# Time blocks
get_dynamic_rate <- function(role_single, current_time, base_ct_rate) {
  # Convert days to hours (0 to 24)
  hr <- (current_time * 24) %% 24
  
  if (hr >= 6 && hr < 11) block <- "morning"
  else if (hr >= 11 && hr < 13) block <- "lunch"
  else if (hr >= 13 && hr < 18) block <- "afternoon"
  else if (hr >= 18 && hr < 20) block <- "dinner"
  else if (hr >= 20 && hr < 22) block <- "evening"
  else block <- "night"
  
  mod <- if(role_single == "staff") mod_staff[block] else mod_vol[block]
  return(base_ct_rate * as.numeric(mod))
}

# Weighted mean of time blocks
block_hours <- c("morning" = 5, "lunch" = 2, "afternoon" = 5, 
                 "dinner" = 2, "evening" = 2, "night" = 8)

weighted_mean_mult <- function(mod_vec) {
  sum(mod_vec * block_hours[names(mod_vec)]) / 24
}

avg_mod_staff <- weighted_mean_mult(mod_staff)
avg_mod_vol   <- weighted_mean_mult(mod_vol)

average_mod <- (avg_mod_staff * 4 + avg_mod_vol * 24) / 28

# =========================================================================
# Population Matrix
# =========================================================================
pop.matrix <- function(n, roles, group_assign = NULL){
  status.matrix <- data.frame(
    infected          = rep(0, n),   # 1:infectee, 0.5:exposed, 0:suceptible, -1: recovered
    time.of.infection = rep(NA, n),
    onset_time        = rep(NA,n),   # time to start the infection
    infector          = rep(NA, n),  # who infected this person
    Recovery          = rep(Inf, n), # when individual recovers from infection
    ct.rates          = rep(0,n),    # ct rate of the individual
    infected_id       = 1:n,
    roles             = roles,
    z_inf             = rep(1, n),   # individual infectivity multiplier
    z_sus             = rep(1, n),   # individual susceptibility multiplier
    group             = if(is.null(group_assign)) rep(1, n) else group_assign
  )
  return(status.matrix)
}

choose_infectee_structured <- function(n, infector, groups, M){
  gi <- groups[infector]              # group of infector
  probs_group <- M[gi, ]              # row = infector group
  
  # Normalize (just in case)
  probs_group <- probs_group / sum(probs_group)
  
  # Sample target group
  g_target <- sample(seq_along(probs_group), size = 1, prob = probs_group)
  
  candidates <- which(groups == g_target & seq_len(n) != infector)
  
  # Fallback if group has only infector (rare edge case)
  if (length(candidates) == 0) candidates <- setdiff(seq_len(n), infector)
  if (length(candidates) > 1) sample(candidates, size = 1) else candidates
}

# =========================================================================
# Main Simulation
# =========================================================================
sim.epi_seir <- function(n, ct.rates, nSeeds, R0, t.stop, roles, average_mod,
                         k = 0.2,     # heterogeneity in infectiousness
                         k_ct = 1.0,  # heterogeneity in contact rates
                         k_sus = 1.0, # heterogeneity in susceptibility
                         mean_e = 5.2, sd_e = 2.8,
                         mean_i = 10.0, sd_i = 3.8,
                         group_assign = NULL, M = NULL,
                         vacc_coverage = 0,
                         vacc_efficacy_sus = 0,
                         vacc_target = "all"){
  
  # Create a matrix in which epidemiological states and other infection
  status.matrix <- pop.matrix(n, roles, group_assign)
  
  # Individual-level heterogeneity in infectiousness:
  # z_inf ~ Gamma(shape=k, scale=1/k) => mean=1, var=1/k
  status.matrix$z_inf <- rgamma(n, shape = k, scale = 1/k)
  
  # Heterogeneity in susceptibility
  status.matrix$z_sus <- rgamma(n, shape = k_sus, scale = 1/k_sus)
  
  # Vaccination
  if (vacc_coverage > 0) {
    if (vacc_target == "vol") {
      eligible_ids <- which(roles == "vol")
    } else if (vacc_target == "staff") {
      eligible_ids <- which(roles == "staff")
    } else {
      eligible_ids <- 1:n  # "all"
    }
      n_vacc <- round(vacc_coverage * length(eligible_ids))
      vacc_ids <- sample(eligible_ids, n_vacc)
      if (vacc_efficacy_sus > 0) {
        status.matrix$z_sus[vacc_ids] <- status.matrix$z_sus[vacc_ids] * (1 - vacc_efficacy_sus)
      }
  }
  
  # Individual-level heterogeneity in contact rates
  status.matrix$ct.rates <- rgamma(n, shape = k_ct, scale = ct.rates / k_ct)
  
  # Calculate infection parameter based on average population contacts
  inf.param <- R0 / (mean(status.matrix$ct.rates) * average_mod)
  
  # Time when symptoms/infectiousness starts
  status.matrix$onset_time <- NA 
  
  # Pre-allocation of variables that will keep track of outbreak characteristics
  max_events <- 5 * n + 200
  time.events <- matrix(NA, nrow = max_events, ncol = 3)
  colnames(time.events) <- c("time", "event", "who")
  event_counter <- 0
  
  # Binary vector for those currently in stage I
  infectives <- rep(0, n) 
  current.time <- 0
  
  # Vector storing the next proposed contact time for every individual
  contact.times <- rep(Inf, n)
  
  # Initial Seeds (S -> E) - adding initial cases
  first.cases <- sample(n, nSeeds)
  
  for (first in first.cases){
    status.matrix$infected[first] <- 0.5 
    status.matrix$time.of.infection[first] <- 0
    
    # Calculate when the person will actually become infectious
    status.matrix$onset_time[first] <- current.time + incubation.period(mean_e, sd_e)
    status.matrix$Recovery[first] <- status.matrix$onset_time[first] + 18 # symptomatic period
    
    event_counter <- event_counter + 1
    time.events[event_counter, ] <- c(current.time, 0.5, first)
  }
  
  while(current.time < t.stop){
    # identify next event times
    exposed_ids <- which(status.matrix$infected == 0.5)
    min_onset <- if(length(exposed_ids) > 0) 
      min(status.matrix$onset_time[exposed_ids]) else Inf
    
    # Events occur in continuous time - selecting and resolving the ones that occur sooner 
    # Only look at active infectives to speed up min()
    active_ids <- which(infectives == 1)
    
    # Check if simulation should end (no more exposed or infectious individuals)
    if(length(exposed_ids) == 0 && length(active_ids) == 0) break
    
    min_ctc_time <- if(length(active_ids) > 0) 
      min(contact.times[active_ids]) else Inf # The earliest time infected made contact
    min_rec_time <- if(length(active_ids) > 0) 
      min(status.matrix$Recovery[active_ids]) else Inf #The earliest time infected has recovered from the illness
    
    # Next event: Contact 
    current.time <- min(min_onset, min_ctc_time, min_rec_time)
    if (current.time > t.stop || current.time == Inf) break
    
    
    if (current.time == min_onset) {
      # Stage E -> I (onset of infectiousness)
      candidates <- exposed_ids[which.min(status.matrix$onset_time[exposed_ids])]
      status.matrix$infected[candidates] <- 1
      infectives[candidates] <- 1
      
      # Contact time using DYNAMIC rate
      d_rate <- get_dynamic_rate(status.matrix$roles[candidates], current.time, status.matrix$ct.rates[candidates])
      contact.times[candidates] <- current.time + rexp(1, d_rate)
      
      event_counter <- event_counter + 1
      time.events[event_counter, ] <- c(current.time, 1, candidates)
      
    } else if (current.time == min_ctc_time) {
      # Contact event
      infector <- active_ids[which.min(contact.times[active_ids])]
      
      # Choose a new contact person based on scenarios
      if (!is.null(group_assign) && !is.null(M)) {
        
        # Scenario 2: Staff-Bridge & SCENARIO 3: Strict-Cohorting
        infectee <- choose_infectee_structured(n, infector, group_assign, M)
      } else {
        # Scenario 1: Free mixing
        infectee <- sample(setdiff(1:n, infector), 1)
      }
      
      # If target is susceptible, try to infect
      if (status.matrix$infected[infectee] == 0) {
        
        # Calculate acceptance rate using the function directly
        t_since_inf <- current.time - status.matrix$onset_time[infector]
        acc.rate <- infectiousness.measure(t_since_inf, mean_i = mean_i, sd_i = sd_i) *
          inf.param * status.matrix$z_inf[infector] * status.matrix$z_sus[infectee]
        
        # acc.rate should always be less or equal to 1 (it's a probability)
        if (runif(1) < min(acc.rate, 1)) { #a contact is successful -> infection
          status.matrix$infected[infectee] <- 0.5
          status.matrix$time.of.infection[infectee] <- current.time
          status.matrix$onset_time[infectee] <- current.time + incubation.period(mean_e, sd_e)
          status.matrix$Recovery[infectee] <- status.matrix$onset_time[infectee] + 18
          status.matrix$infector[infectee] <- infector
          
          event_counter <- event_counter + 1
          time.events[event_counter, ] <- c(current.time, 0.5, infectee)
        }
      }   
      
      # Schedule next contact for the infector
      d_rate <- get_dynamic_rate(status.matrix$roles[infector], current.time, status.matrix$ct.rates[infector])
      contact.times[infector] <- current.time + rexp(1, d_rate)
      
    } else {
      # Stage I -> R (Recovery)
      candidates_rec <- active_ids[which.min(status.matrix$Recovery[active_ids])]
      
      # Update the status matrix for the individual who recovers (not infectious anymore-> so no generating contacts)
      for (recovered in candidates_rec) {
        status.matrix$infected[recovered] <- -1
        infectives[recovered] <- 0
        contact.times[recovered] <- Inf
        event_counter <- event_counter + 1
        time.events[event_counter, ] <- c(current.time, -1, recovered)
      }
    }
  }
  
  # Trim Matrix
  if(event_counter > 0) time.events <- time.events[1:event_counter, , drop=FALSE]
  
# =========================================================================
# Summary Calculation
# =========================================================================
  
  if(event_counter > 0){
    # Create daily counts
    all_days <- floor(time.events[,1] + 0.001)
    max_day <- floor(max(time.events[,1])) + 1
    
    inf_events <- time.events[time.events[,2] == 0.5, 1]
    rec_events <- time.events[time.events[,2] == -1, 1]
    
    day_counts_inf <- table(factor(floor(inf_events + 0.001), levels = 1:max_day))
    day_counts_rec <- table(factor(floor(rec_events + 0.001), levels = 1:max_day))
    
    daily_inc <- as.numeric(day_counts_inf)
    daily_prev <- cumsum(daily_inc) - cumsum(as.numeric(day_counts_rec))
    
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
    
  } else {
    # No events occurred
    sim.info <- data.frame(
      day = 1,
      dailyInc = 0,
      dailyPrev = nSeeds,
      week = 1,
      weeklyInc = 0,
      weeklyPrev = nSeeds)
  }
  return(list(time.events=time.events, status.matrix=status.matrix, sim.info=sim.info))
}

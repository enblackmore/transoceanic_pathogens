
#(1) A function which builds the state vector (x0), propensity vector (a),
#state-change matrix (nu), and model parameters list (parms) needed to run simulations
#using GillespieSSA, given the following paramenters:

#N: total people on ship 
#S: total susceptible people on ship
#e0: number of exposed people at t=0

#bdd: transmission rate under density-dependent transmission
#bfd: transmission rate under frequency-dependent transmission
#q:   index determining relative contribution of density- and 
#frequency-dependent transmission (beta = bfd^(1-q) * (bdd*N)^q)

#mue: mean dwell time in state E
#mui: mean dwell time in state I

#ke:  shape of Erlang-distributed dwell time in state E
#ki:  shape of Erlang-distributed dwell time in state I

#generation_tracking: logical. Determines whether model should track generations
#setting this parameter to TRUE increases run time
#so we default to FALSE for analyses that do not consider the number of transmissin generations

#generation_max: the number of generations the model should track
#adding generations increases run time, so best to set this as low as possible
#.setting too low will underestimate the total number of generations
#because any cases in generations after generation_max will be stored in that final compartment.
#e.g. if generation_max is set to 5, then any cases in generations 6 or higher will be added to the generation 5 totals, 
#so the generations tracked would be 1,2,3,4,5+.

#the function generation_test() (below) tests whether generation_max is appropriate for a given simulation
#by checking how often cases appear in generation compartments close to generation_max 

build.model<- function(
    N=500, 
    S=100, 
    e0=1,
    bdd=1/1000,
    bfd=0,
    q=1,
    mue=5,
    mui=5,
    ke=3,
    ki=3, 
    generation_tracking=FALSE,
    generation_max=10) {
  
  #(1) Check inputs
  
  #Shape of the Erlang distributions for mue and mui must be an integer
  if(as.integer(ke) != ke){
    return(print("Error: ke must be specified as integer"))
  }
  if(as.integer(ki) != ki){
    return(print("Error: ki must be specified as integer"))
  }
  
  #Number of generations to track must be an integer
  if(generation_tracking==TRUE & as.integer(generation_max) != generation_max){
    return(print("Error: generation tracking must be specified as integer"))
  }
  
  #Total people in S and e0 can't be greater than N
  if(S+e0>N){
    return(print("S+e0 greater than N"))
  }
  
  #q must be between 1 and 0
  if(q < 0 | q > 1){
    return("error: q must be a value between 1 and 0")
  }
  
  if(generation_tracking==TRUE & generation_max<2){
    return(print("error: generation_max must be an integer of at least 2"))
  }
  
  #(2) Rewrite some parameters
  #calculate the rate for Erlang distributed E and I dwell times
  gammae <- ke/mue
  gammai <- ki/mui
  
  #calculate beta from bdd, bfd, q, and N
  beta <- (((bdd*N)^q)*bfd^(1-q))/N
  
  #combine in list parms for GillespieSSA
  parms <- c(beta=beta, gammae=gammae, gammai=gammai)
  
  ##Case A: generation_tracking=FALSE
  
  if(generation_tracking==FALSE){
    #(3A) Build state vector x0, with appropriate names for all compartments
    
    #3A.1: calculate the total numbers of compartments required
    n_EI <- ke+ki #total E and I compartments 
    n_compartments_total <- n_EI+2 #total number of model compartments (S+R=2)
    
    #3A.2: write appropriate names for all compartments
    
    names_E <- paste(rep("E", times=ke), letters[1:ke], sep="") 
    names_I <- paste(rep("I", times=ki), letters[1:ki], sep="") 
    
    #add states S and R for a full vector of state names
    names_SEIR <- c("S", names_E, names_I, "R") 
    
    #3A.3 use names_SEIR to write a x0, a vector of states at t=0
    x0 <-  numeric(n_compartments_total)
    names(x0) <- names_SEIR
    
    #Add state values at t=0
    x0[1] <- S 
    x0[length(x0)] <- N-S-e0 #model assumes that anyone who isn't in states S or state E at t=0 is in state R
    
    #NB: initial e0 values will be added to this function's output for each run using a separate function
    #This allows us to randomly distribute the e0 infectious individuals across the ke exposed states for each simulation
    #Without having to re-construct x0, a, nu, and parms each time
    
    
    #(4A) Build state change matrix, nu
    
    #nu_indices: a working matrix
    #each row represents one transition in nu
    #column 1 records vector x0 index of the state from which individuals transition
    #column 2 records vector x0 index of the state into which individuals transition
    
    #total number of transitions is all disease progression transitions + 1 for infection
    nu_indices <- matrix(NA, nrow=(n_EI+1), ncol=2) 
    colnames(nu_indices) <- c("From state", "To state")
    
    #4A.1: indices for disease progression and infection
    
    #every state progresses to the next state 
    #with the exception of state R, which does not progress
    
    indices_SEI <- seq(1,length(x0)-1) # indices of all S, E and I compartments 
    
    nu_indices[c(1:length(indices_SEI)),1] <- indices_SEI #input progression indices into matrix nu_indices: from state i...
    nu_indices[c(1:length(indices_SEI)),2] <- indices_SEI+1  #...to state i+1
    
    #4A.2: use  nu_indices to construct a state transition matrix, nu
    #in the state required by Gillespie SSA
    #where each row represents a state in vector x0
    #and each column represents a transition
    #such that a -1 in row i, column a represents a transition from state i
    #and a 1 in row j, column a represents a transition to state j
    #for any given transition, a
    
    nu <- matrix(0, ncol=(nrow(nu_indices)), nrow=(length(x0))) 
    row.names(nu) <- names(x0); nu
    
    for(i in 1:nrow(nu_indices)){
      nu[nu_indices[i,], i] <- c(-1,1)
    }
    
    #(5A) Write a propensity vector, a, of state transition equations
    #using standard mathematical notation
    #and the state names from vector x0
    #where each index in a corresponds to a state transition column in matrix nu
    
    #5A.1: progression and recovery equations
    
    #write progression equations for all generations:
    #dgEe/dt = gammae*gEe for all g generations and all e substates of E
    #dgIi/dt = gammae*gIi for all g generations and all i substates of I
    
    progression_equations <- paste(c(paste(rep("E", times=ke), letters[1:ke], sep=""),
                                     paste(rep("I", times=ki), letters[1:ki], sep="")), sep="")
    
    
    #add rate parameters 
    progression_equations <- paste(c(rep("gammae*", times=ke), rep("gammai*", times=ki-1)),
                                   progression_equations, sep="")
    
    
    #write the sum of the infectious compartments
    sum_infectious_compartments <- paste(names_I, collapse="+")
    
    #add S and beta
    transmission_equations <- paste("beta*S*(", sum_infectious_compartments, ")", sep="")
    
    #combine all state transition equations in propensity vector, a
    a <- c(transmission_equations, progression_equations)
    
    #return all parameters as a list; this ends the function
    return(list(x0, a, nu, parms, e0, ke))
  }
  
  #Case B: generation_tracking=TRUE
  
  #(3B) Build state vector x0, with appropriate names for all compartments
  
  #3.1: calculate the total numbers of compartments required
  n_EI_per_generation <- ke+ki #total E and I compartments per generation 
  n_EI_all_generations <- generation_max*n_EI_per_generation #total E and I compartments across all generations
  n_compartments_total <- n_EI_all_generations+2 #total number of model compartments (S+R=2)
  
  #3B.2: write appropriate names for all compartments
  
  #general names for E and I compartments without generation indicators: Ea, Eb,....,Eke; Ia, Ib,...,Iki
  names_E_per_generation <- paste(rep("E", times=ke), letters[1:ke], sep="") 
  names_I_per_generation <- paste(rep("I", times=ki), letters[1:ki], sep="") 
  
  #names for E and I compartments with generation indicators: 1Ea,...,1Eke, 2Ea,...,2Eke....
  names_EI <- paste(rep(c(1:generation_max), each=n_EI_per_generation), c(names_E_per_generation,names_I_per_generation), sep="") 
  
  #add states S and R for a full vector of state names
  names_SEIR <- c("S", names_EI, "R") 
  
  #3B.3 use names_SEIR to write a x0, a vector of states at t=0
  x0 <-  numeric(n_compartments_total)
  names(x0) <- names_SEIR
  
  #Add state values at t=0
  x0[1] <- S 
  x0[length(x0)] <- N-S-e0 #model assumes that anyone who isn't in states S or state E at t=0 is in state R
  
  #NB: initial e0 values will be added to this function's output for each run using a separate function
  #This allows us to randomly distribute the e0 infectious individuals across the ke exposed states for each simulation
  #Without having to re-construct x0, a, nu, and parms each time
  
  
  #(4B) Build state change matrix, nu
  
  #nu_indices: a working matrix
  #each row represents one transition in nu
  #column 1 records vector x0 index of the state from which individuals transition
  #column 2 records vector x0 index of the state into which individuals transition
  
  #the total number of transitions is:
  #all disease progression transitions (one for every E and I state)
  #all infection transitions (one per generation, minus the first generation)
  
  nu_indices <- matrix(NA, nrow=(n_EI_all_generations+(generation_max-1)), ncol=2) 
  colnames(nu_indices) <- c("From state", "To state")
  
  #4B.1: indices for disease progression transitions
  
  #for all E and I states
  #every state progresses to the next state 
  # with the exception of states I_ki, which progress to state R
  
  indices_EI_all <- seq(2,length(x0)-2) # indices of all E and I compartments in vector x0.
  indices_I_final <- seq(1+ke+ki, (generation_max*(ke+ki)+1), by=(n_EI_per_generation)) #removes I_ki compartments, which occur every 'n_EI_per_generation'th index
  indices_EI_except_final <- setdiff(indices_EI_all, indices_I_final) #index of all E and I compartments that progress to adjacent compartments
  
  nu_indices[c(1:length(indices_EI_except_final)),1] <- indices_EI_except_final #input progression indices into matrix nu_indices: from state i...
  nu_indices[c(1:length(indices_EI_except_final)),2] <- indices_EI_except_final+1  #...to state i+1
  
  #4B.2: indices for recovery transitions
  
  #every I_ki state progresses to state R
  #(the last state)
  
  n_indices_EI_except_final <- n_EI_all_generations - generation_max #record number of non-recovery EI transitions (one recovery transition per generation)
  n_indices_I_final <- generation_max #record number of recovery EI transitions (one per generation)
  
  #add n_indices_I_final recovery transitions to nu_indices after the n_indices_EI_except_final disease progression transitions
  nu_indices[c((n_indices_EI_except_final+1):(n_indices_EI_except_final+n_indices_I_final)),1] <- indices_I_final #from the (ke+ki)th compartment
  nu_indices[c((n_indices_EI_except_final+1):(n_indices_EI_except_final+n_indices_I_final)),2] <- length(x0) #to state R
  
  #4B.3: indices for infection transitions
  #state S can transition to states 2Ea, 3Ea,..., gEa.
  
  #isolate the indices for gEa for all values of g
  
  infection_indices <- seq(2+n_EI_per_generation, n_compartments_total-n_EI_per_generation, by=n_EI_per_generation)
  
  
  #use this to add to infection transitions to nu_indices
  nu_indices[(n_EI_all_generations+1):nrow(nu_indices),1] <- 1 #from state S
  nu_indices[(n_EI_all_generations+1):nrow(nu_indices),2] <- infection_indices #to state gEa
  
  #4B.4: finally, use the nu_indices matrix to construct a state transition matrix, nu
  #in the state required by Gillespie SSA
  #where each row represents a state in vector x0
  #and each column represents a transition
  #such that a -1 in row i, column a represents a transition from state i
  #and a 1 in row j, column a represents a transition to state j
  #for any given transition, a
  
  nu <- matrix(0, ncol=(nrow(nu_indices)), nrow=(length(x0))) 
  row.names(nu) <- names(x0); nu
  
  for(i in 1:nrow(nu_indices)){
    nu[nu_indices[i,], i] <- c(-1,1)
  }
  
  #(5B) Write a propensity vector, a, of state transition equations
  #using standard mathematical notation
  #and the state names from vector x0
  #where each index in a corresponds to a state transition column in matrix nu
  
  #5B.1: progression and recovery equations
  
  #write progression equations for all generations:
  #dgEe/dt = gammae*gEe for all g generations and all e substates of E
  #dgIi/dt = gammae*gIi for all g generations and all i substates of I
  
  progression_equations <- paste(rep(1:generation_max, each=n_EI_per_generation),
                                 c(paste(rep("E", times=ke), letters[1:ke], sep=""),
                                   paste(rep("I", times=ki), letters[1:ki], sep="")), sep="")
  
  
  #extract the recovery states, gIki, in a separate recovery vector to be included after all progression equations
  #since the order of vector 'a' needs to match column order in matrix 'nu'
  recovery_equations <- progression_equations[seq(n_EI_per_generation, 
                                                  generation_max*n_EI_per_generation, 
                                                  by=n_EI_per_generation)]
  progression_equations <- setdiff(progression_equations, recovery_equations)
  
  #add rate parameters to progression and recovery equations (gammae for state E, gammai for state I)
  progression_equations <- paste(c(rep(c(rep("gammae*", times=ke), rep("gammai*", times=ki-1)), 
                                       times=generation_max)),
                                 progression_equations, sep="")
  recovery_equations <- paste("gammai*", recovery_equations, sep="")
  
  #5B.2: write transmission equations
  #entry into state g+1Ea is a function of the sum of all ki subcompartments of state I in the previous generation, g:
  #d(g+1)Ea/dt = beta*S*sum(gIa+gIb+...+gIki)
  
  #5B.3:  transmission equations
  #entry into state g+1Ea is a function of the sum of all ki subcompartments of state I in the previous generation, g:
  #d(g+1)Ea/dt = beta*S*sum(gIa+gIb+...+gIki)
  
  #start by making make a matrix of infectious compartment names for all tracked generations
  infectious_compartments <- matrix(paste(rep(c(1:generation_max), each=ki), names_I_per_generation, sep=""), ncol=ki, byrow=TRUE)
  
  #collapse using lapply
  sum_infectious_compartments <- apply(infectious_compartments, 1, paste, collapse="+") 
  
  #add S and beta
  transmission_equations <- paste("beta*S*(", sum_infectious_compartments, ")", sep="")
  
  #create a loop so that non-tracked generations (generations>gen) produce generations in generation gen
  transmission_equations[generation_max-1] <- paste(transmission_equations[generation_max], transmission_equations[generation_max-1], sep="+")
  transmission_equations <- transmission_equations[-generation_max]
  
  #combine all state transition equations in propensity vector, a
  a <- c(progression_equations, recovery_equations, transmission_equations)
  
  #(6) return all outputs as a list
  return(list(x0, a, nu, parms, e0, ke))
}

#(2) A function to run a Gillespie Algorithm using the vectors and transmission matrix above
#Using any given 

run.model <- function(model){
  #first, randomly assign e0 exposed individual across the ke exposed states
  #ke and e0 are stored in model[[6]] and model[[5]] respectively
  model[[1]][2:(model[[6]]+1)] <- tabulate(sample(c(1:model[[6]]), length(model[[5]]), replace=TRUE), nbins=model[[6]])
  #next, use modified model to run GillespieSSA
  out <- GillespieSSA::ssa(
    x0 = model[[1]],
    a = model[[2]],
    nu = model[[3]],
    parms = model[[4]],
    tf = 1000000,
    method= GillespieSSA::ssa.d(),
    censusInterval = 0,
    verbose = FALSE
  )
  return(out$data)
}


#(3) A function to run n series of Gillespie Algorithm simulations

run.simulations <- function(
    N=500,
    S=100, 
    e0=1,
    bdd=1/1000, 
    bfd=0, 
    q=1,
    ke=3, 
    mue=5, 
    ki=3, 
    mui=2,
    generation_tracking=FALSE,
    generation_max=10, 
    runs=100){
  
  #Build necessary model components
  model <- build.model(N=N,
                       S=S,
                       e0=e0,
                       bdd=bdd,
                       bfd=bfd,
                       q=q,
                       mue=mue,
                       mui=mui,
                       ke=ke,
                       ki=ki,
                       generation_tracking=generation_tracking,
                       generation_max=generation_max)
  
  #Create list of length sim to store simulation outputs
  simulation_outputs <- as.list(numeric(runs))
  
  for(ii in 1:runs){
    simulation_outputs[[ii]] <- run.model(model)
  }
  
  return(as.list(simulation_outputs))
}

#(4) Three functions for analysing simulation runs
#(and one to check for appropriate generation_max)

#return.end: takes a simulation run (returned as a state by time matrix)
#and extracts the end time (last value in the time column)
return.end <- function(run){
  return(run[nrow(run),1])
}

#return.gen: a function to track the number of transmission generations 
return.gen <- function(run){
  #use colSums to identify all state columns with at least one infection,
  #and extract only the numbers from their names
  #since generation indicators are the only numbers in state names
  #this gives all generations with infections>0
  generations_with_infection <- as.numeric(stringr::str_extract(names((which(colSums(run[,4:(ncol(run)-1)])>0))), "[[:digit:]]+"))
  
  #return the maximum 
  return(max(generations_with_infection))
}

#return.cases: a function to return the total number of cases in an outbreak
return.cases <- function(run){
  #total cases is the number in R at the end minus the number in R at the beginning
  return(run[nrow(run),ncol(run)]-run[1,ncol(run)])
}


#(5) analyse.runs: a function to create summary report from a list of runs
analyse.runs <- function(list){
  #make an output matrix
  summary <- matrix(0, nrow=(length(list)), ncol=3)
  colnames(summary) <- c("Duration", "Generations",
                         "Cases")
  
  #use lapply to analyse each run
  summary[,1] <- unlist(lapply(list, return.end))
  summary[,2] <- unlist(lapply(list, return.gen))
  summary[,3] <- unlist(lapply(list, return.cases))
  return(summary)
}

#(6) run.analysis: a function which takes some inputs as vectors 
#and runs an equal number of simulations for each


run.analysis <- function(
    N=c(100, 200, 500), 
    S=50, 
    e0=1,
    bdd=1/1000,
    bfd=0,
    q=1,
    mue=5,
    mui=5,
    ke=3,
    ki=3, 
    generation_max=10,
    generation_tracking=FALSE,
    runs=10){
  
  #identify which inputs are vectors of length >1
  #and note some values that will be useful later on
  inputs <- list(N, S, e0, bdd, bfd, q, mue, mui, ke, ki)
  
  vector_indices <- which(unlist(lapply(inputs, length))>1)
  vector_lengths <- unlist(lapply(inputs[vector_indices], length))
  
  n_combinations <- prod(vector_lengths)
  n_vectors <- length(vector_lengths)
  
  vector_names <- c("N", "S", "e0", "bdd", "bfd", "q", "mue", "mui", "ke", "ki")[vector_indices]
  
  #list every possible combination of variables
  vector_combinations <- expand.grid(inputs[vector_indices])
  colnames(vector_combinations) <- vector_names
  
  #build an appropriately sized and named list to store raw model output
  results_raw <- vector(mode='list', length=n_combinations)
  list_names <- character(n_vectors)
  for(i in 1:(n_vectors-1)){
    list_names <- paste(list_names, vector_names[i],"=", vector_combinations[,i], "; ", sep="")
  }
  list_names <- paste(list_names, vector_names[n_vectors],"=", vector_combinations[,n_vectors])
  names(results_raw) <- list_names
  
  #build an appropriately sized and named matrix to store  results of analyse.runs
  results_analysis <- matrix(nrow=(runs*n_combinations), ncol=(n_vectors + 3)) #+3 for the output of analyse.runs
  colnames(results_analysis) <- c(vector_names, "Duration", "Generations", "Cases")
  results_analysis <- as.data.frame(results_analysis)
  for(i in 1:ncol(vector_combinations)){
    results_analysis[,i] <- rep(vector_combinations[,i], each=runs)
  }
  
  
  #run simulations for every variable combination and store output
  start <- 1
  end <- runs
  for(i in 1:nrow(vector_combinations)){
    combination_inputs <- inputs
    combination_inputs[vector_indices] <- vector_combinations[i,]
    combination_inputs <- unlist(combination_inputs)
    results_raw[[i]] <- run.simulations(N=combination_inputs[1],
                                        S=combination_inputs[2],
                                        e0=combination_inputs[3],
                                        bdd=combination_inputs[4],
                                        bfd=combination_inputs[5],
                                        q=combination_inputs[6],
                                        mue=combination_inputs[7],
                                        mui=combination_inputs[8],
                                        ke=combination_inputs[9],
                                        ki=combination_inputs[10],
                                        runs=runs,
                                        generation_max=generation_max,
                                        generation_tracking=generation_tracking)
    results_analysis[c(start:end),c((n_vectors+1):(n_vectors+3))] <- analyse.runs(results_raw[[i]])
    print(paste((end*100)/nrow(results_analysis), '%', sep=""))
    start <- start+runs
    end <- end+runs 
  }
  
  #return both raw and analysed results as a list
  return(list(raw = results_raw, analysis=results_analysis))
}

#run.analysis 2: an adapted version of run.analysis 
#for cases where we need to co-vary variables

run.analysis2 <- function(
    N=c(100, 200, 500), 
    S=c(50, 100, 250), 
    e0=1,
    bdd=1/1000,
    bfd=0,
    q=1,
    mue=5,
    mui=5,
    ke=3,
    ki=3, 
    generation_max=10,
    generation_tracking=FALSE,
    runs=100){
  
  #identify which inputs are vectors of length >1
  #and note some values that will be useful later on
  inputs <- list(N, S, e0, bdd, bfd, q, mue, mui, ke, ki)
  
  vector_indices <- which(unlist(lapply(inputs, length))>1)
  vector_lengths <- unlist(lapply(inputs[vector_indices], length))
  if(var(vector_lengths) != 0){
    return(print("All input vectors must be the same length"))
  }
  
  n_combinations <- vector_lengths[1]
  n_vectors <- length(vector_lengths)
  
  vector_names <- c("N", "S", "e0", "bdd", "bfd", "q", "mue", "mui", "ke", "ki")[vector_indices]
  
  #list every possible combination of variables
  vector_combinations <- matrix(unlist(inputs[vector_indices]), ncol=n_vectors, nrow=n_combinations, byrow=FALSE)
  colnames(vector_combinations) <- vector_names
  
  #build an appropriately sized and named list to store raw model output
  results_raw <- vector(mode='list', length=n_combinations)
  list_names <- character(n_vectors)
  for(i in 1:(n_vectors-1)){
    list_names <- paste(list_names, vector_names[i],"=", vector_combinations[,i], "; ", sep="")
  }
  list_names <- paste(list_names, vector_names[n_vectors],"=", vector_combinations[,n_vectors])
  names(results_raw) <- list_names
  
  #build an appropriately sized and named matrix to store  results of analyse.runs
  results_analysis <- matrix(nrow=(runs*n_combinations), ncol=(n_vectors + 3)) #+3 for the output of analyse.runs
  colnames(results_analysis) <- c(vector_names, "Duration", "Generations", "Cases")
  results_analysis <- as.data.frame(results_analysis)
  for(i in 1:ncol(vector_combinations)){
    results_analysis[,i] <- rep(vector_combinations[,i], each=runs)
  }
  
  
  #run simulations for every variable combination and store output
  start <- 1
  end <- runs
  for(i in 1:nrow(vector_combinations)){
    combination_inputs <- inputs
    combination_inputs[vector_indices] <- vector_combinations[i,]
    combination_inputs <- unlist(combination_inputs)
    results_raw[[i]] <- run.simulations(N=combination_inputs[1],
                                        S=combination_inputs[2],
                                        e0=combination_inputs[3],
                                        bdd=combination_inputs[4],
                                        bfd=combination_inputs[5],
                                        q=combination_inputs[6],
                                        mue=combination_inputs[7],
                                        mui=combination_inputs[8],
                                        ke=combination_inputs[9],
                                        ki=combination_inputs[10],
                                        runs=runs,
                                        generation_max=generation_max,
                                        generation_tracking=generation_tracking)
    results_analysis[c(start:end),c((n_vectors+1):(n_vectors+3))] <- analyse.runs(results_raw[[i]])
    print(paste((end*100)/nrow(results_analysis), '%', sep=""))
    start <- start+runs
    end <- end+runs 
  }
  
  #return both raw and analysed results as a list
  return(list(raw = results_raw, analysis=results_analysis))
}

#(7) introduction.risk.time:
#a fuction to bootstrap introduction risk as a function of time

introduction.risk.time <- function(df, resolution=1){
  #identify number variables in data
  variables_index <- c(1:ncol(df))[-which(colnames(df) %in% c("Duration", "Generations", "Cases"))]
  
  variables_names <- colnames(df)[variables_index]
  
  variables_combinations <- as.data.frame(unique(df[,variables_index]))
  colnames(variables_combinations) <- variables_names
  
  #identify longest outbreak duration in dataset
  time_max <- round(max(df$Duration),0)
  
  #list of times to assess introduction risk
  times <- seq(0, time_max, by=resolution)
  
  #create output data frame to store results
  output <- as.data.frame(matrix(nrow=(length(times)*nrow(variables_combinations)),
                                 ncol=(ncol(variables_combinations)+2)))
  colnames(output) <- c(variables_names, 'time', 'p_introduction')
  output$time <- times
  
  count <- 1
  for(i in 1:nrow(variables_combinations)){
    output[c(count:(count+length(times))), c(1:ncol(variables_combinations))] <- variables_combinations[i,]
    count <- count+length(times)
  }
  
  count <- 1
  end <- length(times)
  for(i in 1:nrow(variables_combinations)){
    #isolate the rows in df whose variable columns match the combinations in variables_combinations[i,]
    df_subset <- which(mapply(paste0, df[,variables_index]) == paste0(variables_combinations[i,])) 
    df_subset_duration <- na.omit(df$Duration[df_subset])
    for(j in 1:length(times)){
      output$p_introduction[count] <- (length(which(df_subset_duration > times[j]))/length(df_subset_duration)) 
      count <- count+1
    }
  }
  return(na.omit(output))
}

##(8) run.analysis.ship:
#a function to simulate outnreaks of a given pathogen on ships of different sizes,
#and then bootstrap introduction risk for a given ship journey time

run.analysis.ship <- function(Ship_name, 
                              Journey_time,
                              N_ship,
                              p_susceptible,
                              e0=1,
                              mue=5,
                              mui=5,
                              ke=3,
                              ki=3,
                              q=1,
                              bdd=1/1000,
                              bfd=0,
                              runs){
  
  if(length(unique(N_ship))<length(N_ship)){
    return(print("Function cannot handle duplicate N_ship values; perform separately for identical values of N_ship"))
  }
  
  #calculate S
  S_ship <- round(N_ship*p_susceptible,0)
  while(any(S_ship > N_ship-e0)){
    S_ship[which(S_ship>N_ship-e0)] <- S_ship[which(S>N_ship-e0)]-1
  }
  
  #run analysis
  analysis <- run.analysis2(N=N_ship,
                            S=S_ship,
                            e0=e0,
                            mue=mue,
                            mui=mui,
                            ke=ke,
                            ki=ki,
                            q=q,
                            bdd=bdd,
                            bfd=bfd,
                            runs=runs)
  
  #add ship names and journey times to run.analysis2 output
  analysis$analysis$Ship_name <- rep(Ship_name, times=runs)
  analysis$analysis$Journey_time <- rep(Journey_time, times=runs)
  
  #make output table
  output <- tibble::tibble("Ship_name"=Ship_name, "Journey_time"=Journey_time, "p_introduction"=NA)
  for(i in 1:nrow(output)){
    subset <- dplyr::filter(analysis$analysis, Ship_name==output$Ship_name & Journey_time==output$Journey_time)
    output$p_introduction[i] <- length(which(subset$Duration >= output$Journey_time[i]))/nrow(subset)
  }
  
  return(list(analysis=analysis, ship_risk=output))
}

#(7) check.generation.max(): a function to check whether the value of generation_max is appropriate for a given simulation

check.generation.max <- function(output){
  #(1) extract generation_max from raw data
  generation_max <- max(
    na.omit(
      as.numeric(
        stringr::str_extract(
          colnames(output$raw[[1]][[1]]),
          "[[:digit:]]+")
      )
    )
  )
  print(paste("Generations tracked:", generation_max), quote=FALSE)
  print(paste("Generation distribution in simulation"), quote=FALSE)
  print(table(output$analysis$Generations))
}

#(8) A function to label outbreaks by whether they...
# i. are single-generation
# ii. end before reaching herd immunity
# iii. end after reaching herd immunity

label.outbreaks <- function(df, N, r0){
  df$label <- NA
  df$label[which(df$Generations==1)] <- 1
  df$label[which(df$Generations > 1 & df$Cases/N < max(1-(1/r0),0))] <- 2
  df$label[which(df$Cases/N >= 1-(1/r0))] <- 3
  df$label <- factor(df$label, levels=c(1, 2, 3),
                     labels=c("Single-generation",
                              "Below herd immunity",
                              "At or above herd immunity"))
  return(df)
}

#(9) get.r0: a function to calculate r0 from bdd, bfd, mui, q, and N
get.r0 <- function(bdd, bfd, mui, q, N){
  return(mui*bfd^(1-q)*(bdd*N)^q)
}

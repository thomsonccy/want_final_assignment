# Function: calculate_water_level_reservoir
# Description:
#   Calculates the water level of a reservoir at next time step.
#
# Parameters:
#   state: Numeric. The current water level in the reservoir.
#   area: Numeric. The area of the reservoir.
#   alpha: Numeric. The decay constant affecting the rate of water level change.
#
# Returns:
#   Numeric. The new water level.

calculate_water_level_reservoir <- function(state, area, alpha) {
  return((-alpha * state) / area)
}


# Function: calculate_euler_forward_reservoir
# Description:
#   Calculates the next state of the reservoir using Euler's Forward Method.
#
# Parameters:
#   state: Numeric. The current state of the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   alpha: Numeric. The decay constant.
#
# Returns:
#   Numeric. The new state of the reservoir after applying Euler's Forward Method.

calculate_euler_forward_reservoir <- function(state, dt, area, alpha) {
  forward <- state + dt * calculate_water_level_reservoir(state, area, alpha)
  return(forward)
}


# Function: calculate_heuns_reservoir
# Description:
#   Computes the next state of the reservoir using Heun's Method.
#
# Parameters:
#   state: Numeric. The current state of the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   alpha: Numeric. The decay constant.
#
# Returns:
#   Numeric. The new state of the reservoir after applying Heun's Method.

calculate_heuns_reservoir <- function(state, dt, area, alpha) {
  s_tilde <- state + dt * calculate_water_level_reservoir(state, area, alpha)
  heuns <- state + dt / 2 * (calculate_water_level_reservoir(state, area, alpha) + calculate_water_level_reservoir(s_tilde, area, alpha))
  return(heuns)
}


# Function: calculate_rk4_reservoir
# Description:
#   Calculates the next state of the reservoir using the fourth-order Runge-Kutta method.
#
# Parameters:
#   state: Numeric. The current state of the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   alpha: Numeric. The decay constant.
#
# Returns:
#   Numeric. The new state of the reservoir after applying the fourth-order Runge-Kutta method.

calculate_rk4_reservoir <- function(state, dt, area, alpha) {
  k1 <- dt * calculate_water_level_reservoir(state, area, alpha)
  k2 <- dt * calculate_water_level_reservoir(state + k1 / 5, area, alpha)
  k3 <- dt * calculate_water_level_reservoir(state + 3 / 40 * k1 + 9 / 40 * k2, area, alpha)
  k4 <- dt * calculate_water_level_reservoir(state + 3 / 10 * k1 - 9 / 10 * k2 + 6 / 5 * k3, area, alpha)
  k5 <- dt * calculate_water_level_reservoir(state - 11 / 54 * k1 + 5 / 2 * k2 - 70 / 27 * k3 + 35 / 27 * k4, area, alpha)
  k6 <- dt * calculate_water_level_reservoir(state + 1631 / 55296 * k1 + 175 / 512 * k2 + 575 / 13824 * k3 + 44275 / 110592 * k4 + 253 / 4096 * k5, area, alpha)
  newstate4 <- state + 37 / 378 * k1 + 250 / 621 * k3 + 125 / 594 * k4 + 512 / 1771 * k6
  return(newstate4)
}


# Function: calculate_rk5_reservoir
# Description:
#   Calculates the next state of the reservoir using the fifth-order Runge-Kutta method.
#
# Parameters:
#   state: Numeric. The current state of the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   alpha: Numeric. The decay constant.
#
# Returns:
#   Numeric. The new state of the reservoir after applying the fifth-order Runge-Kutta method.

calculate_rk5_reservoir <- function(state, dt, area, alpha) {
  k1 <- dt * calculate_water_level_reservoir(state, area, alpha)
  k2 <- dt * calculate_water_level_reservoir(state + k1 / 5, area, alpha)
  k3 <- dt * calculate_water_level_reservoir(state + 3 / 40 * k1 + 9 / 40 * k2, area, alpha)
  k4 <- dt * calculate_water_level_reservoir(state + 3 / 10 * k1 - 9 / 10 * k2 + 6 / 5 * k3, area, alpha)
  k5 <- dt * calculate_water_level_reservoir(state - 11 / 54 * k1 + 5 / 2 * k2 - 70 / 27 * k3 + 35 / 27 * k4, area, alpha)
  k6 <- dt * calculate_water_level_reservoir(state + 1631 / 55296 * k1 + 175 / 512 * k2 + 575 / 13824 * k3 + 44275 / 110592 * k4 + 253 / 4096 * k5, area, alpha)
  newstate5 <- state + 2825 / 27648 * k1 + 18575 / 48384 * k3 + 13525 / 55296 * k4 + 277 / 14336 * k5 + 1 / 4 * k6
  return(newstate5)
}


# Function: compare_methods_reservoir
# Description:
#   Compares the results of two numerical methods for calculating the reservoir state.
#
# Parameters:
#   state: Numeric. The current state of the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   alpha: Numeric. The decay constant.
#   method1: Function. The first numerical method for state calculation.
#   method2: Function. The second numerical method for state calculation.
# Returns:
#   Numeric. The absolute difference between the results of the two methods.

compare_methods_reservoir <- function(state, dt, area, alpha, method1, method2) {
  result_method1 <- method1(state, dt, area, alpha)
  result_method2 <- method2(state, dt, area, alpha)
  difference <- abs(result_method1 - result_method2)
  return(difference)
}


# Function: simulate_reservoir
# Description:
#   Simulates the state of a reservoir over time using specified numerical methods and a variable time-stepping approach.
#   It decreases the time step by factor every time the absolute difference of the two methods is larger than the absolute tolerance.
#   It returns the result of method1 if the difference is smaller enough.
#
# Parameters:
#   begin_time: Numeric. The start time of the simulation.
#   end_time: Numeric. The end time of the simulation.
#   dt_start: Numeric. The initial time step size.
#   initial_state: Numeric. The initial state of the reservoir.
#   area: Numeric. The area of the reservoir.
#   alpha: Numeric. The decay constant.
#   factor: Numeric. The factor by which the time step is adjusted.
#   tolerance: Numeric. The tolerance level for the difference between methods.
#   method1: Function. The primary numerical method for state calculation.
#   method2: Function. The secondary numerical method for comparison.
# Returns:
#   List. Contains the simulation time, state, and number of iterations.

simulate_reservoir <- function(begin_time, end_time, dt_start, initial_state, area, alpha, factor, tolerance, method1, method2) {
  # Initialize variables
  time <- begin_time
  result_state <- c(initial_state)
  result_time <- c(time)
  current_state <- initial_state
  iterations <- 0

  # Simulation loop
  while (time < end_time) {
    dt <- dt_start
    while (compare_methods_reservoir(current_state, dt, area, alpha, method1, method2) > tolerance) {
      iterations <- iterations + 1
      dt <- dt * factor
    }
    iterations <- iterations + 1
    current_state <- method1(current_state, dt, area, alpha)
    result_state <- c(result_state, current_state)
    time <- time + dt
    result_time <- c(result_time, time)
  }

  return(list(time = result_time, state = result_state, iterations = iterations))
}

# Function: simulate_reservoir_fixed
#
# Purpose:
#   Simulates the behavior of an emptying reservoir over a specified time period 
#   using fixed time-stepping
#
# Parameters:
#   begin_time    - The starting time of the simulation (numeric).
#   end_time      - The ending time of the simulation (numeric).
#   dt_start      - The initial time step for the simulation (numeric).
#   initial_state - The initial state of the reservoir (numeric).
#   area          - The area parameter of the reservoir (numeric).
#   alpha         - The decay constant (numeric).
#   method1       - Function. The numerical method for state calculation.
#
# Returns:
#   A list containing three elements:
#     - time: A vector of time steps at which the state was calculated (numeric).
#     - state: A vector of states of the reservoir at each time step (numeric).
#     - run: the number of runs 

simulate_reservoir_fixed = function(begin_time, end_time, dt_start, initial_state, area, alpha, factor, tolerance, method1) {
  # Initialize variables
  time = begin_time
  result_state = c(initial_state)
  result_time = c(time)
  current_state = initial_state
  runs = 1
  
  # Simulation loop
  while(time < end_time) {
    dt = dt_start
    current_state = method1(current_state, dt, area, alpha)
    result_state = c(result_state, current_state)
    time = time + dt
    result_time = c(result_time, time)
    runs = runs + 1
  }
  return(list(time = result_time, state = result_state, runs = runs))
}


# Define an error function
error_decrease = function(begin_time, end_time, dt_start, initial_state, A, alpha, factor, tolerance,calculate_heuns_reservoir, calculate_euler_forward_reservoir, simulate_reservoir, simulate_reservoir_fixed){
  
  # -------------------- PART 1 of function: Find RMSE for variable solution ----------------------- #
  # at the defined starting dt value 
  
  ## Run the simulation for the variable solution for the given starting dt value
  simulation_results = simulate_reservoir(begin_time, end_time, dt_start, initial_state, A, alpha, factor, tolerance,calculate_heuns_reservoir, calculate_euler_forward_reservoir)
  
  # Extract results for plotting and further analysis
  result_time = simulation_results$time
  result_state = simulation_results$state
  result_runs = simulation_results$iterations
  
  # Find the number of evaluations
  num_eval_variable = length(result_time)
  
  # Calculate the analytical solution and state at the relevant time-steps
  time_sequence = seq(begin_time,end_time, by=0.1)
  analytic_state=initial_state * exp(-time_sequence*alpha/A)
  analytic_state_at_discrete_time_step=initial_state * exp(-result_time*alpha/A)
  
  # Calculate the error at each time step
  analytical_comparison = abs(analytic_state_at_discrete_time_step - result_state)
  
  # Calculate the root-mean-square deviation (RMSD, or RMSE)
  RMSD_variable = sqrt(sum(analytical_comparison^2)/length(analytical_comparison))
  
  # -------------------- PART 2 of function: Find RMSE for variable solution ----------------------- #
  # at the defined starting dt value 
  
  # Run the simulation for the FIXED time-step - using heun's method
  simulation_results_fixed = simulate_reservoir_fixed(begin_time, end_time, dt_start, initial_state, A, alpha, factor, tolerance,calculate_heuns_reservoir)
  
  # Extract results for plotting and further analysis
  result_time_fixed = simulation_results_fixed$time
  result_state_fixed = simulation_results_fixed$state
  
  # Calculate the analytical solution and state at the relevant time-steps
  analytic_state_at_discrete_time_step_fixed=initial_state * exp(-result_time_fixed*alpha/A)
  analytical_comparison_fixed = abs(analytic_state_at_discrete_time_step_fixed - result_state_fixed)
  
  # Calculate the root-mean-square deviation (RMSD, or RMSE)
  RMSD_fixed = sqrt(sum(analytical_comparison_fixed^2)/length(analytical_comparison_fixed))
  
  
  # --------------- PART 3 of function: Compare RMSD and reduce dt of fixed method ------------------- #
  
  # Empty vectors to store value from the while loop 
  dt_fixed_loop = c()
  evaluations_fixed_loop = c()
  RMSD_fixed_loop = c()
  
  # Set the initial amount of runs to track loops through while loop
  runs_fixed = 1
  
  ### Create while loop to compare RMSE fixed and RMSE variable
  while (RMSD_fixed > RMSD_variable) {
    
    # When condition is fulfilled, reduce the dt_start value 
    dt_start = dt_start * 0.98
    
    # Simulate the fixed stepping with the new dt values 
    simulation_results_fixed = simulate_reservoir_fixed(begin_time, end_time, dt_start, initial_state, A, alpha, factor, tolerance,calculate_heuns_reservoir)
    
    # Extract values for further analysis
    result_time_fixed = simulation_results_fixed$time
    result_state_fixed = simulation_results_fixed$state
    
    # Calculate the RMSE 
    analytic_state_at_discrete_time_step_fixed=initial_state * exp(-result_time_fixed*alpha/A)
    analytical_comparison_fixed = abs(analytic_state_at_discrete_time_step_fixed - result_state_fixed)
    RMSD_fixed = sqrt(sum(analytical_comparison_fixed^2)/length(analytical_comparison_fixed))
    
    # Store values from each iteration into the empty vector
    dt_fixed_loop = cbind(dt_fixed_loop, dt_start)
    evaluations_fixed_loop = cbind(evaluations_fixed_loop,length(result_time_fixed) )
    RMSD_fixed_loop = cbind(RMSD_fixed_loop, RMSD_fixed)
    
    # Keep track of the number of loops through the while loop 
    runs_fixed = runs_fixed + 1
    
  }
  # Store the dt_start that has been adjusted as dt_fixed 
  dt_fixed = dt_start
  
  # Calculate the final number of evaluations 
  evaluations_fixed = length(result_time_fixed)
  
  # Add the last run 
  runs_fixed = runs_fixed + 1
  
  # Return the list of all the relevant variables 
  return(list(RMSD_fixed = RMSD_fixed, RMSD_variable = RMSD_variable, runs_fixed = runs_fixed, dt_fixed = dt_fixed, runs_variable = result_runs, evaluations_variable = num_eval_variable, evaluations_fixed=evaluations_fixed, dt_fixed_loop  = dt_fixed_loop , evaluations_fixed_loop = evaluations_fixed_loop, RMSD_fixed_loop = RMSD_fixed_loop))
  
}
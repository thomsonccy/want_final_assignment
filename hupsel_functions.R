# Function: calculate_water_level_hupsel
# Description:
#   Calculates the next water level in the Hupsel reservoir for a given time, state, and set of reservoir parameters.
#
# Parameters:
#
#   time: Numeric. The current time in the simulation.
#   state: Numeric. The current water level in the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.
# Returns:
#   Numeric. The water level at next time step.

calculate_water_level_hupsel <- function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  if (state <= res_lvl_up) {
    state <- (approximate_Qin(time) / area) - (res_alpha_lw * state) / area
  } else {
    state <- (approximate_Qin(time) - (res_alpha_lw * state) - res_alpha_up * (state - res_lvl_up)) / area
  }
  return(state)
}


# Function: calculate_euler_forward_hupsel
# Description:
#   Calculates the next state of the Hupsel reservoir using Euler's Forward Method.
#
# Parameters:
#   time: Numeric. The current time in the simulation.
#   state: Numeric. The current water level in the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.
#
# Returns:
#   Numeric. The new state of the reservoir after applying Euler's Forward Method.

calculate_euler_forward_hupsel <- function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  forward <- state + dt * calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  return(forward)
}


# Function: calculate_heuns_hupsel
# Description:
#   Computes the next state of the Hupsel reservoir using Heun's Method.
#
# Parameters:
#   time: Numeric. The current time in the simulation.
#   state: Numeric. The current water level in the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.
#
# Returns:
# Numeric. The new state of the reservoir after applying Heun's Method.

calculate_heuns_hupsel <- function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  s_tilde <- state + dt * calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  heuns <- state + dt / 2 * (calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) + calculate_water_level_hupsel(time + dt, s_tilde, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up))
  return(heuns)
}
# Function: calculate_rk4_hupsel
# Description:
#   Calculates the next state of the Hupsel reservoir using the fourth-order Runge-Kutta method.
#
# Parameters:
#   time: Numeric. The current time in the simulation.
#   state: Numeric. The current water level in the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.

# Returns:
#   Numeric. The new state of the reservoir after applying the fourth-order Runge-Kutta method.

calculate_rk4_hupsel <- function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  k1 <- dt * calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k2 <- dt * calculate_water_level_hupsel(time + dt / 5, state + k1 / 5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k3 <- dt * calculate_water_level_hupsel(time + 3 / 10 * dt, state + 3 / 40 * k1 + 9 / 40 * k2, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k4 <- dt * calculate_water_level_hupsel(time + 3 / 5 * dt, state + 3 / 10 * k1 - 9 / 10 * k2 + 6 / 5 * k3, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k5 <- dt * calculate_water_level_hupsel(time + 1 / 1 * dt, state - 11 / 54 * k1 + 5 / 2 * k2 - 70 / 27 * k3 + 35 / 27 * k4, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k6 <- dt * calculate_water_level_hupsel(time + 7 / 8 * dt, state + 1631 / 55296 * k1 + 175 / 512 * k2 + 575 / 13824 * k3 + 44275 / 110592 * k4 + 253 / 4096 * k5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  newstate4 <- state + 37 / 378 * k1 + 250 / 621 * k3 + 125 / 594 * k4 + 512 / 1771 * k6
  return(newstate4)
}


# Function: calculate_rk5_hupsel
# Description:
#   Calculates the next state of the Hupsel reservoir using the fifth-order Runge-Kutta method.
#
# Parameters:
#   time: Numeric. The current time in the simulation.
#   state: Numeric. The current water level in the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.

# Returns:
#   Numeric. The new state of the reservoir after applying the fifth-order Runge-Kutta method.

calculate_rk5_hupsel <- function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  k1 <- dt * calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k2 <- dt * calculate_water_level_hupsel(time + dt / 5, state + k1 / 5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k3 <- dt * calculate_water_level_hupsel(time + 3 / 10 * dt, state + 3 / 40 * k1 + 9 / 40 * k2, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k4 <- dt * calculate_water_level_hupsel(time + 3 / 5 * dt, state + 3 / 10 * k1 - 9 / 10 * k2 + 6 / 5 * k3, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k5 <- dt * calculate_water_level_hupsel(time + 1 / 1 * dt, state - 11 / 54 * k1 + 5 / 2 * k2 - 70 / 27 * k3 + 35 / 27 * k4, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k6 <- dt * calculate_water_level_hupsel(time + 7 / 8 * dt, state + 1631 / 55296 * k1 + 175 / 512 * k2 + 575 / 13824 * k3 + 44275 / 110592 * k4 + 253 / 4096 * k5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  newstate5 <- state + 2825 / 27648 * k1 + 18575 / 48384 * k3 + 13525 / 55296 * k4 + 277 / 14336 * k5 + 1 / 4 * k6
  return(newstate5)
}


# Function: compare_methods_hupsel
# Description:
#   Compares the results of two numerical methods for calculating the state of the Hupsel reservoir.
#
# Parameters:
#   time: Numeric. The current time in the simulation.
#   state: Numeric. The current water level in the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.
#   method1: Function. The primary numerical methods to be compared.
#   method2: Function. The secondary numerical methods to be compared.
#
# Returns:
#   Numeric. The absolute difference between the results of the two methods.

compare_methods_hupsel <- function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up, method1, method2) {
  result_method1 <- method1(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  result_method2 <- method2(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  difference <- abs(result_method1 - result_method2)
  return(difference)
}


# Function: simulate_hupsel
# Description:
#   Simulates the state of the Hupsel reservoir over time using specified numerical methods and a variable time-stepping approach.
#   It decreases the time step by factor every time the absolute difference of the two methods is larger than the absolute tolerance.
#   It returns the result of method1 if the difference is smaller enough.
#
# Parameters:
#   begin_time: Numeric. The start time of the simulation.
#   end_time: Numeric. The end time of the simulation.
#   dt_start: Numeric. The initial time step size.
#   initial_state: Numeric. The initial state of the Hupsel reservoir.
#   area: Numeric. The area of the Hupsel reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.
#   factor: Numeric. The factor by which the time step is adjusted.
#   tolerance: Numeric. The tolerance level for the difference between methods.
#   method1: Function. The primary numerical method for state calculation.
#   method2: Function. The secondary numerical method for comparison.
#
# Returns:
#   List. Contains the simulation time, state, and number of iterations.

simulate_hupsel <- function(begin_time, end_time, dt_start, initial_state, area, res_alpha_lw, res_alpha_up, res_lvl_up, factor, tolerance, method1, method2) {
  # Initialize variables
  time <- begin_time
  result_state <- c(initial_state)
  result_time <- c(time)
  current_state <- initial_state
  result_Qout1 <- c()
  result_Qout2 <- c()
  result_Qin <- c()
  iterations <- 0
  # Simulation loop
  while (time < end_time) {
    dt <- dt_start
    # dt = factor * dt * (tol/abs(sys.fun_CK4(current.state) - sys.fun_CK5(current.state)))^1/5
    while (compare_methods_hupsel(time, current_state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up, method1, method2) > tolerance) {
      iterations <- iterations + 1
      dt <- dt * factor
    }
    iterations <- iterations + 1
    current_state <- method1(time, current_state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
    result_state <- c(result_state, current_state)
    time <- time + dt
    result_time <- c(result_time, time)
  }

  return(list(time = result_time, state = result_state, iterations = iterations))
}


# Function: calculate_water_level_above_upper
# Description:
#   Determines which states are above the upper outlet level in the Hupsel reservoir simulation.
#
# Parameters:
#   res_lvl_up: Numeric. The upper level of the reservoir.
#   simulation_result: List. The results of the reservoir simulation.
#
# Returns:
#   Numeric Vector. Water levels above the upper outlet level.

calculate_water_level_above_upper <- function(res_lvl_up, simulation_result) {
  water_level_above_upper <- c()
  for (i in 1:length(simulation_result$time))
  {
    current_water_level <- simulation_result$state[i]
    if (current_water_level > res_lvl_up) {
      water_level_above_upper <- c(water_level_above_upper, current_water_level)
    } else {
      water_level_above_upper <- c(water_level_above_upper, 0)
    }
  }
  return(water_level_above_upper)
}


# Function: calculate_water_level_fixed_hupsel
# Description:
#   Calculates the next water level in the Hupsel reservoir for a given time, state, and set of reservoir parameters with a fixed time step.
#
# Parameters:
#   time: Numeric. The current time in the simulation.
#   state: Numeric. The current water level in the reservoir.
#   dt: Numeric. Time step size.
#   area: Numeric. The area of the reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.
#
# Returns:
#   Numeric. The water level at next time step.

calculate_water_level_fixed_hupsel <- function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  if (state <= res_lvl_up) {
    state <- state + dt / area * (approximate_Qin(time) - res_alpha_lw * state)
  } else {
    state <- state + dt / area * (approximate_Qin(time) - res_alpha_up * (state - res_lvl_up) - res_alpha_lw * state)
  }
  return(state)
}

# Function: simulate_fixed_hupsel
# Description:
#   Simulates the state of the Hupsel reservoir over time using a fixed time-stepping approach.
#
# Parameters:
#   begin_time: Numeric. The start time of the simulation.
#   end_time: Numeric. The end time of the simulation.
#   dt: Numeric. Time step size.
#   initial_state: Numeric. The initial state of the Hupsel reservoir.
#   area: Numeric. The area of the Hupsel reservoir.
#   res_alpha_lw: Numeric. The lower outlet decay constant.
#   res_alpha_up: Numeric. The upper outlet decay constant.
#   res_lvl_up: Numeric. The height of the upper outlet of the reservoir.
#
# Returns:
#   List. Contains the simulation time and state.

simulate_fixed_hupsel <- function(begin_time, end_time, dt, initial_state, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  # Initialize variables
  time <- begin_time
  result_state <- c(initial_state)
  result_time <- c(time)
  current_state <- initial_state
  result_Qout1 <- c()
  result_Qout2 <- c()
  result_Qin <- c()
  while (time < end_time) {
    current_state <- calculate_water_level_fixed_hupsel(time, current_state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
    result_state <- c(result_state, current_state)
    time <- time + dt
    result_time <- c(result_time, time)
  }
  return(list(time = result_time, state = result_state))
}

calculate_water_level_hupsel = function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up){
  if(state <= res_lvl_up){
    state = (approximate_Qin(time)/area) - (res_alpha_lw*state)/area
  }
  else{
    state = (approximate_Qin(time) - (res_alpha_lw*state) - res_alpha_up*(state-res_lvl_up))/area
  }
  return(state)
}

# Function: calculate_euler_forward_hupsel
#
# Purpose:
#   Calculates the next state of the system using the Euler Forward method.
#
# Parameters:
#   state - The current state of the system (numeric).
#   dt    - The time step to be used for calculation (numeric).
#   A     - The area parameter of the reservoir (numeric).
#   alpha - The decay constant (numeric).
#
# Returns:
#   The next state of the system after applying the Euler Forward method (numeric).

calculate_euler_forward_hupsel = function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  forward = state + dt * calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  return(forward)
}


# Function: calculate_heuns_method
#
# Purpose:
#   Calculates the next state of the system using Heun's method.
#
# Parameters:
#   state - The current state of the system (numeric).
#   dt    - The time step to be used for calculation (numeric).
#   A     - The area parameter of the reservoir (numeric).
#   alpha - The decay constant (numeric).
#
# Returns:
#   The next state of the system after applying Heun's method (numeric).

calculate_heuns_hupsel = function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) {
  s_tilde = state + dt * calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  heuns = state + dt / 2 * (calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) + calculate_water_level_hupsel(time+dt, s_tilde, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)) 
  return(heuns)
}

#system function
calculate_rk4_hupsel = function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
{
  k1 = dt*calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k2 = dt*calculate_water_level_hupsel(time + dt / 5, state+ k1 / 5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) 
  k3 = dt*calculate_water_level_hupsel(time + 3/10*dt,state+ 3/40*k1+9/40*k2, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k4 = dt*calculate_water_level_hupsel(time + 3/5*dt, state+3/10*k1-9/10*k2+6/5*k3, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k5 = dt*calculate_water_level_hupsel(time + 1/1*dt, state-11/54*k1+5/2*k2-70/27*k3+35/27*k4, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k6 = dt*calculate_water_level_hupsel(time + 7/8*dt, state+1631/55296*k1+175/512*k2+575/13824*k3+44275/110592*k4+253/4096*k5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  newstate4 = state + 37/378*k1+250/621*k3+125/594*k4+512/1771*k6
  return(newstate4)
}

calculate_rk5_hupsel = function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
{
  k1 = dt*calculate_water_level_hupsel(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k2 = dt*calculate_water_level_hupsel(time + dt / 5, state+ k1 / 5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up) 
  k3 = dt*calculate_water_level_hupsel(time + 3/10*dt,state+ 3/40*k1+9/40*k2, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k4 = dt*calculate_water_level_hupsel(time + 3/5*dt, state+3/10*k1-9/10*k2+6/5*k3, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k5 = dt*calculate_water_level_hupsel(time + 1/1*dt, state-11/54*k1+5/2*k2-70/27*k3+35/27*k4, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  k6 = dt*calculate_water_level_hupsel(time + 7/8*dt, state+1631/55296*k1+175/512*k2+575/13824*k3+44275/110592*k4+253/4096*k5, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  newstate5 = state + 2825/27648*k1+18575/48384*k3+13525/55296*k4+277/14336*k5+1/4*k6
  return(newstate5)
}


# Function: compare_methods
#
# Purpose:
#   Compares the results of the Euler Forward method and Heun's method to determine
#   the difference between these two approaches for a given state and time step.
#
# Parameters:
#   state - The current state of the system (numeric).
#   dt    - The time step to be used for calculation (numeric).
#   A     - The area parameter of the reservoir (numeric).
#   alpha - The decay constant (numeric).
#
# Returns:
#   The absolute difference between the results of the Euler Forward method
#   and Heun's method (numeric).

# compare_methods = function(state, dt=dt, area=A, alpha=alpha) {
#     euler_forward = calculate_euler_forward_method(state, dt, A, alpha)
#     heuns_method = calculate_heuns_method(state, dt, A, alpha)
#     difference = abs(heuns_method - euler_forward)
#     return(difference)
# }

compare_methods_hupsel = function(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up, method1, method2) {
  result_method1 = method1(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  result_method2 = method2(time, state, dt, area, res_alpha_lw, res_alpha_up, res_lvl_up)
  difference = abs(result_method1 - result_method2)
  return(difference)
}


simulate_hupsel = function(begin_time, end_time, dt_start, initial_state, area, alpha_lower, alpha_upper, level_upper, factor, tolerance, method1, method2) {
  # Initialize variables
  time = begin_time
  result_state = c(initial_state)
  result_time = c(time)
  current_state = initial_state
  result_Qout1    = c()
  result_Qout2    = c()
  result_Qin      = c()
  
  # Simulation loop
  while(time < end_time) {
    dt = dt_start
    #dt = factor * dt * (tol/abs(sys.fun_CK4(current.state) - sys.fun_CK5(current.state)))^1/5
    while (compare_methods_hupsel(time, current_state, dt, area, alpha_lower, alpha_upper, level_upper, method1, method2) > tolerance) {
      iterations = iterations + 1
      dt = dt * factor
    }
    iterations = iterations + 1
    current_state = method1(time, current_state, dt, area, alpha_lower, alpha_upper, level_upper)
    result_state = c(result_state, current_state)
    time = time + dt
    result_time = c(result_time, time)
  }
  
  return(list(time = result_time, state = result_state))
}

##first a function (state.above.up) to determine which states are above the upper outlet
calculate_water_level_above_upper = function(res_lvl_up)
{
  water_level_above_upper = c()
  for (i in 1:length(simulation_hepsel_results$time))
  {
    current_water_level = simulation_hepsel_results$state[i]
    if(current_water_level>res_lvl_up)
    {
      water_level_above_upper=c(water_level_above_upper,current_water_level)
    }else{
      water_level_above_upper=c(water_level_above_upper,0)
    }
  }
  return(water_level_above_upper)
}
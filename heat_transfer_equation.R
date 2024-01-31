# Jan Kalęba
# Równania różniczkowe i różnicowe

# Heat Transfer Equation


# Number of elements
input_value <- readline("Enter n:")
num_elements <- as.numeric(input_value)

# Function to calculate x_i
calc_x_i <- function(i) {
  return(2 * i / num_elements)
}

calc_k <- function(x) {
  if (x > 1) {
    return (2)
  }
  return (1)
}

# Basis function e
basis_func_e <- function(x, i) {
  if (x > calc_x_i(i - 1) && x <= calc_x_i(i)) {
    return (num_elements / 2 * x - i + 1)
  } 
  else if (x > calc_x_i(i) && x < calc_x_i(i + 1)) {
    return (-num_elements / 2 * x + i + 1)
  } 
  
  return (0)
}

# Derivative of the basis function e
basis_func_e_prime <- function(x, i) {
  left_bound = calc_x_i(i - 1)
  right_bound = calc_x_i(i)
  next_bound = calc_x_i(i + 1)
  
  (x > left_bound && x <= right_bound) * (num_elements / 2) +
    (x > right_bound && x < next_bound) * (-num_elements / 2)
}


# Gauss-Legendre Quadrature for integration
# 2 point quadrature - roots of Lagendre polynomial
# of degree 2 are (-1 / sqrt(3)) and (1 / sqrt(3))
gauss_legendre_integration <- function(f, a, b) {
  root1 = (1 / sqrt(3))
  root2 = (-1 / sqrt(3))
  return(
    (b - a) / 2 * (
      # Scaling & shifting for the interval [a, b]
      f((b - a) / 2 * root1 + (a + b) / 2) +
        f((b - a) / 2 * root2 + (a + b) / 2)
    )
  )
}

# Function for u' * v'
u_prime_v_prime <- function(i, j) {
  return (function(x) {
    return(basis_func_e_prime(x, i) * basis_func_e_prime(x, j))
  })
}

# Matrix B calculation
calculate_B_matrix <- function(i, j) {
  
  start_point = max(0, calc_x_i(i - 1), calc_x_i(j - 1))
  end_point = min(calc_x_i(i + 1), calc_x_i(j + 1))
  
  integrand <- function(x) {
    k_value = calc_k(x)
    return(k_value * basis_func_e_prime(x, i) * basis_func_e_prime(x, j))
  }
  
  return(
    basis_func_e(0, i) * basis_func_e(0, j) - gauss_legendre_integration(integrand, start_point, end_point)
  )
}


# Linear form for the right-hand side
linear_form_L <- function(i) {
  return(20 * basis_func_e(0, i))
}

# Constructs the system matrix for the equation by calculating each element based on the basis functions.
build_system_matrix <- function() {
  system_matrix <- matrix(0, nrow = num_elements, ncol = num_elements)
  for (i in 1:num_elements) {
    for (j in 1:num_elements) {
      system_matrix[i, j] <- calculate_B_matrix(j - 1, i - 1)
    }
  }
  return(system_matrix)
}

# Creates the right-hand side vector for the equation system using the linear form function.3
construct_rhs_vector <- function() {
  rhs_vector <- sapply(0:(num_elements - 1), linear_form_L)
  return(rhs_vector)
}

# Solves the linear system of equations using the system matrix and right-hand side vector.
solve_linear_system <- function(system_matrix, rhs_vector) {
  solution <- solve(system_matrix, rhs_vector)
  return(solution)
}

# Builds the combined function representing the solution, using the solution vector and basis functions.
assemble_combined_function <- function(solution_vector) {
  combined_func <- function(x) {
    result <- 0
    for (i in 1:num_elements) {
      result <- result + solution_vector[i] * basis_func_e(x, i - 1)
    }
    return(result)
  }
  return(combined_func)
}


# Main function to find the solution of the equation
find_solution <- function() {
  system_matrix <- build_system_matrix()
  rhs_vector <- construct_rhs_vector()
  
  solution_vector <- solve_linear_system(system_matrix, rhs_vector)
  
  combined_function <- assemble_combined_function(solution_vector)
  
  return(combined_function)
}


# Plotting the basis functions
plot_basis_functions <- function() {
  dev.new()
  
  plot(seq(0, 2, 1 / (100 * num_elements)), mapply(basis_func_e, seq(0, 2, 1 / (100 * num_elements)), 1), 
       main = 'Basis Functions',
       xlab = '',
       ylab = '',
       type = 'l')
  
  for (i in 1:(num_elements - 1)) {
    lines(seq(0, 2, 1 / (10 * num_elements)), mapply(basis_func_e, seq(0, 2, 1 / (10 * num_elements)), i))
  }
}


# Function to plot the solution of the equation
plot_heat_transfer_solution <- function() {
  dev.new()
  
  heat_transfer_solution <- find_solution()
  plot(seq(0, 2, 1 / (100 * num_elements)), 
       mapply(heat_transfer_solution, 
              seq(0, 2, 1 / (100 * num_elements))),
       main = 'Heat Transfer Equation Solution',
       xlab = ' ',
       ylab = ' ',
       type = 'l')
}

# Execute the plot function
plot_basis_functions()
plot_heat_transfer_solution()
    
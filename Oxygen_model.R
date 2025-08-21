###### steady state including prorosity 

# Parameters salp
r <- 0.6  # Particle radius in cm
D_O2 <- 2e-5  # Diffusivity of oxygen in seawater in cm^2/s
Resp <- 0.04  # µmol O2 day-1 (Respiration rate per day)
vol <- 0.0005 # cm3
phi <- 0.9  # Porosity of the particle (fraction)
R <- Resp / vol / 24 / 60 / 60   # Adjust respiration rate for porosity
C_ext <- 120  # Oxygen concentration outside the particle in mol/cm^3
#v <- 0.1  # Sinking velocity in cm/s

# Parameters copepod
r <- 0.01  # Particle radius in cm
D_O2 <- 2e-5  # Diffusivity of oxygen in seawater in cm^2/s
Resp <- 0.01  # µmol O2 day-1 (Respiration rate per day)
vol <- 0.000005 # cm3
phi <- 0.75  # Porosity of the particle (fraction)
R <- Resp / vol / 24 / 60 / 60   # Adjust respiration rate for porosity
C_ext <- 120  # Oxygen concentration outside the particle in mol/cm^3
#v <- 0.1  # Sinking velocity in cm/s

# Spatial grid
dr <- 0.001  # Radial grid spacing in cm
r_grid <- seq(0, r, by=dr)
n <- length(r_grid)

# Function to calculate steady-state concentration with porosity affecting diffusion
steady_state_concentration <- function(D_O2, R, phi, dr, n, C_ext) {
  # Adjust diffusion by porosity
  D_eff <- D_O2 * phi  # Effective diffusion coefficient considering porosity
  
  # Set up matrix for second derivative (finite difference method)
  A <- matrix(0, n, n)
  b <- rep(0, n)
  
  # Boundary conditions
  A[1, 1] <- 1  # C(0) is set to be the internal concentration (zero derivative at the center)
  b[1] <- C_ext  # C(r_max) is set to external concentration
  
  for (i in 2:(n-1)) {
    A[i, i-1] <- D_eff / dr^2
    A[i, i] <- -2 * D_eff / dr^2 - R
    A[i, i+1] <- D_eff / dr^2
  }
  
  # Boundary at the outermost point
  A[n, n] <- 1
  b[n] <- C_ext
  
  # Solve the system of equations
  C_new <- solve(A, b)
  return(C_new)
}

# Find the steady-state concentration profile with porosity affecting both diffusion and respiration
C_steady_state <- steady_state_concentration(D_O2, R, phi, dr, n, C_ext)

# Plot the oxygen concentration gradient inside the particle
plot(r_grid, C_steady_state, type="l", col="blue", xlab="Radius (cm)", ylab="Oxygen Concentration (mol/cm^3)", main="Steady-State Oxygen Concentration with Porosity Impact on Diffusion")



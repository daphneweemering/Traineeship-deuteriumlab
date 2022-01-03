#### Traineeship UMC Utrecht 

#### Daphne Weemering, 3239480, MSBBSS Utrecht University.
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
library(deSolve)
library(boot)

# Model p = d homogeneous [proliferation rate equal to death rate of cells].
model1<- function(t, state, pars){   # A function with as list of arguments: time,
  # state (L) and parameters (p). State indicates how many dL we have (1 for 
  # homogeneous pop, more for heterogeneous pop).
  with(as.list(c(state, pars)),{   
    dL <- ifelse(t <= tau, p-p*L, -p*L)
    return(list(c(dL)))
  })
}

# Set a seed for reproducibility
set.seed(171017)

# Time we want to simulate
time.sim <- c(7, 21, 35, 42)
ntime <- time.sim

# Number of repetitions to be executed
n <- 100
nboot <- 200

# Make some storage 
store <- matrix(0, n, 1)
store_boot <- matrix(0, nboot, 1)
cov <- matrix(0, n, 1)

# Loop
for (i in 1:n){
  
  # 'ode' solves differential equations
  out <- ode(y = c(L = 0), times = ntime, func = model1, parms = c(p = 0.1, tau = 35))
  
  data.simulated <- out[, 2]  
  
  # Add some noise 
  noise <- rnorm(ntime, 0, 0.001) 
  data.sim.noise <- data.simulated + noise
  
  # Estimate the parameter (p) with LS approach
  cost.LS <- function(psi){
    tim <- ntime
    tau2 <- 35
    
    out <- ode(y = c(L = 0), times = tim, func = model1, parms = c(p = exp(psi[1]), tau = tau2))
    data_pred <- out[, 2]

    out1 <- ifelse(data_pred < 0, 0, asin(sqrt(data_pred)))
    out_cost <- (data.fit - out1) 
    
    SSR <- sum(out_cost^2)
    return(SSR)
  }
  
  # Transform the data
  data.fit <- asin(sqrt(data.sim.noise)) 
  
  # nlminb optimizes the cost function. 
  estim <- nlminb(c(-5.11), cost.LS, upper = c(0), lower = c(-15))  
  
  # Estimated p
  prm.LS <- exp(estim$par)  
  
  # Store the p values
  store[i,] <- prm.LS
  
  tau2 <- 35
  
  # Get the residuals 
  pred <- ode(y = c(L = 0), times = ntime, func = model1, parms =c(p = prm.LS[1], tau = tau2))
  res <- data.sim.noise - pred[,2]
  
  ## Bootstrap the residuals
  for (a in 1:nboot){
    
    # Randomly assign the residuals (random order)
    random <- sample(res)
    
    # Add re-assigned residuals to the simulated data
    simdat_res <- data.sim.noise + random
    
    # Estimate the parameters 
    data.fit <- asin(sqrt(simdat_res)) # Transform the data
    estim <- nlminb(c(-5.11), cost.LS, upper = c(0), lower = c(-15)) # Optimize 
                                      # the cost function.
    prm.LS <- exp(estim$par)
    
    # Store the parameters 
    store_boot[a, ] <- prm.LS
  }
  
  # Take the 2.5 and 97.5 percentile to obtain the 95% confidence interval
  CI <- unname(quantile(store_boot, probs = c(0.025, 0.975)))
  
  # Separate them
  CI_lower <- CI[[1]]; CI_upper <- CI[[2]]
  
  # Indicate which p values fall within the confidence interval 
  cov[i,] <- ifelse(store[i] < CI_lower || store[i] > CI_upper, 0, 1)
}

# Calculate the coverage 
coverage <- sum(cov) / n * 100

# Take the mean of all the p (proliferation) values
L_mean <- mean(store_boot[, 1])

# And compute the standard deviation of all the p (proliferation) values
L_sd <- sd(store_boot[, 1])

# Compute Relative Standard Deviation (RSD)
RSD <- (L_sd / L_mean) * 100

# Compute the bias
bias <- store - 0.1  # 0.1 is the value (the real value) for the proliferation rate
bias_mean <- mean(bias)

# Compute Mean Squared Error (MSE)
MSE <- sum(store - 0.1)^2 * 1/n  # 0.1 is the value (the real value) for the 
                                 # proliferation rate


## Make a plot with simulated data and estimated fit
xfit <- seq(0, 130)
tau2 <- 63
out.LS <- ode(y = c(L = 0), times = xfit, func = model1, parms = c(p = prm.LS[1], tau = tau2))
yfit <- out.LS[, 2]
plot(time.sim, data.sim.noise)
lines(xfit, yfit)





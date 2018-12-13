Sys.setenv(RETICULATE_PYTHON="~/anaconda3/bin/python3.6")
library(devtools)
library(gurobi)
library(reticulate)
library(Rglpk)
library(microbenchmark)
load_all("detrendr")
library(osqp)

i <- 10
n <- 5000
simDesign <- "peaks"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
tau <- c(.05, 0.1)

y <- df$y
k <- 3
quant_loss <- function(u, tau) { 0.5 * abs(u) + (tau - 0.5) * u }

theta <- Variable(length(y))
D <- get_Dk(length(y), 3)
lambda <- n
obj <- sum(quant_loss(y - theta, t = tau)) + lambda * cvxr_norm(D%*%theta, 1)
prob <- Problem(Minimize(obj))

prob_data <- get_problem_data(prob, solver = "ECOS")

solver_output <- ECOSolveR::ECOS_csolve(c = prob_data[["c"]],
                                        G = prob_data[["G"]],
                                        h = prob_data[["h"]], 
                                        dims = prob_data[["dims"]])

direct_soln <- unpack_results(prob, "ECOS", solver_output)
thetaOut2 <- direct_soln$getValue(theta)


model <- list()
model$obj <- prob_data$c
model$A <- prob_data$G
model$rhs <- prob_data$h
model$sense <- "<"
model$lb <- rep(-1e3, length(prob_data$c))
params <- list(OutputFlag=0)
result <- gurobi(model, params)
theta_g <- result$x[1:n]

result <- solve(prob, solver = "MOSEK", verbose = TRUE)
thetaOut <- result$getValue(theta)[1:n]

plot(y~x, df)
lines(thetaOut~df$x, col="blue")
lines(theta_g~df$x, col="red")

quant_loss <- function(u, tau) { 0.5 * abs(u) + (tau - 0.5) * u }

get_trend <- function(y, tau, lambda, k){
  theta <- Variable(length(y))
  D <- get_Dk(length(y), k)
  obj <- sum(quant_loss(y - theta, t = tau)) + 
    lambda * cvxr_norm(D%*%theta, 1)
  prob <- Problem(Minimize(obj))
  prob_data <- get_problem_data(prob, solver = "ECOS")
  model <- list()
  model$obj <- prob_data$c
  model$A <- prob_data$G
  model$rhs <- prob_data$h
  model$sense <- rep("<", nrow(model$A))
  model$lb <- rep(-1e2, length(prob_data$c))
  params <- list(OutputFlag=0)
  result <- gurobi(model, params)
  return(result$x[1:n])
}

get_trend_mosek <- function(y, tau, lambda, k){
  theta <- Variable(length(y))
  D <- get_Dk(length(y), k)
  obj <- sum(quant_loss(y - theta, t = tau)) + 
    lambda * cvxr_norm(D%*%theta, 1)
  prob <- Problem(Minimize(obj))
  prob_data <- get_problem_data(prob, solver = "ECOS")
  prob <-  mosek_lptoprob(f = prob_data$c, 
                   A = prob_data$G, 
                   b = prob_data$h, 
                   Aeq = NA, 
                   beq = NA, 
                   lb = rep(-100, length(prob_data$c)), 
                   ub = rep(100, length(prob_data$c)))
  result <- mosek(prob, opts = list(verbose = 0))
  return(result$sol$bas$xx[1:n])
}


gurobi_trend2 <- function(y, tau, lambda, k){
  
  tau <- sort(tau)
  D <- get_Dk(length(y), k)
  n <- length(y)
  m <- nrow(D)
  np <- 2*n + 2*m
  nT <- length(tau)
  missInd <- which(is.na(y))
  y[missInd] <- 0
  model <- list()
  
  for (i in 1:nT){
    model$obj <- c(model$obj,
                   rep(tau[i], n), rep((1-tau[i]), n), rep(lambda[i], 2*m))
    model$obj[missInd] <- 0
    model$obj[missInd + n] <- 0
    model$rhs <- c(model$rhs, as.numeric(D%*%y))
  }
  
  model$rhs <- c(model$rhs, rep(0, n*(nT-1)))
  
  # Constraint Matrix
  if (length(tau) == 1){
    model$A  <- cbind(D, -D, Diagonal(m), -Diagonal(m))
  } else {
    model$A <- Matrix(0, nrow =  m*nT + n*(nT-1), ncol= np*nT, sparse=TRUE)
    for (i in 1:nT){
      
      # D%*%theta = eta constraint
      model$A[(1+m*(i-1)):(m*i), (1+np*(i-1)):(np*i)] <-
        cbind(D, -D, Diagonal(m), -Diagonal(m))
      
      # Non-crossing theta(tau) constrains
      if (i < nT){
        model$A[(m*nT+1+n*(i-1)):(m*nT+n*(i)), (1+(i-1)*np):(2*n + (i-1)*np)] <-
          cbind(Diagonal(n), -Diagonal(n))
        
        model$A[(m*nT+1+n*(i-1)):(m*nT+n*(i)),
                (1+i*np):(2*n + i*np)] <-
          cbind(-Diagonal(n), Diagonal(n))
      }
    }
  }
  
  # Constraint Types
  model$sense <- c(rep('=', m*nT), rep(">", n*(nT-1)))
  model$modelsense <- "min"
  
  # Type of variables (binary)
  params <- list(OutputFlag=0)
  result <- gurobi(model, params)
  
  theta <- matrix(0, nrow=n, ncol = nT)
  for (i in 1:nT){
    theta[,i] <- y - result$x[(1+np*(i-1)):(n+np*(i-1))] +
      result$x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }
  return(theta)
}

theta <- get_trend(df$y, tau, lambda=n, k=3)

model_list <- get_model_data(df$y, tau, lambda=n, k=3)

theta <- solve_model(model_list, df$y, "gurobi")

result <- get_lambda_trend(df$y, tau, k=3,
                         lambdaSeq = n^seq(0, 1.5, length.out=20), 
                         df_tol = 1e-9, 
                         gamma = 1,
                         plot_lambda = TRUE, 
                         solver = NULL, 
                         criteria = "SIC")
result$lambda
plot(result$theta, type="l")


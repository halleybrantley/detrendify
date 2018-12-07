
library(devtools)
load_all("detrendr")
library(gurobi)
library(CVXR)
library(reticulate)
installed_solvers()

load_all("detrendr")
use_python("~/anaconda3/bin/python3.6")
repl_python()


i <- 10
n <- 500
simDesign <- "gaus"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
tau <- .1

y <- df$y
k <- 3
quant_loss <- function(u, tau) { 0.5 * abs(u) + (tau - 0.5) * u }
theta <- Variable(n)
D <- get_Dk(length(y), 3)
lambda <- 1000
obj <- sum(quant_loss(y - theta, t = tau)) + lambda * cvxr_norm(D%*%theta, 1)
prob <- Problem(Minimize(obj))

prob_data <- get_problem_data(prob, solver = "SCS")

solver_output <- scs::scs(obj = prob_data[["c"]],
                                        cone = prob_data[["dims"]],
                                        A = prob_data[["A"]],
                                        b = prob_data[["b"]])

direct_soln <- unpack_results(prob, "SCS", solver_output)

model <- list()
model$obj <- prob_data$c
model$A <- prob_data$A
model$rhs <- prob_data$b
params <- list(OutputFlag=0)
result <- gurobi(model, params)


result <- solve(prob, solver = "GUROBI", verbose = TRUE)
thetaOut <- result$getValue(theta)
plot(y~x, df)
lines(theta~df$x, col="blue")

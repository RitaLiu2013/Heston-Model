# Defining the characteristic function.
characteristicFunctionHeston <- function(u, S=100, v=0.06, r=0.05, kappa=1, theta=0.06, sigma=0.3, rho=-0.5, lambda=0.01, tau=1, j=1){
  a <- kappa*theta
  i <- complex(imaginary =1)
  x <- log(S)
  if (j==1){
    u1 <- 1/2
    b1 <- kappa + lambda - rho*sigma
    d1 <- sqrt((rho*sigma*u*i-b1)^2 - sigma^2*(2*u1*u*i - u^2))
    g1 <- (b1 - rho*sigma*u*i + d1)/(b1-rho*sigma*u*i-d1)
    C1 <- r*u*i*tau + a/sigma^2 * ((b1-rho*sigma*u*i+d1)*tau -
                                     2* log((1-g1*exp(d1*tau))/(1-g1)))
    D1 <- (b1 - rho*sigma*u*i + d1)/(sigma^2) * (1-exp(d1*tau))/(1-g1*exp(d1*tau))
    Psi1 <- exp(C1+D1*v+i*u*x)
    
    return(Psi1)
  } else if (j==2){
    u2 <- -1/2
    b2 <- kappa + lambda
    d2 <- sqrt((rho*sigma*u*i-b2)^2 - sigma^2*(2*u2*u*i - u^2))
    g2 <- (b2- rho*sigma*u*i + d2)/(b2-rho*sigma*u*i-d2)
    C2 <- r*u*i*tau + a/sigma^2 * ((b2-rho*sigma*u*i+d2)*tau -
                                     2* log((1-g2*exp(d2*tau))/(1-g2)))
    D2 <- (b2 - rho*sigma*u*i + d2)/(sigma^2) * (1-exp(d2*tau))/(1-g2*exp(d2*tau))
    Psi2 <- exp(C2+D2*v+i*u*x)
    
    return(Psi2)
  }
  
}

# The function characteristicFunctionHeston is given an abbreviating alias, Psi.
Psi <- characteristicFunctionHeston


#source("characteristicFunctionHeston.R")
# Generating u values within the interval [-20, 20].
u_values <- seq(-20, 20, length.out = 400)

# Computing Psi1 for each generated u value.
Psi1_values <- Psi(u_values)

# Plotting Psi1.
plot(u_values, Re(Psi1_values), type = 'l', col = 'blue', ylim = range(c(Re(Psi1_values),
                                                                         Im(Psi1_values))),xlab = expression(u), ylab = expression(Psi[1](u)),
     main = expression(Psi[1](u) ~ "for" ~ u %in% "[-20, 20]"))
lines(u_values, Im(Psi1_values), col = 'red')
legend("topright", legend = c("Real Part", "Imaginary Part"), col = c("blue", "red"), lty = 1)

# Computing Psi2 for each generated u value.
Psi2_values <- Psi(u_values,j=2)

# Plotting Psi2.
plot(u_values, Re(Psi2_values), type = 'l', col = 'blue', ylim = range(c(Re(Psi2_values),
                                                                         Im(Psi2_values))), xlab = expression(u), ylab = expression(Psi[2](u)),
     main = expression(Psi[2](u) ~ "for" ~ u %in% "[-20, 20]"))
lines(u_values, Im(Psi2_values), col = 'red')
legend("topright", legend = c("Real Part", "Imaginary Part"), col = c("blue", "red"), lty = 1)



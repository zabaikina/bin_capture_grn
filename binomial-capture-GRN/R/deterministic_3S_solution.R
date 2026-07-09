

repressilator <- function(t, y, p) {
  with(as.list(c(y, p)), {
    P1 <- y3^n13          
    P2 <- y1^n21          
    P3 <- y2^n32          
    
    a1 <- (rho_u1*K1 + rho_b1*P1)/(K1 + P1)
    a2 <- (rho_u2*K2 + rho_b2*P2)/(K2 + P2)
    a3 <- (rho_u3*K3 + rho_b3*P3)/(K3 + P3)
    
    dy1 <- beta1 * a1 - k1 * y1
    dy2 <- beta2 * a2 - k2 * y2
    dy3 <- beta3 * a3 - k3 * y3
    list(c(dy1, dy2, dy3))
  })
}

# symmetric parameter set
p <- c(
  beta1=10, beta2=10, beta3=10,
  k1=1, k2=1, k3=1,
  K1=9*1e4, K2=9*1e4, K3=9*1e4,
  rho_u1=140, rho_u2=140, rho_u3=140,
  rho_b1=2.9, rho_b2=2.9, rho_b3=2.9,
  n13=3, n21=3, n32=3
)

y0    <- c(y1=95, y2=75, y3=100)        
times <- seq(0, 1000, by=0.05)

sol <- ode(y=y0, times=times, func=repressilator, parms=p,
           rtol=1e-9, atol=1e-9, method="lsoda")

matplot(sol[, "time"], sol[, c("y1","y2","y3")], type="l", lty=1, lwd=1.1,
        xlab="t", ylab="y", 
        xlim = c(500, 510)
        )


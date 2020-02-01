library(tidyverse)
library(deSolve)
library(ggplot2)
library(graphics)
# parameter values derived from Holt and Polis 1997 in American Naturalist
parameters = c(r = 1, 
               K = 10,
               a1 = 1, a2 = 0.3,
               b1 = 1, b2 = 0.1,
               m1 = 0.1, m2 = 0.5,
               IGP = 1, f12 = 1)

init = c(R = 10, C = 3, P = 3)

times = seq(0, 200, by = 0.01)

IGP_div<-function(t, init, parameters) {
  with (as.list(c(init, parameters)), {
    dR = r*R*(1-(R)/K) - a1*R*C - a2*R*P
    dC = b1*a1*R*C - m1*C - IGP*C*P
    dP = b2*a2*R*P - m2*P + f12*IGP*C*P
    
    list(c(dR, dC, dP))
  }
  )
}

out = ode(y = init, times = times, func = IGP_div, parms = parameters)
realts = as.data.frame(out[seq(1,20000, by=100),])
#head(realts)
realts %>%
  gather(key = "Consumer", value = "Density", -c(time, R)) %>%
  ggplot(aes(x = time, y = Density, group = Consumer, colour = Consumer)) +
  geom_line()


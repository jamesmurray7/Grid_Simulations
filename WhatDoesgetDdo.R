# Working out what joineR:::getD does
joineR:::getD
res <- 3
dt <- 0.001
gridt <- seq(0,5,dt)[-1]
joineR:::getD(3,gridt) %>% dim # 3 x 5000 matrix 

jmD.1 <- function(q, arg){
  D <- matrix(0, q, length(arg))
  D
}
jmD.1(3, gridt) %>% dim # Initialises a 3 x 5000 matrix of zeroes (-> empty matrix)

jmD.2 <- function(q, arg){
  D <- matrix(0,q,length(arg))
  for(i in 1:q){ # Over each row...
    D[i,] <- arg^(i-1) # Give values of time^0/time^1/time^2 as necessary
  }
  D
}
jmD.2(3, gridt)

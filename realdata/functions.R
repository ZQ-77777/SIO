# Helper functions for `SIO.R`.
#  - `polynomial_basis1`--`polynomial_basis2` define the polynomial basis for inputs
#       with 1--2 arguments, respectively, and are used to construct the proposed estimator.
#  - The remaining functions are standard utilities for data generation and standardization.


polynomial_basis2 <- function(x1,x2,x3){
  cbind(1,
        x1,
        x3,
        x1^2, 
        x3^2, 
        x1^3,
        x3^3,
        x2[,1], 
        x2[,2], 
        x2[,3], 
        x2[,4], 
        x2[,5], 
        x2[,6],
        x2[,1]^2, 
        x2[,2]^2, 
        x2[,5]^2, 
        x2[,1]^3, 
        x2[,2]^3, 
        x2[,5]^3,
         deparse.level = 0)
}




polynomial_basis1 <- function(x1,x2){
  cbind(1,
        x2,
        x2^2,
        x2^3,
        x1[,1], 
        x1[,2], 
        x1[,3], 
        x1[,4], 
        x1[,5], 
        x1[,6],
        x1[,1]^2, 
        x1[,2]^2, 
        x1[,5]^2, 
        x1[,1]^3, 
        x1[,2]^3, 
        x1[,5]^3,  deparse.level = 0)
}


ind <- function(x){
  ifelse(x > 0, 1, 0)
}


normalize <- function(x){
  (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))
}

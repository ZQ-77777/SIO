# Helper functions for `1.2Estimation.R`.
#  - `polynomial_basis1`--`polynomial_basis4` define the polynomial basis for inputs
#       with 1--4 arguments, respectively, and are used to construct the proposed estimator.
#  - The remaining functions are standard utilities for data generation and standardization.

polynomial_basis4 <- function(x1,x2,x3,x4){
  cbind(1, 
        x1,
        x2, 
        x3[,1], 
        x3[,2],
        x4[,1], 
        x4[,2],
        x1^2, 
        x2^2, 
        x3[,1]^2, 
        x3[,2]^2,
        x4[,1]^2, 
        x1^3, 
        x2^3,
        x3[,1]^3,
        x3[,2]^3,
        x4[,1]^3,
        deparse.level = 0)
}


polynomial_basis3 <- function(x1,x2,x3){
  cbind(1, 
        x1, 
        x2, 
        x3[,1], 
        x3[,2], 
        x3[,3],
        x1^2,
        x2^2, 
        x3[,1]^2, 
        x3[,2]^2, 
        x1^3, 
        x2^3, 
        x3[,1]^3,
        x3[,2]^3,
        deparse.level = 0)
}



polynomial_basis2 <- function(x1,x3){
  cbind(1, x1,
        x3[,1],
        x3[,2],
        x3[,3],   
        x1^2, 
        x3[,1]^2, 
        x3[,2]^2, 
        x1^3, 
        x3[,1]^3,
        x3[,2]^3,
        deparse.level = 0)
}


polynomial_basis1 <- function(x3){
  cbind(1,x3[,1], 
        x3[,2], 
        x3[,3],  
        x3[,1]^2, 
        x3[,2]^2,
        x3[,1]^3,
        x3[,2]^3,
        deparse.level = 0)
}

qr_basis1 <- function(x3){
  cbind(1,x3[,1], x3[,1]^2,
        x3[,2], x3[,2]^2,
        x3[,3],  
        deparse.level = 0)
}


ind <- function(x){
  ifelse(x > 0, 1, 0)
}


expit <- function(x){exp(x)/(1+exp(x))}
iexpit <- function(x){1+exp(-x)}
iplus <- function(x){ifelse(x>=0,1,0)} 
normalize <- function(x){(x-min(x))/(max(x)-min(x))}
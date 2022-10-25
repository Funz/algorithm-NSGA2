## This file should provide following objects, when loaded:
# f : function
# input.f : list of input dimensions, contains list of properties like lower & upper bounds of each dimensions
# output.f : list of output dimensions
# *.f : list of math properties. To be compared with algorithm results
# [print.f] : method to print/plot the function for information

f <- function(x) {
	x1 <- x[1]*15-5   
	x2 <- x[2]*15     
    c(
        (x2 - 5/(4*pi^2)*(x1^2) + 5/pi*x1 - 6)^2 + 10*(1 - 1/(8*pi))*cos(x1) + 10
    , x[2]-x[1]) # : constraint x2>x1
}
input.f = list(
    x1=list(min=0,max=1),
    x2=list(min=0,max=1)
)

output.f = c("branin","x2_sup_x1")
argmin.f_constr = c(0.1238946, 0.8166644)
min.f = 0.3978874

library(testthat)
if (!isTRUE(test_that("f(armgin.f) == f.min",{expect_equal(f(argmin.f_constr)[1],min.f,tolerance = .0001)}))) quit(status=1)

test = function(algorithm_file) {
    results = run.algorithm(algorithm_file, options=list(generations='10', constraints_index='2'),fun=list(input=input.f,output=output.f,fun=f))
    if (!isTRUE(test_that("branin min",{expect_equal(min(as.numeric(results$min)),min.f,tolerance = .01)}))) quit(status=1)
}


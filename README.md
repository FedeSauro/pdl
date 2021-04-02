# pdl
Code polynomial distribuited lag model with constraints

`pdl` is an R package implementing functionalities to estimate a polynomial distribuited lag model with constraints.

You can estimate three classes of models: 
  - polynomial distribuited lag model (parameter: degree of the polynomial - pdl_L)
  - polynomial distribuited lag model subjected to a a bound end (parameter: degree 1 of the polynomial - pdl_L1)
  - polynomial distribuited lag model subjected to two bound extremes (parameter: degree 2 of the polynomial - pdl_L2)

Each class has a function to:
  - create the polynomial matrix
  - estimate the model
  - plot the distribuited lags
  - estimate the confident intervals

----------------------------------------------------------------------------------------------------------------------

R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `pdl` package. R can be downloaded from https://www.r-project.org/.

To install the `pdl` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("FedeSauro/pdl")
```

For any request or feedback, please write to <federicasg94@gmail.com> (Federica Sauro Graziano)

Below, you find some examples of use of the package.
----------------------------------------------------------------------------------------------------------------------

Load simulated data
```
library(pdl)

# load simulated data
data(datasim)
summary(datasim)
```

```
# 2nd order polynomial lag with endpoint constraints for X1
m1 <- pdl_L2(data=datasim, y.name="y", x.name="X1", a=5, b=15)
summary(m1)     ## summary of coefficients beta
summary.lm(m1)  ## summary of parameters gamma
plot(m1, main=expression(paste(X[1])))

# 2nd order polynomial lag with endpoint constraints for X2
m2 <- pdl_L2(data=datasim, y.name="y", x.name="X2", a=5, b=15)
summary(m2)     ## summary of coefficients beta
summary.lm(m2)  ## summary of parameters gamma
plot(m2, main=expression(paste(X[2])))
```


# This clears everything in the environment!
rm(list = ls())
# Clear the R Console
cat("\014")
cat("\n------------------\n")

source("ols.r")
source("print_ols_results.r")
source("hwhite.r")
source("print_hwhite_results.r")
source("nwest.r")
source("print_nwest_results.r")

set.seed(10)

n <- 100
k <- 3

xtmp <- matrix(rnorm((n*k)-n, mean = 0, sd = 1)*10, ncol = k-1)

tt <- 1:n

e <- rnorm(n, mean = 0, sd = sqrt(tt)) # heteroscedastic error term

u <- rep(0, n) # serial heteroscedastic error term
u[1] <- e[1]
for (i in 2:n){
  u[i] <- 0.8*u[i-1] + e[i]
}

b <- rep(1, k)

iota <- rep(1, n)

x <- cbind(iota, xtmp)

# generate y-data
y <- x %*% b + u

vnames <- c('yvar', 'iota', 'x1', 'x2')

# Load the package
library(readxl)

# Specify the path to your Excel file
filepath <- './nwest_data.xlsx'

# Read the sheets
x <- read_excel(filepath, sheet = "x")
y <- read_excel(filepath, sheet = "y")

# Convert the dataframe to a matrix
x <- data.matrix(x)
y <- data.matrix(y)

################################################################################
# OLS
# do ols regression
res_ols1 <- lm(y ~ x - 1) # '-1' to remove the intercept that lm() adds by default
summary(res_ols1)

res_ols2 <- ols(y, x)
print_ols_results(res_ols2, c("Y", "iota", "x1", "x2"))

################################################################################
# HC
# compare to HC regression
library(car)
res_hwhite1 <- hccm(lm(y ~ x - 1))
print(res_hwhite1)

res_hwhite2 <- hwhite(y, x)
print_hwhite_results(res_hwhite2, c("Y", "iota", "x1", "x2"))

################################################################################
# NWEST
# compare to Newey-West regression
library(sandwich)
library(lmtest)
nlag <- 2
res_nwest1 <- coeftest(lm(y ~ x - 1), vcov. = NeweyWest(lm(y ~ x - 1), lag = nlag, prewhite = FALSE))
print(res_nwest1)

res_nwest2 <- nwest(y, x, nlag)
print_nwest_results(res_nwest2, c("Y", "iota", "x1", "x2"))

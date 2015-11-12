# Code to compute asymptotic p-values & critical values for the Johansen et al. (2000) modified Trace tests for cointegration in the presence of structural breaks.

crit.val.ca.jomoni <- function (z, pr.max = 2, T1 = NULL , T2 = NULL, trace.value = NULL, type = c("Hl(r)", "Hc(r)")){
# Written by Ryan Godwin & David Giles (Dept. of Economics, Univesity of Victoria, Canada)
# Additional code to be called from a function by Johannes Lips
# Last Update: 10 November 2015
# z - object of class ca.jomoni
# T1 - observation of the first structural break
# T2 - observation of the second structural break
# trace.value - value of the trace test-statistic
# type - type of Johansen et al (2000) test used
#======================================
if (!(class(z) == "ca.jomoni")) {
    stop("\nPlease, provide object of class 'ca.jomoni' as 'z'.\n")
}
z <- z
n.breaks <- ncol(z@break.matrix)
if (is.null(n.breaks)){
  n.breaks <- 0
} else {
  n.breaks <- as.integer(ncol(z@break.matrix))
}

q <- n.breaks + 1        # q=1 implies no breaks & the values of v1 & v2 are ignored          
#Position the breaks happened in the sample
totobs <- nrow(z@x)
if(!is.null(T2))
{
v1 <- round(T1/totobs,2)               
v2 <- round(T2/totobs,2)
} else {
  v1 <- round(T1/totobs,2)
}
#======================================
# If user requires p-value for one or both Trace statistics, 
# alter one or both of the next 2 lines
if ( type == "Hl(r)")
{
traceL <- trace.value      # Value of Hl(r) statistic
traceC <- NULL 
} else if (type == "Hc(r)")
{
  traceL <- NULL
  traceC<- trace.value            # Value of Hc(r) statistic
}
#=========================================
if (pr.max < 2) {
  stop("\nPlease set pr.max at least to 2.\n")
} else{
  pr_max<- pr.max    # Do NOT make (p-r) greater than 10 - see Johansen et al. (2000)
}
# The values of  "a" and "b" depend on the number (q-1) of structural breaks.
# When q =1, set a=b=0, & this is the case of no structural breaks 
# When q=2, set a=0 &  b = min[ V1 , (1-V1)] 
# When q=3, set a=min[V1, (V2-V1), (1-V2)] & b = min[remaining two V expressions ]

a = c(0, 0, min(v1, v2-v1, 1-v2))[q]
b = c(0, min(v1, 1-v1), median(c(v1,v2-v1,1-v2)))[q]

# lm denotes the logarithm of the mean of the asymptotic distribution
# lv  denotes the logarithm of the variance of the asymptotic distribution
# Add L or C to the names to reflect the H(L) test or the H(c) test.
# See Table 4 of Johansen et al. (2000). 
# First construct critical values for the Hl(r) test; and then for the Hc(r) test

pr<- c(1:pr_max)

lmL<- 3.06+0.456*pr+1.47*a+0.993*b-0.0269*pr^2-0.0363*a*pr-0.0195*b*pr-4.21*a^2-2.35*b^2+0.000840*pr^3+6.01*a^3-1.33*a^2*b+2.04*b^3-2.05/pr-0.304*a/pr+1.06*b/pr+9.35*a^2/pr+3.82*a*b/pr+2.12*b^2/pr-22.8*a^3/pr-7.15*a*b^2/pr-4.95*b^3/pr+0.681/pr^2-0.828*b/pr^2-5.43*a^2/pr^2+13.1*a^3/pr^2+1.5*b^3/pr^2
lvL<- 3.97+0.314*pr+1.79*a+0.256*b-0.00898*pr^2-0.0688*a*pr-4.08*a^2+4.75*a^3-0.587*b^3-2.47/pr+1.62*a/pr+3.13*b/pr-4.52*a^2/pr-1.21*a*b/pr-5.87*b^2/pr+4.89*b^3/pr+0.874/pr^2-0.865*b/pr^2
meanL<- exp(lmL)-(3-q)*pr
varL<- exp(lvL)-2*(3-q)*pr
# Use the asymptotic mean and variance to obtain the shape and scale parameters 
# of the gamma distribution to be used to approximate the asymptotic distribution 
# of the test statistics, and hence obtain the desired quantiles under the null:
thetaL<- varL/meanL
kL<- meanL^2/varL

lmC<- 2.80+0.501*pr+1.43*a+0.399*b-0.0309*pr^2-0.0600*a*pr-5.72*a^2-1.12*a*b-1.70*b^2+0.000974*pr^3+0.168*a^2*pr+6.34*a^3+1.89*a*b^2+1.85*b^3-2.19/pr-0.438*a/pr+1.79*b/pr+6.03*a^2/pr+3.08*a*b/pr-1.97*b^2/pr-8.08*a^3/pr-5.79*a*b^2/pr+0.717/pr^2-1.29*b/pr^2-1.52*a^2/pr^2+2.87*b^2/pr^2-2.03*b^3/pr^2
lvC<- 3.78+0.346*pr+0.859*a-0.0106*pr^2-0.0339*a*pr-2.35*a^2+3.95*a^3-0.282*b^3-2.73/pr+0.874*a/pr+2.36*b/pr-2.88*a^2/pr-4.44*b^2/pr+4.31*b^3/pr+1.02/pr^2-0.807*b/pr^2
meanC<- exp(lmC)-(3-q)*pr
varC<- exp(lvC)-2*(3-q)*pr

# Use the asymptotic mean and variance to obtain the shape and scale parameters of the 
# gamma distribution to be used to approximate the asymptotic distribution of the test 
# statistics:

thetaC<- varC/meanC
kC<- meanC^2/varC

# Now create a list of critical values, and p-values (if either traceL or traceC is non-zero):

crit <- cbind(sapply(c(.9,.95,.99) , function(x) round(c(qgamma(x, shape=kL,scale=thetaL)),2)),
	sapply(c(.9,.95,.99) , function(x) round(c(qgamma(x, shape=kC,scale=thetaC)),2)))
colnames(crit) <- rep(c(0.9,0.95,0.99),2)
rownames(crit) <- paste("(p-r) =", pr, sep=" ")
crit.list <- list("Asymptotic Critical Values" = c(
                  list("Hc(r)" = crit[ , 4:6]),
                  list("Hl(r)" = crit[ , 1:3]))
                  )
if(!is.null(traceL)){
    pvalL <- matrix(round(1 - pgamma(traceL, shape=kL, scale = thetaL),3))
    colnames(pvalL) <- c(paste("Pr((H_l)  ",traceL,")", sep = ""))
    rownames(pvalL) <-  paste("(p-r) =", pr, sep=" ")
    pval.list <- list("Asymptotic p-Values" = pvalL)
}
if(!is.null(traceC)){
    pvalL <- matrix(round(1 - pgamma(traceC, shape=kC, scale = thetaC),3))
    colnames(pvalL) <- c(paste("Pr((H_c)  ",traceC,")", sep = ""))
    rownames(pvalL) <-  paste("(p-r) =", pr, sep=" ")
    pval.list <- list("Asymptotic p-Values" = pvalL)
}
  return(crit.list)
  return(pval.list)
}

crit.val.ca.jomoni(vecm.EEXPeak.jomoni.nc , r = 1, T1 = 320, T2 = 832, trace.value= 13.34, type = "Hl(r)")



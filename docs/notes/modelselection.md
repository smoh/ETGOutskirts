# model selection

Model selection is difficult.

Feigelson+12
"The common procedure in astronomy based on the reduced chi-squared χ2ν ≃ 1 is a primitive technique not used by statisticians or researchers in other fields."

### F-test
### Levene's test
### Bartlett's test
### Brown-Forsythe test

### Bayesian evidence
costly to compute since we need to map the entire paramter space.

### Cross-validation
* Leave-one-out

### Akaike Information Criterion (AIC, Akaike ????)
AIC = -2*lnL_max + 2k = chi_min^2 + 2k
L_max: maximum likelihood
k: number of free parameters
Jeffreys' scale - delta(AIC) > 5: strong, delta(AIC) > 10: decisive

maximizing likelihood is equivalent to minimizing chi^2 in the case of Gaussian noise.

### Bayesian Information Criterion (BIC, Schwarz 1978)

- assumes the data points are independent and identically distributed
- as an easier way to estimate Bayesian evidence?

BIC = -2*lnL_max + k*lnN
N: number of data points

I am uncertain if all these criteria will prove to be useful in the end... except that it makes dealing with a large data easier.

## References
Liddle 2007
Feigelson & Babu 2012

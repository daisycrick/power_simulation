rm(list=ls())

set.seed(12345)

Nrep = 1000	#Number of simulation replicates
N = 500000*0.3	#Number of individuals who smoked during pregnancy
p = 0.3		#Decreaser allele frequency
q <- 1-p 	#Increaser allele frequency
vg <- 0.03	#Variance in CPD explained by smoking SNP
beta_zx <- sqrt(vg)
beta_xy <- 0.2	#Causal effect of CPD on liability to left handedness. Measured in SDs.
prevalence <- 0.1	#Prevalence of left handedness
threshold <- qnorm(p=prevalence, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)	#Threshold value for declaring an individual affected
alpha <- 0.05	#Significance level
vex <- 1 - beta_zx^2	#Residual variance in cigarettes per day
vey <- 1 - beta_xy^2	#Residual variance in liability to left handedness

pval <- vector(length = Nrep)

a <- sqrt(1/(2*p*q)) #Create genetic variable of variance one. Assume no dominance.

for(j in 1:Nrep) {

	#Sample mothers' genotypes
	Zm <- sample(x = c(-a,0,a), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
	#Sample fathers' genotypes
	Zf <- sample(x = c(-a,0,a), size = N, replace = TRUE, prob = c(p^2, 2*p*q, q^2))

	Zo <- vector(length = N)

	#Simulate offspring genotype
	r <- runif(N)
	for (i in 1:N) { #This is not efficient...

	  if((Zm[i]==-a) && (Zf[i]==-a)) {Zo[i] = -a}
	  if((Zm[i]==-a) && (Zf[i]==0))  {if(r[i] <= 0.5) {Zo[i] = -a} else {Zo[i] = 0}}
	  if((Zm[i]==-a) && (Zf[i]==a))  {Zo[i] = 0}
	  if((Zm[i]==0) &&  (Zf[i]==-a)) {if(r[i] <= 0.5) {Zo[i] = -a} else {Zo[i] = 0}}
	  if((Zm[i]==0) &&  (Zf[i]==0))  {
 	    if(r[i] <= 0.25) {Zo[i] = -a}
	    if(r[i] > 0.25 && r[i] <= 0.75) {Zo[i] = 0}
	    if(r[i] > 0.75) {Zo[i] = a}
	  }
	  if((Zm[i]==0) &&  (Zf[i]==a))  {if(r[i] <= 0.5) {Zo[i] = a} else {Zo[i] = 0}}
	  if((Zm[i]==a) &&  (Zf[i]==-a)) {Zo[i] = 0}
	  if((Zm[i]==a) &&  (Zf[i]==0))  {if(r[i] <= 0.5) {Zo[i] = a} else {Zo[i] = 0}}
	  if((Zm[i]==a) &&  (Zf[i]==a))  {Zo[i] = a}
	}

	#Simulate maternal number of cigarettes
	X_m <- beta_zx*Zm + rnorm(N,0,sqrt(vex))

	#Simulate underlying distribution of liability to left handedness
	Y_l <- beta_xy*X_m + rnorm(N,0,sqrt(vey))

	Y_o <- as.integer(Y_l > threshold)

	result <- glm(Y_o ~ Zo, family='binomial')
	pval[j] <- summary(result)$coefficients[2,4]
}

power <- sum(as.integer(pval < alpha))/Nrep

power

li_power <- power - 1.96*sqrt((power/(1-power))/Nrep)
li_power
hi_power <- power + 1.96*sqrt((power/(1-power))/Nrep)
hi_power

#Compare logistic beta approximation to simulated logistic coefficient
beta_xy*(pi/sqrt(3))
summary(glm(Y_o ~ X_m, family='binomial'))

#Approximate odds ratio per 1SD CPD
OR<-exp(beta_xy*(pi/sqrt(3)))
OR

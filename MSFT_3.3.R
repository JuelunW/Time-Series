library(quantmod)
library(PerformanceAnalytics)
library(fUnitRoots)#adf.test
library(TSA)#eacf, gBox,(conflict with function acf)
#detach("package:TSA", unload = TRUE)
library(tseries)#GARCH
library(forecast)
library(MTS)
library(urca)#cointegration test

setwd("O:/FIN 620 - Financial Econometrics/01 Project")

getSymbols(Symbol="MSFT",from = "2006-01-01", to = "2020-12-31")

# Price
P.daily <- MSFT$MSFT.Adjusted
colnames(P.daily) <- c("daily adjusted closing prices")
P.monthly <- apply.monthly(MSFT$MSFT.Adjusted, function(x)tail(x, 1))
colnames(P.monthly) <- c("monthly adjusted closing prices")
P.log <- log(MSFT$MSFT.Adjusted)
colnames(P.log) <- c("log daily adjusted closing prices")

# Return
R.daily <- periodReturn(MSFT, period = "daily", type = "arithmetic")
colnames(R.daily) <- c("simple daily returns")
R.monthly <- periodReturn(MSFT, period = "monthly", type = "arithmetic")
colnames(R.monthly) <- c("simple monthly returns")
R.log <- na.omit(log(MSFT$MSFT.Adjusted/lag(MSFT$MSFT.Adjusted)))
colnames(R.log) <- c("log returns")
R.logm <- P.monthly/lag(P.monthly)
colnames(R.logm) <- c("monthly log returns")

# Volume
V.daily <- MSFT$MSFT.Volume
colnames(V.daily) <- c("daily volumes")
V.monthly <- apply.monthly(MSFT$MSFT.Volume, sum)
colnames(V.monthly) <- c("monthly volumes")


#### Descriptive Statistics #### 
##### Statistics table ####
nn <- list(P.daily, P.monthly, R.daily, R.monthly, R.log, R.logm, V.daily, V.monthly)
library(fBasics)
bb <- function(x)lapply(x, basicStats)
bb(nn)

#mm <- function(x)lapply(x, summary)
#mm(nn)
#f <- function(x)lapply(x, function(x)c(sd(x), skewness(x), kurtosis(x)))
#f(nn)

##### Histogram ####
j <- c("Price", "Price", "Return", "Return", "Return","Return", "Volume", "Volume")
for (i in c(1:length(nn))) {
  png(filename = paste("Hist of", colnames(nn[[i]]), ".png"))
  hist(nn[[i]], main = paste("Histogram of", colnames(nn[[i]])), xlab = j[i], col = "gray")
  #a <- coredata(nn[[i]])
  #a <- a[order(a)]
  #y <- lines(a, dnorm(a, mean = mean(nn[[i]]), sd = sd(nn[[i]])), col = "red")
  #plot(density(nn[[i]]), col = "red")
  dev.off()
}


##### Box plot ####
for (i in c(1:length(nn))) {
  png(filename = paste("Boxplot of", colnames(nn[[i]]), ".png"))
  boxplot(nn[[i]], main = paste("Boxplot of", colnames(nn[[i]])), col = "gray")
  dev.off()
}

#### Visualization ####
par(mfrow = c(2,1))

plot(P.daily)
plot(V.daily)
plot(P.monthly)
plot(V.monthly)
plot(R.daily)
chart.CumReturns(R.daily, main = "The Cumulative Daily Return")
plot(R.monthly)
chart.CumReturns(R.monthly, main = "The Cumulative Monthly Return")
plot(R.log)
chart.CumReturns(R.log, main = "The Cumulative Log Return")
plot(R.logm)
chart.CumReturns(R.logm, main = "The Cumulative monthly Log Return")

dev.off()

#### Prep ####
##### Test for normal dist ####
###### Perform test for mean return being zero ####
t.test(as.vector(R.monthly))

###### Perform normality test using the Jaque-Bera method ####
normalTest(as.vector(R.monthly),method="jb")
t.test(as.vector(R.monthly)/(6/length(R.monthly))^0.5)#Skewness
t.test((as.vector(R.monthly) - 3)/(24/length(R.monthly))^0.5)#Kurtosis

###### qq plot ####
qqnorm(R.monthly)
qqline(R.monthly, col = "red")

###### Perform skewness and kurtosis tests ####
s3=skewness(R.monthly)
T <- length(R.monthly)
tst = abs(s3/sqrt(6/T))  # t-ratio for testing skewness
pv <- 2*(1-pnorm(tst)) # calculate p value for 2-tailed value 2(1-cdf(|t|))
pv

k4 <- kurtosis(R.monthly)  # already subtracts 3 from kurtosis 
tst <-abs( k4/sqrt(24/T))  # test excess kurtosis
pv <- 2*(1-pnorm(tst)) # calculate p value for 2-tailed value 2(1-cdf(|t|))
pv


##### 1.Test for stationary/unit-root/seasonality ####
#Null: non-stationary, unit-root
#p-value should <5% 
R.nm <- na.omit(diff(R.monthly))
adf.test(R.nm)

P.nm <- na.omit(diff(diff(P.monthly, lag = 4)))
adf.test(P.nm)


##### 2.Test for autocorrelation ####
#Null: no autocorrelation
#p-value should <5% 
log(length(R.nm))  # check what should be m --> m ~ log(T) #freedom = lag = m

Box.test(as.numeric(R.nm), lag = 5, type = "Ljung")
log(length(P.nm))
Box.test(as.numeric(P.nm), lag = 5, type = "Ljung")
Box.test(as.numeric(P.nm1), lag = 5, type = "Ljung")

#### I. ARIMA ####
##### Step 1: identify seasonality and trends for differencing ####

p4 <- function(P){
  d1 <- diff(P)
  d4 <- diff(P, lag = 4)
  d14 <- diff(diff(P, lag = 4))
  par(mfrow=c(3, 4))
  ts.plot(P, main = c("Series P"))
  ts.plot(d1, main = c("Series diff(P)"))
  ts.plot(d4, main = c("Series diff(P, lag = 4)"))
  ts.plot(d14, main = c("Series diff(diff(P, lag = 4))"))
  acf(P, lag.max = 20, main = NA)
  acf(d1, lag.max = 20, main = NA)
  acf(d4, lag.max = 20, main = NA)
  acf(d14, lag.max = 20, main = NA)
  pacf(P, lag.max = 20, main = NA)
  pacf(d1, lag.max = 20, main = NA)
  pacf(d4, lag.max = 20, main = NA)
  pacf(d14, lag.max = 20, main = NA)
}



P.daily <- as.ts(P.daily)
P.log <- as.ts(P.log)
R.daily <- as.ts(R.daily)
R.log <- as.ts(R.log)
P.monthly <- as.ts(P.monthly)
R.monthly <- as.ts(R.monthly)

p4(as.ts(P.daily))
p4(as.ts(P.log))
p4(as.ts(R.daily))
p4(as.ts(R.log))
p4(as.ts(P.monthly))
p4(as.ts(R.monthly))


dev.off()

###### ACF & PACF ####
par(mfrow = c(2, 1))

acf(R.nm)
#acf(R.log,lag=20) # specify the number of ACF to compute  
pacf(R.nm) #compute partial ACF of the GNP growth rate 

acf(P.nm)
pacf(P.nm)

acf(P.nm1)
pacf(P.nm1)

dev.off()


##### Step 2: identifying order ####

eacf(R.monthly, ar.max = 10, ma.max = 10)
eacf(R.nm, ar.max = 10, ma.max = 10)
m1 = arima(R.monthly, order = c(0, 1, 1))
m2 = arima(R.monthly, order = c(2, 1, 0))
m3 = arima(R.monthly, order = c(5, 1, 0))
m4 = arima(R.monthly, order = c(6, 1, 0))
m5 = arima(R.monthly, order = c(1, 1, 1))
m6 = arima(R.monthly, order = c(2, 1, 2))
m7 = arima(R.monthly, order = c(4, 1, 3))
m8 = arima(R.monthly, order = c(4, 1, 4))
m9 = arima(R.monthly, order = c(0, 0, 1))
m1$aic
m2$aic
m3$aic
m4$aic
m5$aic
m6$aic
m7$aic
m8$aic
m9$aic

eacf(P.monthly, ar.max = 10, ma.max = 10)
eacf(P.nm, ar.max = 20, ma.max = 20)
m11 = arima(P.monthly, order = c(1, 1, 7), seasonal = c(0, 1, 1))
m12 = arima(P.monthly, order = c(1, 1, 7), seasonal = c(1, 1, 0))
m13 = arima(P.monthly, order = c(2, 1, 7), seasonal = c(0, 1, 1))
m14 = arima(P.monthly, order = c(2, 1, 7), seasonal = c(1, 1, 0))
m15 = arima(P.monthly, order = c(3, 1, 7), seasonal = c(0, 1, 1))
m16 = arima(P.monthly, order = c(3, 1, 7), seasonal = c(1, 1, 0))
m17 = arima(P.monthly, order = c(2, 1, 2), seasonal = c(1, 1, 0))
m18 = arima(P.monthly, order = c(3, 1, 2), seasonal = c(1, 1, 0))
m19 = arima(P.monthly, order = c(3, 1, 3), seasonal = c(1, 1, 0))

m11$aic
m12$aic
m13$aic
m14$aic
m15$aic
m16$aic
m17$aic
m18$aic
m19$aic



m111 <- auto.arima(R.log)
m222 <- auto.arima(P.daily)
m333 <- auto.arima(P.nm)

##### Step 3: Estimation ARIMA ####
m1   # see parameter estimates 
names(m1) # See the output

sqrt(m1$sigma2) # Residual standard error


##### Step 4: Model checking ####

tsdiag(m1) # Model checking 
tsdiag(m1, gof=20) ## increasing the number of residual ACFs used in checking. 
Box.test(m1$residuals, lag = 12, type = 'Ljung')
#1-pchisq(8.4823, 9) # Compute p-value using 12-3 = 9 degrees of freedom

tsdiag(m13, gof=20)

library(fpp2)
checkresiduals(m1)
checkresiduals(m4)
checkresiduals(m13)

qqnorm(residuals(m1))
qqline(residuals(m1), col = "red")

qqnorm(residuals(m13))
qqline(residuals(m13), col = "red")

##### Step 5: Prediction/forecasting ####
# Confidence Interval
pp <- function(m){
  predict(m, 10) # Prediction 
  p <- predict(m, 10) 
  names(p) 
  lcl <- p$pred - 1.96*p$se ## calculate lower 95% interval
  ucl <- p$pred + 1.96*p$se ## calculate upper 95% interval
  ci <- cbind(lcl, p$pred, ucl) 
  round(ci,2)
  print(ci)
}

pp(m1)
pp(m4)
pp(m13)

# plot
autoplot(forecast(m1))
autoplot(forecast(m13))
#au <- auto.arima(R.monthly)
#forward <- forecast(au, h=10,level=c(90,95))
#plot(forward)
###### Comparison between real and predicted value ####
rnp <- function(m){
  predicted <- data.frame(forecast(m, 10))$Point.Forecast
  MM <- get(getSymbols(Symbol = "MSFT",from = "2021-01-01", to = "2021-10-30"))
  cc <- periodReturn(MM, period = "monthly", type = "arithmetic")
  cc$predicted <- predicted
  return(cc)
}
rnp1 <- function(m){
  predicted <- data.frame(forecast(m, 10))$Point.Forecast
  MM <- get(getSymbols(Symbol = "MSFT",from = "2021-01-01", to = "2021-10-30"))
  cc <- apply.monthly(MM$MSFT.Adjusted, function(x)tail(x, 1))
  cc$predicted <- predicted
  return(cc)
}

MTSplot(rnp(m1))
plot(rnp(m1), legend.loc = "top")
MTSplot(rnp(m4))
MTSplot(rnp1(m13))
plot(rnp1(m13), legend.loc = "top")

#### II. VAR ####
###### Second variate ####
getSymbols(Symbol="SPY",from = "2006-01-01", to = "2020-12-31")
# Price
PS.daily <- SPY$SPY.Adjusted
colnames(PS.daily) <- c("daily adjusted closing prices")
PS.monthly <- apply.monthly(SPY$SPY.Adjusted, function(x)tail(x, 1))
colnames(PS.monthly) <- c("monthly adjusted closing prices")
PS.log <- log(SPY$SPY.Adjusted)
colnames(PS.log) <- c("log daily adjusted closing prices")

# Return
RS.daily <- periodReturn(SPY, period = "daily", type = "arithmetic")
colnames(RS.daily) <- c("simple daily returns")
RS.monthly <- periodReturn(SPY, period = "monthly", type = "arithmetic")
colnames(RS.monthly) <- c("simple monthly returns")
RS.log <- na.omit(log(SPY$SPY.Adjusted/lag(SPY$SPY.Adjusted)))
colnames(RS.log) <- c("log returns")
RS.logm <- PS.monthly/lag(PS.monthly)
colnames(RS.logm) <- c("monthly log returns")


PS.monthly <- as.ts(PS.monthly)
RS.monthly <- as.ts(RS.monthly)

##### Step 1: Test for ? ####

regre <- lm(R.monthly~RS.monthly) 
summary(regre)
wt = regre$residuals
plot(wt, type = 'l', ylab = "residuals of regression")
abline(h = mean(wt))
abline(h = mean(wt) + sd(wt), lty = 2)
abline(h = mean(wt)- sd(wt), lty = 2)

P.nnm <- diff(log(P.monthly))
PS.nnm <- diff(log(PS.monthly))

regre1 <- lm(P.nnm~PS.nnm) 
summary(regre1)
wt1 = regre1$residuals
plot(wt1, type = 'l', ylab = "residuals of regression")
abline(h = mean(wt1))
abline(h = mean(wt1) + sd(wt1), lty = 2)
abline(h = mean(wt1)- sd(wt1), lty = 2)

###### plots with two variates ####
par(mfrow = c(2, 1))

plot(R.monthly)
plot(RS.monthly)

plot(P.monthly)
plot(PS.monthly)

dev.off()

##### Step 1: Test for ?, CCF plot ####
zt <- cbind(as.numeric(R.monthly), as.numeric(RS.monthly))
colnames(zt) <- c("MSFT monthly return", "SPY monthly return")
ccm(zt)
gt <- diffM(zt)
ccm(gt)
MTSplot(gt)

zt1 <- cbind(as.numeric(P.monthly), as.numeric(PS.monthly))
colnames(zt1) <- c("MSFT monthly price", "SPY monthly price")
ccm(zt1)
gt1 <- diffM(zt1)
ccm(gt1)
MTSplot(gt1)
gt2 <- diffM(log(zt1))
ccm(gt2)
MTSplot(gt2)

##### Step 2: Identifying order and Cointegration test ####
VARorder(zt)
m111 <- VAR(zt,p = 1) ### Fit a VAR(3) model
m111$coef
VARorder(gt)
m112 <- VAR(gt,p = 6) ### Fit a VAR(3) model
m112$coef

m122=ca.jo(gt, K=6)
summary(m122)
# null hypothesis is of cointegration 
# In our case, the test fails to reject 
# the null hypothesis for the first time when r = 1 
# Therefore, we have one cointegration relationship.

VARorder(gt1)
m211 <- VAR(gt1,p = 10) ### Fit a VAR(10) model
m211$coef

m221=ca.jo(gt1, K=10)
summary(m221)


VARorder(gt2)
m212 <- VAR(gt2,p = 1)
m212$coef

m222=ca.jo(gt2, K=1)
summary(m222)


##### Step 3: Model checking ####
MTSdiag(m111)
MTSdiag(m112)
MTSdiag(m211)
MTSdiag(m212)
#MTSdiag(m111, gof = 24, adj = 0, level = F)

# Multivariate Ljung-Box Q Statistics
mb <- function(x){
  res <- x$residuals
  #mq(res1,adj=18)
  mq(res,lag = 24,adj=18)
}

mb(m111)
mb(m112)
mb(m211)
mb(m212)

#m112 <- refVAR(m111, thres = 1.96) # threshhold is the t -stat
#m112$coef

##### Step 4: Forecast ####
library(vars)
VARpred(m112,8)
VARpred(m211,8)
VARpred(m212,8)
###### Plot #### 
vp <- function(x){
  forecasts <- predict(x)
  plot(forecasts)
}
detach("package:vars", unload = TRUE)
vp(m111)
vp(m112)
vp(m211)
vp(m212)


#### III. GARCH ####
##### Step 2: Test for ARCH effects, Ljung-Box test ####
r <- resid(m1)

rr <- r^2
r_abs <- abs(r)

archTest(r)

Box.test(R.nm, lag = 12, type = 'Ljung')
at <- R.nm - mean(R.nm)
Box.test(at^2, lag = 12, type = 'Ljung')#should <5%

Box.test(P.nm, lag = 12, type = 'Ljung')
at1 <- P.nm - mean(P.nm)
Box.test(at1^2,lag=12,type='Ljung')#should <5%

P.nnm <- as.ts(na.omit(diff(diff(log(P.monthly), lag = 4))))
Box.test(P.nnm, lag = 12, type = 'Ljung')
at2 <- P.nnm - mean(P.nnm)
Box.test(at2^2,lag=12,type='Ljung')#should <5%

##### Step 3: Identifying order ####
###### ACF & PACF ####
par(mfrow=c(2,1))

acf(abs(R.nm))
pacf(abs(R.nm))

acf(R.nm^2)
pacf(R.nm^2)

acf(abs(P.nm))
pacf(abs(P.nm))

acf(P.nm^2)
pacf(P.nm^2)

par(mfrow=c(1,2))
acf(as.vector(rr),main="ACF of Squared Residuals"); 
pacf(as.vector(rr),main="PACF of Squared Residuals") # homoscedasticity check

par(mfrow=c(1,2))
acf(as.vector(r_abs),main="ACF of Absolute Residuals"); 
pacf(as.vector(r_abs),main="PACF of Absolute Residuals") # homoscedasticity check
dev.off()

###### Construction ####
eacf(R.nm^2)
eacf(abs(R.nm))

eacf(P.monthly^2)
eacf(abs(P.monthly))

g01 <- garch(R.nm, order = c(1, 1))
g02 <- garch(R.nm, order = c(2, 1))
g03 <- garch(R.nm, order = c(1, 2))
g04 <- garch(R.nm, order = c(1, 0))
g05 <- garch(R.nm, order = c(2, 0))

summary(g01)
summary(g02)

AIC(g01)
AIC(g02)
AIC(g03)
AIC(g04)
AIC(g05)


g11 <- garch(P.nm, order = c(1, 1))
g12 <- garch(P.nm, order = c(2, 1))
g13 <- garch(P.nm, order = c(1, 2))
g14 <- garch(P.nm, order = c(1, 0))
g15 <- garch(P.nm, order = c(2, 0))

summary(g11)
summary(g12)

AIC(g11)
AIC(g12)
AIC(g13)
AIC(g14)
AIC(g15)



g21 <- garch(P.nnm, order = c(1, 1))
g22 <- garch(P.nnm, order = c(2, 1))
g23 <- garch(P.nnm, order = c(1, 2))
g24 <- garch(P.nnm, order = c(1, 0))
g25 <- garch(P.nnm, order = c(2, 0))

summary(g21)
summary(g22)

AIC(g21)
AIC(g22)
AIC(g23)
AIC(g24)
AIC(g25)

##### Step 3: MODEL Diagnostics ####
plot(residuals(g01), type = 'h', ylab = 'Standardized Residuals')
plot(residuals(g12), type = 'h', ylab = 'Standardized Residuals')
plot(residuals(g21), type = 'h', ylab = 'Standardized Residuals')


# checking normality assumption 
qqnorm(residuals(g01))
qqline(residuals(g01), col = "red")

qqnorm(residuals(g12))
qqline(residuals(g12), col = "red")

qqnorm(residuals(g21))
qqline(residuals(g21), col = "red")

# ACF and Generalized Portmanteau Test of residual  
par(mfrow=c(2,1))

acf(residuals(g01), na.action = na.omit)
Box.test(residuals(g01), lag = 12, type = 'Ljung')
gBox(g01)

acf(residuals(g12), na.action = na.omit)
Box.test(residuals(g12), lag = 12, type = 'Ljung')
gBox(g12)

acf(residuals(g21), na.action = na.omit)
Box.test(residuals(g21), lag = 12, type = 'Ljung')
gBox(g21)

# ACF and Generalized Portmanteau Test of residual square 
acf(residuals(g01)^2, na.action = na.omit)
gBox(g01, method = 'squared')

acf(residuals(g12)^2, na.action = na.omit)
gBox(g12, method = 'squared')

acf(residuals(g21)^2, na.action = na.omit)
gBox(g21, method = 'squared')

# ACF and Generalized Portmanteau Test of absolute residual square 
acf(abs(residuals(g01)), na.action = na.omit)
gBox(g01, method = 'absolute')

acf(abs(residuals(g12)), na.action = na.omit)
gBox(g12, method = 'absolute')

acf(abs(residuals(g21)), na.action = na.omit)
gBox(g21, method = 'absolute')

dev.off()


##### Step 4: Forecasting ####

plot(fitted(g01)[, 1]^2, type = 'l', ylab = 'Conditional Var - GARCH(1,1)', xlab = 't')
plot(fitted(g12)[, 1]^2, type = 'l', ylab = 'Conditional Var - GARCH(2,1)', xlab = 't')
plot(fitted(g21)[, 1]^2, type = 'l', ylab = 'Conditional Var - GARCH(1,1)', xlab = 't')

pm1 <- predict(g01,10)

#### IV. VaR ####

source("RMfit.R")
source("RMeasure.R")

### RiskMetrics 
LossVal <- -log(R.monthly+1) # Loss value
RMfit(LossVal,estim = F)

### Econometric modeling
gfit1=garchFit(~garch(1,1), data = LossVal,trace=FALSE) 
summary(gfit1)
pm1=predict(gfit1,10)
pm1
RMeasure(pm1[1,1], pm1[1,2])


#### 10-day VaR
names(pm1)
v1=sqrt(sum(pm1$standardDeviation^2))
RMeasure(-pm1[1,1],v1)
gfit2=garchFit(~garch(1,1),data=LossVal,trace=F,cond.dist="std")
summary(gfit2)
pm2=predict(gfit2,1)
pm2
RMeasure(-0.01657841,0.06199958,cond.dist="std",df=6.0583797)

### Empirical quantile and quantile regression


### Extreme value theory
require("evir")
m4=gev(LossVal,block=21)
m4
source("evtVaR.R")
evtVaR(-0.37233314, 0.04747955, 0.09309143 )


### codes for other methods in the book
### Peaks over threshold
m4a=pot(LossVal,thres=0.01)
plot(m4a)
riskmeasures(m4a,c(0.95,0.99,0.999))

### generalized Pareto distribution
m5=gpd(LossVal,0.01)
m5
plot(m5)
riskmeasures(m5,c(0.95,0.99,0.999))








#### ‘other GARCH ####
library(rugarch)
# default model.
spec = ugarchspec()
print(spec)

def.fit = ugarchfit(spec = spec, data = R.nm)
print(def.fit)

# own model
spec = ugarchspec(variance.model=list(model="fGARCH",submodel = "TGARCH", garchOrder=c(1,1)), 
                  mean.model=list(armaOrder=c(1,0),include.mean=TRUE),  
                  distribution.model="norm" ) 
def.fit = ugarchfit(spec = spec, data = R.nm)
print(def.fit)
plot(def.fit,which=3)
plot(def.fit,which=9)


### egarch 
spec = ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1)), 
                  mean.model=list(armaOrder=c(1,2), include.mean=TRUE),  
                  distribution.model="norm" )
def.fit = ugarchfit(spec = spec, data = R.nm)
print(def.fit)
plot(def.fit,which=3)
plot(def.fit,which=9)


### igarch 
spec = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)), 
                  mean.model=list(armaOrder=c(1,0), include.mean=TRUE),  
                  distribution.model="norm" ) 

def.fit = ugarchfit(spec = spec, data = R.nm)
print(def.fit)
plot(def.fit,which=3)
plot(def.fit,which=9)









 


 



 




rm(list=ls())

data(macrodata)
str(macrodata)
colnames(macrodata)

# subset
data <- macrodata[,c("RGDP_Q_G","INT_KTB3Y","USDKRW_G","INT_CALL")]
plot.xts(data,multi.panel=TRUE,yaxis.same = FALSE)


# obtain VAR coefficients ans Covariance by OLS
mod1 <- VAR_OLS(data,p=3)
class(mod1)
plot(mod1)

newdata <- xts(data.frame(INT_CALL=c(0.81,0.91,1.0,1.2)),order.by=as.yearqtr(seq(2021.5,2022.25,by=0.25)))
Ypred <- predict(mod1, newdata=newdata)
Ypred

IR <- impulse_response(mod1,variable="INT_KTB3Y", period=6, p=3, type="structural",ncol.fig = 1)



# Bayesian VAR
mod2 <- VAR_bayes(data,p=2, N=120, warmup = 50, exos="INT_CALL",
                   minnesota_par=list(rho=0.8, lambdas = c(0.1, 0.5, 2, 100)))
class(mod2)
plot(mod2, ncol.fig = 2)
mod2$fitted

# unconditional forecasting
pred1 <- predict(mod2, newdata=newdata)
plot(pred1,ncol.fig=2)
pred1$interval %>% arrange(var,time)

# conditional forecasting
condition <- xts( data.frame(RGDP_Q_G = c(4.3,0.5,-4.5,-2.3),
                            INT_KTB3Y = c(1.2,1.3,1.6,1.7)),
                  order.by=as.yearqtr(seq(2021.5,2022.25,by=0.25)))

pred2 <- predict(mod2, newdata=newdata, condition=condition)
plot(pred2,ncol.fig=2)
pred2$fitted



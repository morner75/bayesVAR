packages <- c("tidyverse","zoo","xts")
sapply(packages, require, character.only=TRUE, warn.conflicts = FALSE)

rm(list=ls())

data(macrodata)
str(macrodata)
colnames(macrodata)

# subset
data <- macrodata[,c("RGDP_Q_G","INT_KTB3Y","USDKRW_G","INT_CALL")]
plot(data,multi.panel=TRUE,yaxis.same = FALSE)


# obtain VAR coefficients ans Covariance by OLS
Res <- VAR_OLS(data,p=3)
Coef_mat <- Res$Coef_mat
Sigma <- Res$Sigma
Y <- Res$model$Y

newdata <- xts(data.frame(INT_CALL=c(0.81,0.91,1.0,1.2)),order.by=as.yearqtr(seq(2021.5,2022.25,by=0.25)))
Ypred <- predict_simple(Y=Y, newdata=newdata, Coef_mat=Coef_mat)
Ypred

Res2 <- impulse_response(Res,"INT_KTB3Y", period=6, p=3, type="structural",ncol.fig = 1)
Res2


# Bayesian VAR
model <- VAR_bayes(data,p=2, N=1200, warmup = 200, exos="INT_CALL", 
                   minnesota_par=list(rho=0.8, lambdas = c(0.1, 0.5, 2, 100)))
plot(model,ncol.fig = 1, type="fitting")
model$fitted


# conditional forcasting
condition <- xts( data.frame(RGDP_Q_G = c(4.3,0.5,-4.5,-2.3),
                            INT_KTB3Y = c(1.2,1.3,1.6,1.7)), 
                  order.by=as.yearqtr(seq(2021.5,2022.25,by=0.25)))

pred <- predict(model, newdata=newdata, condition=condition)
plot(pred,ncol.fig=2,type="forecast")
pred$fitted



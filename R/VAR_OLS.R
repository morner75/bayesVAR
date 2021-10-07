#' Ordinary Least Square (OLS) estimation for VAR model
#
#' provide OLS estimates of VAR models
#' @param data time-series data with a xts format
#' @param p a integer, time lags of endogenous variables
#' @param exos a character vector, exogenous variables
#' @return Coef_mat: a coefficient matrix, Sigma: covariance matrix,
#' model: a list of two matrices: endogenous variables and explanatory variables
#' @export
VAR_OLS <- function(data, p, exos=colnames(data)[ncol(data)]){
  if(!is.xts(data)) stop("data must be xts class")
  endo <- data[, setdiff(colnames(data),exos)]
  k=ncol(endo)
  m=length(exos)+1
  N = nrow(endo)
  Y <- endo[(p+1):N,]
  Ylag <- map(seq_len(p), function(i) coredata(endo)[(p+1-i):(N-i),]) %>% do.call(cbind,.)
  colnames(Ylag) <- outer(str_c(colnames(endo),"_L"),1:p,str_c) %>% as.vector()
  x <- data[(p+1):N,exos] %>% merge(Intercept=1,.)
  X <- cbind(Ylag,x)
  Beta <- solve(t(X)%*%X)%*%t(X)%*%Y
  Resid <- Y-X%*%Beta
  Sigma <- t(Resid)%*%Resid/(N-k*p-m-1)
  res <- list(Coef_mat=Beta,Sigma=Sigma,model=list(Y=Y,X=X))
  class(res) <- "OLS"
  return(res)
}

#' Impulse-Response function of estimated VAR model
#
#' provide an impulse-response function and its plot
#' @param model OLS object from \code{VAR_OLS} function
#' @param variable a character, a variable which gives a shock
#' @param period a integer, time period
#' @param p a integer, time lags of endogenous variables
#' @param type shock to be used
#' @param ncol.fig a integer, number of figures plotted in a row
#' @return impulse-response function and data to be used in the plot
#' @export
impulse_response <- function(model, p, variable, period,
                                 type=c("origin","structural"), ncol.fig=2)
  UseMethod("impulse_response")

#' @export
impulse_response.OLS <- function(model, p, variable, period,
                                    type=c("origin","structural"), ncol.fig=2){
  Coef_mat <- model$Coef_mat
  Sigma <- model$Sigma
  if(type=="structural") type <- "triangle"

  res <-  IR_mat_generator(period, p, Coef_mat,Sigma,type) %>%
    unnest(cols=data) %>%
    select(term, response, one_of(variable)) %>%
    rename(value=one_of(variable))

  p <- ggplot(res,aes(term,value))+
    geom_line(aes(col=response),size=2)+
    facet_wrap(facet=vars(response),ncol=ncol.fig,scales = "free_y")+
    theme_bw()+
    theme(legend.position = "bottom")
  print(p)
  return(res)
}



#' Simple point predicition with OLS parameter estimates
#'
#' @param Y a xts object, endogeneous variables
#' @param newdata a xts object, exogenous variables for the prediction period
#' @param Coef_mat a matrix, coefficients matrix
#' @export
predict_simple <- function(Y, newdata, Coef_mat){
  if(!is.xts(newdata)) stop("external variables must be xts class")
  period <- nrow(newdata)
  p <- (nrow(Coef_mat)-ncol(newdata)-1)/ncol(Y)
  Ynew <- tail(Y,n=p)
  x <- cbind(1,newdata)
  for(i in seq_len(period)){
    Y_update <- c(tail(Ynew,n=p) %>% rev.zoo() %>% t() %>% as.vector(),x[i,]) %*% Coef_mat %>%
      xts(order.by=index(newdata[i,]))
    Ynew <- rbind(Ynew,Y_update)
  }
  Ynew
}


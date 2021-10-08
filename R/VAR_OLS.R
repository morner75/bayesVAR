#' Ordinary Least Square (OLS) estimation for VAR model
#
#' provide OLS estimates of VAR models
#' @param data time-series data with a xts format
#' @param p a integer, time lags of endogenous variables
#' @param exos a character vector, exogenous variables
#' @param confint level of confidence interval
#' @return OLS object, a coefficient matrix, covariance matrix, fitted values etc.
#' @export
VAR_OLS <- function(data, p, exos=colnames(data)[ncol(data)], confint=0.95){
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
  fitted <- X%*%Beta %>% xts(order.by=index(X))
  Resid <- Y-fitted
  Sigma <- t(Resid)%*%Resid/(N-k*p-m-1)
  interval <- merge(fitted,fitted,fitted) %>%
    setNames(outer(colnames(fitted),c("lwr","med","upr"),str_c,sep=".") %>% as.vector()) +
    outer(diag(Sigma), c(qnorm((1-confint)/2),0,-qnorm((1-confint)/2)), "*") %>%
    as.vector() %>% matrix(nrow=3*k,ncol=nrow(fitted)) %>% t()
  interval <- as.data.frame(interval) %>%
                rownames_to_column(var="time") %>%
                mutate(time=as.yearqtr(time)) %>%
                pivot_longer(-time,names_to=c("var",".value"),names_pattern="(.*)\\.(.*)")
  res <- list(formula = as.formula(str_c(str_c(colnames(Y), collapse=" + "),
                        " ~ " ,str_c(c(colnames(Ylag), names(exos)), collapse=" + "))))
  class(res) <- "OLS"
  res$parameters <- list(Coef_mat=Beta, Sigma=Sigma)
  res$fitted <- fitted
  res$interval <- interval
  res$model_var <- list(Y=Y, X=X)
  return(res)
}



#' Impulse-Response function of estimated VAR model
#
#' provide an impulse-response function and its plot
#' @param model a object
#' @return impulse-response function and data to be used in the plot
#' @export
impulse_response <- function(model,...)
  UseMethod("impulse_response")

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
impulse_response.OLS <- function(model, p, variable, period,
                                    type=c("origin","structural"), ncol.fig=2){
  Coef_mat <- model$parameter$Coef_mat
  Sigma <- model$parameters$Sigma
  if(type=="structural") type <- "triangle"

  res <-  IR_mat_generator(period=period, p=p, Coef_mat=Coef_mat,Sigma=Sigma,type=type) %>%
    unnest(cols=data) %>%
    dplyr::select(term, response, one_of(variable)) %>%
    dplyr::rename(value=one_of(variable))

  p <- ggplot(res,aes(term,value))+
    geom_line(aes(col=response),size=2)+
    facet_wrap(facet=vars(response),ncol=ncol.fig,scales = "free_y")+
    theme_bw()+
    theme(legend.position = "bottom")
  print(p)
  return(res)
}

#' Figures from OLS VAR models
#
#' provide figures from OLS VAR methods
#' @param model OLS VAR model
#' @param ncol.fig a integer, number of figures plotted in a row
#' @export
plot.OLS <- function(model, ncol.fig=2){

  data <- model$interval %>%
    mutate(time=as.yearqtr(time))

  if(!is.null(model$model_var)){
    Y <- model$model_var$Y
    actual <- as.data.frame(Y) %>%
      rownames_to_column("time") %>%
      as_tibble() %>%
      mutate(time=as.yearqtr(time)) %>%
      pivot_longer(-time,"var",values_to="actual")
    data <- full_join(data,actual,by=c("time","var"))
  }

  ggplot(data=data,aes(time,col=var,group=1))+
    geom_line(aes(y=med))+
    geom_ribbon(aes(ymin=lwr,ymax=upr,fill=var),alpha=0.3)+
    facet_wrap(facets=vars(var),ncol=ncol.fig, scales="free_y")+
    theme_bw()+labs(y="")+
    theme(legend.position = "bottom")+
    {if(!is.null(model$model_var)) geom_point(aes(y=actual))}

}




#' Simple point prediction with OLS parameter estimates
#'
#' @param model a OLS model object
#' @param newdata a xts object, exogenous variables for the prediction period
#' @export
predict.OLS <- function(model, newdata){
  if(!is.xts(newdata)) stop("external variables must be xts class")
  period <- nrow(newdata)
  Y <- model$model_var$Y
  Coef_mat <- model$parameters$Coef_mat
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


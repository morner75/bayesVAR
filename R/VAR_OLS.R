#' Ordinary Least Square (OLS) estimation for VAR model
#'
#' provide OLS estimates of VAR models
#' @param data time-series data with a xts format
#' @param p a integer, time lags of endogenous variables
#' @param exos a character vector, exogenous variables
#' @param confint a level of confidence interval
#' @return OLS object, a coefficient matrix, covariance matrix, fitted values etc.
#' @export
VAR_OLS <- function(data, p, exos=colnames(data)[ncol(data)], confint=0.95){
  if(!is.xts(data)) stop("data must be xts class")
  endo <- data[, setdiff(colnames(data),exos)]
  k <- ncol(endo)
  m <- length(exos)+1
  N <-  nrow(endo)
  Y <- endo[(p+1):N,]
  Ylag <- do.call(cbind,lapply(seq_len(p), function(i) coredata(endo)[(p+1-i):(N-i),]))
  colnames(Ylag) <- as.vector(outer(paste0(colnames(endo),"_L"),1:p,paste0))
  X <- cbind(Ylag,merge(Intercept=1, data[(p+1):N,exos]))
  qrdecom <- qr(X)
  Beta <- qr.coef(qrdecom,Y)
  fitted <- qr.fitted(qrdecom,Y)
  Sigma <- t(qr.resid(qrdecom,Y))%*%qr.resid(qrdecom,Y)/(N-k*p-m-1)
  interval <- (matrix(rep(1,3),nrow=1) %x% coredata(fitted)) +
    (matrix(outer(diag(Sigma), c(qnorm((1-confint)/2),0,-qnorm((1-confint)/2)), "*"),nrow=1) %x% matrix(1,nrow=nrow(fitted)))
  colnames(interval) <- as.vector(outer(colnames(fitted),c("lwr","med","upr"),paste,sep="."))
  ans <- list(terms = list( Y =colnames(Y),
                            X= as.formula(paste0(" ~ " ,paste(c(colnames(Ylag), exos), collapse=" + ")))))
  ans$parameters <- list(Coef_mat=Beta, Sigma=Sigma)
  ans$fitted <- fitted
  ans$interval <- as.data.frame(interval) |>
                    dplyr::mutate(time=zoo::index(fitted)) |>
                    tidyr::pivot_longer(-time,names_to=c("var",".value"),names_pattern="(.*)\\.(.*)")
  ans$model_var <- list(Y=Y, X=X)
  class(ans) <- "OLS"
  return(ans)
}


#' @describeIn VAR_OLS fitted values for response variables
#' @param object a OLS VAR model
#' @export
fitted.OLS <- function(object) return (object$fitted)

#' @describeIn VAR_OLS a coefficient matrix of a VAR model
#' @param object a OLS VAR model
#' @export
coef.OLS <- function(object) return (object$parameters$Coef_mat)

#' @describeIn VAR_OLS model Summary
#' @param object a OLS VAR model
#' @export
summary.OLS <- function(object){
  Coef_mat <- coef(object)
  Sigma <- object$parameters$Sigma

  tb2 <-data.frame(response=factor(rownames(Sigma),levels=rownames(Sigma)),
                   `S.E`=sqrt(diag(Sigma))) |>
    tidyr::expand(response,explanatory=rownames(Coef_mat),`S.E`)

  ans <- bayesVAR:::mat2Longtb(Coef_mat,c("explanatory","response","coef.value"))[c("response","explanatory","coef.value")] |>
    dplyr::right_join(tb2,by=c("response","explanatory")) |>
    dplyr::mutate(`Z.value`=`coef.value`/`S.E`,
                  `p.value`=format(1-pnorm(`Z.value`),digits=5,justify="right"))
    print.data.frame(ans)
}



#' Impulse-Response function of estimated VAR model
#
#' provide an impulse-response function and its plot
#' @param object  a model object
#' @param variable a character, a variable which gives a shock
#' @param period a integer, time period
#' @param p a integer, time lags of endogenous variables
#' @param type shock to be used, \code{structural} stands for independent shocks.
#' @param ncol.fig a integer, number of figures plotted in a row
#' @param ... extra arguments
impulse_response <- function(object, p, variable, period,
                             type=c("origin","structural"), ncol.fig=2,...)   UseMethod("impulse_response")

#' @describeIn impulse_response default method
impulse_response.default <- function(object, p, variable, period,
                                     type=c("origin","structural"), ncol.fig=2){
  bayesVAR:::impulse_response.OLS(object, p, variable, period,
                                  type=c("origin","structural"), ncol.fig=2)
}

#' @describeIn impulse_response impulse-response function of a OLS VAR model
#' @export
impulse_response.OLS <- function(object, p, variable, period,
                                    type=c("origin","structural"), ncol.fig=2){
  Coef_mat <- object$parameter$Coef_mat
  Sigma <- object$parameters$Sigma
  if(type=="structural") type <- "triangle"

  res <-  IR_mat_generator(names=colnames(Coef_mat),period=period, p=p, Coef_mat=Coef_mat,Sigma=Sigma,type=type) %>%
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


#' @describeIn VAR_OLS Summary plot
#' @param x OLS VAR model
#' @param ncol.fig a integer, number of figures plotted in a row
#' @export
plot.OLS <- function(x,y=NULL, ncol.fig=2){

  data <- dplyr::mutate(x$interval, time=as.yearqtr(time))

  if(!is.null(x$model_var)){
    Y <- x$model_var$Y
    actual <- as.data.frame(Y) |>
      tibble::rownames_to_column("time") |>
      tibble::as_tibble() |>
      dplyr::mutate(time=as.yearqtr(time)) |>
      tidyr::pivot_longer(-time,"var",values_to="actual")
    data <- dplyr::full_join(data,actual,by=c("time","var"))
  }

  ggplot2::ggplot(data=data,aes(time,col=var,group=1))+
            geom_line(aes(y=med))+
            geom_ribbon(aes(ymin=lwr,ymax=upr,fill=var),alpha=0.3)+
            facet_wrap(facets=vars(var),ncol=ncol.fig, scales="free_y")+
            theme_bw()+labs(y="")+
            theme(legend.position = "bottom")+
            {if(!is.null(x$model_var)) geom_point(aes(y=actual))}

}



#' Simple point prediction with OLS parameter estimates
#'
#' @param object a OLS model object
#' @param newdata a xts object, exogenous variables for the prediction period
#' @export
predict.OLS <- function(object, newdata){
  if(!is.xts(newdata)) stop("external variables must be xts class")
  period <- nrow(newdata)
  Y <- object$model_var$Y
  Coef_mat <- object$parameters$Coef_mat
  p <- (nrow(Coef_mat)-ncol(newdata)-1)/ncol(Y)
  Ynew <- tail(Y,n=p)
  x <- cbind(1,newdata)
  for(i in seq_len(period)){
    Y_update <- c(as.vector(t(zoo::rev.zoo(tail(Ynew,n=p)))),x[i,]) %*% Coef_mat |>
      xts(order.by=index(newdata[i,]))
    Ynew <- rbind(Ynew,Y_update)
  }
  Ynew
}


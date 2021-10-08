#' Bayesian VAR model
#
#' implement Bayesian VAR methods
#' @param data time-series data with a xts format
#' @param p a integer, time lags of endogenous variables
#' @param exos a character vector, exogenous variables
#' @param N  a integer, the number of iterations
#' @param warmup a integer, the number of iterations to discard before convergence
#' @param minnesota_par a list, parameters in Minnesota prior
#' @param prior  a character,  priors and methods to be used
#' @return posterior of parameters, credible bands of responses in Bayesian VAR model
#' @export
VAR_bayes <- function(data, p, exos=colnames(data)[ncol(data)], N=1500, warmup=500, prior="Minnetota-Wishart",
                      minnesota_par=list(rho=1,lambdas=c(0.1, 0.5, 2, 100))){

  # initial values from OLS
  Init <- VAR_OLS(data, p, exos)
  Y <- Init$model_var$Y
  X <- Init$model_var$X
  beta <- as.vector(Init$parameters$Coef_mat)
  Sigma <- Init$parameters$Sigma

  k <- ncol(Y) ; m <- length(exos)+1
  rho <- minnesota_par$rho
  lambdas <- minnesota_par$lambdas

  par.samples <- vector("list", length = N)

  for(i in seq_len(N)) {
    beta <- NW_beta(Sigma, X, Y, rho=rho, lambdas=lambdas)
    Sigma <- NW_Sigma(beta, X, Y)
    par.samples[[i]] <- list(Beta=Coef_vec2mat(beta,k,p,m), Sigma=Sigma)
  }

  par.samples <- par.samples[(warmup+1):N]
  fitted.samples <- imap(par.samples, ~ (X%*% .x$Beta) %>%
                           as_tibble(.name_repair = "minimal") %>%
                           set_names(colnames(Y)) %>%
                           mutate(sample=.y,
                                  time=index(Y))) %>%
    do.call(rbind,.)

  interval <- fitted.samples %>%
    pivot_longer(cols=!one_of("time","sample"), names_to="var") %>%
    group_by(time,var) %>% nest() %>%
    mutate(lwr=map_dbl(data, ~quantile(.x[["value"]], probs=0.025)),
           med=map_dbl(data, ~quantile(.x[["value"]], probs=0.5)),
           upr=map_dbl(data, ~quantile(.x[["value"]], probs=0.975))) %>%
    dplyr::select(-data) %>% ungroup()

  fitted <- interval %>%
    dplyr::select(time,var,med) %>%
    pivot_wider(names_from=var,values_from=med) %>%
    column_to_rownames("time") %>%
    as.xts(order.by=as.yearqtr(rownames(.)))

  model <- list(fitted=fitted,
                posterior=list(parameter=par.samples,predictive=fitted.samples),
                interval=interval,
                data=data,
                model_var=list(Y=Y,X=X))
  class(model) <- "bayesVAR"
  return(model)
}


#' Figures from Bayesian VAR model
#
#' provide figures from Bayesian VAR methods
#' @param model Bayesian VAR model
#' @param ncol.fig a integer, number of figures plotted in a row
#' @export
plot.bayesVAR <- function(model, ncol.fig=2){

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

#' Posterior predictive distributions from Bayesian VAR model
#
#' sampling posterior predictive distributions from the Bayesian VAR model
#' @param model Bayesian VAR model
#' @param newdata a data of xts class, exogenous variables for the prediction period
#' @param condition a data of xts class, endogenous variables to be conditioned
#' @param type plotting for goodness-of-fitted or forecasting
#' @export
predict.bayesVAR <- function(model, newdata, condition, type=c("unconditional","conditional")){

  if(missing(condition) & missing(type)){ condition <- xts() ; type <- "unconditional" }
  if(!missing(condition) & missing(type)){ type <- "conditional" }
  if(!is.xts(newdata) | !is.xts(condition)) stop("input data format must be a xts class")

  Y <- model$model_var$Y
  X <- model$model_var$X
  k <- ncol(Y)
  p <- (ncol(X)-ncol(newdata)-1)/ncol(Y)

  par <- model$posterior$parameter

  if(type=="unconditional"){
    predictive <- imap(par, ~(predict_draw(Y, newdata=newdata, Coef_mat=.x$Beta, Sigma=.x$Sigma)) %>%
                         as.data.frame() %>%
                         rownames_to_column("time") %>%
                         mutate(sample=.y)) %>%
      do.call(rbind,.)

  } else if(type=="conditional"){
    predictive <- imap(par, ~(predict_cond(Y, newdata, condition, Coef_mat=.x$Beta, Sigma=.x$Sigma)) %>%
                         as.data.frame() %>%
                         rownames_to_column("time") %>%
                         mutate(sample=.y)) %>%
      do.call(rbind,.)
  }

  # predict_cond(Y, newdata, condition, Coef_mat=par[[1]]$Beta, Sigma=par[[1]]$Sigma)
  # Coef_mat=par[[1]]$Beta ;  Sigma=par[[1]]$Sigma
  interval <- predictive %>%
    pivot_longer(cols=!one_of("time","sample"), names_to="var") %>%
    group_by(time,var) %>% nest() %>%
    mutate(lwr=map_dbl(data, ~quantile(.x[["value"]], probs=0.025)),
           med=map_dbl(data, ~quantile(.x[["value"]], probs=0.5)),
           upr=map_dbl(data, ~quantile(.x[["value"]], probs=0.975))) %>%
    dplyr::select(-data) %>% ungroup()

  fitted <- interval %>%
    dplyr::select(time,var,med) %>%
    pivot_wider(names_from=var,values_from=med) %>%
    column_to_rownames("time") %>%
    as.xts(order.by=as.yearqtr(rownames(.)))

  Fred <- list(fitted=fitted,
               posterior=predictive,
               interval=interval)
  class(Fred) <- "bayesVAR"
  return(Fred)
}



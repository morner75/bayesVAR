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

  par.samples <- list(Coef_mat=vector("list", length=N), Sigma=vector("list", length=N))
  for(i in seq_len(N)) {
    Coef_mat <- NW_beta(Sigma, X, Y, rho=rho, lambdas=lambdas,p=p,m=m) |> Coef_vec2mat(k,p,m)
    Sigma <- NW_Sigma(beta, X, Y)
    par.samples$Coef_mat[[i]] <- Coef_mat
    par.samples$Sigma[[i]] <- Sigma
  }

  par.samples <- lapply(par.samples, function(lst) lst[(warmup+1):N])
  fitted.samples <- lapply(seq_len(N-warmup), function(i) X%*% par.samples$Coef_mat[[i]] |>
                                                          as.xts(order.by=index(Y)) |> setNames(colnames(Y)))

  fitted <- Reduce(`+`, fitted.samples)/(N-warmup)

  interval <- array(unlist(fitted.samples),dim=c(nrow(fitted.samples[[1]]),ncol(fitted.samples[[1]]),N-warmup)) |>
                apply(c(1,2),quantile,probs=c(0.025,0.5,0.975)) |> purrr::array_branch(margin=1)
  interval <- lapply(interval, function(lst) xts(lst, order.by=index(Y)) |> setNames(colnames(Y))) |>
                setNames(c("lower","median","upper"))

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
#' @param x Bayesian VAR model
#' @param y NULL
#' @param ncol.fig a integer, number of figures plotted in a row
#' @export
plot.bayesVAR <- function(x, y=NULL, ncol.fig=2){

  data <- purrr::map_dfr(x$interval,
          ~(dplyr::mutate(tibble::as_tibble(.x),time=as.yearqtr(index(x$interval[[1]])),.before=1)),.id="type") |>
          tidyr::pivot_longer(!one_of("time","type"),names_to="var") |>
          tidyr::pivot_wider(names_from="type",values_from="value")


  if(!is.null(x$model_var)){
    Y <- x$model_var$Y
    actual <- tibble::rownames_to_column(as.data.frame(Y), "time") |>
      dplyr::mutate(time=as.yearqtr(time)) |>
      tidyr::pivot_longer(-time,"var",values_to="actual")
    data <- dplyr::left_join(data,actual,by=c("time","var"))
  }

  ggplot(data=data,aes(time,col=var,group=1))+
    geom_line(aes(y=median))+
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=var),alpha=0.3)+
    facet_wrap(facets=vars(var),ncol=ncol.fig, scales="free_y")+
    theme_bw()+labs(y="")+
    theme(legend.position = "bottom")+
    {if(!is.null(x$model_var)) geom_point(aes(y=actual))}

}

#' Posterior predictive distributions from Bayesian VAR model
#
#' sampling posterior predictive distributions from the Bayesian VAR model
#' @param object Bayesian VAR model
#' @param newdata a data of xts class, exogenous variables for the prediction period
#' @param condition a data of xts class, endogenous variables to be conditioned
#' @param type plotting for goodness-of-fitted or forecasting
#' @export
predict.bayesVAR <- function(object, newdata, condition, type=c("unconditional","conditional")){

  if(missing(condition) & missing(type)){ condition <- xts() ; type <- "unconditional" }
  if(!missing(condition) & missing(type)){ type <- "conditional" }
  if(!is.xts(newdata) | !is.xts(condition)) stop("input data format must be a xts class")

  Y <- object$model_var$Y
  X <- object$model_var$X
  k <- ncol(Y)
  p <- (ncol(X)-ncol(newdata)-1)/ncol(Y)

  par <- object$posterior$parameter

  if(type=="unconditional"){
    predictive <- lapply(seq_along(par[[1]]), function(i)
                          predict_draw(Y, newdata=newdata, Coef_mat=par$Coef_mat[[i]], Sigma=par$Sigma[[i]]))


  } else if(type=="conditional"){
    predictive <- lapply(seq_along(par[[1]]), function(i)
                        predict_cond(Y, newdata, condition, Coef_mat=par$Coef_mat[[i]], Sigma=par$Sigma[[i]]))

  }

  fitted <- Reduce(`+`, predictive)/length(par$Coef_mat)

  interval <- array(unlist(predictive),dim=c(nrow(predictive[[1]]),ncol(predictive[[1]]),length(par$Coef_mat))) |>
    apply(c(1,2),quantile,probs=c(0.025,0.5,0.975)) |> purrr::array_branch(margin=1)
  interval <- lapply(interval, function(lst) xts(lst, order.by=index(fitted)) |> setNames(colnames(Y))) |>
    setNames(c("lower","median","upper"))


  Fred <- list(fitted=fitted,
               posterior=predictive,
               interval=interval)
  class(Fred) <- "bayesVAR"
  return(Fred)
}



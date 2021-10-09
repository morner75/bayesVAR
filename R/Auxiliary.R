# S3 generic functions -------------------------------------------------------------------

plot <- function(model,...) UseMethod("plot")
plot.default <- function(model,...) base::plot(model)

predict <- function(model,...) UseMethod("predict")
predict.default <- function(model,...) stats::predict.lm(model)


# helper functions between different forms of coefficient values -------------------------

Coef_mat2list <- function(Coef_mat,m){
  k=ncol(Coef_mat) ; p <- (nrow(Coef_mat)-m)/k

  res <- vector(mode="list",length=p+1)
  for(i in 1:p) res[[i]] <- t(Coef_mat)[1:k,(i-1)*k+1:k]
  res[[p+1]] <- t(Coef_mat)[1:k,p*k+1:m,drop=FALSE]
  res
}

Coef_list2mat <- function(Coef_list) do.call(cbind,Coef_list) %>% t()

Coef_vec2mat <- function(Coef_vec,k,p,m){
  n=k*p+m
  matrix(Coef_vec,nrow=n,ncol=k)
}

Coef_vec2list<- function(Coef_vec,k,p,m){
  Coef_mat2list(Coef_vec2mat(Coef_vec,k,p,m),m)
}


# Impulse Response matrix generator
IR_mat_generator <- function(period, p, Coef_mat, Sigma, type=c("origin","chol","triangle")){
  # type="origin"
  if(is.null(colnames(Coef_mat))) colnames(Coef_mat) <- str_c("var",seq_len(ncol(Coef_mat)))
  k <- ncol(Coef_mat)
  res <- vector("list",k)

  for(i in seq_len(k)){
    Ynew <- matrix(0,ncol=k,nrow=p)
    Ynew[p,i] <- 1
    for(j in seq_len(period)){
      Y_update <- matrix(t(tail(Ynew,n=p)[rev(seq_len(p)),]),nrow=1,byrow=TRUE) %*% Coef_mat[1:(k*p),]
      Ynew <- rbind(Ynew,Y_update)
    }

    res[[i]] <- as_tibble(Ynew[p:(p+period),]) %>% mutate(term=0:period,.before=1)
  }
  names(res) <-  colnames(Coef_mat)
  res1 <- imap(res, ~mutate(.x, impulse=.y,.before=1) %>%
                 pivot_longer(!one_of("term","impulse"),names_to = "response"))  %>%
    do.call(rbind,.) %>%
    arrange(term,factor(response,levels = colnames(Coef_mat))) %>%
    pivot_wider(names_from="impulse",values_from="value") %>%
    group_by(term) %>% nest() %>% ungroup()

  if(type=="origin") return(res1)

  res2 <- res1 %>% transmute(term=term, data=map(data, ~(as.matrix(.x[,-1])%*%t(chol(Sigma)) %>%
                                                           as_tibble()) %>%
                                                   set_names(colnames(Coef_mat)) %>%
                                                   mutate(response=colnames(Coef_mat),.before=1)))
  if(type=="chol") return(res2)

  res3 <- res1 %>% transmute(term=term, data=map(data, ~(as.matrix(.x[,-1])%*%t(chol(Sigma))%*%diag(1/sqrt(diag(Sigma))) %>%
                                                           as_tibble(.name_repair="minimal")) %>%
                                                   set_names(colnames(Coef_mat)) %>%
                                                   mutate(response=colnames(Coef_mat),.before=1)))
  if(type=="triangle") return(res3)
}



# forecasting -----------------------------------------------------------------

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

predict_draw <- function(Y, newdata,Coef_mat, Sigma){
  if(!is.xts(newdata)) stop("external variables must be xts class")
  period <- nrow(newdata)
  p <- (nrow(Coef_mat)-ncol(newdata)-1)/ncol(Y)
  Ynew <- tail(Y,n=p)
  x <- cbind(1,newdata)
  for(i in seq_len(period)){
    Y_update <- matrix(c(tail(Ynew,n=p) %>% rev.zoo() %>% t() %>% as.vector(),x[i,]) %*% Coef_mat +
                         MASS::mvrnorm(n=1, mu=rep(0,ncol(Y)), Sigma=Sigma),nrow=1) %>%
      xts(order.by=index(newdata[i,]))
    Ynew <- rbind(Ynew,Y_update)
  }
  Ynew
}



predict_cond <- function(Y, newdata, condition, Coef_mat, Sigma){

  period <- nrow(newdata)
  p <- (nrow(Coef_mat) - ncol(newdata)-1)/ncol(Coef_mat)
  cond.var <- colnames(condition)
  pred_wo_shock <- predict_simple(Y, newdata, Coef_mat) %>%  tail(n=period)
  r <-  map(seq_len(period), ~coredata(condition - pred_wo_shock[,cond.var])[.x,]) %>% do.call(c,.)

  IRmat <- IR_mat_generator(period=(period-1), p=(nrow(Coef_mat)-ncol(newdata)-1)/ncol(Y), Coef_mat=Coef_mat,Sigma=Sigma,type="chol")
  R <- transmute(IRmat, R0=map(data, ~(as_tibble(.x,.name_repair="minimal") %>%
                                         set_names(c("response",colnames(Y))) %>%
                                         mutate(response=colnames(Y))   %>%
                                         filter(response%in% cond.var) %>%
                                         dplyr::select(-response))))
  names <- str_c("R",1:(period-1))
  for(i in seq_along(names)) R <- R %>% mutate(!!sym(names[i]) := lag(R0,n=i))
  R <- unnest(R, cols = everything(), names_repair = "minimal") %>% zoo::na.fill(0)

  eta_bar <- t(R)%*%solve(R%*%t(R))%*%r
  Gamma_bar <- diag(dim(R)[2]) - t(R)%*%solve(R%*%t(R))%*%R
  eta_vec <- MASS::mvrnorm(1,eta_bar, Gamma_bar)

  # P_inv <- solve(diag(svd(R)$d))
  # V1 <- t(svd(R,nv=ncol(R))$v)[,1:nrow(R)]
  # V2 <- t(svd(R,nv=ncol(R))$v)[,(nrow(R)+1):ncol(R)]
  # U <- t(svd(R,nu=nrow(R))$u)
  # eta_vec <- V1%*%P_inv%*%U%*%r + V2%*%MASS::mvrnorm(1,rep(0,diff(dim(R))),diag(diff(dim(R))))
  eta_list <- map(seq_len(period), ~eta_vec[(1+(.x-1)*ncol(Y)):(.x*ncol(Y))])
  cond_shocks <- mutate(IRmat, data=map(data, ~dplyr::select(.x,-response) %>% as.matrix(),.name_repair="minimal"),
                        eta=eta_list,
                        eta_accum=accumulate(eta, function(x,y) append(x,y)),
                        IRmat_accum=accumulate(data, function(x,y) cbind(y,x)),
                        shocks= map2(IRmat_accum,eta_accum,~(.x %*%.y) %>% t())) %>%
    dplyr::select(shocks) %>%
    unnest(cols=c(shocks),names_repair = "minimal") %>%
    as.xts(order.by=index(condition)) %>% setNames(colnames(Y))

  rbind(tail(Y,n=p), cond_shocks + pred_wo_shock)
}


# helper functions ----------------------------------------------------------------------

rmmvnorm <- function(M,Sigma,Phi) {
  M+chol(Phi)%*%matrix(rnorm(nrow(Sigma)*nrow(Phi)),nrow(Phi),nrow(Sigma))%*%chol(Sigma)
}

ind_finder <- function(k,p,m,val1,val2){
  ind1 <- ind2 <- numeric()
  cp <- 0
  for(i in 1:k){
    new <- (0:(p-1))*k+i+cp
    new <- c(new,cp+k*p+1:m)
    ind1 <- c(ind1,new)
    ind2 <- c(ind2,tail(ind1,n=2))
    cp <- tail(ind2,1)
  }
  vec1 <- vec2 <- rep(1,k*(k*p+m))
  vec1[ind1] <- val1
  vec2[ind2] <- val2
  vec1*vec2
}

mat2Longtb <- function(mat, name= c("var1","var2","value")){
  as.data.frame(mat) |>
    tibble::rownames_to_column(var="var1") |>
    tidyr::pivot_longer(-var1,names_to="var2",values_to = "value") |>
    rlang::set_names(name)
}






# Minnesota Priors ----------------------------------------------------------------------

# beta0 generator for Minnesota prior

Minnesota_b0 <- function(k,p,m,rho=1){
  n <- k*p+m
  b0 <- rep(0,k*n)
  for(i in seq_len(k)) b0[i+n*(i-1)] <- rho
  b0
}

# Omega0 generator for Minnesota prior
Minnesota_Omega0 <- function(k,p,m,Sigma, lambdas=c(0.1, 0.5, 2, 100)){
  sigmas <- diag(Sigma)
  rep(sigmas,each=k*p+m)*rep(c(rep(1/sigmas,p),rep(1,m)),k)*lambdas[1]^2*lambdas[2]^2/
    (rep(rep(c(1:p,1) ,times=c(rep(k,p),m)),k)^lambdas[3])^2*
    ind_finder(k,p,m,1/lambdas[2]^2,lambdas[4]^2)
}



# sampler -------------------------------------------------------------------------------

# beta sampler
NW_beta <- function(Sigma, X, Y, rho=1, lambdas=c(0.1, 0.5, 2, 100)){
  k <- ncol(Y) ; p <- ncol(X)%/%ncol(Y) ; m <- ncol(X)%%ncol(Y)
  b0 <- Minnesota_b0(k,p,m,rho=rho)
  Omega0 <- Minnesota_Omega0(k,p,m,Sigma,lambdas=lambdas)
  Omega.bar <- solve(diag(1/Omega0)+kronecker(solve(Sigma),t(X)%*%X))
  beta.bar <- Omega.bar%*%(diag(1/Omega0)%*%b0+ kronecker(solve(Sigma),t(X))%*%as.vector(Y))
  MASS::mvrnorm(n=1,mu=beta.bar,Sigma=Omega.bar)
}

# Sigma sampler
NW_Sigma <- function(beta, X, Y){
  k <- ncol(Y) ; p <- ncol(X)%/%ncol(Y) ; m <- ncol(X)%%ncol(Y)
  S0 <-  diag(k)
  alpha0 <- k+2
  Beta <- Coef_vec2mat(beta,k,p,m)
  S.hat.inv <- solve(t(Y-X%*%Beta)%*%(Y-X%*%Beta)+S0)
  alpha.hat <- alpha0 + nrow(Y)
  rWishart(1,alpha.hat,S.hat.inv)[,,.drop=TRUE] %>% solve()
}



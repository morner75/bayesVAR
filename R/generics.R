# generic functions ---------------------------------------

fitted <- function(object, ...) UseMethod("fitted")
fitted.default <- function(object, ...) stats::fitted(object, ...)

coef <- function(object, ...) UseMethod("coef")
coef.default <- function(object, ...) stats::coef(object, ...)


plot <- function(x, y=NULL, ...) UseMethod("plot")
plot.default <- function(x, y=NULL, ...) plot.bayesVAR(x, y=NULL, ...)


summary <- function(object, ...) UseMethod("summary")
summary.default <- function(object, ...) base::summary(object, ...)

predict <- function(model,...) UseMethod("predict")
predict.default <- function(model,...) predict.bayes(model, ...)

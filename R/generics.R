# generic functions ---------------------------------------

fitted <- function(object, ...) UseMethod("fitted")
fitted.default <- function(object, ...) stats::fitted(object, ...)

coef <- function(object, ...) UseMethod("coef")
coef.default <- function(object, ...) stats::coef(object, ...)

plot <- function(x, y, ...) UseMethod("plot")
plot.default <- function(x, y, ...) graphics::plot(x, y, ...)

summary <- function(object, ...) UseMethod("summary")
summary.default <- function(object, ...) base::summary(object, ...)

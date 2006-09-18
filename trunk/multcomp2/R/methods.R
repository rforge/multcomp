
coef.glht <- function(object, null = FALSE, ...) 
{
    if (null) return(object$m)
    betahat <- coef(object$model)
    if (inherits(object$model, "lmer"))
        betahat <- coeflmer(object$model)
    drop(object$K %*% betahat)
}

vcov.glht <- function(object, ...) 
    object$K %*% tcrossprod(as.matrix(vcov(object$model)), object$K)

summary.glht <- function(object, test = adjusted(), ...) 
{
    pq <- test(object)
    if (class(pq) == "summary.glht.global") {
        object$gtest <- pq
        class(object) <- c("summary.glht.global", "glht")
        return(object)
    }
    object$mtests <- cbind(pq$coefficients, coef(object, null = TRUE), 
                           pq$sigma, pq$tstat, pq$pvalues)
    attr(object$mtests, "error") <- attr(pq$pvalues, "error")
    colnames(object$mtests) <- c("Estimate", "Hypothesis", "Std. Error",
        ifelse(object$df == 0, "z value", "t value"), "p value")
    attr(object$mtests, "type") <- pq$type
    class(object) <- c("summary.glht", "glht")
    return(object)
}

confint.glht <- function(object, parm, level = 0.95, ...) 
{
    pq <- pqglht(object)
    ci <- pq$qfunction(conf.level = level, ...)
    object$confint <- cbind(pq$coefficients, ci)
    colnames(object$confint) <- c("Estimate", "lwr", "upr")
    attr(object$confint, "conf.level") <- level
    attr(object$confint, "calpha") <- attr(ci, "calpha")
    attr(object$confint, "error") <- attr(ci, "error")
    class(object) <- c("confint.glht", "glht")
    return(object)
}

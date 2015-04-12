#' @include Symbol.R

#' @title Represent an effect symbol as defined by the model, extends Symbol
#'
#' @description
#'        The instance tracks the name and the coefficient of a
#'        single occurence of an effect symbol.
#'        In \eqn{ 3*(x2 - x1) - 4*(x3 -x1) == 0} the same effect
#'        x1 gets referenced multiple times. This can only be supported
#'        if the two instances do not share the same memory region,
#'        as it is the case with symbol nodes returned from \code{parse()}
#'
#' @method  is.Effect return true if applied to an instance of Effect
#' @rdname internal/Effect
setClass( Class = 'Effect',  contains = 'Symbol')


# Establish the 'is.Effect' method:

#' @title  Generic method is.Effect
#' @param  object  any object
#' @return True if the instance is of class Effect
#' @rdname internal/Effect
setGeneric('is.Effect', function(object) standardGeneric('is.Effect'))

#' @title  Tests if an object is an instance of class Effect
#' @param  object  Any object (except an instance of class Effect)
#' @return Always false
#' @rdname internal/Effect
setMethod('is.Effect', signature = 'ANY', definition = function(object) F )

#' @title  Tests if an object is an instance of class Effect
#' @param  object  Instance of class Effect
#' @return Always true
#' @rdname internal/Effect
setMethod('is.Effect', signature = 'Effect', definition = function(object) T )


#' @title  Effect constructor
#' @param  name   Effect name
#' @param  coef   Effect coefficient 
#' @return A new instance of class Effect
#' @rdname internal/Effect
Effect <- function(name, coef=1) new("Effect", name, coef)


# ^ Done with S4 methodic hassles for class Effect ...




#' @include Symbol.R

#' @title  Represent a numeric literal by a symbol, extends Symbol
#' @method is.Literal  
setClass( Class = 'Literal', contains = 'Symbol')


#' @title   Initialize method for instances of class Literal
#' @param  .Object Instance derived from class Literal
#' @param   value Numeric value, also serving as symbol name 
#' @return  object instance 
setMethod(f = 'initialize', signature = 'Literal', definition =  function( .Object, value ) {
   if ( ! is.numeric(value) || length(value) != 1 )
          stop('Literal::initialize',': parameter \'value\' = ', sQuote(value), ' ',
                                       'should of type numeric and of length 1')
   .Object@name <- as.character(value)
   .Object@coef <- value
    if ( validObject(.Object) )
         return(.Object) 
    NULL
 })


# Establish the 'is.Literal' methods:

#' @title  Generic method Is.Effect
#' @param  object  Any object
#' @return True if object is an instance of class Effect
setGeneric('is.Literal', function(object) standardGeneric('is.Literal'))

#' @title  Tests if object is an instance of class Literal
#' @param  object  Any object (except an instance derived from Literal)
#' @return Always false
setMethod('is.Literal', signature = 'ANY', definition = function(object) F )

#' @title   Tests if object is an instance of class Literal
#' @param   object  An instance derived from Literal
#' @return  Always true
setMethod('is.Literal', signature = 'Literal', definition = function(object) T )

#' @title   Literal constructor
#' @param   value Numeric value, also serving as symbol name 
#' @return  A new instance of class Literal
Literal <- function(value) new("Literal", value)

# ^ Done with S4 methodic hassles for class Literal ...


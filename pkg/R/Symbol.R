#' @title   Represent an instance of an external symbol, e.g. `pi`
#' @slot    name         Name of the symbol
#' @slot    coef         Coefficient of the symbol
#' @method  Symbol       Create the instance
#' @method  coef         Return the coefficient associated with the symbol
#' @method  coef<-       Set coefficient, return instance
#' @method  name         Return symbol name
#' @method  name<-       Set symbol name, return instance
#' @method  ==           Return true if  names of two instances are equal or the instance name matches a string
#' @method  is.Symbol    Return true if applied to an instance of \class{Symbol}
#' @method  as.character Same as name
#' @method  quoted       Return the symbol name enclosed in backticks
setClass( Class     = 'Symbol',
          slots     =  c( name = 'character', coef = 'numeric'),
          prototype =  list( name = character(0), coef = 1),
          validity  =  function(object) {
                       if  ( nchar(gsub('[[:blank:]]+','',object@name)) < 1 )
                             return('Instance should carry a non-blank name')
                       if  ( ! is.finite(object@coef) )
                             return('Instance should carry a finite coefficient')
                       return(TRUE)

          })

#' @title  Initialize method for instances of class Symbol
#' @param .Object instance of class Symbol
#' @param  name Symbol name
#' @param  coef Symbol coefficient
#' @return An instance of Symbol or NULL if initialization fails
setMethod(f = 'initialize', signature = 'Symbol', definition = function( .Object, name, coef ) {
   .Object@name <- name
   .Object@coef <- coef
    if ( validObject(.Object) ) 
         return(.Object)
    NULL
})


# Establish the 'is.Symbol' methods:

#' @title  Generic method is.Symbol
#' @param  object  Any object
#' @return True if object is an instance derived from class Symbol
setGeneric('is.Symbol', function(object) standardGeneric('is.Symbol'))

#' @title   Test if object is an instance derived from  class Symbol
#' @param   object Any object (except an instance derived from class Symbol) 
#' @return  Always false
setMethod('is.Symbol', signature = 'ANY', definition = function(object) F )

#' @title   Test if object is an instance derived from  class Symbol
#' @param   object An instance derived from class Symbol
#' @return  Always true
setMethod('is.Symbol', signature = 'Symbol', definition = function(object) T )


# Establish the 'coef' methods:

#' @title   Generic method coef
#' @param   object Any object
#' @param   ...  Other arguments prescribed by \code{\link[stats]{coef}}, ignored
#' @return  Coefficient value
setGeneric('coef', function(object) standardGeneric('coef'))

#' @title   Get the value of the coefficient associated with the symbol
#' @param   object  An instance derived from class Symbol
#' @param   ...  Other arguments prescribed by \code{\link[stats]{coef}}, ignored
#' @return  Value of the coefficient
setMethod('coef', signature = 'Symbol', definition = function(object, ...) object@coef )

#' @title   Generic method  coef<-
#' @param   object  Any object
#' @param   value   A scalar
#' @return  Instance of object
setGeneric('coef<-', function(object, value) standardGeneric('coef<-'))

#' @title   Sets the coefficient of the Symbol instance to a new value
#' @param   object   Instance of class Symbol
#' @param   value    A scalar
#' @return  Instance of object
setReplaceMethod('coef',  signature = 'Symbol', definition = function(object, value) { 
   object@coef <- value
   validObject(object)
   object 
})



# Establish the 'name' method:

#' @title  Generic method 'name'
#' @param   object Any object
#' @return  Name of object as a character string
setGeneric('name', function(object) standardGeneric('name'))

#' @title   Get the name of the Symbol instance
#' @param   object An instance derived from class Symbol
#' @return  Name of object as a character string
setMethod('name', signature = 'Symbol', definition = function(object) object@name )

#' @title   Get the name of the Symbol instance
#' @param   x  An instance derived from the Symbol
#' @return  Name of object as a character string (same as name method)
setMethod('as.character', signature = 'Symbol', definition = function(x) name(x) )

# Establish the 'quoted' method:

#' @title  Quote the name of the Symbol instance using backticks ('`')
#' @param  object Any object
#' @return Name of object as a quoted character string
setGeneric('quoted', function(object) standardGeneric('quoted'))

#' @title   Quote the name of the Symbol instance using backticks ('`')
#' @param   object An instance derived from class Symbol
#' @return  Name of object as a quoted character string
setMethod('quoted', signature = 'Symbol', definition = function(object) paste0('`',object@name,'`' ))


# Establish the '==' methods for Symbol:

# A generic method called '==' is already known

#' @title  Compare the name of a Symbol instance against the name of another instance of class Symbol
#' @param  e1  Instance of Symbol
#' @param  e2  Instance of Symbol
#' @return true if both names are equal
setMethod('==', signature = c(e1 = 'Symbol', e2 = 'Symbol'), definition = function(e1, e2) {
   name(e1) == name(e2)  
})

#' @title  Compares the name of an instance derived from class Symbol against a character string
#' @param  e1  Instance of Symbol
#' @param  e2  A character string
#' @return True if name for instance is equal to string e2
setMethod('==', signature = c(e1 = 'Symbol', e2 = 'character'), definition = function(e1, e2) {
   name(e1) == e2
})

#' @title  Symbol constructor
#' @param  name  Name of the Symbol
#' @param  coef  Symbol coefficient 
#' @return A new instance of class Symbol
Symbol <- function(name, coef=1) new("Symbol", name, coef)


# ^ Done with S4 methodic hassles for class Symbol ...


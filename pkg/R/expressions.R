
# $Id$

#' @title   \code{expression2coef} is an interpreter for a single null hypothesis
#' @description
#'         Extract the effect coefficients from an equation or inequation describing a single null hypothesis and
#'         also provide the rhs value.
#
#' @param  ex           R expression describing the null hypothesis
#' @param  vars         Character vector of effect names
#' @param  sh0          Label of h0 as in ex = c('A greater B' = ' A - B >= 0 ' )
#' @param  verbose      Trace parser states and actions, if true
#' @param  zerocheck    Check for effects being near zero in terms of machine precision
#'
#' @return A list with the following elements:
#'         \itemize{
#'              \item  coef          Vector of effect coefficients
#'              \item  names         Vector of effect names
#'              \item  m             Right-hand side value (a numeric value of length 1)
#'              \item  alternative   One of: 'less', 'greater' or 'two.sided'
#'              \item  lhs           Left-hand side terms
#'              \item  rhs           Right-hand side terms
#'         }
#'
#' @examples 
#' e1 <- parse(text = ' (x2 - x1) >= 1 + sin(pi/4) * ( 2*x3 - x4 ) ')
#' e2 <- parse(text = ' (x2 - x1) - sin(pi/4) * ( 2*x3 - x4 ) -1 >= 0 ')
#' 
#' r1 <- multcomp:::expression2coef( e1, c('x1', 'x2', 'x3', 'x4'), 'h0.1a' )
#' r2 <- multcomp:::expression2coef( e2, c('x1', 'x2', 'x3', 'x4'), 'h0.1b' )
#' 
#' stopifnot( identical( r1$coef, r2$coef ) )
#' stopifnot( identical( r1$m, r2$m ) )
#'
#' @note    Implementation needs package \pkg{codetools} and also depends on S4 classes defined in Symbol.R, Effect.R and Literal.R
#' @seealso chrlinfct2matrix
#' @note
#'          Starting from the results of a call to \code{parse(text='some equation')}, 
#'          the function \code{expression2coef} gets called with an object of type expression passed 
#'          in its parameter `ex`. 
#'
#'          Expressions in R are internally based on a tree structure, which can be mapped to a list of lists.
#' 
#'          Each entry of these recursive lists is then going to be inspected by the interaction 
#'          of a pair of very concisely written functions:  \code{makeCodeWalker(..., handler, call, leaf)} 
#'          and \code{walkCode(e, w = makeCodeWalker())}  with `e` representing the  symbol node 
#'          to be processed. Both can be found in package \pkg{codetools}. The walkCode function 
#'          itself spans only a few lines:
#'          \preformatted{
#'            walkCode <- 
#'            function ( e, w = makeCodeWalker() ) {
#'               if ( typeof(e) == 'language' ) {
#'                    if ( typeof([[e1]]) %in% c( 'symbol', 'character' )) 
#'                    {
#
# (roxygen2 bug: need to move the curly 
#  brace from the end of the previous line 
#  to avoid an error saying unmatched 
#  braces ... within a block tagged as 'preformatted' ...)
#
#'                       h <- w$handler(as.character(e[[1]]), w)
#'                       if (!is.null(h) ) 
#'                            h(e, w) 
#'                       else w$call(e, w)
#'                    } else  w$call(e, w) 
#'               } else w$leaf(e, w)
#'            }
#'          }
#'          There is also the \code{\link[codetools]{showTree}} utility which can print 
#'          out something very close to a polish notation of any expression by 
#'          simply using makeCodeWalker's default implementation.
#'
#'          Standing on the fundament laid by the codetools framework, a call to the 
#'          core of an interpreter for arbitrary expressions will become 
#'          straightforward: 
#'          \preformatted{
#'          result <-
#'          walkCode( e = parse( text = ' (a - x ) * 2 >= 0 ' )[[1]], 
#'                    w = makeCodeWalker( some helper functions ..., 
#'                                        handler = function (v, w) { (...)}, 
#'                                        call    = function (v, w) { (...)}, 
#'                                        leaf    = function (v, w) { (...)}))
#'          }
#'          The symbol tree  for \eqn{ ( a - x ) * 2 >= 0} would  look like this:
#'          \preformatted{
#'                      >=        <- node is addressed by as.list(exp[[1]])[1]
#'                     /  \
#'                    *    0      <- nodes are addressed by as.list(exp[[1]])[2:3]
#'                   / \
#'                  -   2
#'                 / \
#'                a   x
#'          }          
#'          So any node can be seen as the head of a recursive list. exp[[1]] would
#'          address the top level node named `>=` as well as all nodes attached to it.
#'          Symbol names provide the node labels. Terminal symbols like `a`, `x`, 2 and 0 
#'          are called leafs. 
#'
#'          With the top node being represented by the symbol `>=`, \code{walkCode(`>=`, w)} 
#'          would recursively visit all descendants of this node, with the first of the two subtrees 
#'          being the left-hand and the right-hand side of the inequation.
#'
#'          A node with leafs attached represents an operator. If \code{walkCode} sees such a node,
#'          it calls the top level function \code{handler(v, w)} with `v` being the symbol 
#'          represented by the node, say `*`.
#'
#'          \code{handler(`*`, w)} then dispatches a function specialised on that particular 
#'          operator, which is \code{mul(v, w)} in our case. So \code{walkCode} calls 
#'          \code{h(e, w)} with `h` actually being the function \code{mul(v , w)} and 
#'          `e` or `v`  being  the symbol node `*` having  attached descendant nodes represented by 
#'          the symbol nodes `-` and  the literal 2.
#'          
#'          From within \code{sub(`-`, w)}, \code{walkCode} would then recursively get
#'          called on the two leafs `a` and  `x`. Since both are terminal symbols, 
#'          \code{walkCode}  would immediately call \code{leaf(e, w)} for each of them.
#'        
#'          \code{leaf(e , w)} would return e, but associated with a coefficient of 1.
#'
#'          If \code{handler(v, w)} returns null for a particular operator or function name, 
#'          the node and its subtree gets passed to \code{call(v, w)}, with the idea to transfer
#'          control to some external functions able to handle that node.
#'
#'          Now, the other cornerstone of this implementation consists in binding a coefficient 
#'          to each node of the parse tree. For  \eqn{(a - x) * 2  >= 0}, this would look like 
#'          a decorated christmas tree:
#'          \preformatted{
#'                             >= ~ coef 1
#'                            /  \
#'                  coef 1 ~ *    0 ~ coef 1
#'                          / \
#'                coef 1 ~ -   2 ~ coef 1
#'                        / \
#'              coef 1 ~ a   x ~ coef  1
#'          }
#'          Again, \code{walkCode(top, w )} will start by calling the \code{handler(v, w)}, 
#'          which will return a call to \code{eqn(`>=`, w)}
#'
#'          \code{eqn(`>=`, w)}  in turn calls \code{walkCode} for the lhs and rhs side. 
#'          On lhs, the next symbol will be `*`, so \code{walkCode} calls \code{handler(`*`, w)}, 
#'          which returns \code{mul(v, w)}
#'
#'          When called, this function first runs \code{walkCode} again on all descendants of 
#'          the `*` node . So the next operation will be \code{sub(`-`, w)}. Itself recursively 
#'          calling \code{walkCode} on her own descandant will result in two calls to 
#'          \code{leaf(v, w)} with `v` being `a` and `x`
#'
#'          \code{leaf(v, w)} in turn is defined as a function which somehow associates 
#'          a coefficient of 1 to each symbol `v` before returning it.
#'
#'          Thus the descendants of the `*` node will be equivalent to the term \eqn{(a - x) * 2},
#'          representing the whole left-hand side of the inequation in our example.
#'
#'          \code{mul(`*`, w)} will then be able to compute the lhs effect coefficents `a` and 
#'          `x` by running two multiplications:
#'          \preformatted{
#'              coef(a) <- coef(a) * coef(2)
#'                         |         |
#'                         1         2
#'
#'              coef(x) <- coef(x) * coef(2)
#'                         |         |
#'                        -1         2
#'          }
#'          The coefficents of `a`  and `x` had been computed before when \code{mul(v, w)} walked its two
#'          subtrees.  \code{walkCode} saw the `-` node and thus got the instruction from
#'          \code{handler(`-`, w)} to call \code{sub(`-`, w)}, which in turn got the two 
#'          nodes `a` and `x` and set their coefficients that way:
#'          \preformatted{
#'              coef(a) <- coef(a)
#'                         |
#'                         1
#'
#'              coef(x) <- coef(x) * - 1
#'                         |
#'                         1
#'          }
#'          So far the concept came along as an elegant solution. But there was a hidden problem 
#'          in the implementation which didn't manifest itself simply because previously constructs 
#'          like \eqn{ x - 2 * x == 0} were recognized and rejected. Therein, the effect `x` 
#'          appears twice in different contexts. Supporting this by simply adding up the coefs 
#'          of all occurences of  `x` as the last processing step appeared to be straightforward, 
#'          except that it really wasn't that simple.
#'
#'          For keeping track of coefs, coefficients were previously linked as an additional 
#'          attribute to the symbol nodes of the tree returned by \code{\link[base]{parse}}. This 
#'          is not very different from instances of S4 classes carrying around their slot values as 
#'          attributes named after the slot names. 
#'
#'          Linking the coefficients as attributes to the variables and expressions they actually 
#'          belong to appeared to be attractive as this spared the effort of inventing 
#'          an explicit mechanism for keeping them in sync.
#'
#'          However, a symbol like `x` is only a name refering operations on it to a certain 
#'          memory region where its value gets stored. This region does not change when the 
#'          reference associated with that region gets copied along by program actions. 
#'          But within the same environment, all references carrying the same name also map 
#'          to the same address, much like copied pointers or references in \verb{C} and 
#'          \verb{C++} point or refer to the very same memory address.
#'
#'          In the example above, the information hold by the \code{coef} attribute would be +1
#'          for the first `x1`. As soon as the next occurence of `x1` gets processed its 
#'          coefficient would be set to -2. But then both occurences of `x1` are being associated 
#'          with a coefficent of -2 and hence the sum of over all coefficients of both occurences of 
#'          `x1` would become -4 instead of -1.
#' 
#'          The magic behind had already been told: while the attribute operation was performed on 
#'          two different instances of the reference, both redirected the attribute operation to 
#'          the very same memory area.
#'
#'          The solution was to create a copy for every instance of all symbols by means of 
#'          creating an instance of the \linkS4class{Effect} and \linkS4class{Symbol} classes
#'          in order to support all these corner cases correctly.
#'
#'          Having solved all this, a quickly written test case \eqn{ 2*(x - 1) +3 == 0} uncovered
#'          yet another problem. With `x` being an effect, \eqn{ 2*(x - 1)} the
#'          \code{add(v, w)} and \code{sub(v, w)} function had to return a list of symbol 
#'          nodes, containing both effects and literal elements at the same time.  Hence the class
#'          \linkS4class{Literal} was introduced to provide placeholder also for 
#'          additive constants as they can not be folded together with effect symbols. 
#'
#'          \code{transform} and \code{merge} utilities running below the equation handler
#'          \code{eqn(v, w)} are thus finally able to add up all these coefficients, transfer the 
#'          result to unique instances of the effect variable representing the left-hand side, and also 
#'          sum up all literals and return the result as right-hand side value.
#' 
#'          Todo: Instances of class 'Literal' aren't really needed any more since coefficients of 
#'          numeric literals are immediately folded. Yet leaving these instances in place for some time
#'          will make any related problem easily recognizable from trace output.
#'
#' @rdname internal/expression2coef
expression2coef <- function(ex, vars, sh0, verbose = F, zerocheck = T) {

   if ( verbose ) {
        message('multcomp::expression2coef(', sh0, '): ',
                'effect names are ', paste0(sQuote(vars), collapse = ', '))
        message('multcomp::expression2coef(', sh0, '): length of ex is ', length(ex))
   }

   # used in verbose mode to indent trace
   indent <- 0


   # used in context of try(...) to retrieve the exception message
   condition <- function(x) attr(x, 'condition')$message

   w <- codetools::makeCodeWalker(  #' @title       determine if a variable name matches an effect name as defined by `vars`
                                    #' @param x     an expression which can be coerced to a string of length 1
                                    #' @param w     environment created by \code{makeCodeWalker}
                                    #' @return      true or false
                                    #' @note        `vars` is a list of effect names passed as a parameter to the outmost function.
                                    is.effect = function(x, w) {

                                       if ( length(x)  != 1 )  
                                             w$bug('is.effect', sQuote(x) , ' has length != 1. Length is ', length(x) ) 

                                        w$msg(w,'is.effect', 'symbol name tested for naming an effect is ', sQuote(as.character(x)))

                                        as.character(x) %in% vars
                                    },

                                    #' @title    Dispatch a suitable handler for registered operators
                                    #' @param v  Operator as character string
                                    #' @param w  Environment created by \code{makeCodeWalker}
                                    #' @return   A list of lists
                                    #' @note     Called directly from within the \code{walkCode} framework
                                    handler = function(v, w) {
                                              w$enter('handler', v, w)

                                              res <- switch( v,
                                                            '-'   = w$sub,
                                                            '+'   = w$add,
                                                            '*'   = w$mul,
                                                            '/'   = w$div,
                                                            '('   = w$exp,
                                                            ':'   = w$ita,
                                                           '<='   = w$leq,
                                                           '>='   = w$geq,
                                                           '=='   = w$eqn,
                                                           '='    = w$eqn,
                                                            NULL )  # -> use default action w$call
                                                            
                                              w$leave('handler', v, v, w)
                                              res
                                    }, # end handler

                                    ## begin of function definitions called directly from within \code{w$handler(v, w)}

                                    #' @title  Support an equation 
                                    #' @param  v  Symbol node  '==' lhs rhs or '=' lhs rhs
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return A list of lists
                                    eqn  = function(v, w) {
                                          w$enter('eqn', v, w)

                                          if ( length(as.list(v)) != 3 ) {
                                               w$fatal('eqn', 'expression ', sQuote(deparse(v)), ' is incomplete')
                                          }

                                          operator <- as.character(v[[1]])
                                          eqn.lhs  <- paste0(deparse(v[[2]], 500), collapse='')
                                          eqn.rhs  <- paste0(deparse(v[[3]], 500), collapse='')
                                          eqn.bth  <- paste(eqn.lhs, operator, eqn.rhs)

                                          w$msg(w, 'eqn', 'equation is ', sQuote(eqn.lhs), ' ', sQuote(operator), ' ', sQuote(eqn.rhs))

                                          sym.lhs <- c()
                                          for ( e in as.list(v)[2]) {

                                                w$msg(w, 'eqn', 'walking lhs term ', sQuote(deparse(e)))

                                                r <- walkCode(e, w)

                                                # catch the corner case of an external symbol standing alone on lhs
                                                if ( is.Symbol(r) && !is.Effect(r) && !is.Literal(r) )
                                                     r <- w$eval(r, w)

                                                sym.lhs <- c( sym.lhs, r )
                                          }


                                          sym.rhs <- c()
                                          for ( e in as.list(v)[3]) {

                                                w$msg(w, 'eqn', 'walking rhs term ', sQuote(deparse(e)))

                                                r <- walkCode(e, w)

                                                # catch the corner case of an external symbol standing alone on rhs
                                                if ( is.Symbol(r) && !is.Effect(r) && !is.Literal(r) )
                                                     r <- w$eval(r, w)

                                                sym.rhs <- c( sym.rhs, r )

                                          }


                                          # shift all effects to lhs and the sum of all offsets to rhs
                                          trx <- w$transform(v, w, sym.lhs, sym.rhs)

                                          if ( length(trx$lhs) == 0 || !any( unlist(lapply( trx$lhs,is.Effect ) ) ) ) {
                                               w$fatal('eqn', 'equation ', sQuote(eqn.bth), ' does not name any effect')
                                          }


                                          if ( !all(unlist(lapply(trx$lhs, is.Effect)) ) ) { 
                                                w$bug('eqn', 'the normalisation of ', sQuote(eqn.bth), ' was not successful. ',
                                                             'left-hand side is [', w$dump(trx$lhs, w), '], ',
                                                             'while right-hand side is ', w$dump(trx$rhs, w))
                                          }

                                 
                                          zero <- c()
                                          for ( e in trx$lhs ) {
                                                if ( zerocheck && abs(coef(e)) <= .Machine$double.eps)  {
                                                     zero <- c( zero, e)
                                                }
                                          }

                                          if ( length(zero) > 0 ) {
                                               w$fatal('eqn', 'within expression ', sQuote(eqn.bth), ' ',
                                                              'the coefficients of the following effects ',
                                                              'evaluate to zero (or near zero): ',
                                                               w$enum(zero), '. Values are ',
                                                               paste(sprintf("%g",lapply(zero, w$getCoef)), collapse = ', '), '. ',
                                                              'Calling glht with the argument zerocheck=F will disable this check')
                                          }

                                          w$leave('eqn', list( symbols  = list( lhs = trx$lhs,
                                                                                rhs = trx$rhs),
                                                               equation = list( opr = operator,
                                                                                alt = w$as.side(operator, w),
                                                                                lhs = eqn.lhs,
                                                                                rhs = eqn.rhs,
                                                                                bth = eqn.bth)), v, w)
                                    }, # end eqn

                                    #' @title  Support inequation less or equal
                                    #' @param  v  Symbol node  '<=' lhs rhs
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return A list of lists
                                    leq = function(v, w) {
                                          w$enter('leq', v, w)
                                          w$leave('leq', w$eqn(v, w), v, w)
                                    },# end leq
                                    
                                    #' @title  Support inequation greater or equal
                                    #' @param  v  Symbol node  '=>' lhs rhs
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return A list of lists
                                    geq = function(v, w) {
                                          w$enter('geq', v, w)
                                          w$leave('geq', w$eqn(v, w), v, w)
                                    },# end geq

                                    #' @title  Support the addition of constants and effects 
                                    #' @param  v  Symbol node  '+' a b  or '+' a
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return A scalar or a list of symbol instances derived from \linkS4class{Symbol} 
                                    add = function(v, w) {
                                          w$enter('add', v, w)

                                          res <- c()
                                          for ( e in as.list(v)[-1] ) {
                                                w$msg(w, 'add', 'walking term ', sQuote(deparse(e)) )
                                                res <- c(res, walkCode(e, w))
                                          }

                                          w$msg(w, 'sub', 'expanded literals  ', w$dump(res, w))

                                          symbols <- c()
                                          sum     <- 0
                                          for ( e in res ) {
                                                w$msg(w,'sub', 'handling term ', w$dump(e, w), ' from ', sQuote(deparse(e)))
                                                switch(class(e),
                                                       'Effect'  = symbols <- c(symbols, e),
                                                       'Literal' = sum     <- sum + coef(e),
                                                       'Symbol'  = sum     <- sum + w$eval(e,w),
                                                       'numeric' = sum     <- sum + e,
                                                       w$bug('add', 'don\'t know how to handle ', w$enum(e), ' ',
                                                                    'within expression ', sQuote(deparse(v))))
                                          }

                                          w$msg(w, 'add', 'after sum ', sQuote(sum) )

                                          # just return constants as a single number if there are no effects
                                          if ( length(symbols) == 0 ) {
                                               return( w$leave('add', sum * w$getCoef(v), v, w))
                                          }

                                          # if we have symbols and a non-zero sum turn the sum into a symbol
                                          if ( sum ) {
                                               symbols <- c(symbols, Literal(sum))
                                          }

                                          w$leave('add', lapply(symbols,function(x,v,w) w$setCoef(x, w$getCoef(x)*w$getCoef(v)), v=v, w=w),v,w)
                                    }, # end add a b

                                    #' @title  Support negation or subtraction of effects and constants
                                    #' @param  v  Symbol node  '-' a b  |  '-' a NULL
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return A scalar or a list of symbol instances derived from \linkS4class{Symbol} 
                                    sub = function(v, w) {
                                          w$enter('sub', v, w)

                                          minuend    <- as.list(v)[2]
                                          subtrahend <- as.list(v)[3]
                                          uminus     <- length(as.list(v)) == 2 # one operator + one operand

                                          if ( uminus )
                                               subtrahend <- NULL

                                          w$msg(w, 'sub', 'minuend is ', w$enum(minuend))
                                          w$msg(w, 'sub', 'subtrahend is ', w$enum(subtrahend))
                                          w$msg(w, 'sub', 'uminus is ', uminus)


                                          minuends <- c()
                                          for ( e in minuend ) {
                                                w$msg(w, 'sub', 'walking minuend term ', sQuote(deparse(e)))
                                                minuends <- c(minuends, walkCode(e, w))
                                          }

                                          subtrahends <- c()
                                          for ( e in subtrahend ) {
                                                w$msg(w, 'sub', 'walking subtrahend term ', sQuote(deparse(e)))
                                                subtrahends <- c(subtrahends, walkCode(e, w))
                                          }

                                          
                                          sign <- ifelse(uminus, -1, 1)
                                          minuends.signed <- c()
                                          for ( e in minuends ) {
                                                minuends.signed <- c( minuends.signed, w$setCoef(e, w$getCoef(e) * sign * w$getCoef(v)))
                                          }


                                          subtrahends.signed <- c()
                                          for ( e in subtrahends ) {
                                                subtrahends.signed <- c( subtrahends.signed, w$setCoef(e, w$getCoef(e) * -w$getCoef(v)))
                                          }

                                          minuends    <- minuends.signed
                                          subtrahends <- subtrahends.signed


                                          sum     <- 0
                                          symbols <- c()

                                          # fold constants and separate them from symbols
                                          for ( e in c(minuends, subtrahends) ) {
                                                w$msg(w,'sub', 'handling term ', w$dump(e, w), ' from ', sQuote(deparse(e)))
                                                switch( class(e),
                                                        'Effect'  = symbols <- c(symbols,e),
                                                        'Symbol'  = sum     <- sum + w$eval(e,w),
                                                        'Literal' = sum     <- sum + coef(e),
                                                        'numeric' = sum     <- sum + e,
                                                         w$bug('sub', 'don\'t know how to handle ', w$enum(e), ' ',
                                                                      'within expression ', sQuote(deparse(v))))
                                          }

                                          # return a single number (coef(v) was already factored in)
                                          if ( length(symbols) == 0 ) {
                                               return(w$leave('sub', sum, v, w))
                                          }

                                          # if we have symbols and also a non-zero sum turn the sum into \linkS4class{Literal}
                                          if ( sum  ) {
                                               symbols <- c(symbols, Literal(sum))
                                          }

                                          # coef(v) was already factored in 
                                          w$leave('sub', symbols, v,w)
                                    }, # end sub a b

                                    #' @title  Support multiplication of effects or constants by constants
                                    #' @param  v  Symbol node  '*' a b 
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return A scalar or a list of symbol instances derived from \linkS4class{Symbol} 
                                    mul = function(v, w) {
                                          w$enter('mul', v, w)

                                          if ( ( lv <- length(v)) != 3 ) {
                                                 w$bug('mul', 'internal error: length of language object ', sQuote(v), ' ',
                                                              'is not 3, but ', lv, '. Expression is ', sQuote(deparse(v)))
                                          }

                                          a <- as.list(v)[2]
                                          b <- as.list(v)[3]

                                          w$msg(w, 'mul', 'a terms are ', w$enum(a))
                                          w$msg(w, 'mul', 'b terms are ', w$enum(b))


                                          res.a    <- c()
                                          effect.a <- 0
                                          for ( e in a ) {
                                                w$msg(w, 'mul', 'walking a term ', sQuote(deparse(e)))
                                                res.a <- c(res.a, r <- walkCode(e, w))
                                                for ( x in c(r) ) {
                                                      if ( is.Effect(x) ) effect.a <- effect.a + 1
                                                }
                                          }

                                          res.b    <- c()
                                          effect.b <- 0
                                          for ( e in b ) {
                                                w$msg(w, 'mul', 'walking b term ', sQuote(deparse(e)))
                                                res.b <- c(res.b,  r <- walkCode(e, w))
                                                for ( x in c(r) ) {
                                                      if ( is.Effect(x) ) effect.b <- effect.b + 1
                                                }
                                          }

                                          w$msg(w, 'mul', 'expanded terms of a  ', w$dump(res.a, w))
                                          w$msg(w, 'mul', 'expanded terms of b ',  w$dump(res.b, w))

                                          if ( effect.a && effect.b ) {
                                               w$fatal('mul', 'multiplication of effects as intended in ', sQuote(deparse(v)), ' ',
                                                              'is not implemented')
                                          }

                                          # Fold constants, evaluate symbolic constants, factor out effect symbols and 
                                          # treat like effect symbols also the instances of \linkS4class{Literal}. These
                                          # may have been created  by \code{w$sum(v,w)} or \code{w$sub(v, w) when encountering 
                                          # a summation of an effect and a numeric literal 

                                          product  <- 1
                                          symbols  <- c()
                                          for ( e in c(res.a, res.b) ) {
                                                w$msg(w, 'mul', 'handling term ', w$dump(e, w), ' from ', sQuote(deparse(e)))
                                                switch( class(e),
                                                       'Effect'  = symbols <- c(symbols, e),
                                                       'Literal' = symbols <- c(symbols, e),
                                                       'Symbol'  = product <- product * w$eval(e, w),
                                                       'numeric' = product <- product * e,
                                                        w$bug('mul', 'don\'t know how to handle ', w$enum(e), ' ',
                                                                     'within expression ', sQuote(deparse(v))))
                                          }


                                          # In case there were only literals just return a single number
                                          if ( length(symbols) == 0 ) {
                                               return(w$leave('mul', w$setCoef( product, w$getCoef(v)), v, w))
                                          }

                                          symbols <- lapply(symbols, function(x,v,w,p) w$setCoef(x, w$getCoef(x)*p*w$getCoef(v)), v, w, product )

                                          w$leave('mul', symbols, v, w)
                                    }, # end mul a b

                                    #' @title  Support division of an expression by a constant expression
                                    #' @param  v  Symbol node  '/' a b 
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return A scalar or a list of symbol instances derived from \linkS4class{Symbol} 
                                    div = function(v, w) {
                                          w$enter('div', v, w)

                                          # const / const  is allowed
                                          # fixed / const  is allowed:  -> coef(fixed) <- 1/const
                                          # const / effect is forbidden

                                          if ( ( lv <- length(v)) != 3 ) {
                                                 w$bug('div', 'length of language object ', sQuote(v), ' is not 3, but ', lv, '. ',
                                                              'Expression is ', sQuote(deparse(v)))
                                          }

                                          dividend <- c()
                                          for ( e in as.list(v)[2]) {

                                                w$msg(w, 'div', 'walking dividend term ', sQuote(deparse(e)))

                                                for ( r in c( walkCode(e, w) ) ) {
                                                      w$msg(w, 'div', 'handle dividend symbol ', w$dump(r, w))
                                                      switch(class(r),
                                                             'Effect'  = dividend <- c(dividend, r),
                                                             'Literal' = dividend <- c(dividend, r),
                                                             'Symbol'  = dividend <- c(dividend, w$eval(r, w)), # `pi`, for example
                                                             'numeric' = dividend <- c(dividend, r),
                                                                         w$bug('div', 'don\'t know how to handle dividend ', w$enum(r), ' ',
                                                                                      'within expression ', sQuote(deparse(v))))
                                                }
                                          }

                                          divisor <- 1
                                          for ( e in as.list(v)[3]) {
                                                w$msg(w, 'div', 'walking divisor term ', sQuote(deparse(e)))
                                                for ( r in c( walkCode(e, w) )) {
                                                      w$msg(w, 'div', 'handle divisor symbol ', w$dump(r, w))
                                                      switch(class(r),
                                                             'Symbol'  = divisor <- divisor * w$eval(r, w), # a symbolic constant like `pi`
                                                             'numeric' = divisor <- divisor * r,
                                                             'Effect'  = w$fatal('div', 'division by effect ', sQuote(r), ' ',
                                                                                        'as in ', sQuote(deparse(v)), ' ',
                                                                                        'is currently not implementable'),

                                                              w$bug('div', 'don\'t know how to handle divisor ', w$enum(r), ' ',
                                                                           'within expression ', sQuote(deparse(v))))
                                             }
                                          }

                                          if ( !is.finite(divisor) || divisor == 0) {
                                               w$fatal('div', 'won\'t divide by ', divisor, ' in ', sQuote(deparse(v)))
                                          }

                                          res <- lapply( dividend, 
                                                         function(x, d)  w$setCoef(x, w$getCoef(x) * w$getCoef(v) / d), d = divisor)

                                          w$msg(w, 'div', 'dividend is ', w$dump(dividend, w))
                                          w$msg(w, 'div', 'divisor is ',  w$dump(divisor, w))

                                          w$leave('div', res, v, w)
                                    }, # end div a b

                                    #' @title  Support for interaction of effects as in A:B:C:D
                                    #' @param  v  Symbol node  ':' a b 
                                    #' @param  w  Environment created by \code{makeCodeWalker}
                                    #' @return An effect symbol of \linkS4class{Effect} representing the interaction
                                    ita = function(v, w) {
                                          w$enter('ita', v, w)

                                          res <- c()
                                          for ( e in as.list(v)[-1] )  {

                                                w$msg(w, 'ita', 'walking interaction ', sQuote(deparse(e)))

                                                res <- c( res, walkCode(e, w))
                                          }

                                          res <- Effect( name = paste0(lapply(res, as.character), collapse=':' ) )

                                          # the interaction effect name must be part of the effect list
                                          if ( !w$is.effect(name(res), w)) {
                                                w$fatal('ita', 'within expression ', sQuote(deparse(v)), ', ',
                                                               'the interaction named ', sQuote(name(res)), ' ',
                                                               'was not fond in the list of effect names')
                                          }
                                          coef(res) <- coef(res) * w$getCoef(v) 

                                          w$leave('ita', res, v, w)
                                    }, # end ita a b

                                    #' @title  Support parenthesized terms '(' expr ')'  
                                    #' @param  v Symbol node
                                    #' @param  w Environment created by \code{makeCodeWalker}
                                    #' @return A list of symbols
                                    exp = function(v, w) {
                                          w$enter('exp', v, w)

                                          res <- c()
                                          for ( e in as.list(v)[-1] ) {
                                                w$msg(w, 'exp', 'walking term ', sQuote(deparse(e)))
                                                res <- c( res, walkCode(e, w))
                                          }

                                          symbols <- c()
                                          for ( e in  res ) {
                                                symbols <- c( symbols, w$setCoef(e, w$getCoef(e) * w$getCoef(v)))
                                          }

                                          w$leave('exp', symbols, v, w)
                                    },# end exp ression

                                    #' @title  Evaluate otherwise unhandled expressions ( e.g. a function call or operators like '^')
                                    #' @param  v Symbol node
                                    #' @param  w Environment created by \code{makeCodeWalker}
                                    #' @return A scalar, else call w$fatal    
                                    #' @note   Called from the walkCode() whenever w$handler(v, w) returns null
                                    call = function(v, w) {
                                           w$enter('call', v, w)

                                           operator <- v[[1]]

                                           w$msg(w, 'call', 'requested function is ', sQuote(operator))

                                           parms <- c()
                                           for ( e in as.list(v)[-1] )  {
                                                 w$msg(w, 'call', 'walking term ', sQuote(deparse(e)))

                                                 for ( p in c( walkCode(e, w) ) ) {
                                                       if ( is.Effect(p)  ) {
                                                            w$fatal('call', 'within expression ', sQuote(deparse(v)), ', ' ,
                                                                             w$enum(name(p)), ' ',
                                                                            'must not denote an effect because its estimate ',
                                                                            'is not yet known here')
                                                       }

                                                       if ( is.Symbol(p) && coef(p) != 1 )
                                                            # coef != 1 so form a product of the symbol and its coefficient
                                                            parms <- c(parms, w$parse('call',text=paste0('(',quoted(p),'*',coef(p),')', w)))
                                                       else if (is.Symbol(p) )
                                                            # coef == 1 so just avoid any surprising error messages
                                                            parms <- c(parms, w$parse('call',text=quoted(p),w))
                                                       else
                                                            parms <- c(parms, p)
                                                 }
                                           }

                                           w$msg(w, 'call', 'running ', w$as.docall(operator,parms), ' ',
                                                            'from expression ', sQuote(deparse(v)))

                                           res <- try ( do.call(as.character(operator), as.list(parms)), silent=T )

                                           if ( inherits(res, 'try-error') ) {
                                                w$fatal('call', w$as.docall(operator, parms), ' ',
                                                               'said ', sQuote(condition(res)), '. ',
                                                               'Expression was ', sQuote(deparse(v)))
                                           }

                                           if ( ! length(res) ) {
                                                  w$fatal('call', 'the call ', sQuote(deparse(v)), ' did not return anything')
                                           }

                                           if ( !all(is.numeric(res)))  {
                                                w$fatal('call',  sQuote(deparse(v)), ' returned  ', w$enum(res), ' ',
                                                                 'which is not numeric. This is currently not supported')
                                           }

                                           if ( !all(is.finite(res))) {
                                                 w$fatal('call',  sQuote(deparse(v)), ' returned non-finite result ', w$enum(res))
                                           }

                                           if ( ( l <- length(res)) > 1 ) {
                                                  w$fatal('call', 'the expression ', sQuote(deparse(v)), ' ',
                                                                  'returned ', l, ' numbers: ' , w$enum(res), '. ',
                                                                  'This is currently not supported')
                                           }

                                           w$leave('call', w$setCoef(res * w$getCoef(v), 1), v, w)
                                    },# end call


                                    #' @title  Handle terminal symbols
                                    #' @param  e Is either a named symbol or a numeric literal. 
                                    #' @param  w Is the environment created by \code{makeCodeWalker}
                                    #' @return If `e` represents
                                    #'         \itemize{
                                    #'               \item {an effect name}    return a new instance of \linkS4class{Effect} 
                                    #'               \item {a symbol name}     return a new instance of \linkS4class{Symbol}.
                                    #'               \item {a numeric literal} return it unchanged
                                    #'         }
                                    leaf = function(e, w) {
                                           w$enter('leaf', e, w)
                                           e <- switch( class(e),
                                                       'numeric' = e,
                                                       'name'    = if ( w$is.effect(e, w) )
                                                                        Effect( name = as.character(e) )
                                                                   else
                                                                        Symbol( name = as.character(e) ),
                                                        w$bug('leaf', 'don`t know how to handle a symbol of class ', sQuote(class(e))))
                                           w$leave('leaf', e, e, w)
                                    },# end leaf


                                    #### end of functions called from within w$handler

                                    #### begin of utility functions




                                    #' @title  Evaluate symbolic constants not naming an effect
                                    #' @param  e An instance of \linkS4class{Symbol}
                                    #' @param  w Is the environment created by \code{makeCodeWalker}
                                    #' @return A scalar numeric
                                    #' @note   eval is used as a utility only.
                                    eval = function(e, w) {
                                           w$enter('eval', e, w)

                                           if ( length(e) !=  1 || ! is.Symbol(e)  || is.Effect(e)) {
                                                w$bug('eval', 'the expression ', sQuote(deparse(e)), ' ',
                                                              'must represent a symbol and not an efffect')
                                           }

                                           expr <- ifelse( coef(e) == 1, 
                                                           quoted(e), 
                                                           paste0(quoted(e), '*', coef(e)))

                                           w$msg(w, 'eval', 'trying to evaluate ', sQuote(expr)) 

                                           pexpr <- w$parse('eval', text = expr, w)

                                           res   <- try ( eval(pexpr), silent = T )

                                           if ( inherits(res, 'try-error' ) ) {
                                                w$fatal('eval', 'evaluation of expression ', sQuote(expr), ' ',
                                                                'failed with ', dQuote(condition(res)))
                                           }

                                           if ( length(res) != 1 || !is.numeric(res) || !is.finite(res)  ) {
                                                w$fatal('eval', 'expression ', sQuote(expr), ' ',
                                                                'did not evaluate to a numeric constant of length 1 ',
                                                                'Result is ', sQuote(res))
                                           }

                                           w$leave('eval', res, e, w) 
                                    },# end eval


                                    #' @title  parse input and handle error
                                    #' @param  name  Name of the calling function
                                    #' @param  text  Text to be parsed
                                    #' @param  w     Environment created by \code{makeCodeWalker}
                                    #' @return resulting from \code{[base]parse()}
                                    #' @note   stops on error
                                    parse = function(name, text, w ) {
                                            res <- try ( parse(text = text ) )
                                            if ( inherits(res, 'try-error') ) {
                                                 w$fatal('name','parsing of ', sQuote(text), ' ',
                                                                'failed with ', sQuote(condition(res))) 
                                            } 
                                            res
                                    },# end  parse


                                    #' @title  Summarize and reorganize both lhs and rhs in order to match the glht internal interface
                                    #'
                                    #' @param  v    Symbol node
                                    #' @param  w    Environment created by \code{makeCodeWalker}
                                    #' @param  lhs  A list of instances of \linkS4class{Effect} or \linkS4class{Literal}
                                    #' @param  rhs  A list of instances of \linkS4class{Effect} or \linkS4class{Literal}
                                    #'
                                    #' @return A list with two named entries:
                                    #'         \itemize{
                                    #'               \item lhs A unique list of \linkS4class{Effect} instances
                                    #'               \item rhs A single numeric value
                                    #'         }
                                    #'
                                    #' @note   Exclusively called from within the topmost handler \code{w$equ(v, w)}
                                    #'
                                    #' @description 
                                    #'         The reorganization works as follows
                                    #'         \itemize{
                                    #'               \item a Shift all effect symbols from `rhs` to `lhs` and flip their sign.
                                    #'               \item b Shift all literal symbols from `lhs` to `rhs` and flip their sign.
                                    #'               \item c Collapse all duplicated effect symbols on  `lhs` into a set of 
                                    #'                       unique symbols by summing up the coefficients from all namesakes.
                                    #'               \item d Sum all coefs from the literal symbols on rhs into a single number. 
                                    #'         }
                                    transform = function(v, w, lhs, rhs) {
                                       w$enter('transform', v, w)

                                       # move symbols from source to destination side of the equation
                                       flip <- function(dst, src, is.match) {
                                          tmp <- c()
                                          for ( s in src ) {
                                                if ( is.match(s) ) {
                                                     dst <- c(dst, w$setCoef(s, -w$getCoef(s)))
                                                } else {
                                                     tmp <- c(tmp, s)
                                                }
                                          }
                                          list(dst = dst, src = tmp)
                                       }


                                       # shift all effects over to the left-hand side
                                       tmp <- flip(dst = lhs, src = rhs, is.match = is.Effect )
                                       lhs <- tmp$dst
                                       rhs <- tmp$src

                                       # switch all literals over to the right-hand side
                                       tmp <- flip(dst = rhs, src = lhs, is.match = function(x) is.numeric(x) || is.Literal(x) )
                                       lhs <- tmp$src
                                       rhs <- tmp$dst

                                       # collect all symbolic placeholders for numbers (created within w$add(v,w) and w$sub(v,w))
                                       rhs <- lapply(rhs, function(x) if( is.Literal(x) ) coef(x) else x )


                                       # combine the placeholders for additive terms into a real number
                                       sum.rhs <- sum(unlist(rhs))

                                       w$msg(w,'transform', 'sum of ', length(rhs), ' additive terms on rhs is ', sum.rhs, ', ',
                                               'additive terms are ', w$dump(rhs, w))

                                       w$leave('transform', list(lhs = w$merge(v, w, lhs), rhs = sum.rhs), v, w)
                                    },# end transform


                                    #' @title  Add up the coefficients of effects which may have been reference multiple times
                                    #'
                                    #' @param  v    Symbol node of \code{w$eqn(v, w)}
                                    #' @param  w    Environment as created by \code{makeCodeWalker}
                                    #' @param  lhs  A list of possibly duplicated instances of \linkS4class{Effect} 
                                    #'
                                    #' @return A list of unique of \linkS4class{Effect} instances with updated coefficients
                                    #'
                                    #' @note   The very same effect may appear several times within an equation. Consider, e.g, 
                                    #'         the effect symbol `x1` in  
                                    #'
                                    #'              \deqn{ 3 * ( x2 - x1 ) - 2 * ( x3 - x1 ) == 0}
                                    #'
                                    #'         where coefficient of `x1` should evaluate to -1 after summation.
                                    #'
                                    #' @description 
                                    #'         The function will return a unique list of \linkS4class{Effect} instances with 
                                    #'         their coef attribute set to the sum of the coefficients of all instances having the
                                    #'         same effect name. These copies came into being when \code{w$leaf(v, w)} was
                                    #'         called on the respective symbol.
                                    merge = function(v, w, lhs) {
                                       w$enter('merge', v, w)

                                       if ( any( idx <- unlist(lapply(lhs, function(x) !is.Effect(x) ))) ) {
                                            w$bug('merge', 'the following symbols must not reside on lhs any more: ', w$enum(lhs[idx]))
                                       }

                                       merged  <- c()


                                       names <- lapply(lhs[unlist(lapply(lhs, is.Effect))], name)

                                       w$msg(w, 'merge', 'unique effect symbols are ', w$enum(names))

                                       for ( name in unique(names) ) {

                                             # find all instances of s in lhs
                                             instances <- lhs[ unlist(lapply(lhs, function(x, s) x  == name)) ]
                                             w$msg(w, 'merge', 'summing up coefficients of ', length(instances), ' instances of ', sQuote(name))

                                             sum <- 0
                                             for ( instance in instances ) {
                                                   w$msg(w, 'merge', 'adding coefficient ', w$dump(instance, w) )
                                                   sum <- sum + coef(instance)
                                             }

                                             w$msg(w, 'merge', 'sum of the coefficients over all instances of ', sQuote(name), ' is ', sum)

                                             merged <- c(merged, Effect(name, sum))
                                       }

                                       w$leave('merge', merged, v, w)
                                    },# end merge


                                    #' @title  Translate the comparison operator into traditional text
                                    #' @param  v  Symbol node of \code{w$eqn(v, w)}
                                    #' @param  w  Environment as created by \code{makeCodeWalker}
                                    #' @return One of: 'less', 'greater' or 'two.sided'
                                    as.side = function(v, w) {
                                              switch( v,
                                                      '<=' = 'greater',
                                                      '>=' = 'less',
                                                      '==' = 'two.sided',
                                                      '='  = 'two.sided',
                                                      w$fatal('side', 'operator ', sQuote(v), ' ',
                                                                      'is not in ', w$enum(c('<=', '>=', '==', '='))))
                                    },# end as.side

                                    #' @title  Wrapper to access the coefficient of symbol, an expression or a constant
                                    #' @param  e   Can be anything, also a node itself representing a complex expression
                                    #' @return The value of the associated coefficient, or 1 if it was not set before
                                    getCoef = function(e) {
                                              switch( class(e),
                                                      'Effect'  = coef(e),
                                                      'Literal' = coef(e),
                                                      'Symbol'  = coef(e),
                                                       ifelse( is.null(a<-attr(e,'coef')), 1, a ))

                                    },# end getCoef

                                    #' @title  Wrapper function to set coefficient of something to a new value
                                    #' @param  e     Can be anything, perhaps also a node itself representing an expression
                                    #' @param  coef  A scalar
                                    #' @return The object `e` with its coefficient set to `coef`
                                    setCoef = function(e, coef) {
                                              if ( length(e) != 1 || length(coef) !=1 ) {
                                                   w$bug('setCoef', 'Don\'t know how to handle a list of values. ',
                                                                    'Parameter e is ',    paste(sQuote(e),collapse=', '), ',  ',
                                                                    'Parameter coef is ', paste(sQuote(coef),collapse=', ')) 
                                              }
                                              if ( class(e) == 'numeric' ) {
                                                    e <- e * coef
                                              } else if ( is.Symbol(e) ) {
                                                    coef(e) <- coef
                                              } else {
                                                    attr(e,'coef') <- coef;
                                              }
                                              e
                                    },# end setCoef

                                    #' @title  Format a comma separated list of entities for printing them along with messages 
                                    #' @param  x   A list of anything which can be coerced to character
                                    #' @return A   A character string containig a dump of the single-quoted entities, 
                                    #'             separated by comma or the string ' and ', if there are only two entries
                                    #' @note   Function isn't always working as expected.
                                    enum = function(x) {
                                           paste(sQuote(sapply(c(x), function(x) as.character(x))), collapse=ifelse( length(x) == 2, ' and ', ', '))
                                    },

                                    #' @title  Print a list of arguments and their associated coefficients
                                    #' @param  x   A list of enties having coefficients 
                                    #' @return A character string representing the list of entities and their 
                                    #          associated coefficients as in  "( 'x1' x -1 ) and ( 'x2' x 0.3 )" 
                                    dump = function(x, w) {
                                           w$enum(lapply(c(x), FUN=function(y,w) paste('(',y,'x', w$getCoef(y),')'), w))
                                    },

                                    #' @title  Prepare a list for tracing a call to do.call()
                                    #' @param  op    An operator name, e.g. 'exp' or '^'
                                    #' @param  args  A list of arguments to operator call (operators are treated like functions in R)
                                    #' @return A string representing the call
                                    as.docall = function(op, args) {
                                      paste0('do.call(',paste(c(op,args),collapse=', '),')')
                                    },

                                    #' @title  Output a trace line describing the interpreter state 
                                    #' @param  w     Environment as created by \code{makeCodeWalker}
                                    #' @param  name  The name of the function from where the trace originates
                                    #' @param  ...   Parameters to be printed by \code{message}
                                    #' @return NULL
                                    msg = function(w, name, ...) {
                                          # R does lazy evaluation of function arguments, so overhead is minimal if verbose is false
                                          if (!verbose) return(NULL) 

                                          # `sh0` is defined as a parameter to the outmost function
                                          message(rep('.', indent), name, '(', sh0, '): ', ...)
                                          NULL
                                    },

                                    #' @title  Print a message describing the problem and stop the program
                                    #' @param  name  The name of the function where the error gets detected
                                    #' @param  ...   Parameters to be printed by \code{stop}
                                    #' @return The function does not return
                                    #' @note   Called if inconsistent user input gets detected
                                    fatal = function(name, ...) {
                                            # `sh0`  is defined as a parameter to the outmost function
                                            stop(paste0('multcomp::expression2coef::walkCode::', name, '(', sh0, ')' ), ': ', ...,'.', 
                                                 call. = F)
                                    },

                                    #' @title  Stops the program with a description of the problem and request to send a bug report.
                                    #' @param  name  The name of the function where the error gets detected
                                    #' @param  ...   Parameters to be printed by \code{stop}
                                    #' @return The function does not return
                                    #' @note   Called on if an inconsistent state of the interpreter itself gets detected
                                    bug = function(name, ...) {
                                          # `sh0` is defined as a parameter to the outmost function
                                          stop(paste0('multcomp::expression2coef::walkCode::', name, '(', sh0, ')' ), ': ', ..., '. ',
                                                      'Please file a bug report!')
                                    },

                                    #' @title  Deparse a variable if she is of type language and contains a complex expression
                                    #' @param  x  The expression to deparse
                                    #' @return An empty string or the deparsed expression `x`
                                    deparsed = function(x) {
                                          ifelse(is.language(x) && length(x) > 1, paste0(' (or ',deparse(x),')'),'')
                                    },

                                    #' @title  Trace the entry to a function if `verbose` is true
                                    #' @param  name  The name of the function (mostly a handler like mul, div, sub, add)
                                    #' @param  v     Symbol node
                                    #' @param  w     Environment created by \code{makeCodeWalker}
                                    #' @return NULL
                                    #' @note   Due to R's lazy evaluation scheme, function arguments are only evaluated if
                                    #          the parameter `verbose` of the outmost function is true.
                                    enter = function(name, v, w) {

                                            if (!verbose) return(NULL) 

                                            # `indent` is defined at top of the outmost function
                                            indent <<- indent + 4 # indent is defined within enclosing function
                                            w$msg(w, name, '(v) x ', w$getCoef(v), ' is ' , w$enum(v),  w$deparsed(v),
                                                           ', mode(v) = ', mode(v),
                                                           ', typeof(v) = ', typeof(v),
                                                           ', length(v) = ', length(v))
                                            NULL
                                    },

                                    #' @title  Trace the result of a function immediately before leaving if `verbose` is true
                                    #' @param  name  The name of the function (mostly a handler like mul, div, sub, add)
                                    #' @param  e     The result of the function 
                                    #' @param  v     Symbol node
                                    #' @param  w     Environment created by \code{makeCodeWalker}
                                    #' @return The unchanged value of parameter `e`
                                    #' @note   Due to R's lazy evaluation scheme, function arguments are only evaluated if
                                    #          the parameter `verbose` of the outmost function is true.
                                    leave = function(name, e, v, w) {

                                            if ( !verbose ) return(e) 

                                            w$msg(w, name, 'returns ( ', w$dump(e, w), ' ) ')

                                            # `indent` is defined at top of the outmost function
                                            indent <<- indent - 4

                                            e
                                    }) # end of walkCode(lhs, makeCodeWalker( ... ))

   hypothesis <-
   codetools::walkCode(as.list(ex)[[1]], w = w)

   if ( is.null(hypothesis$symbols)) {
        stop('multcomp::expression2coef(', sh0, '): ',
             'no valid comparison found in ', sQuote(deparse(ex[[1]], 500)))
   }

   effect.names <- c()
   effect.coefs <- c()

   # create interface to chrlinfct2matrix
   for ( effect in c(hypothesis$symbols$lhs)) {
         effect.names  <- c( effect.names, as.character(effect))
         effect.coefs  <- c( effect.coefs, attr(effect, 'coef'))
   }

   if ( verbose ) {
        message('multcomp::expression2coef(', sh0, '): returning names ', paste0(sQuote(effect.names), collapse=', '))
        message('multcomp::expression2coef(', sh0, '): returning coefs ', paste0(sQuote(effect.coefs), collapse=', '))
   }

   list( coef        = effect.coefs,
         names       = effect.names,
         m           = hypothesis$symbols$rhs,
         alternative = hypothesis$equation$alt,
         lhs         = paste0(hypothesis$equation$lhs, collapse = ''),
         rhs         = paste0(hypothesis$equation$rhs, collapse = ''))
}

#' @title   Interpreter for character representations of linear functions
#'
#' @description
#'          Extract the effect coefficients from a potentially named list of null hypotheses 
#'
#' @param   ex        An (optionally named) character vector of one or more null hypotheses
#' @param   var       Character vector of effect names
#' @param   verbose   Trace the interpreter actions if true
#' @param   zerocheck Enable checking for zero or near zero coefficients.
#' @param   ...       passed down from ghlt, unused
#' @return  list(K, m, alternative, h0 = list( labels = list( original, generated ))) with
#'          \itemize{
#'               \item K            Matrix of numeric effect coefficients, with \code{rownames(K)} 
#'                                  set to the character string representing the left-hand side
#'               \item m            Vector of right hand side values
#'               \item alternative  One of: 'less', 'greater' or 'two.sided'
#'               \item h0           A list describing the hypotheses passed in argument `ex`
#'                   \itemize{
#'                        \item lhs          Array of character strings representing the left hand side
#'                        \item rhs          Array of character strings representing the right hand side
#'                        \item labels       Labels for the list of equations argument passed in argument ex
#'                             \itemize{
#'                                  \item original   Names passed along with argument ex (may be null)
#'                                  \item generated  Names like 'h0.1', 'h0.2', ...      (null if original was not null )
#'                   }
#'               }
#'         }
#' @seealso expression2coef
#' @examples 
#' a <- '1/2*( x1 -   x2 ) +1 = 1/3 * x1:x2'
#' b <- '1/3*( x2 - 2*x1 ) -3 = 2 + exp((-2)^3) * x1:x2'
#' multcomp:::chrlinfct2matrix(ex = c('test x' = a,  'test y' = b), var = c('x1', 'x2', 'x1:x2'))
#' @rdname internal/chrlinfct2matrix
chrlinfct2matrix <- function(ex, var, verbose=F, zerocheck=T, ...) {

    enum <- function(x) {
      paste(paste0(sQuote(x), collapse=', '))
    }

    fatal <- function(...) {
       stop('multcomp::chrlinfct2matrix',': ', ..., '.', call. = F)
    }

    if ( !is.character(ex) )
         fatal('Argument ex is not of type character but ', typeof(ex))

    if ( !is.character(var) )
         fatal('Argument var is not of type character but ', typeof(var))

    # enforce a consistent naming scheme for the equations listed in argument ex
    if ( any( dups<-duplicated(names(ex) )) )
         fatal('Within argument ex, the following linear function names are duplicated: ', enum(names(ex)[dups]) )

    # enforce unique effect names 
    if ( any( dups<-duplicated(var)) )
         fatal('There are duplicate effect names in argument var: ', enum(var[dups]) )

    generatedNames <- NULL
    # if the hypotheses are not named, just set default names for reference
    if ( is.null( lfNames <- names(ex) ) ) {
         names(ex) <- generatedNames <- sprintf('h0.%02d', seq_along(ex))
    }

    K <- matrix(0, nrow = length(ex), ncol = length(var))

    colnames(K) <- var
    rownames(K) <- 1:length(ex)
    m <- rep(0, length(ex))

    lhs <- c()
    rhs <- c()

    for ( i in seq_along(ex) ) {

          expr <- try( parse(text = ex[i]), silent = T)
          if ( inherits(expr, 'try-error')  || length(expr[[1]]) != 3 ) {
               fatal('argument ', sQuote(ex[i]), ' cannot be interpreted as expression')
          }

          tmp <- expression2coef( expr, vars=var, sh0 = names(ex[i]),
                                  verbose=verbose, zerocheck=zerocheck)

          if (!all(tmp$names %in% var)) {
               missing <- tmp$names[!tmp$names %in% var]
               fatal('variable(s) ', paste(sQuote(missing), collapse=', '), ' not found')
          }

          for ( n in tmp$names )
                K[i, var == n] <- tmp$coef[tmp$names == n]

          m[i] <- tmp$m

          if ( i == 1 ) {
              alternative <- tmp$alternative
          } else if (tmp$alternative != alternative) {
              # limitation to a single alternative appears to be outdated!
              fatal("mix of alternatives is currently not implemented")
          }

          rownames(K)[i] <- tmp$lhs

          lhs <- c(lhs, tmp$lhs)
          rhs <- c(rhs, tmp$rhs)

    }

    list( K = K,
          m = m,
          alternative = alternative,
          h0 = list( lhs    =  lhs,
                     rhs    =  rhs,
                     labels =  list( original = lfNames, generated = generatedNames ) ))
}


# $Id$

# extract effect    coefficients from a single null hypothesis
# @param  ex        already parsed string containing the null
# @param  vars      list of known effect names
# @param  sh0       short name of H0: the label as in ex = c('A greater B' = ' A - B >= 0 ' )
# @param  verbose   trace parser states and actions
# @param  zerocheck check for effects being near zero if true
# @return list(coef, names, alternative, lhs, rhs )
expression2coef <- function(ex, vars, sh0, verbose, zerocheck) {
   if ( verbose ) {
        message('multcomp::expression2coef(', sh0, '): ',
                'effect names are ', paste0(sQuote(vars), collapse = ', '))
        message('multcomp::expression2coef(', sh0, '): length of ex is ', length(ex))
   }


   # used in verbose mode to indent trace
   indent <- 0

   # used with as.copy(x) ( see there )
   copy.id <- 0

   # used with as.addend(x) ( see there )
   addend.id <- 0

   # used in context of try(...) to retrieve error message
   condition <- function(x) attr(x, 'condition')$message

   hypothesis <-
   walkCode(ex[[1]],
            makeCodeWalker(
                            # determine if a variable name matches an effect name as defined by `vars`, which is a list
                            # of effects passed as a parameter to the outmost function.
                            is.effect = function(x) {
                                        as.character(x) %in% vars
                            },

                            # called directly from the codetools::walkCode() framework and serves as a operator dispatcher
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
                                                   '=='   = w$equ,
                                                   '='    = w$equ,
                                                    NULL )  # -> use default action which is w$call
                                      w$leave('handler', v, v, w)
                                      res
                            }, # end of handler ( called directly from codetools::codeWalker )

     #### begin of the functions called from within w$handler (see above)

                           # == lhs rhs : called from within w$equ()
                           # >= lhs rhs : called from within w$geq()
                           # <= lhs rhs : called from within w$leq()
                           equ  = function(v, w) {
                                  w$enter('equ', v, w)

                                  if ( length(as.list(v)) != 3 ) {
                                       w$fatal('equ', 'expression ', sQuote(deparse(v)), ' is incomplete')
                                  }

                                  operator <- as.character(v[[1]])
                                  equ.lhs  <- paste0(deparse(v[[2]], 500), collapse='')
                                  equ.rhs  <- paste0(deparse(v[[3]], 500), collapse='')
                                  equ.bth  <- paste(equ.lhs, operator, equ.rhs)

                                  w$msg(w, 'equ', sQuote(equ.lhs), ' ', sQuote(operator), ' ', sQuote(equ.rhs))

                                  sym.lhs <- c()
                                  for ( e in as.list(v)[2]) {
                                        w$msg(w, 'equ', 'walking lhs term ', sQuote(deparse(e)), ' x 1')
                                        sym.lhs <- c( sym.lhs, walkCode(w$setCoef(e, 1), w))
                                  }


                                  sym.rhs <- c()
                                  for ( e in as.list(v)[3]) {
                                        w$msg(w, 'equ', 'walking rhs term ', sQuote(deparse(e)), ' x 1')
                                        sym.rhs <- c( sym.rhs, walkCode(w$setCoef(e, 1), w))
                                  }


                                  # shift all effects to lhs and collapse any offsets into rhs
                                  trx <- w$transform(v, w, sym.lhs, sym.rhs)

                                  if ( any(  unlist(lapply(trx$lhs, is.numeric)) )
                                     | any(  unlist(lapply(trx$rhs, is.symbol))  )
                                     | any( !unlist(lapply(trx$rhs, is.finite))  )) {
                                       w$bug('equ', 'the normalisation of ', sQuote(equ.bth), ' was not successful. ',
                                                    'Symbols on left-hand side are [', w$dump(w, trx$lhs), '], ',
                                                    'right-hand side is ', w$enum(trx$rhs))
                                  }

                                  if ( length(trx$lhs) == 0 ) {
                                       w$fatal('equ', 'no effects found in ', sQuote(equ.bth))
                                  }


                                  zero <- c()
                                  for ( e in trx$lhs ) {
                                        if ( max( abs(w$getCoef(e)))  <= 2 * .Machine$double.eps)  {
                                             zero <- c( zero, e)
                                        }
                                  }

                                  if ( zerocheck && length(zero) > 0 ) {
                                       w$fatal('equ', 'within expression ', sQuote(equ.bth), ' ',
                                                      'the coefficients of the following effects ',
                                                      'evaluate to zero (or near zero, entailing numeric instability): ',
                                                       w$enum(zero), '. Values are ',
                                                       paste(sprintf("%g",lapply(zero, w$getCoef)), collapse = ', '), '. ',
                                                      'Calling glht with the argument zerocheck=F will disable this test.')
                                  }

                                  w$leave('equ', list( symbols  = list( lhs = trx$lhs,
                                                                        rhs = trx$rhs),
                                                       equation = list( opr = operator,
                                                                        alt = w$side(w, operator),
                                                                        lhs = equ.lhs,
                                                                        rhs = equ.rhs,
                                                                        bth = equ.bth)), v, w)
                            },

                            # '<=' lhs rhs
                            leq = function(v, w) {
                                  w$enter('leq', v, w)
                                  w$leave('leq', w$equ(v, w), v, w)
                            },



                            # '>=' lhs rhs
                            geq = function(v, w) {
                                  w$enter('geq', v, w)
                                  w$leave('geq', w$equ(v, w), v, w)
                            },



                            # '-' a b : support for subtraction of effects or constants (but not both)
                            # '-' a   : negation of effect coefficients, constants or expressions
                            sub = function(v, w) {
                                  w$enter('sub', v, w)

                                  minuend    <- as.list(v)[2]
                                  subtrahend <- as.list(v)[3]
                                  uminus     <- length(as.list(v)) == 2 # one operator + one operand

                                  if ( uminus )
                                       subtrahend <- NULL

                                  w$msg(w, 'sub', 'minuend is ', sQuote(minuend))
                                  w$msg(w, 'sub', 'subtrahend is ', sQuote(subtrahend))
                                  w$msg(w, 'sub', 'uminus is ', uminus)

                                  exp.coef <- ifelse( uminus, -w$getCoef(v), w$getCoef(v))
                                  res      <- c()
                                  for ( e in minuend ) {
                                        w$msg(w, 'sub', 'walking minuend ', sQuote(deparse(e)), ' x ', exp.coef)
                                        res <- c( res, walkCode(w$setCoef(e, exp.coef ), w))
                                  }

                                  for ( e in subtrahend ) {
                                        w$msg(w, 'sub', 'walking subtrahend ', sQuote(deparse(e)), ' x ', -exp.coef)
                                        res  <- c( res, walkCode(w$setCoef(e, -exp.coef), w))
                                  }

                                  sum     <- 0
                                  symbols <- c()

                                  # split up result set into constants and symbols
                                  for ( e in res ) {

                                        if ( w$is.effect(e) || w$is.copy(e)) {
                                             if ( any( w$lookup(e, symbols) ) ) {
                                                  symbols <- c(symbols, w$as.copy(w, e) )
                                             } else {
                                                  symbols <- c(symbols, e)
                                             }
                                        } else if ( is.numeric(e)) {
                                             sum <- sum + e
                                        } else if ( w$is.addend(e) ) {
                                             w$msg(w, 'sub', 'sum before adding ', Quote(e), ' is ', sum)
                                             sum <- sum + w$getCoef(e)
                                             w$msg(w, 'sub', 'sum after adding  ', Quote(e), ' is ', sum)
                                        } else if ( is.symbol(e)) {
                                             # evaluate the symbol possibly representing a numeric constant like `pi`
                                             sum <- sum + w$eval(e, w)
                                        } else {
                                             w$bug('sub', "don't know how to handle ", sQuote(e), " ",
                                                          "within expression ", sQuote(deparse(v)))
                                        }
                                  }

                                  # fold constants into single number
                                  if ( length(symbols) == 0 ) {
                                       return(w$leave('sub', sum, v, w))
                                  }

                                  if ( sum  ) {
                                       symbols <- c(symbols, w$as.addend(w, sum))
                                  }

                                  w$leave('sub', symbols, v, w)
                            }, # end sub

                            # ':' a b : support for the interaction of effects as in A:B:C:D
                            ita = function(v, w) {
                                  w$enter('ita', v, w)

                                  res <- c()
                                  for ( e in as.list(v)[-1] )  {

                                        w$msg(w, 'ita', 'walking interaction ', sQuote(deparse(e)), ' x 1')

                                        res <- c( res, r <- walkCode(w$setCoef(e, 1), w))

                                        # because cell means models provide estimates only for the interactions,
                                        # just test for `e` being a name here, but not for `e` being an effect.
                                        # The estimate may not have been defined by the model
                                        if ( ! is.symbol(r) ) {
                                               w$fatal('ita', 'within expression ', sQuote(deparse(v)), ', ',
                                                              'the term ', sQuote(e), ' does not represent a symbol')
                                        }
                                  }

                                  res <- w$setCoef(as.name(paste0(lapply(res, as.character), collapse=':' )), w$getCoef(v))

                                  # however, an estimate for the interaction must be present anyway
                                  if ( !w$is.effect(res)) {
                                        w$fatal('ita', 'within expression ', sQuote(deparse(v)), ', ',
                                                       'the term ', sQuote(res), ' does not name an effect')
                                  }

                                  w$leave('ita', res, v, w)
                            }, # end interaction (ita)

                            # '(' expr ')' : support parenthesized terms
                            exp = function(v, w) {
                                  w$enter('exp', v, w)

                                  res <- c()
                                  for ( e in as.list(v)[-1] ) {
                                        w$msg(w, 'exp', 'walking term ', sQuote(deparse(e)), ' x 1')
                                        res <- c( res, walkCode(w$setCoef(e, 1), w))
                                  }

                                  symbols <- c()
                                  for ( e in  res ) {
                                        symbols <- c( symbols, w$setCoef(e, w$getCoef(e) * w$getCoef(v)))
                                  }

                                  w$leave('exp', symbols, v, w)
                            }, # end '(' expr ')'

                            # '+' a b : support for the addition of constants or effects (but not both)
                            # '+' a   : mostly a no-op
                            add = function(v, w) {
                                  w$enter('add', v, w)

                                  res <- c()
                                  for ( e in as.list(v)[-1] ) {
                                        w$msg(w, 'add', 'walking term ', sQuote(deparse(e)), ' x 1' )
                                        res <- c( res, walkCode( w$setCoef(e, 1), w))
                                  }

                                  symbols <- c()
                                  sum     <- 0
                                  for ( e in res ) {
                                        w$msg(w, 'add', 'handle symbol ', w$dump(w, e) )

                                        if ( w$is.effect(e) || w$is.copy(e) ) {
                                             if ( any( w$lookup(e, symbols) ) ) {
                                                  symbols <- c(symbols, w$as.copy(w, e) )
                                             } else {
                                                  symbols <- c(symbols, e)
                                             }
                                        } else if ( is.numeric(e)) {
                                             # fold constant
                                             sum <- sum + e
                                        } else if ( w$is.addend(e) ) {
                                             sum <- sum + w$getCoef(e)
                                        } else if ( is.symbol(e)) {
                                             # evaluate and fold non-effect numeric constants like `pi`
                                             sum <- sum + w$eval(e, w)
                                        } else {
                                             w$bug('add', "don't know how to handle ", sQuote(e), ' ',
                                                          'within expression ', sQuote(deparse(v)))
                                        }
                                  }

                                  sum <- w$setCoef( sum * w$getCoef(v), 1)
                                  w$msg(w, 'add', 'after sum ', sQuote(sum) )

                                  # just return constants as a single number if there are no effects
                                  if ( length(symbols) == 0 ) {
                                       return( w$leave('add', sum, v, w))
                                  }

                                  # turn the summand into a special variable undergoing transformations like
                                  # an effect variable
                                  if ( sum ) {
                                       symbols <- c(symbols, w$as.addend(w, sum))
                                  }

                                  # associate the expression coefficient with all leafs
                                  res <- c()
                                  for ( e in symbols ) {
                                        res <- c( res, w$setCoef(e, w$getCoef(e) * w$getCoef(v)))
                                  }

                                  w$leave('add', res, v, w)
                            }, # end '+' a b

                            # '*' a b : support multiplication of effects or constants by constants
                            mul = function(v, w) {
                                  w$enter('mul', v, w)

                                  if ( ( lv <- length(v)) != 3 ) {
                                         w$bug('mul', 'internal error: length of language object ', sQuote(v), ' ',
                                                      'is not 3, but ', lv, '. Expression is ', sQuote(deparse(v)))
                                  }

                                  a <- as.list(v)[2]
                                  b <- as.list(v)[3]

                                  res.a  <- c()
                                  for ( e in a ) {
                                        w$msg(w, 'mul', 'walking term a ', sQuote(deparse(e)), ' x 1')
                                        res.a <- c(res.a, walkCode(w$setCoef(e, 1), w))
                                  }

                                  res.b <- c()
                                  for ( e in b ) {
                                        w$msg(w, 'mul', 'walking term b ', sQuote(deparse(e)), ' x 1')
                                        res.b <- c(res.b, walkCode(w$setCoef(e, 1), w))
                                  }


                                  # prevent multiplication of fixed effects as in 'A * B' but still
                                  # allow the multiplication effects by a constant
                                  for ( e.a in res.a ) {
                                        for ( e.b in res.b ) {
                                              if (   ( w$is.effect(e.a) || w$is.copy(e.a) )
                                                  && ( w$is.effect(e.b) || w$is.copy(e.b) ) ) {
                                                  w$fatal('mul', 'multiplication of effects ', w$enum(c(e.a, e.b)), ' ',
                                                                 'as in ', sQuote(deparse(v)), ' ', 'is not implemented')
                                              }
                                        }
                                  }

                                  res      <- c(res.a, res.b)
                                  product  <- 1
                                  symbols  <- c()
                                  addends  <- c()
                                  for ( e in res ) {
                                        w$msg(w, 'mul', 'handle symbol ', w$dump(w, e))
                                        if ( w$is.effect(e) || w$is.copy(e) ) {
                                             symbols <- c(symbols, e)
                                        } else if ( is.numeric(e)) {
                                             # fold constant
                                             product <- product * e
                                        } else if ( w$is.addend(e) ) {
                                             # treat a summand like a symbol in this context
                                             addends <- c(addends, e)
                                        } else if ( is.symbol(e)) {
                                             # fold non-effect numeric constants like `pi`
                                             product <- product * w$eval(e, w)
                                        } else {
                                             w$bug('mul', "don't know how to handle ", sQuote(e), ' ',
                                                          'within expression ', sQuote(deparse(v)))
                                        }
                                  }

                                  # also take the expression coefficient into account
                                  product <- product * w$getCoef(v)

                                  # if there are only literals return a single number
                                  if ( length(symbols) == 0 ) {
                                       return(w$leave('mul', product, v, w))
                                  }

                                  # associate the product of all literals as a coefficient
                                  # with all symbols, treating additive terms as in 2*(a + b + 1)
                                  # like symbols
                                  res <- c()
                                  for ( e in c(symbols, addends) ) {
                                        res <- c(res, w$setCoef(e, w$getCoef(e) * product))
                                  }

                                  w$leave('mul', res, v, w)
                            }, # end mul a b

                            # short hand function
                            classify = function(x, fun) {
                                  unlist(lapply(x, fun))
                            }, # end classify

                            # '/' a b : division of an expression by a constant expression
                            div = function(v, w) {
                                  w$enter('div', v, w)

                                  # const / const  is allowed
                                  # fixed / const  is allowed:  -> coef(fixed) <- 1/const
                                  # const / effect is forbidden

                                  # collect all leafs, including constants

                                  if ( ( lv <- length(v)) != 3 ) {
                                         w$bug('div', 'internal error: length of language object ', sQuote(v), ' ',
                                                      'is not 3, but ', lv, '. Expression is ', sQuote(deparse(v)))
                                  }

                                  dividend <- c()
                                  for ( e in as.list(v)[2]) {

                                        w$msg(w, 'div', 'walking dividend term ', sQuote(deparse(e)), ' x 1')

                                        r <- walkCode(w$setCoef(e, 1), w)

                                        if ( is.numeric(r) || w$is.effect(r) || w$is.addend(r) ) {
                                             dividend <- c(dividend, r)
                                        } else if ( is.symbol(r)) {
                                             # non-effect numeric constants like `pi`
                                             dividend <- c(dividend, w$eval(r, w))
                                        } else {
                                             w$bug('div', "don't know how to handle the dividend ", sQuote(e), ' ',
                                                          'within expression ', sQuote(deparse(v)))
                                        }
                                  }

                                  divisor <- 1
                                  for ( e in as.list(v)[3]) {

                                        w$msg(w, 'div', 'walking divisor term ', sQuote(deparse(e)), ' x 1')
                                        r <- walkCode(w$setCoef(e, 1), w)

                                        w$msg(w, 'div', 'handle divisor symbol ', w$dump(w, r))
                                        if ( w$is.effect(r) || w$is.copy(r) ) {
                                             w$fatal('div', 'division by effect ', sQuote(r), ' ',
                                                            'as in ', sQuote(deparse(v)), ' ', 'is not implemented')
                                        } else if ( w$is.addend(r) ) {
                                             w$warn(w, 'div', 'the result from using the summation term placeholder ', sQuote(r), ' ',
                                                             'may be correct, or it maybe not. ',
                                                             'Please check the coefficients of expression ', sQuote(deparse(v)))
                                             divisor <- divisor * w$getCoef(r)
                                        } else if ( is.numeric(r) ) {
                                             divisor <- divisor * r
                                        } else if ( is.symbol(r) ) {
                                             # non-effect numeric constants like `pi`
                                             divisor <- divisor * w$eval(r, w)
                                        } else {
                                             w$bug('div', "don't know how to handle divisor ", sQuote(e), ' ',
                                                          'within expression ', sQuote(deparse(v)))
                                        }
                                  }

                                  if ( !is.finite(divisor) || divisor == 0) {
                                       w$fatal('div', "won't divide by ", divisor, ' in ', sQuote(deparse(v)))
                                  }

                                  res <- c()
                                  for ( e in dividend ) {
                                        if ( is.numeric(e)) {
                                             res <- c(res, e * w$getCoef(v) / divisor)
                                        } else {
                                             res <- c(res, w$setCoef(e, w$getCoef(e) * w$getCoef(v) / divisor))
                                        }
                                  }

                                  w$msg(w, 'div', 'dividend is ', w$dump(w, dividend))
                                  w$msg(w, 'div', 'divisor is ',  w$dump(w, divisor))

                                  w$leave('div', res, v, w)
                            }, # end div a b

                            # The function gets called from within the codetools::codeWalker() framework with `e`
                            # being either a named symbol or a numeric literal
                            leaf = function(e, w) {
                                   w$enter('leaf', e, w)
                                   if ( !( is.symbol(e) || is.numeric(e) ) ) {
                                           w$bug('leaf',  sQuote(e), ' is not a symbol or a numeric constant')
                                   }
                                   w$leave('leaf', e, e, w)
                            }, # end leaf

                            # called from the codetools:walkCode() framework whenever w$handler() returns NULL and the symbol
                            # does not represent a terminal symbol (a leaf, see above)
                            call = function(v, w) {
                                   w$enter('call', v, w)

                                   operator <- v[[1]]

                                   if ( is.numeric(operator) ) {
                                        w$fatal('call', 'an operator appears to be missing somewhere within expression ', sQuote(deparse(v)))
                                   }

                                   w$msg(w, 'call', 'requested function is ', sQuote(operator))

                                   parms <- c()
                                   for ( e in as.list(v)[-1] )  {
                                         w$msg(w, 'call', 'walking term ', sQuote(deparse(e)), ' x 1')

                                         parms <- c(parms, p <- walkCode(w$setCoef(e, 1), w))

                                         if ( any( w$is.effect(p) ) ) {
                                              w$fatal('call', 'within expression ', sQuote(deparse(v)), ', ' ,
                                                              'the term(s) ', w$enum(p), ' must not denote an effect')
                                         }
                                   }

                                   cparms <- c()
                                   for ( e in parms ) {
                                         cparms <- c(cparms, parse(text = paste(e, '*', w$getCoef(e))))
                                   }

                                   w$msg(w, 'call', 'is ', w$enum(c(operator, cparms)), '). Expression is ', deparse(v))

                                   res <- try ( do.call(as.character(operator), as.list(cparms)), silent=T )

                                   if ( class(res) == 'try-error' ) {
                                        w$fatal('call', 'the evaluation of the expression ', sQuote(deparse(v)), ' ',
                                                        'failed with ', dQuote(condition(res)))
                                   }

                                   if ( length(res) != 1 || !is.numeric(res) || !is.finite(res)) {
                                        w$fatal('call', 'the expression ', sQuote(deparse(v)), ' ',
                                                        'did not evaluate to a real valued numeric constant. ',
                                                        'Result is ', sQuote(res))
                                   }

                                   w$leave('call', w$setCoef(res * w$getCoef(v), 1), v, w)
                            }, # end call. 

     #### end of functions called from within w$handler 

     #### begin of utility functions

                            # exclusively used to evaluate symbolic constants not naming an effect
                            # eval is a utility, not a handler
                            eval = function(e, w) {
                                   w$enter('eval', e, w)

                                   if ( length(e) !=  1 || ! is.symbol(e)  || w$is.effect(e)) {
                                        w$bug('eval', 'reaching here, the expression ', sQuote(deparse(e)), ' ',
                                                      'must represent a single symbol and not an efffect')
                                   }

                                   res <- try ( eval(e), silent=T )

                                   if ( class(res) == 'try-error' ) {
                                        w$fatal('eval', 'the evaluation of the expression ', sQuote(deparse(e)), ' ',
                                                        'failed with ', dQuote(condition(res)))
                                   }

                                   if ( length(res) != 1 || !is.numeric(res) || !is.finite(res)  ) {
                                        w$fatal('eval', 'the expression ', sQuote(deparse(e)), ' ',
                                                        'did not evaluate to a real valued numeric constant. ',
                                                        'Result is ', sQuote(res))
                                   }

                                   w$leave('eval', w$setCoef(res * w$getCoef(e), 1), e, w)
                            }, # end eval

                            # transform summarize and reorganize both left-hand and right-hand side, it gets
                            # called from within the topmost handler equ(). 
                            #
                            # The processing scheme goes as follows:
                            #
                            #  a.) shift all effects from `rhs` to `lhs` and flip their sign.
                            #  b.) shift all constants from `lhs` to `rhs` and flip their sign.
                            #  c.) collapse all duplicated symbols on  `lhs` into a set of unique symbols with their
                            #      coefs set to the sum of the coefficients of their respective namesakes
                            #  d.) sum up all constants with the offsets on rhs.
                            #
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

                               w$msg(w, 'transform', 'lhs before shift is ', w$enum(lhs))
                               w$msg(w, 'transform', 'rhs before shift is ', w$enum(rhs))

                               # shift all effects over to the left-hand side
                               tmp <- flip(dst = lhs, src = rhs, is.match = function(x) w$is.effect(x) || w$is.copy(x) )
                               lhs <- tmp$dst
                               rhs <- tmp$src

                               # switch all summands over to the right-hand side
                               tmp <- flip(dst = rhs, src = lhs, is.match = function(x) is.numeric(x) ||  w$is.addend(x) )
                               lhs <- tmp$src
                               rhs <- tmp$dst

                               # turn all symbolic placeholders for numbers (created within add and sub) back into real numbers
                               rhs <- unlist(lapply(rhs, function(x, w) ifelse( w$is.addend(x), w$getCoef(x), x ), w))

                               w$msg(w, 'transform', 'lhs after shift is ', w$enum(lhs))
                               w$msg(w, 'transform', 'rhs after shift is ', w$enum(rhs))


                               w$leave('transform', list(lhs = w$merge(v, w, lhs), rhs = sum(rhs)), v, w)
                            },

                            # merge() gets called from within transform() to add up the coefficients 
                            # of all referenced effects, which may have appeared several times in 
                            # the equation.
                            #
                            # The function returns a unique list of the effect symbols. Their coef 
                            # attribute will be equal to the sum of the coefficients of their 
                            # placeholders previously created by as.copy(). 
                            #
                            # See the documentation for as.copy() below in order to understand why this 
                            # operation is needed.
                            merge = function(v, w, lhs) {
                               w$enter('merge', v, w)

                               if ( !all( idx <- w$classify(lhs, function(x) w$is.effect(x) || w$is.copy(x) ))  ) {
                                     w$bug('merge', 'there are still some non-effect symbols left on lhs:', w$enum(lhs[!idx]))
                               }

                               merged <- c()
                               for ( s in lhs[ w$classify(lhs, w$is.effect ) ] ) {
                                     nc <- length(zoo <- lhs[ w$getCopies(s, lhs) ])
                                     w$msg(w, 'merge', 'adding up the coefs of ', nc, ' copies of ', sQuote(s))

                                     sum <- 0
                                     # sum up the coefs of all symbols being a copy of as `s`, including `s` itself
                                     for ( ds in  zoo ) {
                                           w$msg(w, 'merge', 'adding coef of ', w$dump(w, ds))
                                           sum <- sum + w$getCoef(ds)
                                     }
                                     w$msg(w, 'merge', 'sum over all copies of ', sQuote(s), ' is ', sum)
                                     merged <- c(merged, w$setCoef(s, sum))
                               }

                               w$leave('merge', merged, v, w)
                           }, # end of merge


                            # as.copy(w, s): create a separate copy of symbol `s` using a different symbol
                            # name
                            #
                            # Problem: With 'x1 - 2 * x1 = 0', the effect 'x1' should be associated with a
                            # coefficient of '-1'. This could be solved by adding up the coefs of all occurences
                            # of `x1`
                            #
                            # In R, however, a symbol is a reference to a certain memory region. This region
                            # does not change when the reference associated with that region gets copied along
                            # by various program actions, while the reference itself may be copied.
                            # Yet all of these instances of a reference in turn map to the same address, much
                            # like copied pointers and/or references in "C" and "C++" always refer to the
                            # very same memory address.
                            #
                            # Now, for keeping track of coefs, the implementation decision was to bind them as an
                            # attribute of the effect symbols. This was done in order to have the coefs around
                            # whenever needed.
                            #
                            # In the example above, the information hold by this attribute would be '+1'
                            # for the first 'x1'. Yet as soon as the next `x1` would get processed, this value
                            # would change to  '-2' for both since there there is only one instance of this
                            # variable.

                            # R symbols themselves can not be copied because there is no copy method for 
                            # symbols. In addition the symbol class is "sealed"; hence no S4 class can be
                            # derived from a symbol.
                            #
                            # The solution implemented here is to create a copy under a different name
                            # and associate the copy with the original symbol. The 'merge' utility running
                            # below the 'equ' handler will then sum up all these coefficients and associate
                            # the result with the original effect variable.

                            as.copy = function(w, s) {
                                  if ( ! w$is.effect(s) )
                                         w$fatal('as.copy', sQuote(s), ' does not represent an effect')

                                  # make sure the name is unique: `copy.id` is defined at top of the outmost function
                                  copy.id <<- copy.id + 1
                                  e <- as.name(paste(s, '~', copy.id))
                                  attr(e, 'copy of') <- as.character(s)
                                  attr(e, 'coef')    <- 1
                                  w$msg(w,'as.copy', 'created the symbol ', sQuote(e), ' ',
                                                     'as a true copy of the effect symbol ', sQuote(s))
                                  e
                            }, # end as.copy

                            # see if the symbol is a copy of another effect symbol
                            is.copy = function(s) {
                                 ! is.null(attr(s, 'copy of'))
                            }, # end is.copy

                            # return a list of copies of the effect symbol `s` contained
                            # in list of `symbols`. The returned list also contains `s`
                            getCopies = function(s,  symbols ) {
                               unlist(lapply(symbols, function(x, s)
                                      s == x || !is.null(attr(x, 'copy of')) && attr(x, 'copy of') == s, s ))
                            }, # end of getCopies

                            # as.adddend create a placeholder for a literal ( e.g., the '3' in  'a - (b + 3) == 0' ).
                            # This avoids a technical problem when trying to perform an addition of a constant
                            # term and an effect of unknown numeric value. 
                            as.addend = function(w, value) {
                                  if ( ! is.numeric(value) )
                                         w$fatal('as.addend', sQuote(value), ' is not numeric')

                                  # `addend.id` is defined at top of the outmost function
                                  addend.id <<- addend.id + 1

                                  # make sure the name is unique.
                                  e <- as.name(paste(value, '~', addend.id))
                                  attr(e, 'addend') <- T
                                  attr(e, 'coef')   <- value

                                  w$msg(w, 'as.addend', 'created the symbol ', sQuote(e), ' ',
                                                        'in place of the the numeric literal ', sQuote(value))
                                  e
                            }, # end as.addend

                            # see if the symbol is a generated placeholder for a literal
                            is.addend = function(e) {
                                 !is.null(attr(e, 'addend'))
                            }, # end is.addend

                            # translate comparison operator to traditional text
                            side = function(w, opr) {
                                   switch( opr,
                                           '<=' = 'greater',
                                           '>=' = 'less',
                                           '==' = 'two.sided',
                                           '='  = 'two.sided',
                                           w$fatal('side', 'operator ', sQuote(opr), ' ',
                                                           'is not in ', w$enum(c('<=', '>=', '==', '='))))
                            }, # end side

                            # return coefficient of `e`, or else 1 if coefficient was not set before
                            getCoef = function(e) {
                                      a <- attr(e, 'coef')
                                      ifelse( is.null(a), 1, a )
                            }, # end getCoef

                            # set coefficient of `e` to `coef` and return `e`
                            # setCoef is a utility, not a handler.
                            setCoef = function(e, coef) {

                                      if ( is.numeric(e) ) {
                                           # fold coefficient directly into the number and return.
                                           return(e * coef) 
                                      }

                                      # associate the coefficient with the symbol
                                      attr(e, 'coef') <- coef 
                                      e
                            }, # end setCoef

                            # lookup symbol `e` in list `symbols`
                            # lookup is a utility, not a handler.
                            lookup = function(e, symbols) {
                                     unlist(lapply(symbols, function(x, e) e == x, e ))
                            }, # end lookup

                            enum = function(x) {
                                   paste0(sQuote(x), collapse=ifelse(length(x) == 2, ' and ', ', '))
                            }, # end enum

                            dump = function(w, x) {
                                   w$enum(lapply(x, function(z, w) paste(as.character(z), 'x', w$getCoef(z)), w))
                            }, #end dump

                            warn = function(w, name, ...) {
                                   # `sh0` is defined as a parameter to the outmost function
                                   warning(name, '(', sh0, '): ', ...)
                            }, # end warn

                            msg = function(w, name, ...) {
                                  if (!verbose) return(NULL)

                                  # `sh0` is defined as a parameter to the outmost function
                                  message(rep('.', indent), name, '(', sh0, '): ', ...)
                            }, #end msg

                            # error is a consequence of user input
                            fatal = function(name, ...) {
                                    # `sh0`  is defined as a parameter to the outmost function
                                    stop(paste0('multcomp::expression2coef::walkCode::', name, '(', sh0, ')' ), ': ', ...)
                            }, # end fatal

                            # error is internal to the program
                            bug = function(name, ...) {
                                  # `sh0` is defined as a parameter to the outmost function
                                  stop(paste0('multcomp::expression2coef::walkCode::', name, '(', sh0, ')' ), ': ', ..., '. ',
                                              'Please file a bug report!')
                            }, # end bug

                            # called on each handler entry in order to trace the course of actions
                            enter = function(name, v, w) {
                                    if (!verbose) return(NULL)

                                    # `indent` is defined at top of the outmost function
                                    indent <<- indent + 4 # indent is defined within enclosing function
                                    w$msg(w, name, 'v is ' , w$enum(v),
                                                   ifelse(is.language(v), paste0(' (or deparsed: ', sQuote(deparse(v)), ' )'), ''),
                                                   ', mode(v) = ', mode(v),
                                                   ', typeof(v) = ', typeof(v),
                                                   ', length(v) = ', length(v),
                                                   ', coef(v) = ', w$getCoef(v))
                            },

                            # called upon each return from handler in order to trace the course of actions
                            leave  = function(name, e, v, w) {
                                     if ( !verbose ) return(e)

                                     w$msg(w, name, 'returns [ ', w$dump(w, e), ' ] ',
                                                    'from expr (', deparse(v), ') x ', w$getCoef(v))

                                     # `indent` is defined at top of the outmost function
                                     indent <<- indent - 4
                                     e
                            })) # end of walkCode(lhs, makeCodeWalker( ... ))

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
         lhs         = hypothesis$equation$lhs,
         rhs         = hypothesis$equation$rhs)
}

# interpret character representations of linear functions
# @param   ex        list of null hypotheses
# @param   var       list of known effect names
# @param   verbose   activate parser trace if true
# @param   zerocheck enable checking for zero or near zero coefficients (default: TRUE)
# @param   ...       passed down from ghlt, unused
# @return  list( K,             # effect coefficients
#                m,             # right hand side constant ( a single value )
#                alternative,   # 'less'|'greater'|'two.sided'
#                h0 = list( labels = list( original,   # names passed along with @see ex (may be null)
#                                          generated   # names like 'H0.1', 'H0.2', ...  (or null if the former entry was not null )
#              ))
chrlinfct2matrix <- function(ex, var, verbose=F, zerocheck=T, ...) {

    enum <- function(arg, x) {
      paste(arg, '=', '[ ', paste0(sQuote(x), collapse=', '), ' ]')
    }

    if ( !is.character(ex) )
         stop("multcomp::chrlinfct2matrix", ": argument ", enum('ex', ex),
                                            "is not of type character")

    if ( !is.character(var) )
         stop("multcomp::chrlinfct2matrix", ": argument ", enum('var', var),
                                            "is not of type character")
    # enforce a consistent naming scheme
    if ( any( dups<-duplicated(names(ex) )) )
         stop("multcomp::chrlinfct2matrix", ": duplicated names in argument ",
                                               enum('ex', names(ex)[dups]))

    # set a default name for each hypothesis
    if ( is.null( exNames <- names(ex) ) ) {
         names(ex) <- paste0('H0.',1:length(ex))
    }

    K <- matrix(0, nrow = length(ex), ncol = length(var))
    colnames(K) <- var
    rownames(K) <- 1:length(ex)
    m <- rep(0, length(ex))

    for ( i in seq_along(ex) ) {

          expr <- parse(text = ex[i])

          if (length(expr[[1]]) != 3)
              stop("argument ", sQuote(ex[i]),
                   " cannot be interpreted as expression")

          tmp <- expression2coef( expr, vars=var, sh0 = names(ex[i]),
                                  verbose=verbose, zerocheck=zerocheck)

          if (!all(tmp$names %in% var))
              stop("variable(s) ", sQuote(tmp$names[!tmp$names %in% var]),
                   " not found")

          for ( n in tmp$names )
                K[i, var == n] <- tmp$coef[tmp$names == n]

          m[i] <- tmp$m

          if (i == 1) {
              alternative <- tmp$alternative
          } else if (tmp$alternative != alternative) {
              # this limitation appears to be outdated!
              stop("mix of alternatives currently not implemented")
          }

          rownames(K)[i] <- paste0(tmp$lhs, collapse = "")
    }

    list( K = K,
          m = m,
          alternative = alternative,
          h0 = list( labels = list( original = exNames, generated = names(ex) ) ))
}

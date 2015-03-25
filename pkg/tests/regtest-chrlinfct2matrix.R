### test the extended interpreter for both sides of linear hypotheses
### Features:
### - fully recursive expression parser built upon a small code fragment copied over from base::codetools
### - left and right-hand side are treated equally: effect names can appear on both sides
### - the parser stops if any of the following conditions is not met:
###   . the expression can be parsed as an equation
###   . the comparison operators are in '==', '=', '>=', '<='
###   . all operators and functions must finally evaluate to a constant values
###   . function parameters must not denote an effect name ( because the estimate is unknown when the parser runs )
###   . effects can not be multiplied or divided by another effect ( for the same reason )
###   . no effect coefficient must evaluate to zero ( because this is likely an error )
###
### - Examples:
###   x1 * x2             == 0                -> not accepted
###   x1 / x2             == 0                -> not accepted
###   x1 - x2             == x3               -> accepted because x3 will be shifted to the lhs, inverting sign
###   x1 + x2 -1          == 0                -> accepted because the literal gets shifted to rhs, inverting sign
###   f(x1)               == 0                -> not accepted if x1 denotes an effect
###   2*3                 == 6                -> not accepted because no effect was named ( no non-zero effect coefficient )
###   x1 - x1             == 0                -> not accepted because coef(x1) == 0
###   x1 + x2*0           == 0                -> not accepted because coef(x2) == 0
###   x1 + 3*(4-5+1)*x2   == 0                -> not accepted because coef(x2) == 0
###   x1*3/0              == 0                -> not accepted because coef(x1) is infinite
###   log(-1)*x1          == 0                -> not accepted because coef(x1) is infinite
###   x1 - x2 + 0         == 0                -> accepted because all offsets on lhs are added up on rhs
###   sin(pi/2) * x1      == 0                -> accepted if 'pi' is not an effect
###   sin(Pi/2) * x1      == 0                -> accepted if 'Pi' is not an effect and `Pi` evaluates to a numeric constant
###   x1 - x2             == -sqrt(sqrt(2^2)) -> accepted because the rhs evaluates to a numeric constant


expectSucc<- function(testname, x, expected ) {
 if ( inherits( tx <- try(x, silent = T), 'try-error' ) ) {
      stop(testname, ' unexpectedly failed. Result is: ', attr(tx,'condition')$message, '\n')
 }

 dK <- tx$K - expected$K
 dm <- tx$m - expected$m

 stopifnot(!any(is.null(dK) | !is.finite(dK)))
 stopifnot(!any(is.null(dm) | !is.finite(dm)))

 rownames(dK) <- names(dm) <- x$h0$labels$original

 e.dm <- error.dm <- 0 # because of short evaluation path with boolean expr
 if (    ( error.dK <- ( e.dK <- max(abs(dK)) ) > sqrt(.Machine$double.eps) )
      || ( error.dm <- ( e.dm <- max(abs(dm)) ) > sqrt(.Machine$double.eps) ) ) {

      rownames(dK) <- names(dm) <- tx$h0$labels$original
      message(testname,': ', ifelse(error.dK, '  ','no'), ' problem with dK: delta is ', e.dK)
      message(testname,': ', ifelse(error.dm, '  ','no'), ' problem with dm: delta is ', e.dm)
      print(dK)
      print(dm)
      stop(testname, ': unexpectedly failed\n')
 } else {
      message(testname, ': expectedly succeeded\n')
 }
}

expectFail <- function(testname, x) {
 if ( ! inherits( tx <- try(x, silent = T), 'try-error') ) {
        stop(testname, ': unexpectedly succeeded. Result is ', paste(tx, collapse = ', '),'\n')
 }
 message(testname, ': expectedly failed. Message is ', attr(tx,'condition')$message, '\n')
}


# external definitions for evaluation within equation 'l13'
l13.f <- function(x) exp(x)
PI     <- pi
common <-  c('x1','x2','x3','x4','x1:x2','x1:x2:x3','GenotypeWT:ReagentsCMH','GenotypeWT:ReagentsCMH+pSOD')

expectSucc('test 0.l06', multcomp:::chrlinfct2matrix(c(l06 = '-(x1 - x2)*-2  - (1/3+2)*( -x3 - 2*x4 )*7/-10 = -3'),c('x1','x2','x3','x4'), verbose = F), expected = list( K = c(2, -2, (1/3+2)*-7/10, 2*(1/3+2)*-7/10), m = -3))


expectSucc('test 0',
           x = multcomp:::chrlinfct2matrix( c( l01 = '  x1 - x2 = 2'
                                             , l02 = '  x2 + 3 * x3 = 1'
                                             , l03 = '  (x1 - x2) - (x3 - x4) =  0'
                                             , l04 = ' +(x1 - x2)*-2 - (1/3+2)*( +x3 - 2*x4 ) = -1'
                                             , l05 = ' -(x1 - x2)*-2 - (1/3+2)*( -x3 - 2*x4 ) = -2'
                                             , l06 = ' -(x1 - x2)*-2  - (1/3+2)*( -x3 - 2*x4 )*7/-10 = -3'
                                             , l07 = ' -1*(x1:x2 - x1:x2:x3) - x3 = -4'
                                             , l08 = ' -(x1:x2 - x1:x2:x3) - x3 = -4'
                                             , l09 = ' -(x1:x2 - 3*x1:x2:x3)*-2 - x3 -5/3*-x4= -5'
                                             , l10 = ' --cos(pi/4)*x1 - 10*(log(10^-3)+1)*-x2 -10^-3*x3 + -exp(-2)*x4= -6'
                                             , l11 = '  x1 + x2 + 0 = -7'
                                             , l12 = ' 1/pi*(1/3*`GenotypeWT:ReagentsCMH` - 2/3*`GenotypeWT:ReagentsCMH+pSOD`) == -8'
                                             , l13 = ' -pi*l13.f(3*-PI)/sqrt(2*PI)*(`x1`*exp(-pi) - `x2`*sqrt(sqrt(4*PI-1))) == -9'
                                             , l14 = ' log(110,12*exp(-pi^2))*(x1-x2) == -sqrt(sqrt(100^2))'
                                             , l15 = ' 1/2*(x1 - x2 -5)  ==  -3*x3 -17'
                                             , l16 = ' 1/2*(x1 - 1/2*x2 -3/5)  ==  -3*x3 -17'
                                             , l17 = ' x1+x1+x2 = x3+1'
                                             , l18 = ' 14 + x1 + 3 + x1 + 3*(2 + x2 + 3/2)/3 = 4*x3+1'
                                             ),
                                            c('x1','x2','x3','x4','x1:x2','x1:x2:x3','GenotypeWT:ReagentsCMH','GenotypeWT:ReagentsCMH+pSOD'),
                                            verbose=F
                                         ), # chrlinfct ...
            expected = list( K = rbind( c(           1,                 -1,                      0,                 0,       0,     0,         0,           0)
                                      , c(           0,                  1,                      3,                 0,       0,     0,         0,           0)
                                      , c(           1,                 -1,                     -1,                 1,       0,     0,         0,           0)
                                      , c(          -2,                  2,               -(1/3+2),         2*(1/3+2),       0,     0,         0,           0)
                                      , c(           2,                 -2,                (1/3+2),         2*(1/3+2),       0,     0,         0,           0)
                                      , c(           2,                 -2,          (1/3+2)*-7/10,   2*(1/3+2)*-7/10,       0,     0,         0,           0)
                                      , c(           0,                  0,                     -1,                 0,      -1,     1,         0,           0)
                                      , c(           0,                  0,                     -1,                 0,      -1,     1,         0,           0)
                                      , c(           0,                  0,                     -1,           -5/3*-1,       2,    -6,         0,           0)
                                      , c( --cos(pi/4),  10*(log(10^-3)+1),                 -10^-3,          -exp(-2),       0,     0,         0,           0)
                                      , c(           1,                  1,                      0,                 0,       0,     0,         0,           0)
                                      , c(           0,                  0,                      0,                 0,       0,     0,   1/(3*pi),  -2/(3*pi))
                                      , c( -pi*l13.f(3*-PI)/sqrt(2*PI)*+1*exp(-pi),
                                           -pi*l13.f(3*-PI)/sqrt(2*PI)*-1*sqrt(sqrt(4*PI-1)),
                                                                                                 0,                 0,       0,     0,         0,           0)
                                      , c( log(110,12*exp(-pi^2)),
                                          -log(110,12*exp(-pi^2)),                               0,                 0,       0,     0,         0,           0)
                                      , c(         1/2,               -1/2,                      3,                 0,       0,     0,         0,           0)
                                      , c(         1/2,               -1/4,                      3,                 0,       0,     0,         0,           0)
                                      , c(           2,                  1,                     -1,                 0,       0,     0,         0,           0)
                                      , c(           2,              3*1/3,                     -4,                 0,       0,     0,         0,           0)
                                      ), # end rbind
                             m = c(   2
                                  ,   1
                                  ,   0
                                  ,  -1
                                  ,  -2
                                  ,  -3
                                  ,  -4
                                  ,  -4
                                  ,  -5
                                  ,  -6
                                  ,  -7
                                  ,  -8
                                  ,  -9
                                  , -sqrt(sqrt(100^2))
                                  , -17+(1/2)*5
                                  , -17+(1/2)*(3/5)
                                  ,   1
                                  , -(14 +  3 +  3*(2 + 3/2)/3)+1
                                  ) ) #end expected =
            )  #end expectSucc


#test 01 previously failed for addressing x1 twice. It still fails, but this time for another reason: the sum of the coefs of all x1 is zero.
#A more apt test is in test 01a
expectFail('test 01',  multcomp:::chrlinfct2matrix(c(l1 = 'x1 - x1  = 0'), c('x1','x2'), verbose=F))
expectSucc('test 01a', multcomp:::chrlinfct2matrix(c(l1 = '3/5*(x1 - 2*x1)  = 0'), c('x1','x2'),verbose=F), expected = list(K = c(3/5*(1-2), 0), m = 0))

expectFail('test 02',  multcomp:::chrlinfct2matrix(c(l1 = 'x1 - X2  = 0'), c('x1','x2')))

#test 03 previously failed, but works now by shifting the constant terms to the right hand side, inverting sign
expectSucc('test 03',  multcomp:::chrlinfct2matrix(c(l1 = 'x1 - x2 -1 = 0'), c('x1','x2'),verbose=F), expected = list(K = c(1, -1), m = 1))

expectFail('test 04',  multcomp:::chrlinfct2matrix(c(l1 = 'x1 * x2  = 0'), c('x1','x2'),verbose=F))

expectFail('test 05',  multcomp:::chrlinfct2matrix(c(l1 = 'x1 / x2  = 0'), c('x1','x2')))

expectFail('test 06',  multcomp:::chrlinfct2matrix(c(l1 = 'x1 - exp(x2)  = 0'), c('x1','x2')))

expectFail('test 07',  multcomp:::chrlinfct2matrix(c(l1 = 'sin(Pi)*x1   = 0'), c('x1','x2')))

expectFail('test 08',  multcomp:::chrlinfct2matrix(c(l1 = '3*4 = 0'), c('x1','x2')))

expectFail('test 09',  multcomp:::chrlinfct2matrix(c(l1 = 'x1 + 3*(4-5+1)*x2 = 0'), c('x1','x2')))

expectFail('test 10',  multcomp:::chrlinfct2matrix(c(l1 = 'x1*3/0 = 0'), c('x1','x2')))

expectFail('test 11',  multcomp:::chrlinfct2matrix(c(l1 = 'log(-1)*x1 = 0'), c('x1','x2')))

expectFail('test 12',  multcomp:::chrlinfct2matrix(c(l1 = '1/2(x1-x2) = x3'), c('x1','x2','x3')))
expectFail('test 12b', multcomp:::chrlinfct2matrix(c(l1 = '1/2(1-2)*(x1-x2) = x3 '), c('x1','x2','x3')))
expectFail('test 12c', multcomp:::chrlinfct2matrix(c(l1 = '1/2(a-2)*(x1-x2) = x3 '), c('x1','x2','x3')))
expectFail('test 12d', multcomp:::chrlinfct2matrix(c(l1 = '1/a(1-2)*(x1-x2) = x3 '), c('x1','x2','x3')))
expectFail('test 12e', multcomp:::chrlinfct2matrix(c(l1 = '1/diff(c(1,2))*(x1-x2) = x3 '), c('x1','x2','x3')))

expectFail('test 13',  multcomp:::chrlinfct2matrix(c(l1 = '1/2(x1-x2) = '), c('x1','x2','x3')))
expectFail('test 13b', multcomp:::chrlinfct2matrix(c(l1 = '1/2*(x1-x2) = '), c('x1','x2','x3')))
expectFail('test 13c', multcomp:::chrlinfct2matrix(c(l1 = '1/2*(x1-x2) '), c('x1','x2','x3')))

expectFail('test 14',  multcomp:::chrlinfct2matrix(c(l1 = '1/2(x1-x2) = x3', l1 = '1/2(x1-x2) = x3'), c('x1','x2','x3')))

expectFail('test 15',  multcomp:::chrlinfct2matrix(c(l1 = '10^-30 *(x1 - x2) = 0'), c('x1','x2')))
expectSucc('test 15a', multcomp:::chrlinfct2matrix(c(l1 = '10^-30 *(x1 - x2) = 0'), c('x1','x2'), zerocheck=F), expected = list(K = c(10^-30, -10^-30), m = 0))

expectSucc('test 16',  multcomp:::chrlinfct2matrix(c('claim 1'= 'A - B + 3*B +10 +10 == 0'),c('A','B'),verbose=F),expected = list(K = c(1,2),m=-20))
expectFail('test 16b', multcomp:::chrlinfct2matrix(c('claim 1'= 'A - B + 3*B +10 +10 == 0'),c('A','B','B'),verbose=F))


expectSucc('test 17',  multcomp:::chrlinfct2matrix(c('claim 1'= 'B - A == pi'), c('A','B'), verbose = F), expected = list(K = c(-1,+1),m = +pi))
expectSucc('test 17b', multcomp:::chrlinfct2matrix(c('claim 1'= 'pi == B - A'), c('A','B'), verbose = F), expected = list(K = c(+1,-1),m = -pi))


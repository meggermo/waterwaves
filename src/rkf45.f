      SUBROUTINE RKF45
      INCLUDE "knd_params.inc"
C
      REAL   (KIND=RK) A_FEHLBERG (5)
      REAL   (KIND=RK) B_FEHLBERG (5, 5)
      REAL   (KIND=RK) C_FEHLBERG (6, 2)
C
      REAL   (KIND=RK) A_CASHKARP (5)
      REAL   (KIND=RK) B_CASHKARP (5, 5)
      REAL   (KIND=RK) C_CASHKARP (6, 2)
C
      REAL   (KIND=RK) A_DORMANDPRINCE (6)
      REAL   (KIND=RK) B_DORMANDPRINCE (6, 6)
      REAL   (KIND=RK) C_DORMANDPRINCE (7, 2)

C     Fehlberg:
      A_FEHLBERG (1) =  1.0 /  4.0
      A_FEHLBERG (2) =  3.0 /  8.0
      A_FEHLBERG (3) = 12.0 / 13.0
      A_FEHLBERG (4) =  1.0
      A_FEHLBERG (5) =  1.0 /  2.0

      B_FEHLBERG = 0.0
      B_FEHLBERG (1, 1) =       1.0 /    4.0
      B_FEHLBERG (1, 2) =       3.0 /   32.0
      B_FEHLBERG (2, 2) =       9.0 /   32.0
      B_FEHLBERG (1, 3) =    1932.0 / 2197.0
      B_FEHLBERG (2, 3) =   -7200.0 / 2197.0
      B_FEHLBERG (3, 3) =    7296.0 / 2197.0
      B_FEHLBERG (1, 4) =     439.0 /  216.0
      B_FEHLBERG (2, 4) =      -8.0
      B_FEHLBERG (3, 4) =    3680.0 /  513.0
      B_FEHLBERG (4, 4) =    -845.0 / 4104.0
      B_FEHLBERG (1, 5) =      -8.0 /   27.0
      B_FEHLBERG (2, 5) =       2.0
      B_FEHLBERG (3, 5) =   -3544.0 / 2565.0
      B_FEHLBERG (4, 5) =    1859.0 / 4104.0
      B_FEHLBERG (5, 5) =     -11.0 / 40.0

      C_FEHLBERG (1, 1) =    25.0 /   216.0
      C_FEHLBERG (2, 1) =     0.0
      C_FEHLBERG (3, 1) =  1408.0 /  2565.0
      C_FEHLBERG (4, 1) =  2197.0 /  4104.0
      C_FEHLBERG (5, 1) =    -1.0 /     5.0
      C_FEHLBERG (6, 1) =     0.0
      C_FEHLBERG (1, 2) =    16.0 /   135.0
      C_FEHLBERG (2, 2) =     0.0
      C_FEHLBERG (3, 2) =  6656.0 / 12825.0
      C_FEHLBERG (4, 2) = 28561.0 / 56430.0
      C_FEHLBERG (5, 2) =    -9.0 /    50.0
      C_FEHLBERG (6, 2) =     2.0 /    55.0

C Cash-Karp 4/5
C The Cash-Karp method was developed to satisfy different constraints,
C namely to deal with non-smooth problems better. They chose the cici,
C the percentage of the timestep in the iith step (i.e. t+ciΔtt+ciΔt
C is the time the iith step is calculated at) to be as uniform as
C possible, yet still achieve order 5. Then it also was derived to
C have embedded 1st, 2nd, 3rd, and 4th order methods with this uniformity
C of the cici. They are spaced in such a manner that you can find out
C where a stiff part starts by which difference is large. Moreover, note
C that the more stiff the equation, the worse a higher order method does
C (because it needs bounds on higher derivatives). So they develop a strategy
C which uses the 5 embedded methods to "quit early": i.e. if you detect
C stiffness, stop at stage i<6i<6 to decrease the number of function
C calls and save time. So in the end, this "pair" was developed with a lot
C of other constraints in mind, and so there's no reason to expect it would
C be "more accurate", at least as a 4/5 pair. If you add all of this other
C achinery then, on (semi-)stiff problems, it will be more accurate (but in
C that case you may want to use a different method like a W-Rosenbrock method).
C This is one reason why this pair hasn't become standard over the DP5 pair,
C but it still can be useful (maybe it would be good for a hybrid method which
C switches to a stiff solver when stiffness is encountered?).
      A_CASHKARP (1) = 1.0 /  5.0
      A_CASHKARP (2) = 3.0 / 10.0
      A_CASHKARP (3) = 3.0 /  5.0
      A_CASHKARP (4) = 1.0
      A_CASHKARP (5) = 7.0 /  8.0

      B_CASHKARP = 0.0
      B_CASHKARP (1, 1) =     1.0 /      5.0
      B_CASHKARP (1, 2) =     3.0 /     40.0
      B_CASHKARP (2, 2) =     9.0 /     40.0
      B_CASHKARP (1, 3) =     3.0 /     10.0
      B_CASHKARP (2, 3) =    -9.0 /     10.0
      B_CASHKARP (3, 3) =     6.0 /      5.0
      B_CASHKARP (1, 4) =   -11.0 /     54.0
      B_CASHKARP (2, 4) =     5.0 /      2.0
      B_CASHKARP (3, 4) =   -70.0 /     27.0
      B_CASHKARP (4, 4) =    35.0 /     27.0
      B_CASHKARP (1, 5) =  1631.0 /  55296.0
      B_CASHKARP (2, 5) =   175.0 /    512.0
      B_CASHKARP (3, 5) =   575.0 /  13828.0
      B_CASHKARP (4, 5) = 44275.0 / 110592.0
      B_CASHKARP (5, 5) =   253.0 /   4096.0

      C_CASHKARP (1, 1) =  37.0 /  378.0
      C_CASHKARP (2, 1) =   0.0
      C_CASHKARP (3, 1) = 250.0 /  621.0
      C_CASHKARP (4, 1) = 125.0 /  594.0
      C_CASHKARP (5, 1) =   0.0
      C_CASHKARP (6, 1) = 512.0 / 1771.0

      C_CASHKARP (1, 2) =  2825.0 / 27648.0
      C_CASHKARP (2, 2) =     0.0
      C_CASHKARP (3, 2) = 18575.0 / 48384.0
      C_CASHKARP (4, 2) = 13525.0 / 55296.0
      C_CASHKARP (5, 2) =   277.0 / 14336.0
      C_CASHKARP (6, 2) =     1.0 /    40.0

C Dormand-Prince 4/5
C The Dormand-Prince method was developed to be accurate as a 4/5 pair with
C local extrapolation usage (i.e. step with the order 5 pair. This is because
C it was designed to be close to optimal (i.e. minimal) principle truncation
C error coefficient (under the restraint of also having the minimal number of
C steps to achieve order 5). It has an order 4 interpolation which is free,
C but needs extra steps for an order 5 interpolation.

      A_DORMANDPRINCE (1) = 1.0 /  5.0
      A_DORMANDPRINCE (2) = 3.0 / 10.0
      A_DORMANDPRINCE (3) = 4.0 /  5.0
      A_DORMANDPRINCE (4) = 8.0 /  9.0
      A_DORMANDPRINCE (5) = 1.0
      A_DORMANDPRINCE (6) = 1.0

      B_DORMANDPRINCE (1, 1) =      1.0 /     5.0
      B_DORMANDPRINCE (1, 2) =      3.0 /    40.0
      B_DORMANDPRINCE (2, 2) =      9.0 /    40.0
      B_DORMANDPRINCE (1, 3) =     44.0 /    45.0
      B_DORMANDPRINCE (2, 3) =    -56.0 /    15.0
      B_DORMANDPRINCE (3, 4) =     32.0 /     9.0
      B_DORMANDPRINCE (1, 4) =  19372.0 /  6561.0
      B_DORMANDPRINCE (2, 4) = -25360.0 /  2187.0
      B_DORMANDPRINCE (3, 4) =  64448.0 /  6561.0
      B_DORMANDPRINCE (4, 4) =   -212.0 /   729.0
      B_DORMANDPRINCE (1, 5) =   9017.0 /  3168.0
      B_DORMANDPRINCE (2, 5) =   -355.0 /    33.0
      B_DORMANDPRINCE (3, 5) =  46732.0 /  5247.0
      B_DORMANDPRINCE (4, 5) =     49.0 /   176.0
      B_DORMANDPRINCE (5, 5) =  -5103.0 / 18656.0
      B_DORMANDPRINCE (1, 6) =     35.0 /   384.0
      B_DORMANDPRINCE (2, 6) =      0.0
      B_DORMANDPRINCE (3, 6) =     500.0 / 1113.0
      B_DORMANDPRINCE (4, 6) =     125.0 /  192.0
      B_DORMANDPRINCE (5, 6) =   -2187.0 / 6784.0
      B_DORMANDPRINCE (6, 6) =      11.0 /   84.0

      C_DORMANDPRINCE (1, 1) =     35.0 /    384.0
      C_DORMANDPRINCE (2, 1) =      0.0
      C_DORMANDPRINCE (3, 1) =    500.0 /   1113.0
      C_DORMANDPRINCE (4, 1) =    125.0 /    192.0
      C_DORMANDPRINCE (5, 1) =  -2187.0 /   6784.0
      C_DORMANDPRINCE (6, 1) =     11.0 /     84.0
      C_DORMANDPRINCE (7, 1) =      0.0
      C_DORMANDPRINCE (1, 2) =   5179.0 /  57600.0
      C_DORMANDPRINCE (2, 2) =      0.0
      C_DORMANDPRINCE (3, 2) =   7571.0 /  16695.0
      C_DORMANDPRINCE (4, 2) =    393.0 /    640.0
      C_DORMANDPRINCE (5, 2) = -92097.0 / 339200.0
      C_DORMANDPRINCE (6, 2) =    187.0 /   2100.0
      C_DORMANDPRINCE (7, 2) =      1.0 /      40.0

      END SUBROUTINE RKF45

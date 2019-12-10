C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++ Copyright S. Zhang and J. M. Jing                               ++
C++ Computation of Special Functions, John Wiley & Sons, New York,  ++
C++ 1996 and http://jin.ece.illinois.edu .                          ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
        SUBROUTINE EIX(X,EI)
C
C       ============================================
C       Purpose: Compute exponential integral Ei(x)
C       Input :  x  --- Argument of Ei(x)
C       Output:  EI --- Ei(x) ( x > 0 )
C       ============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0) THEN
           EI=-1.0D+300
        ELSE IF (X.LE.40.0) THEN
           EI=1.0D0
           R=1.0D0
           DO 15 K=1,100
              R=R*K*X/(K+1.0D0)**2
              EI=EI+R
              IF (DABS(R/EI).LE.1.0D-15) GO TO 20
15         CONTINUE
20         GA=0.5772156649015328D0
           EI=GA+DLOG(X)+X*EI
        ELSE
           EI=1.0D0
           R=1.0D0
           DO 25 K=1,20
              R=R*K/X
25            EI=EI+R
           EI=DEXP(X)/X*EI
        ENDIF
        RETURN
        END

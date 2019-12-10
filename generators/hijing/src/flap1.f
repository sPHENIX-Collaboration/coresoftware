C
C
        FUNCTION FLAP1(X)
        COMMON/PACT/BB,B1,PHI,Z1
        SAVE  /PACT/
        R=SQRT(BB**2+X**2)
        FLAP1=WDSAX1(R)
        RETURN
        END

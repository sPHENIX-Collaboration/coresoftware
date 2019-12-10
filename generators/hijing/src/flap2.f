C
C
        FUNCTION FLAP2(X)
        COMMON/PACT/BB,B1,PHI,Z1
        SAVE  /PACT/
        R=SQRT(BB**2+X**2)
        FLAP2=WDSAX2(R)
        RETURN
        END

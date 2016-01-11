C
        FUNCTION FGP3(X)
        COMMON/PACT/BB,B1,PHI,Z1
        SAVE  /PACT/
        R1=SQRT(B1**2+Z1**2)
        R2=SQRT(BB**2+B1**2-2.0*B1*BB*COS(PHI)+X**2)
        FGP3=B1*WDSAX1(R1)*WDSAX2(R2)
        RETURN
        END

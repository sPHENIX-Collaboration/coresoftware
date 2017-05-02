         
      SUBROUTINE AIZ(IFUN,IFAC,X0,Y0,GAIR,GAII,IERRO)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C COMPUTATION OF THE AIRY FUNCTION AI(Z) OR ITS DERIVATIVE AI'(Z)
C THE CODE USES:
C      1. MACLAURIN SERIES FOR |Y|<3 AND -2.5<X<1.3 (Z=X+I*Y)
C      2. GAUSS-LAGUERRE QUADRATURE  FOR |Z|<15 AND  WHEN 
C         MACLAURIN SERIES ARE NOT USED.
C      3. ASYMPTOTIC EXPANSION FOR |Z|>15.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  INPUTS: 
C    IFUN:
C         * IFUN=1, THE CODE COMPUTES AI(Z)
C         * IFUN=2, THE CODE COMPUTES AI'(Z)
C    IFAC:
C         * IFAC=1, THE CODE COMPUTES  AI(Z) OR AI'(Z)
C         * IFAC=2, THE CODE COMPUTES NORMALIZED AI(Z) OR AI'(Z)
C    X0:   REAL PART OF THE ARGUMENT Z
C    Y0:   IMAGINARY PART OF THE ARGUMENT  Z
C
C  OUTPUTS:
C    GAIR: REAL PART OF AI(Z) OR AI'(Z)
C    GAII: IMAGINARY PART OF AI(Z) OR AI'(Z)
C    
C    IERRO: ERROR FLAG
C          * IERRO=0, SUCCESSFUL COMPUTATION       
C          * IERRO=1, COMPUTATION OUT OF RANGE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   ACCURACY:                                        
C
C     1) SCALED AIRY FUNCTIONS:
C        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO
C        THE ZEROS, WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
C        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z|>1000 
C        (REACHING 10**(-8) ABSOLUTE ACCURACY FOR |Z| CLOSE 
C        TO 10**(6)) IN THE CASE OF PHASE(Z) CLOSE TO PI.
C     2) UNSCALED AIRY FUNCTIONS:
C        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR 
C        3/2*|Z|**(3/2)>LOG(OVER).
C        FOR |Z|<30:
C        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE
C           ZEROS) BETTER THAN 10**(-13).
C        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER
C           THAN 10**(-13), WHERE R(Z)=REAL(AI)/IMAG(AI) 
C           OR R(Z)=REAL(AI')/IMAG(AI').
C        FOR |Z|>30, GRADUAL LOSS OF PRECISION TAKES PLACE
C        AS |Z| INCREASES.    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     AUTHORS:                                               
C        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN). 
C                      E-MAIL: AMPARO.GIL@UAM.ES
C        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
C                      E-MAIL: JSEGURA@MATH.UC3M.ES
C        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
C                      E-MAIL: NICO.TEMME@CWI.NL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    REFERENCES:                                                
C         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
C         NUMERICAL ALGORITHMS (2001).
C         A. GIL, J. SEGURA, N.M. TEMME
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION X0,Y0,GAIR,GAII,X,W,XD,WD 
      DOUBLE PRECISION OVER,UNDER,DL1,DL2,COVER,D1MACH
      DOUBLE PRECISION PI,PIHAL,PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,
     * R32,FACTO,TH025,S3,C3,F23,PI23,SQRT3,XA,YA,F23R,DF1,DF2,
     * S11,C11,DEX,DRE,DIMA,GAR,GAI,C,S,U,V,V0,AR,AI,AR1,AI1,
     * RO,COE1,COE2,REX,DFR,DFI,AR11,AI11,PHASE
      INTEGER IFUN,IFAC,IERRO,IEXPF,IEXPF2,N
      DIMENSION X(25),W(25)
      DIMENSION XD(25),WD(25)
      COMMON/PARAM1/PI,PIHAL
      COMMON/PARAM2/PIH3,PISR,A,ALF
      COMMON/PARAM3/THET,R,TH15,S1,C1,R32
      COMMON/PARAM4/FACTO,TH025,S3,C3
      SAVE X,W
      SAVE XD,WD
      DATA X,W/.283891417994567679D-1,.170985378860034935D0, 
     *.435871678341770460D0,.823518257913030858D0,1.33452543254227372D0,
     *1.96968293206435071D0,2.72998134002859938D0,3.61662161916100897D0,  
     *4.63102611052654146D0,5.77485171830547694D0,7.05000568630218682D0,
     *8.45866437513237792D0,10.0032955242749393D0,11.6866845947722423D0,  
     *13.5119659344693551D0,15.4826596959377140D0,17.6027156808069112D0,
     *19.8765656022785451D0,22.3091856773962780D0,24.9061720212974207D0,  
     *27.6738320739497190D0,30.6192963295084111D0,33.7506560850239946D0,
     *37.0771349708391198D0,40.6093049694341322D0,.143720408803313866D0, 
     *.230407559241880881D0,.242253045521327626D0,.203636639103440807D0,
     *.143760630622921410D0,.869128834706078120D-1,.4541750018329
     * 15883D-1,.206118031206069497D-1,.814278821268606972D-2,.280266
     *075663377634D-2,.840337441621719716D-3,.219303732907765020D-3,
     *.497401659009257760D-4,.978508095920717661D-5,.166542824603725
     *563D-5,.244502736801316287D-6,.308537034236207072D-7,.3332960
     *72940112245D-8,.306781892316295828D-9,.239331309885375719D-10,
     *.157294707710054952D-11,.864936011664392267D-13,.394819815
     *638647111D-14,.148271173082850884D-15,.453390377327054458D-17/  
      DATA XD,WD/.435079659953445D-1,.205779160144678D0,
     *.489916161318751D0,.896390483211727D0,1.42582496737580D0,
     *2.07903190767599D0,2.85702335104978D0,3.76102058198275D0,
     *4.79246521225895D0,5.95303247470003D0,7.24464710774066D0,
     *8.66950223642504D0,10.2300817341775D0,11.9291866622602D0,
     *13.7699665302828D0,15.7559563095946D0,17.8911203751898D0,
     *20.1799048700978D0,22.6273004064466D0,25.2389175786164D0,
     *28.0210785229929D0,30.9809287996116D0,34.1265753192057D0,
     *37.4672580871163D0,41.0135664833476D0,.576354557898966D-1, 
     *.139560003272262D0,.187792315011311D0,.187446935256946D0,
     *.150716717316301D0,.101069904453380D0,.575274105486025D-1,  
     *.280625783448681D-1,.117972164134041D-1,.428701743297432D-2,
     *.134857915232883D-2,.367337337105948D-3,.865882267841931D-4,  
     *.176391622890609D-4,.309929190938078D-5,.468479653648208D-6,
     *.607273267228907D-7,.672514812555074D-8,.633469931761606D-9,  
     *.504938861248542D-10,.338602527895834D-11,.189738532450555D-12,
     *.881618802142698D-14,.336676636121976D-15,.104594827170761D-16/     
CC CONSTANTS CCCCCCCCCCCCCCCCCCCCCCC    
      PI=3.1415926535897932385D0
      PIHAL=1.5707963267948966192D0      
      PIH3=4.71238898038469D0
      F23=.6666666666666666D0
      PI23=2.09439510239320D0
      PISR=1.77245385090552D0
      SQRT3=1.7320508075688772935D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      YA=Y0
      XA=X0
      IERRO=0
      IEXPF=0
      IEXPF2=0
      IF (YA.LT.0.D0) YA=-YA 
      R=SQRT(XA*XA+YA*YA)
      R32=R*SQRT(R)
      THET=PHASE(XA,YA)
      COVER=2.D0/3.D0*R32*ABS(COS(1.5D0*THET))
CCC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
      OVER=D1MACH(2)*1.D-3 
CCC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
      UNDER=D1MACH(1)*1.D+3 
      DL1=LOG(OVER)
      DL2=-LOG(UNDER)
      IF (DL1.GT.DL2) OVER=1/UNDER
      IF (IFAC.EQ.1) THEN
        IF (COVER.GE.LOG(OVER)) THEN
CCC OVERFLOW/UNDERFLOW PROBLEMS. 
CCC   CALCULATION ABORTED
          IERRO=1
          GAIR=0
          GAII=0
        ENDIF
        IF (COVER.GE.(LOG(OVER)*0.2)) IEXPF2=1
      ELSE
        IF (COVER.GE.(LOG(OVER)*0.2)) IEXPF=1
      ENDIF  
      IF (IERRO.EQ.0) THEN
        IF (IFUN.EQ.1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC CALCULATION OF AI(Z) CCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCC SERIES, INTEGRALS OR EXPANSIONS CCCCCCCCCCCCCCCCCCCCCCCC
          IF ((YA.LT.3.D0).AND.(XA.LT.1.3D0).AND.(XA.GT.-2.5D0)) THEN
CCC SERIES CCC
            CALL SERAI(XA,YA,GAR,GAI)
            IF (IFAC.EQ.2) THEN 
              THET=PHASE(XA,YA)         
              TH15=1.5D0*THET
              S1=SIN(TH15)
              C1=COS(TH15)
              F23R=F23*R32
              DF1=F23R*C1
              DF2=F23R*S1
              S11=SIN(DF2)
              C11=COS(DF2)
              DEX=EXP(DF1)
              DRE=DEX*C11
              DIMA=DEX*S11
              GAIR=DRE*GAR-DIMA*GAI
              GAII=DRE*GAI+DIMA*GAR
            ELSE
              GAIR=GAR
              GAII=GAI
              IF (Y0.EQ.0.) GAII=0.D0 
            ENDIF
          ELSE
            IF (R.GT.15.D0) THEN
CCC ASYMPTOTIC EXPANSIONS CCC 
              THET=PHASE(XA,YA)  
              FACTO=0.5D0/PISR*R**(-0.25D0)
              IF (THET.GT.PI23) THEN    
CCCCCCCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCC
CCC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
                N=1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL EXPAI(AR1,AI1)
                IF (V0.LT.0.D0) AI1=-AI1
                AR=-(C*AR1-S*AI1)
                AI=-(S*AR1+C*AI1)
                IF (IEXPF.EQ.0) THEN
                  IF (IEXPF2.EQ.0) THEN
CCC     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
                    N=-1
                    C=-0.5D0
                    S=N*0.5*SQRT3
                    U=XA*C-YA*S
                    V=XA*S+YA*C
                    V0=V
                    IF (V.LT.0.D0) V=-V
                    THET=PHASE(U,V)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    TH025=THET*0.25D0
                    S3=SIN(TH025)
                    C3=COS(TH025)
                    CALL EXPAI(AR1,AI1)
                    IF (V0.LT.0.D0) AI1=-AI1
                    THET=PHASE(XA,YA)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    RO=1.333333333333333D0*R32
                    COE1=RO*C1
                    COE2=RO*S1
                    REX=EXP(COE1)
                    DFR=REX*COS(COE2)
                    DFI=REX*SIN(COE2)
                    AR11=DFR*AR1-DFI*AI1
                    AI11=DFR*AI1+DFI*AR1
                    GAIR=AR-(C*AR11-S*AI11)
                    GAII=AI-(S*AR11+C*AI11) 
                  ELSE
                    THET=PHASE(XA,YA)
                    TH15=1.5D0*THET
                    S1=SIN(TH15)
                    C1=COS(TH15)
                    GAIR=AR
                    GAII=AI
                  ENDIF
                ELSE
                  GAIR=AR
                  GAII=AI
                ENDIF
              ELSE
CCCCCCC  ASYMPTOTIC EXPANSION CCCCCCCCCCCCCCC    	      
                THET=PHASE(XA,YA) 
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL EXPAI(GAIR,GAII)
              ENDIF
            ELSE
CCC INTEGRALS
              A=0.1666666666666666D0    
              ALF=-A      
              FACTO=0.280514117723058D0*R**(-0.25D0)
              THET=PHASE(XA,YA) 
              IF (THET.LE.PIHAL) THEN
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL AIRY1(X,W,GAIR,GAII)
              ENDIF
              IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)
                CALL AIRY2(X,W,GAIR,GAII)
              ENDIF   
              IF (THET.GT.PI23) THEN
                N=1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                IF (THET.LE.PIHAL) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY1(X,W,AR1,AI1)
                ENDIF
                IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY2(X,W,AR1,AI1)
                ENDIF
                IF (V0.LT.0.D0) AI1=-AI1
                AR=-(C*AR1-S*AI1)
                AI=-(S*AR1+C*AI1)
                N=-1               
                C=-0.5D0          
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                IF (THET.LE.PIHAL) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY1(X,W,AR1,AI1)
                ENDIF
                IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                  TH15=1.5D0*THET
                  S1=SIN(TH15)
                  C1=COS(TH15)
                  TH025=THET*0.25D0
                  S3=SIN(TH025)
                  C3=COS(TH025)
                  CALL AIRY2(X,W,AR1,AI1)
                ENDIF
                IF (V0.LT.0.D0) AI1=-AI1
                THET=PHASE(XA,YA)
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                RO=1.333333333333333D0*R32
                COE1=RO*C1
                COE2=RO*S1
                REX=EXP(COE1)
                DFR=REX*COS(COE2)
                DFI=REX*SIN(COE2)
                AR11=DFR*AR1-DFI*AI1
                AI11=DFR*AI1+DFI*AR1
                GAIR=AR-(C*AR11-S*AI11)
                GAII=AI-(S*AR11+C*AI11) 
              ENDIF
            ENDIF
            IF (IFAC.EQ.1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC CALCULATION OF THE UNSCALED AI(Z) CCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              F23R=F23*R32
              DF1=F23R*C1
              DF2=F23R*S1
              S11=SIN(DF2)
              C11=COS(DF2)
              DEX=EXP(-DF1)
              DRE=DEX*C11
              DIMA=-DEX*S11
              GAR=DRE*GAIR-DIMA*GAII
              GAI=DRE*GAII+DIMA*GAIR
              GAIR=GAR
              GAII=GAI
              IF (Y0.EQ.0.) GAII=0.D0
            ENDIF
          ENDIF
        ELSE
CCCC CALCULATION OF AI´(Z) CCCCCCCCCCC 
          ALF=0.1666666666666666D0        
          FACTO=-0.270898621247918D0*R**0.25D0   
CCCCCCCCCCCCCCC SERIES OR INTEGRALS CCCCCCCCCCCCCCCCCCCCCCCCCC
          IF ((YA.LT.3.D0).AND.(XA.LT.1.3D0).AND.(XA.GT.-2.5D0)) THEN
CCC SERIES
            CALL SERAID(XA,YA,GAR,GAI)
            IF (IFAC.EQ.2) THEN 
              THET=PHASE(XA,YA)         
              TH15=1.5D0*THET
              S1=SIN(TH15)
              C1=COS(TH15)
              F23R=F23*R32
              DF1=F23R*C1
              DF2=F23R*S1
              S11=SIN(DF2)
              C11=COS(DF2)
              DEX=EXP(DF1)
              DRE=DEX*C11
              DIMA=DEX*S11
              GAIR=DRE*GAR-DIMA*GAI
              GAII=DRE*GAI+DIMA*GAR
            ELSE
              GAIR=GAR
              GAII=GAI
              IF (Y0.EQ.0.) GAII=0.D0
            ENDIF
          ELSE
            IF (R.GT.15.D0) THEN
CCC  ASYMPTOTIC EXPANSIONS CCCCCCCCCCCCC 
              THET=PHASE(XA,YA) 
              FACTO=0.5D0/PISR*R**0.25D0
              IF (THET.GT.PI23) THEN
CCCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
                N=1
                C=-0.5D0
                S=N*0.5*SQRT3
                U=XA*C-YA*S
                V=XA*S+YA*C
                V0=V
                IF (V.LT.0.D0) V=-V
                THET=PHASE(U,V)
                TH15=1.5D0*THET
                S1=SIN(TH15)
                C1=COS(TH15)
                TH025=THET*0.25D0
                S3=SIN(TH025)
                C3=COS(TH025)              
                CALL EXPAID(AR1,AI1)         
                IF (V0.LT.0.D0) AI1=-AI1         
                AR=-(C*AR1+S*AI1)            
                AI=-(-S*AR1+C*AI1)
                IF (IEXPF.EQ.0) THEN
                  IF (IEXPF2.EQ.0) THEN
CCC     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
                    N=-1
                    C=-0.5D0         
                    S=N*0.5*SQRT3           
                    U=XA*C-YA*S           
                    V=XA*S+YA*C             
                    V0=V           
                    IF (V.LT.0.D0) V=-V             
                    THET=PHASE(U,V)             
                    TH15=1.5D0*THET             
                    S1=SIN(TH15)              
                    C1=COS(TH15)               
                    TH025=THET*0.25D0            
                    S3=SIN(TH025)               
                    C3=COS(TH025)            
                    CALL EXPAID(AR1,AI1)              
                    IF (V0.LT.0.D0) AI1=-AI1               
                    THET=PHASE(XA,YA)            
                    TH15=1.5D0*THET            
                    S1=SIN(TH15)           
                    C1=COS(TH15)          
                    RO=1.333333333333333D0*R32           
                    COE1=RO*C1            
                    COE2=RO*S1              
                    REX=EXP(COE1)              
                    DFR=REX*COS(COE2)             
                    DFI=REX*SIN(COE2)             
                    AR11=DFR*AR1-DFI*AI1              
                    AI11=DFR*AI1+DFI*AR1               
                    GAIR=AR-(C*AR11+S*AI11)               
                    GAII=AI-(-S*AR11+C*AI11)                
                  ELSE                 
                    THET=PHASE(XA,YA)            
                    TH15=1.5D0*THET             
                    S1=SIN(TH15)              
                    C1=COS(TH15)           
                    GAIR=AR          
                    GAII=AI                 
                  ENDIF                
                ELSE                   
                  GAIR=AR               
                  GAII=AI               
                ENDIF                 
              ELSE          
                TH15=1.5D0*THET            
                S1=SIN(TH15)                
                C1=COS(TH15)                 
                TH025=THET*0.25D0               
                S3=SIN(TH025)             
                C3=COS(TH025)                  
                CALL EXPAID(GAIR,GAII)               
              ENDIF                  
            ELSE               
CCC INTEGRALS CCCCCCCCCCCCCCCC
              THET=PHASE(XA,YA)                 
              IF (THET.LE.PIHAL) THEN                 
                TH15=1.5D0*THET               
                S1=SIN(TH15)                
                C1=COS(TH15)            
                TH025=THET*0.25D0             
                S3=SIN(TH025)                 
                C3=COS(TH025)                 
                CALL AIRY1D(XD,WD,GAIR,GAII)                
              ENDIF              
              IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN             
                TH15=1.5D0*THET              
                S1=SIN(TH15)                 
                C1=COS(TH15)                     
                TH025=THET*0.25D0                 
                S3=SIN(TH025)                
                C3=COS(TH025)             
                CALL AIRY2D(XD,WD,GAIR,GAII)               
              ENDIF              
              IF (THET.GT.PI23) THEN
                N=1
                C=-0.5D0               
                S=N*0.5*SQRT3            
                U=XA*C-YA*S         
                V=XA*S+YA*C              
                V0=V               
                IF (V.LT.0.D0) V=-V             
                THET=PHASE(U,V)                 
                IF (THET.LE.PIHAL) THEN
                  TH15=1.5D0*THET                 
                  S1=SIN(TH15)                
                  C1=COS(TH15)              
                  TH025=THET*0.25D0             
                  S3=SIN(TH025)            
                  C3=COS(TH025)            
                  CALL AIRY1D(XD,WD,AR1,AI1)          
                ENDIF             
                IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                  TH15=1.5D0*THET         
                  S1=SIN(TH15)            
                  C1=COS(TH15)               
                  TH025=THET*0.25D0           
                  S3=SIN(TH025)             
                  C3=COS(TH025)              
                  CALL AIRY2D(XD,WD,AR1,AI1)              
                ENDIF                    
                IF (V0.LT.0.D0) AI1=-AI1              
                  AR=-(C*AR1+S*AI1)                
                  AI=-(-S*AR1+C*AI1)
                  N=-1                 
                  C=-0.5D0               
                  S=N*0.5*SQRT3                 
                  U=XA*C-YA*S               
                  V=XA*S+YA*C              
                  V0=V                
                  IF (V.LT.0.D0) V=-V           
                  THET=PHASE(U,V)                
                  IF (THET.LE.PIHAL) THEN                                    
                    TH15=1.5D0*THET             
                    S1=SIN(TH15)               
                    C1=COS(TH15)               
                    TH025=THET*0.25D0             
                    S3=SIN(TH025)               
                    C3=COS(TH025)             
                    CALL AIRY1D(XD,WD,AR1,AI1)            
                  ENDIF                                        
                  IF ((THET.GT.PIHAL).AND.(THET.LE.PI23)) THEN
                    TH15=1.5D0*THET              
                    S1=SIN(TH15)        
                    C1=COS(TH15)            
                    TH025=THET*0.25D0         
                    S3=SIN(TH025)        
                    C3=COS(TH025)              
                    CALL AIRY2D(XD,WD,AR1,AI1)          
                  ENDIF                   
                  IF (V0.LT.0.D0) AI1=-AI1                 
                  THET=PHASE(XA,YA)              
                  TH15=1.5D0*THET           
                  S1=SIN(TH15)             
                  C1=COS(TH15)              
                  RO=1.333333333333333D0*R32            
                  COE1=RO*C1           
                  COE2=RO*S1              
                  REX=EXP(COE1)           
                  DFR=REX*COS(COE2)            
                  DFI=REX*SIN(COE2)               
                  AR11=DFR*AR1-DFI*AI1            
                  AI11=DFR*AI1+DFI*AR1             
                  GAIR=AR-(C*AR11+S*AI11)          
                  GAII=AI-(-S*AR11+C*AI11)               
                ENDIF              
              ENDIF              
              IF (IFAC.EQ.1) THEN                                        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC CALCULATION OF THE UNSCALED AI'(z) CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                F23R=F23*R32
                DF1=F23R*C1
                DF2=F23R*S1
                S11=SIN(DF2)
                C11=COS(DF2)
                DEX=EXP(-DF1)
                DRE=DEX*C11
                DIMA=-DEX*S11
                GAR=DRE*GAIR-DIMA*GAII
                GAI=DRE*GAII+DIMA*GAIR
                GAIR=GAR
                GAII=GAI
                IF (Y0.EQ.0) GAII=0.D0
              ENDIF                                           
            ENDIF
          ENDIF
        ENDIF
        IF (Y0.LT.0.D0) GAII=-GAII
        RETURN
        END
        SUBROUTINE  AIRY1(X,W,GAIR,GAII)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
CCC              0 <= PHASE(Z) <= PI/2                        C
CCC                                                           C
CCC INPUTS:                                                   C
CCC      X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN         C
CCC                QUADRATURE                                 C
CCC OUTPUTS:                                                  C
CCC      GAIR, GAII,  REAL AND IMAGINARY PARTS OF AI(Z)       C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION X,W,GAIR,GAII
        DOUBLE PRECISION PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,
     *  R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,
     *  S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2,PHASE
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SUMAR=0.D0
        SUMAI=0.D0
        DO 1 I=1,25
          DF1=1.5D0*X(I)/R32
          DF1C1=DF1*C1
          PHI=PHASE(2.D0+DF1C1,DF1*S1)
          PHI6=PHI/6.D0
          S2=SIN(PHI6)
          C2=COS(PHI6)
          DMODU=SQRT(4.D0+DF1*DF1+4.D0*DF1C1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMAR=SUMAR+W(I)*FUNR
          SUMAI=SUMAI+W(I)*FUNI
 1      CONTINUE
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR+FAC2*SUMAI
        GAII=FAC1*SUMAI-FAC2*SUMAR
        RETURN 
        END
        SUBROUTINE  AIRY2(X,W,GAIR,GAII)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
CCC              PI/2 < PHASE(Z) <= 2PI/3                     C
CCC                                                           C
CCC INPUTS:                                                   C
CCC      X,W,        NODES AND WEIGHTS FOR THE GAUSSIAN       C
CCC                  QUADRATURE                               C
CCC OUTPUTS:                                                  C
CCC      GAIR, GAII, REAL AND IMAGINARY PARTS OF AI(Z)        C       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION X,W,GAIR,GAII
        DOUBLE PRECISION PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,
     *  R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,
     *  S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2,PHASE
        DOUBLE PRECISION SQR2,SQR2I,TAU,TGTAU,B,ANG,CTAU,CFAC,CT,ST,
     *  SUMR,SUMI,TTAU,BETA
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SQR2=1.41421356237310D0
        SQR2I=0.707106781186548D0
        TAU=TH15-PIH3*0.5D0
        TGTAU=DTAN(TAU)
        B=5.D0*A
        ANG=TAU*B
        CTAU=COS(TAU)
        CFAC=CTAU**(-B)
        CT=COS(ANG)
        ST=SIN(ANG)
        SUMR=0.D0
        SUMI=0.D0
        DO 2 I=1,25
          DF1=3.D0*X(I)/(CTAU*R32)
          DF1C1=DF1*SQR2I*0.5D0
          PHI=PHASE(2.D0-DF1C1,DF1C1)
          PHI6=PHI/6.D0
          TTAU=X(I)*TGTAU
          BETA=PHI6-TTAU
          S2=SIN(BETA)
          C2=COS(BETA)
          DMODU=SQRT(4.D0+DF1*DF1*0.25D0-SQR2*DF1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMR=SUMR+W(I)*FUNR
          SUMI=SUMI+W(I)*FUNI
 2      CONTINUE
        SUMAR=CFAC*(CT*SUMR-ST*SUMI)
        SUMAI=CFAC*(CT*SUMI+ST*SUMR)   
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR+FAC2*SUMAI
        GAII=FAC1*SUMAI-FAC2*SUMAR
        RETURN 
        END
        SUBROUTINE AIRY1D(X,W,GAIR,GAII)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
CCC              0 <= PHASE(Z) <= PI/2                         C
CCC                                                            C
CCC INPUTS:                                                    C
CCC       X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN         C
CCC                 QUADRATURE                                 C
CCC OUTPUTS:                                                   C
CCC       GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)        C  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION X,W,GAIR,GAII
        DOUBLE PRECISION PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,
     *  R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,
     *  S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2,PHASE
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SUMAR=0.D0
        SUMAI=0.D0
        DO 3 I=1,25
          DF1=1.5D0*X(I)/R32
          DF1C1=DF1*C1
          PHI=PHASE(2.D0+DF1C1,DF1*S1)
          PHI6=-PHI*ALF
          S2=SIN(PHI6)
          C2=COS(PHI6)
          DMODU=SQRT(4.D0+DF1*DF1+4.D0*DF1C1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMAR=SUMAR+W(I)*FUNR
          SUMAI=SUMAI+W(I)*FUNI
 3      CONTINUE
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR-FAC2*SUMAI
        GAII=FAC1*SUMAI+FAC2*SUMAR
        RETURN 
        END  
        SUBROUTINE  AIRY2D(X,W,GAIR,GAII)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
CCC              PI/2 < PHASE(Z) <= 3PI/2                      C
CCC                                                            C
CCC INPUTS:                                                    C
CCC      X,W,   NODES AND WEIGHTS FOR THE GAUSSIAN             C
CCC                   QUADRATURE                               C
CCC OUTPUTS:                                                   C
CCC      GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION X,W,GAIR,GAII
        DOUBLE PRECISION PIH3,PISR,A,ALF,THET,R,TH15,S1,C1,
     *  R32,FACTO,TH025,S3,C3,SUMAR,SUMAI,DF1,DF1C1,PHI,PHI6,
     *  S2,C2,DMODU,DMODU2,FUNR,FUNI,FAC1,FAC2,PHASE
        DOUBLE PRECISION SQR2,SQR2I,TAU,TGTAU,B,ANG,CTAU,CFAC,CT,ST,
     *  SUMR,SUMI,TTAU,BETA
        INTEGER I
        DIMENSION X(25),W(25)
        COMMON/PARAM2/PIH3,PISR,A,ALF
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SQR2=1.41421356237310D0
        SQR2I=0.707106781186548D0
        TAU=TH15-PIH3*0.5D0
        TGTAU=DTAN(TAU)
        B=7.D0*ALF
        ANG=TAU*B
        CTAU=COS(TAU)
        CFAC=CTAU**(-B)
        CT=COS(ANG)
        ST=SIN(ANG)
        SUMR=0.D0
        SUMI=0.D0
        DO 4 I=1,25
          DF1=3.D0*X(I)/(CTAU*R32)
          DF1C1=DF1*SQR2I*0.5D0
          PHI=PHASE(2.D0-DF1C1,DF1C1)
          PHI6=-PHI/6.D0
          TTAU=X(I)*TGTAU
          BETA=PHI6-TTAU
          S2=SIN(BETA)
          C2=COS(BETA)
          DMODU=SQRT(4.D0+DF1*DF1*0.25D0-SQR2*DF1)
          DMODU2=DMODU**ALF
          FUNR=DMODU2*C2
          FUNI=DMODU2*S2
          SUMR=SUMR+W(I)*FUNR
          SUMI=SUMI+W(I)*FUNI
 4      CONTINUE
        SUMAR=CFAC*(CT*SUMR-ST*SUMI)
        SUMAI=CFAC*(CT*SUMI+ST*SUMR)   
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR-FAC2*SUMAI
        GAII=FAC1*SUMAI+FAC2*SUMAR
        RETURN 
        END
        DOUBLE PRECISION FUNCTION PHASE(X,Y)
        DOUBLE PRECISION PI,PIHAL,X,Y,AY,P
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC  COMPUTES THE PHASE OF Z = X + IY, IN (-PI,PI]
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMMON/PARAM1/PI,PIHAL
        IF ((X.EQ.0).AND.(Y.EQ.0)) THEN
          P=0.D0
        ELSE
          AY=ABS(Y) 
          IF (X.GE.AY) THEN
            P=ATAN(AY/X)
          ELSEIF ((X+AY).GE.0.D0) THEN
            P=PIHAL-ATAN(X/AY)
          ELSE
            P=PI+ATAN(AY/X)
          ENDIF
          IF (Y.LT.0.D0) P=-P
        ENDIF
        PHASE=P
        END        
        SUBROUTINE FGP(X,Y,EPS,FR,FI,GR,GI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    COMPUTES THE FUNCTIONS F AND G FOR THE SERIES  C
C    OF AI'(Z).                                     C 
C    THIS ROUTINE IS CALLED BY SERAID.              C
C                                                   C
C    INPUTS:                                        C
C         X,Y,  REAL AND IMAGINARY PARTS OF Z       C
C         EPS,  PRECISION FOR THE COMPUTATION OF    C
C               THE SERIES                          C
C    OUTPUTS:                                       C
C         FR,FI, REAL AND IMAGINARY PARTS OF F      C
C         GR,GI, REAL AND IMAGINARY PARTS OF G      C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION X,Y,EPS,FR,FI,GR,GI
        INTEGER A,B,K3
        DOUBLE PRECISION X2,Y2,U,V,P,Q,CR,CI,DR,DI
        X2=X*X
        Y2=Y*Y
        K3=0
        U=X*(X2-3*Y2)
        V=Y*(3*X2-Y2)
        CR=0.5D0
        CI=0.D0
        DR=1.D0
        DI=0.D0
        FR=0.5D0
        FI=0.D0
        GR=1.D0
        GI=0.D0                
 70     A=(K3+5)*(K3+3)
        B=(K3+1)*(K3+3)
        P=(U*CR-V*CI)/A
        Q=(V*CR+U*CI)/A
        CR=P
        CI=Q
        P=(U*DR-V*DI)/B
        Q=(V*DR+U*DI)/B
        DR=P
        DI=Q
        FR=FR+CR
        FI=FI+CI
        GR=GR+DR
        GI=GI+DI
        K3=K3+3  
        IF ((ABS(CR)+ABS(DR)+ABS(CI)+ABS(DI)).GE.EPS) GOTO 70          
        U=X2-Y2
        V=2.D0*X*Y
        P=U*FR-V*FI
        Q=U*FI+V*FR
        FR=P
        FI=Q
        RETURN
        END
        SUBROUTINE FG(X,Y,EPS,FR,FI,GR,GI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    COMPUTES THE FUNCTIONS F AND G IN EXPRESSION   C
C    10.4.2 OF ABRAMOWITZ & STEGUN FOR THE SERIES   C
C    OF AI(Z).                                      C 
C    THIS ROUTINE IS CALLED BY SERAI.               C
C                                                   C
C    INPUTS:                                        C
C          X,Y,  REAL AND IMAGINARY PARTS OF Z      C
C          EPS,  PRECISION FOR THE COMPUTATION      C
C                OF THE SERIES.                     C
C    OUTPUTS:                                       C
C          FR,FI, REAL AND IMAGINARY PARTS OF F     C
C          GR,GI, REAL AND IMAGINARY PARTS OF G     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        INTEGER A,B,K3
        DOUBLE PRECISION X2,Y2,U,V,P,Q,CR,CI,DR,DI
        DOUBLE PRECISION X,Y,EPS,FR,FI,GR,GI
        X2=X*X
        Y2=Y*Y
        K3=0
        U=X*(X2-3.D0*Y2)
        V=Y*(3.D0*X2-Y2)
        CR=1.D0
        CI=0.D0
        DR=1.D0
        DI=0.D0
        FR=1.D0
        FI=0.D0
        GR=1.D0
        GI=0.D0
 71     A=(K3+2)*(K3+3)
        B=(K3+4)*(K3+3)
        P=(U*CR-V*CI)/A
        Q=(V*CR+U*CI)/A
        CR=P
        CI=Q
        P=(U*DR-V*DI)/B
        Q=(V*DR+U*DI)/B
        DR=P
        DI=Q
        FR=FR+CR
        FI=FI+CI
        GR=GR+DR
        GI=GI+DI
        K3=K3+3 
        IF ((ABS(CR)+ABS(DR)+ABS(CI)+ABS(DI)).GE.EPS) GOTO 71
        P=X*GR-Y*GI
        Q=X*GI+Y*GR
        GR=P
        GI=Q
        RETURN
        END
        SUBROUTINE SERAI(X,Y,AIR,AII)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   AIRY AI(Z), TAYLOR, COMPLEX Z      CCC
CCC                                      CCC
CCC   INPUTS:                            CCC
CCC        X,Y,    REAL AND IMAGINARY    CCC
CCC                PARTS OF Z            CCC 
CCC   OUTPUTS:                           CCC
CCC        AIR,AII, REAL AND IMAGINARY   CCC
CCC                 PARTS OF AI(Z)       CCC               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION X,Y,EPS,AIR,AII
        DOUBLE PRECISION FZR,FZI,GZR,GZI,CONS1,CONS2
        DOUBLE PRECISION D1MACH
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15
        CONS1=0.355028053887817239260D0
        CONS2=0.258819403792806798405D0
        CALL FG(X,Y,EPS,FZR,FZI,GZR,GZI)
        AIR=CONS1*FZR-CONS2*GZR
        AII=CONS1*FZI-CONS2*GZI
        RETURN
        END
        SUBROUTINE SERAID(X,Y,AIR,AII)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   AIRY AI'(Z), TAYLOR, COMPLEX Z     CCC
CCC                                      CCC
CCC   INPUTS:                            CCC
CCC        X,Y,   REAL AND IMAGINARY     CCC
CCC               PARTS OF Z             CCC
CCC   OUTPUTS:                           CCC
CCC        AIR,AII, REAL AND IMAGINARY   CCC
CCC                 PARTS OF AI'(Z)      CCC                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
        DOUBLE PRECISION X,Y,EPS,AIR,AII
        DOUBLE PRECISION FZR,FZI,GZR,GZI,CONS1,CONS2
        DOUBLE PRECISION D1MACH  
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15
        CONS1=0.355028053887817239260D0
        CONS2=0.258819403792806798405D0
        CALL FGP(X,Y,EPS,FZR,FZI,GZR,GZI)
        AIR=CONS1*FZR-CONS2*GZR
        AII=CONS1*FZI-CONS2*GZI 
        RETURN
        END
        SUBROUTINE EXPAI(GAIR,GAII)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   AIRY AI(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
CCC                                               CCC
CCC   OUTPUTS:                                    CCC
CCC        GAIR, GAII,  REAL AND IMAGINARY        CCC
CCC                     PARTS OF AI(Z)            CCC   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION EPS,GAIR,GAII
        DOUBLE PRECISION THET,R,TH15,S1,
     *  C1,R32,FACTO,TH025,S3,C3
        DOUBLE PRECISION DF1,PSIIR,PSIII,CK,DFRR,DFII,SUMAR,SUMAI,
     *  DFR,DFI,DELTAR,DELTAI,FAC1,FAC2
        DOUBLE PRECISION CO,DF
        DOUBLE PRECISION D1MACH
        INTEGER K
        DIMENSION CO(20)
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3         
        SAVE CO  
        DATA CO/-6.944444444444445D-2,3.713348765432099D-2,
     *  -3.799305912780064D-2,5.764919041266972D-2,-0.116099064025515D0,    
     *  0.291591399230751D0,-0.877666969510017D0,3.07945303017317D0,   
     *  -12.3415733323452D0,55.6227853659171D0,-278.465080777603D0,
     *  1533.16943201280D0,-9207.20659972641D0,59892.5135658791D0,     
     *  -419524.875116551D0,3148257.41786683D0,-25198919.8716024D0,    
     *  214288036.963680D0,-1929375549.18249D0,18335766937.8906D0/
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15  
        DF1=1.5D0/R32
        PSIIR=C1
        PSIII=-S1
        K=0
        CK=1.D0
        DF=1.D0
        DFRR=1.D0
        DFII=0.D0
        SUMAR=1.D0
        SUMAI=0.D0
80      DF=DF*DF1
        CK=CO(K+1)*DF
        DFR=DFRR
        DFI=DFII    
        DFRR=DFR*PSIIR-DFI*PSIII
        DFII=DFR*PSIII+DFI*PSIIR
        DELTAR=DFRR*CK
        DELTAI=DFII*CK
        SUMAR=SUMAR+DELTAR
        SUMAI=SUMAI+DELTAI
        K=K+1
        IF (SUMAR.NE.0) THEN
          IF (ABS(DELTAR/SUMAR).GT.EPS)  GOTO 80
        ENDIF
        IF (SUMAI.NE.0) THEN
          IF (ABS(DELTAI/SUMAI).GT.EPS) GOTO 80
        ENDIF 
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=FAC1*SUMAR+FAC2*SUMAI
        GAII=FAC1*SUMAI-FAC2*SUMAR
        RETURN                                                 
        END
        SUBROUTINE EXPAID(GAIR,GAII)                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   AIRY AI'(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
CCC                                                CCC
CCC   OUTPUTS:                                     CCC   
CCC        GAIR, GAII,  REAL AND IMAGINARY         CCC
CCC                     PARTS OF AI'(Z)            CCC  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION EPS,GAIR,GAII
        DOUBLE PRECISION THET,R,TH15,S1,
     *  C1,R32,FACTO,TH025,S3,C3
        DOUBLE PRECISION DF1,PSIIR,PSIII,VK,DFRR,DFII,SUMAR,SUMAI,
     *  DFR,DFI,DELTAR,DELTAI,FAC1,FAC2
        DOUBLE PRECISION CO,DF
        DOUBLE PRECISION D1MACH
        INTEGER K
        DIMENSION CO(20)
        COMMON/PARAM3/THET,R,TH15,S1,C1,R32
        COMMON/PARAM4/FACTO,TH025,S3,C3
        SAVE CO   
        DATA CO/9.722222222222222D-2,-4.388503086419753D-2,
     *  4.246283078989484D-2,-6.266216349203230D-2,
     *  0.124105896027275D0,-0.308253764901079D0,     
     *  0.920479992412945D0,-3.21049358464862D0,     
     *  12.8072930807356D0,-57.5083035139143D0,     
     *  287.033237109221D0,-1576.35730333710D0,     
     *  9446.35482309593D0,-61335.7066638521D0,     
     *  428952.400400069D0,-3214536.52140086D0,     
     *  25697908.3839113D0,-218293420.832160D0,     
     *  1963523788.99103D0,-18643931088.1072D0/
        EPS=D1MACH(3)
        IF (EPS.LT.1.D-15) EPS=1.D-15                                 
        DF1=1.5D0/R32                                    
        PSIIR=C1
        PSIII=-S1
        K=0
        DF=1.D0
        DFRR=1.D0
        DFII=0.D0
        SUMAR=1.D0
        SUMAI=0.D0
 81     DF=DF*DF1
        VK=CO(K+1)*DF
        DFR=DFRR
        DFI=DFII    
        DFRR=DFR*PSIIR-DFI*PSIII
        DFII=DFR*PSIII+DFI*PSIIR
        DELTAR=DFRR*VK
        DELTAI=DFII*VK
        SUMAR=SUMAR+DELTAR
        SUMAI=SUMAI+DELTAI
        K=K+1
        IF (SUMAR.NE.0) THEN
          IF (ABS(DELTAR/SUMAR).GT.EPS)  GOTO 81
        ENDIF
        IF (SUMAI.NE.0) THEN
          IF (ABS(DELTAI/SUMAI).GT.EPS) GOTO 81
        ENDIF 
        FAC1=FACTO*C3
        FAC2=FACTO*S3
        GAIR=-(FAC1*SUMAR-FAC2*SUMAI)
        GAII=-(FAC1*SUMAI+FAC2*SUMAR)
        RETURN
        END

       SUBROUTINE BIZ(IFUN,IFAC,X0,Y0,GBIR,GBII,IERRO)   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C COMPUTATION OF THE AIRY FUNCTION BI(Z) OR ITS DERIVATIVE BI'(Z)
C THE CODE USES THE CONNECTION OF BI(Z) WITH AI(Z). 
C                BIZ CALLS THE ROUTINE AIZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  INPUTS: 
C    IFUN:
C         * IFUN=1, THE CODE COMPUTES BI(Z)
C         * IFUN=2, THE CODE COMPUTES BI'(Z)
C    IFAC:
C         * IFAC=1, THE CODE COMPUTES  BI(Z) OR BI'(Z)
C         * IFAC=2, THE CODE COMPUTES NORMALIZED BI(Z) OR BI'(Z)
C    X0:   REAL PART OF THE ARGUMENT Z
C    Y0:   IMAGINARY PART OF THE ARGUMENT  Z
C
C  OUTPUTS:
C    GBIR: REAL PART OF BI(Z) OR BI'(Z)
C    GBII: IMAGINARY PART OF BI(Z) OR BI'(Z)  
C   
C    IERRO: ERROR FLAG
C          * IERRO=0, SUCCESSFUL COMPUTATION       
C          * IERRO=1, COMPUTATION OUT OF RANGE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   ACCURACY:                                        
C
C     1) SCALED AIRY FUNCTIONS:
C        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO
C        THE ZEROS, WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
C        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z|>1000 
C        IN THE CASE OF PHASE(Z) CLOSE TO +3*PI/2 OR -3*PI/2.
C     2) UNSCALED AIRY FUNCTIONS:
C        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR 
C        3/2*|Z|**(3/2)>LOG(OVER).
C        FOR |Z|<30:
C        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE
C           ZEROS) BETTER THAN 10**(-13).
C        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER
C           THAN 10**(-13), WHERE R(Z)=REAL(BI)/IMAG(BI) 
C           OR R(Z)=REAL(BI')/IMAG(BI').
C        FOR |Z|>30, GRADUAL LOSS OF PRECISION TAKES PLACE
C        AS |Z| INCREASES.    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     AUTHORS:                                               
C        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN). 
C                      E-MAIL: AMPARO.GIL@UAM.ES
C        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
C                      E-MAIL: JSEGURA@MATH.UC3M.ES
C        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
C                      E-MAIL: NICO.TEMME@CWI.NL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    REFERENCES:                                                
C         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
C         NUMERICAL ALGORITHMS (2001).
C         A. GIL, J. SEGURA, N.M. TEMME
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DOUBLE PRECISION X0,Y0,GBIR,GBII 
       DOUBLE PRECISION OVER,UNDER,DL1,DL2,COVER,D1MACH     
       DOUBLE PRECISION PI,PI3,PI23,SQRT3,C,S,C1,S1,U,V,X,Y,AR,AI,
     * APR,API,BR,BI,BPR,BPI,BBR,BBI,BBPR,BBPI,PHASE
       DOUBLE PRECISION THET,R,R32,THET32,A1,B1,DF1,EXPO,EXPOI,
     * ETAR,ETAI,ETAGR,ETAGI,PIHAL  
       INTEGER IFUN,IFAC,IEXPF,IERR,IERRO
       COMMON/PARAM1/PI,PIHAL
       SQRT3=1.7320508075688772935D0
       PI=3.1415926535897932385D0
       PIHAL=1.5707963267948966192D0
       PI3=PI/3.D0
       PI23=2.D0*PI3
       X=X0
       C=0.5D0*SQRT3
       S=0.5D0       
       IERRO=0
       IEXPF=0
       IF (Y0.LT.0.D0) THEN
         Y=-Y0
       ELSE
         Y=Y0
       ENDIF
       R=SQRT(X*X+Y*Y)
       R32=R*SQRT(R)
       THET=PHASE(X,Y)  
       COVER=2.D0/3.D0*R32*ABS(COS(1.5D0*THET))
CCC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
       OVER=D1MACH(2)*1.D-3
CCC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
       UNDER=D1MACH(1)*1.D+3 
       DL1=LOG(OVER)
       DL2=-LOG(UNDER)
       IF (DL1.GT.DL2) OVER=1/UNDER     
       IF (IFAC.EQ.1) THEN
         IF (COVER.GE.LOG(OVER)) THEN
CCC OVERFLOW/UNDERFLOW PROBLEMS. 
CCC   CALCULATION ABORTED
           IERRO=1
           GBIR=0
           GBII=0
         ENDIF
       ELSE
         IF (COVER.GE.(LOG(OVER)*0.2)) IEXPF=1
       ENDIF 
       IF (IERRO.EQ.0) THEN
         IF (IFAC.EQ.1) THEN
           IF (Y.EQ.0.D0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC  TAKE TWICE THE REAL PART OF EXP(-PI I/6) AI_(1)(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             C1=-0.5D0
             S1=-0.5D0*SQRT3
             U=X*C1-Y*S1
             V=X*S1+Y*C1
             CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
             IF (IFUN.EQ.1) THEN
               BR=SQRT3*AR+AI
               BI=0.D0
             ELSE
               U=AR*C1-AI*S1
               V=AR*S1+AI*C1
               APR=U
               API=V
               BPR=SQRT3*APR+API
               BPI=0.D0
             ENDIF
           ELSE
             IF ((X.LT.0.D0).AND.(Y.LT.-X*SQRT3)) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC      2 PI/3  < PHASE(Z) < PI
CCC      BI(Z)=EXP(I PI/6) AI_(-1)(Z) + EXP(-I PI/6) AI_(1)(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               C1=-0.5D0
               S1=0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1 
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BR=C*AR-S*AI
                 BI=C*AI+S*AR
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 BPR=C*APR-S*API
                 BPI=C*API+S*APR
               ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   WE NEED ALSO AI_(1)(Z) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               C1=-0.5D0
               S1=-0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1  
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 S=-S
                 BR=BR+C*AR-S*AI
                 BI=BI+C*AI+S*AR
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 S=-S
                 BPR=BPR+C*APR-S*API
                 BPI=BPI+C*API+S*APR
               ENDIF
             ELSE   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   BI(Z) = I AI(Z) + 2 EXP(-I PI/6) AI_(1)(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               C1=-0.5D0
               S1=-0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1 
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)  
               IF (IFUN.EQ.1) THEN
                 BR=SQRT3*AR+AI
                 BI=-AR+SQRT3*AI
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 BPR=SQRT3*APR+API
                 BPI=-APR+SQRT3*API
               ENDIF
               CALL AIZ(IFUN,IFAC,X,Y,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BR=BR-AI
                 BI=BI+AR
               ELSE
                 BPR=BPR-AI
                 BPI=BPI+AR
               ENDIF
             ENDIF
           ENDIF
         ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC         SCALED BI AIRY FUNCTIONS         C
CCC   WE USE THE FOLLOWING NORMALIZATION:    C
CCC   LET ARGZ=ARG(Z), THEN:                 C
CCC   A) IF  0 <= ARGZ <= PI/3               C
CCC      BI=EXP(-2/3Z^3/2)BI                 C 
CCC   B) IF  PI/3 <= ARGZ <= PI              C
CCC      BI=EXP(2/3Z^3/2)BI                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           THET=PHASE(X,Y)
           IF (THET.LE.PI3) THEN
             C1=-0.5D0
             S1=-0.5D0*SQRT3
             U=X*C1-Y*S1
             V=X*S1+Y*C1 
             CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)  
             IF (IFUN.EQ.1) THEN
               BR=SQRT3*AR+AI
               BI=-AR+SQRT3*AI
             ELSE
               U=AR*C1-AI*S1
               V=AR*S1+AI*C1
               APR=U
               API=V
               BPR=SQRT3*APR+API
               BPI=-APR+SQRT3*API
             ENDIF
             IF (IEXPF.EQ.0) THEN   
               R=SQRT(X*X+Y*Y)
               R32=R*SQRT(R)
               THET32=THET*1.5D0
               A1=COS(THET32)
               B1=SIN(THET32)
               DF1=4.D0/3.D0*R32         
               EXPO=EXP(DF1*A1)
               EXPOI=1.D0/EXPO
               ETAR=EXPO*COS(DF1*B1)
               ETAI=EXPO*SIN(DF1*B1)
               ETAGR=EXPOI*COS(-DF1*B1)
               ETAGI=EXPOI*SIN(-DF1*B1)
               CALL AIZ(IFUN,IFAC,X,Y,AR,AI,IERR)
               IF (IFUN.EQ.1) THEN
                 BR=BR-AR*ETAGI-ETAGR*AI
                 BI=BI+AR*ETAGR-ETAGI*AI
               ELSE
                 BPR=BPR-AR*ETAGI-ETAGR*AI
                 BPI=BPI+AR*ETAGR-ETAGI*AI
               ENDIF
             ENDIF 
           ENDIF
           IF ((THET.GT.PI3).AND.(THET.LE.PI23)) THEN
             IF (IEXPF.EQ.0) THEN
               R=SQRT(X*X+Y*Y)
               R32=R*SQRT(R)
               THET32=THET*1.5D0
               A1=COS(THET32)
               B1=SIN(THET32)
               DF1=4.D0/3.D0*R32
               EXPO=EXP(DF1*A1)
               EXPOI=1.D0/EXPO
               ETAR=EXPO*COS(DF1*B1)
               ETAI=EXPO*SIN(DF1*B1)
               ETAGR=EXPOI*COS(-DF1*B1)
               ETAGI=EXPOI*SIN(-DF1*B1)
               C1=-0.5D0
               S1=-0.5D0*SQRT3
               U=X*C1-Y*S1
               V=X*S1+Y*C1 
               CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)  
               IF (IFUN.EQ.1) THEN
                 BBR=SQRT3*AR+AI
                 BBI=-AR+SQRT3*AI
                 BR=BBR*ETAR-BBI*ETAI
                 BI=BBI*ETAR+BBR*ETAI
               ELSE
                 U=AR*C1-AI*S1
                 V=AR*S1+AI*C1
                 APR=U
                 API=V
                 BBPR=SQRT3*APR+API
                 BBPI=-APR+SQRT3*API
                 BPR=BBPR*ETAR-BBPI*ETAI
                 BPI=BBPI*ETAR+BBPR*ETAI
               ENDIF
             ELSE
               IF (IFUN.EQ.1) THEN
                 BR=0.D0
                 BI=0.D0
               ELSE
                 BPR=0.D0
                 BPI=0.D0
               ENDIF
             ENDIF  
             CALL AIZ(IFUN,IFAC,X,Y,AR,AI,IERR)
             IF (IFUN.EQ.1) THEN
               BR=BR-AI
               BI=BI+AR
             ELSE
               BPR=BPR-AI
               BPI=BPI+AR
             ENDIF
           ENDIF
           IF (THET.GT.PI23) THEN
             C1=-0.5D0
             S1=0.5D0*SQRT3
             U=X*C1-Y*S1
             V=X*S1+Y*C1 
             CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
             IF (IFUN.EQ.1) THEN
               BR=C*AR-S*AI                
               BI=C*AI+S*AR
             ELSE
               U=AR*C1-AI*S1
               V=AR*S1+AI*C1
               APR=U
               API=V
               BPR=C*APR-S*API
               BPI=C*API+S*APR
             ENDIF
             IF (IEXPF.EQ.0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   WE NEED ALSO AI_(1)(Z) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                R=SQRT(X*X+Y*Y)
                R32=R*SQRT(R)
                THET32=THET*1.5D0
                A1=COS(THET32)
                B1=SIN(THET32)
                DF1=4.D0/3.D0*R32
                EXPO=EXP(DF1*A1)
                EXPOI=1.D0/EXPO
                ETAR=EXPO*COS(DF1*B1)
                ETAI=EXPO*SIN(DF1*B1)
                ETAGR=EXPOI*COS(-DF1*B1)
                ETAGI=EXPOI*SIN(-DF1*B1)
                C1=-0.5D0
                S1=-0.5D0*SQRT3
                U=X*C1-Y*S1
                V=X*S1+Y*C1  
                CALL AIZ(IFUN,IFAC,U,V,AR,AI,IERR)
                IF (IFUN.EQ.1) THEN
                  S=-S
                  BBR=C*AR-S*AI
                  BBI=C*AI+S*AR
                  BR=BR+ETAR*BBR-ETAI*BBI
                  BI=BI+BBI*ETAR+ETAI*BBR
                ELSE
                  U=AR*C1-AI*S1
                  V=AR*S1+AI*C1
                  APR=U
                  API=V
                  S=-S
                  BBPR=C*APR-S*API
                  BBPI=C*API+S*APR
                  BPR=BPR+ETAR*BBPR-ETAI*BBPI
                  BPI=BPI+BBPI*ETAR+ETAI*BBPR
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          IF (Y0.LT.0) THEN
            BI=-BI 
            BPI=-BPI
          ENDIF
          IF (IFUN.EQ.1) THEN
            GBIR=BR
            GBII=BI
          ELSE
            GBIR=BPR
            GBII=BPI
          ENDIF
        ENDIF
        RETURN
        END

 





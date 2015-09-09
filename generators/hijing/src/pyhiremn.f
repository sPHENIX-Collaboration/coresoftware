    
C*********************************************************************  
    
      SUBROUTINE PYHIREMN(IPU1,IPU2)  
    
C...Adds on target remnants (one or two from each side) and 
C...includes primordial kT. 
      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/
      COMMON/HISTRNG/NFP(300,15),PPHI(300,15),NFT(300,15),PTHI(300,15)
      SAVE  /HISTRNG/
C...COMMON BLOCK FROM HIJING
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYHIPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYHIPARS/ 
      COMMON/PYHIINT1/MINT(400),VINT(400) 
      SAVE /PYHIINT1/ 
      DIMENSION KFLCH(2),KFLSP(2),CHI(2),PMS(6),IS(2),ROBO(5)   
    
C...Special case for lepton-lepton interaction. 
      IF(MINT(43).EQ.1) THEN    
        DO 100 JT=1,2   
        I=MINT(83)+JT+2 
        K(I,1)=21   
        K(I,2)=K(I-2,2) 
        K(I,3)=I-2  
        DO 100 J=1,5    
  100   P(I,J)=P(I-2,J) 
      ENDIF 
    
C...Find event type, set pointers.  
      IF(IPU1.EQ.0.AND.IPU2.EQ.0) RETURN    
      ISUB=MINT(1)  
      ILEP=0    
      IF(IPU1.EQ.0) ILEP=1  
      IF(IPU2.EQ.0) ILEP=2  
      IF(ISUB.EQ.95) ILEP=-1    
      IF(ILEP.EQ.1) IQ=MINT(84)+1   
      IF(ILEP.EQ.2) IQ=MINT(84)+2   
      IP=MAX(IPU1,IPU2) 
      ILEPR=MINT(83)+5-ILEP 
      NS=N  
    
C...Define initial partons, including primordial kT.    
  110 DO 130 JT=1,2 
      I=MINT(83)+JT+2   
      IF(JT.EQ.1) IPU=IPU1  
      IF(JT.EQ.2) IPU=IPU2  
      K(I,1)=21 
      K(I,3)=I-2    
      IF(ISUB.EQ.95) THEN   
        K(I,2)=21   
        SHS=0.  
      ELSEIF(MINT(40+JT).EQ.1.AND.IPU.NE.0) THEN    
        K(I,2)=K(IPU,2) 
        P(I,5)=P(IPU,5) 
        P(I,1)=0.   
        P(I,2)=0.   
        PMS(JT)=P(I,5)**2   
      ELSEIF(IPU.NE.0) THEN 
        K(I,2)=K(IPU,2) 
        P(I,5)=P(IPU,5) 
C...No primordial kT or chosen according to truncated Gaussian or   
C...exponential.
C
c     X.N. Wang (7.22.97)
c
        RPT1=0.0
        RPT2=0.0
        SS_W2=(PPHI(IHNT2(11),4)+PTHI(IHNT2(12),4))**2
     &       -(PPHI(IHNT2(11),1)+PTHI(IHNT2(12),1))**2
     &       -(PPHI(IHNT2(11),2)+PTHI(IHNT2(12),2))**2
     &       -(PPHI(IHNT2(11),3)+PTHI(IHNT2(12),3))**2
C
C********this is s of the current NN collision
        IF(SS_W2.LE.4.0*PARP(93)**2) GOTO 1211
c
        IF(IHPR2(5).LE.0) THEN
120	     IF(MSTP(91).LE.0) THEN
               PT=0. 
             ELSEIF(MSTP(91).EQ.1) THEN
               PT=PARP(91)*SQRT(-LOG(RLU(0)))
             ELSE    
               RPT1=RLU(0)   
               RPT2=RLU(0)   
               PT=-PARP(92)*LOG(RPT1*RPT2)   
             ENDIF   
             IF(PT.GT.PARP(93)) GOTO 120 
	     PHI=PARU(2)*RLU(0)  
	     RPT1=PT*COS(PHI)  
	     RPT2=PT*SIN(PHI)
	ELSE IF(IHPR2(5).EQ.1) THEN
	     IF(JT.EQ.1) JPT=NFP(IHNT2(11),11)
	     IF(JT.EQ.2) JPT=NFT(IHNT2(12),11)
1205	     PTGS=PARP(91)*SQRT(-LOG(RLU(0)))
	     IF(PTGS.GT.PARP(93)) GO TO 1205
	     PHI=2.0*HIPR1(40)*RLU(0)
	     RPT1=PTGS*COS(PHI)
	     RPT2=PTGS*SIN(PHI)
	     DO 1210 I_INT=1,JPT-1
		PKCSQ=PARP(91)*SQRT(-LOG(RLU(0)))
		PHI=2.0*HIPR1(40)*RLU(0)
		RPT1=RPT1+PKCSQ*COS(PHI)
		RPT2=RPT2+PKCSQ*SIN(PHI)
1210	     CONTINUE
             IF(RPT1**2+RPT2**2.GE.SS_W2/4.0) GO TO 1205
	ENDIF
C     X.N. Wang
C			********When initial interaction among soft partons is
C				assumed the primordial pt comes from the sum of
C				pt of JPT-1 number of initial interaction, JPT
C				is the number of interaction including present
C				one that nucleon hassuffered 
1211    P(I,1)=RPT1
	P(I,2)=RPT2  
        PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2   
      ELSE  
        K(I,2)=K(IQ,2)  
        Q2=VINT(52) 
        P(I,5)=-SQRT(Q2)    
        PMS(JT)=-Q2 
        SHS=(1.-VINT(43-JT))*Q2/VINT(43-JT)+VINT(5-JT)**2   
      ENDIF 
  130 CONTINUE  
    
C...Kinematics construction for initial partons.    
      I1=MINT(83)+3 
      I2=MINT(83)+4 
      IF(ILEP.EQ.0) SHS=VINT(141)*VINT(142)*VINT(2)+    
     &(P(I1,1)+P(I2,1))**2+(P(I1,2)+P(I2,2))**2 
      SHR=SQRT(MAX(0.,SHS)) 
      IF(ILEP.EQ.0) THEN    
        IF((SHS-PMS(1)-PMS(2))**2-4.*PMS(1)*PMS(2).LE.0.) GOTO 110  
        P(I1,4)=0.5*(SHR+(PMS(1)-PMS(2))/SHR)   
        P(I1,3)=SQRT(MAX(0.,P(I1,4)**2-PMS(1))) 
        P(I2,4)=SHR-P(I1,4) 
        P(I2,3)=-P(I1,3)    
      ELSEIF(ILEP.EQ.1) THEN    
        P(I1,4)=P(IQ,4) 
        P(I1,3)=P(IQ,3) 
        P(I2,4)=P(IP,4) 
        P(I2,3)=P(IP,3) 
      ELSEIF(ILEP.EQ.2) THEN    
        P(I1,4)=P(IP,4) 
        P(I1,3)=P(IP,3) 
        P(I2,4)=P(IQ,4) 
        P(I2,3)=P(IQ,3) 
      ENDIF 
      IF(MINT(43).EQ.1) RETURN  
    
C...Transform partons to overall CM-frame (not for leptoproduction).    
      IF(ILEP.EQ.0) THEN    
        ROBO(3)=(P(I1,1)+P(I2,1))/SHR   
        ROBO(4)=(P(I1,2)+P(I2,2))/SHR   
        CALL LUDBRB(I1,I2,0.,0.,-DBLE(ROBO(3)),-DBLE(ROBO(4)),0D0)  
        ROBO(2)=ULANGL(P(I1,1),P(I1,2)) 
        CALL LUDBRB(I1,I2,0.,-ROBO(2),0D0,0D0,0D0)  
        ROBO(1)=ULANGL(P(I1,3),P(I1,1)) 
        CALL LUDBRB(I1,I2,-ROBO(1),0.,0D0,0D0,0D0)  
        NMAX=MAX(MINT(52),IPU1,IPU2)    
        CALL LUDBRB(I1,NMAX,ROBO(1),ROBO(2),DBLE(ROBO(3)),DBLE(ROBO(4)),    
     &  0D0)    
        ROBO(5)=MAX(-0.999999,MIN(0.999999,(VINT(141)-VINT(142))/   
     &  (VINT(141)+VINT(142)))) 
        CALL LUDBRB(I1,NMAX,0.,0.,0D0,0D0,DBLE(ROBO(5)))    
      ENDIF 
    
C...Check invariant mass of remnant system: 
C...hadronic events or leptoproduction. 
      IF(ILEP.LE.0) THEN    
        IF(MSTP(81).LE.0.OR.MSTP(82).LE.0.OR.ISUB.EQ.95) THEN   
          VINT(151)=0.  
          VINT(152)=0.  
        ENDIF   
        PEH=P(I1,4)+P(I2,4)+0.5*VINT(1)*(VINT(151)+VINT(152))   
        PZH=P(I1,3)+P(I2,3)+0.5*VINT(1)*(VINT(151)-VINT(152))   
        SHH=(VINT(1)-PEH)**2-(P(I1,1)+P(I2,1))**2-(P(I1,2)+P(I2,2))**2- 
     &  PZH**2  
        PMMIN=P(MINT(83)+1,5)+P(MINT(83)+2,5)+ULMASS(K(I1,2))+  
     &  ULMASS(K(I2,2)) 
        IF(SHR.GE.VINT(1).OR.SHH.LE.(PMMIN+PARP(111))**2) THEN  
          MINT(51)=1    
          RETURN    
        ENDIF   
        SHR=SQRT(SHH+(P(I1,1)+P(I2,1))**2+(P(I1,2)+P(I2,2))**2) 
      ELSE  
        PEI=P(IQ,4)+P(IP,4) 
        PZI=P(IQ,3)+P(IP,3) 
        PMS(ILEP)=MAX(0.,PEI**2-PZI**2) 
        PMMIN=P(ILEPR-2,5)+ULMASS(K(ILEPR,2))+SQRT(PMS(ILEP))   
        IF(SHR.LE.PMMIN+PARP(111)) THEN 
          MINT(51)=1    
          RETURN    
        ENDIF   
      ENDIF 
    
C...Subdivide remnant if necessary, store first parton. 
  140 I=NS  
      DO 190 JT=1,2 
      IF(JT.EQ.ILEP) GOTO 190   
      IF(JT.EQ.1) IPU=IPU1  
      IF(JT.EQ.2) IPU=IPU2  
      CALL PYHISPLI(MINT(10+JT),MINT(12+JT),KFLCH(JT),KFLSP(JT))  
      I=I+1 
      IS(JT)=I  
      DO 150 J=1,5  
      K(I,J)=0  
      P(I,J)=0. 
  150 V(I,J)=0. 
      K(I,1)=3  
      K(I,2)=KFLSP(JT)  
      K(I,3)=MINT(83)+JT    
      P(I,5)=ULMASS(K(I,2)) 
    
C...First parton colour connections and transverse mass.    
      KFLS=(3-KCHG(LUCOMP(KFLSP(JT)),2)*ISIGN(1,KFLSP(JT)))/2   
      K(I,KFLS+3)=IPU   
      K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I    
      IF(KFLCH(JT).EQ.0) THEN   
        P(I,1)=-P(MINT(83)+JT+2,1)  
        P(I,2)=-P(MINT(83)+JT+2,2)  
        PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2   
    
C...When extra remnant parton or hadron: find relative pT, store.   
      ELSE  
        CALL LUPTDI(1,P(I,1),P(I,2))    
        PMS(JT+2)=P(I,5)**2+P(I,1)**2+P(I,2)**2 
        I=I+1   
        DO 160 J=1,5    
        K(I,J)=0    
        P(I,J)=0.   
  160   V(I,J)=0.   
        K(I,1)=1    
        K(I,2)=KFLCH(JT)    
        K(I,3)=MINT(83)+JT  
        P(I,5)=ULMASS(K(I,2))   
        P(I,1)=-P(MINT(83)+JT+2,1)-P(I-1,1) 
        P(I,2)=-P(MINT(83)+JT+2,2)-P(I-1,2) 
        PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2 
C...Relative distribution of energy for particle into two jets. 
        IMB=1   
        IF(MOD(MINT(10+JT)/1000,10).NE.0) IMB=2 
        IF(IABS(KFLCH(JT)).LE.10.OR.KFLCH(JT).EQ.21) THEN   
          CHIK=PARP(92+2*IMB)   
          IF(MSTP(92).LE.1) THEN    
            IF(IMB.EQ.1) CHI(JT)=RLU(0) 
            IF(IMB.EQ.2) CHI(JT)=1.-SQRT(RLU(0))    
          ELSEIF(MSTP(92).EQ.2) THEN    
            CHI(JT)=1.-RLU(0)**(1./(1.+CHIK))   
          ELSEIF(MSTP(92).EQ.3) THEN    
            CUT=2.*0.3/VINT(1)  
  170       CHI(JT)=RLU(0)**2   
            IF((CHI(JT)**2/(CHI(JT)**2+CUT**2))**0.25*(1.-CHI(JT))**CHIK    
     &      .LT.RLU(0)) GOTO 170    
          ELSE  
            CUT=2.*0.3/VINT(1)  
            CUTR=(1.+SQRT(1.+CUT**2))/CUT   
  180       CHIR=CUT*CUTR**RLU(0)   
            CHI(JT)=(CHIR**2-CUT**2)/(2.*CHIR)  
            IF((1.-CHI(JT))**CHIK.LT.RLU(0)) GOTO 180   
          ENDIF 
C...Relative distribution of energy for particle into jet plus particle.    
        ELSE    
          IF(MSTP(92).LE.1) THEN    
            IF(IMB.EQ.1) CHI(JT)=RLU(0) 
            IF(IMB.EQ.2) CHI(JT)=1.-SQRT(RLU(0))    
          ELSE  
            CHI(JT)=1.-RLU(0)**(1./(1.+PARP(93+2*IMB))) 
          ENDIF 
          IF(MOD(KFLCH(JT)/1000,10).NE.0) CHI(JT)=1.-CHI(JT)    
        ENDIF   
        PMS(JT)=PMS(JT+4)/CHI(JT)+PMS(JT+2)/(1.-CHI(JT))    
        KFLS=KCHG(LUCOMP(KFLCH(JT)),2)*ISIGN(1,KFLCH(JT))   
        IF(KFLS.NE.0) THEN  
          K(I,1)=3  
          KFLS=(3-KFLS)/2   
          K(I,KFLS+3)=IPU   
          K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I    
        ENDIF   
      ENDIF 
  190 CONTINUE  
      IF(SHR.LE.SQRT(PMS(1))+SQRT(PMS(2))) GOTO 140 
      N=I   
    
C...Reconstruct kinematics of remnants. 
      DO 200 JT=1,2 
      IF(JT.EQ.ILEP) GOTO 200   
      PE=0.5*(SHR+(PMS(JT)-PMS(3-JT))/SHR)  
      PZ=SQRT(PE**2-PMS(JT))    
      IF(KFLCH(JT).EQ.0) THEN   
        P(IS(JT),4)=PE  
        P(IS(JT),3)=PZ*(-1)**(JT-1) 
      ELSE  
        PW1=CHI(JT)*(PE+PZ) 
        P(IS(JT)+1,4)=0.5*(PW1+PMS(JT+4)/PW1)   
        P(IS(JT)+1,3)=0.5*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)  
        P(IS(JT),4)=PE-P(IS(JT)+1,4)    
        P(IS(JT),3)=PZ*(-1)**(JT-1)-P(IS(JT)+1,3)   
      ENDIF 
  200 CONTINUE  
    
C...Hadronic events: boost remnants to correct longitudinal frame.  
      IF(ILEP.LE.0) THEN    
        CALL LUDBRB(NS+1,N,0.,0.,0D0,0D0,-DBLE(PZH/(VINT(1)-PEH)))  
C...Leptoproduction events: boost colliding subsystem.  
      ELSE  
        NMAX=MAX(IP,MINT(52))   
        PEF=SHR-PE  
        PZF=PZ*(-1)**(ILEP-1)   
        PT2=P(ILEPR,1)**2+P(ILEPR,2)**2 
        PHIPT=ULANGL(P(ILEPR,1),P(ILEPR,2)) 
        CALL LUDBRB(MINT(84)+1,NMAX,0.,-PHIPT,0D0,0D0,0D0)  
        RQP=P(IQ,3)*(PT2+PEI**2)-P(IQ,4)*PEI*PZI    
        SINTH=P(IQ,4)*SQRT(PT2*(PT2+PEI**2)/(RQP**2+PT2*    
     &  P(IQ,4)**2*PZI**2))*SIGN(1.,-RQP)   
        CALL LUDBRB(MINT(84)+1,NMAX,ASIN(SINTH),0.,0D0,0D0,0D0) 
        BETAX=(-PEI*PZI*SINTH+SQRT(PT2*(PT2+PEI**2-(PZI*SINTH)**2)))/   
     &  (PT2+PEI**2)    
        CALL LUDBRB(MINT(84)+1,NMAX,0.,0.,DBLE(BETAX),0D0,0D0)  
        CALL LUDBRB(MINT(84)+1,NMAX,0.,PHIPT,0D0,0D0,0D0)   
        PEM=P(IQ,4)+P(IP,4) 
        PZM=P(IQ,3)+P(IP,3) 
        BETAZ=(-PEM*PZM+PZF*SQRT(PZF**2+PEM**2-PZM**2))/(PZF**2+PEM**2) 
        CALL LUDBRB(MINT(84)+1,NMAX,0.,0.,0D0,0D0,DBLE(BETAZ))  
        CALL LUDBRB(I1,I2,ASIN(SINTH),0.,DBLE(BETAX),0D0,0D0)   
        CALL LUDBRB(I1,I2,0.,PHIPT,0D0,0D0,DBLE(BETAZ)) 
      ENDIF 
    
      RETURN    
      END   

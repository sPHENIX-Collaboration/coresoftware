      SUBROUTINE PYHIINIT(FRAME,BEAM,TARGET,WIN)  
    
C...Initializes the generation procedure; finds maxima of the   
C...differential cross-sections to be used for weighting.   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/LUDAT4/CHAF(500)   
      CHARACTER CHAF*8  
      SAVE /LUDAT4/ 
      COMMON/PYHISUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYHISUBS/ 
      COMMON/PYHIPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYHIPARS/ 
      COMMON/PYHIINT1/MINT(400),VINT(400) 
      SAVE /PYHIINT1/ 
      COMMON/PYHIINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYHIINT2/ 
      COMMON/PYHIINT5/NGEN(0:200,3),XSEC(0:200,3) 
      SAVE /PYHIINT5/ 
      CHARACTER*(*) FRAME,BEAM,TARGET   
      CHARACTER CHFRAM*8,CHBEAM*8,CHTARG*8,CHMO(12)*3,CHLH(2)*6 
      DATA CHMO/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',  
     &'Oct','Nov','Dec'/, CHLH/'lepton','hadron'/   
    
C...Write headers.  
C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1000) MSTP(181),MSTP(182),  
C     &MSTP(185),CHMO(MSTP(184)),MSTP(183)   
      CALL LULIST(0)
C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1100)  
    
C...Identify beam and target particles and initialize kinematics.   
      CHFRAM=FRAME//' ' 
      CHBEAM=BEAM//' '  
      CHTARG=TARGET//' '    
      CALL PYHIINKI(CHFRAM,CHBEAM,CHTARG,WIN) 
    
C...Select partonic subprocesses to be included in the simulation.  
      IF(MSEL.NE.0) THEN    
        DO 100 I=1,200  
  100   MSUB(I)=0   
      ENDIF 
      IF(MINT(43).EQ.1.AND.(MSEL.EQ.1.OR.MSEL.EQ.2)) THEN   
C...Lepton+lepton -> gamma/Z0 or W. 
        IF(MINT(11)+MINT(12).EQ.0) MSUB(1)=1    
        IF(MINT(11)+MINT(12).NE.0) MSUB(2)=1    
      ELSEIF(MSEL.EQ.1) THEN    
C...High-pT QCD processes:  
        MSUB(11)=1  
        MSUB(12)=1  
        MSUB(13)=1  
        MSUB(28)=1  
        MSUB(53)=1  
        MSUB(68)=1  
        IF(MSTP(82).LE.1.AND.CKIN(3).LT.PARP(81)) MSUB(95)=1    
        IF(MSTP(82).GE.2.AND.CKIN(3).LT.PARP(82)) MSUB(95)=1    
      ELSEIF(MSEL.EQ.2) THEN    
C...All QCD processes:  
        MSUB(11)=1  
        MSUB(12)=1  
        MSUB(13)=1  
        MSUB(28)=1  
        MSUB(53)=1  
        MSUB(68)=1  
        MSUB(91)=1  
        MSUB(92)=1  
        MSUB(93)=1  
        MSUB(95)=1  
      ELSEIF(MSEL.GE.4.AND.MSEL.LE.8) THEN  
C...Heavy quark production. 
        MSUB(81)=1  
        MSUB(82)=1  
        DO 110 J=1,MIN(8,MDCY(21,3))    
  110   MDME(MDCY(21,2)+J-1,1)=0    
        MDME(MDCY(21,2)+MSEL-1,1)=1 
      ELSEIF(MSEL.EQ.10) THEN   
C...Prompt photon production:   
        MSUB(14)=1  
        MSUB(18)=1  
        MSUB(29)=1  
      ELSEIF(MSEL.EQ.11) THEN   
C...Z0/gamma* production:   
        MSUB(1)=1   
      ELSEIF(MSEL.EQ.12) THEN   
C...W+/- production:    
        MSUB(2)=1   
      ELSEIF(MSEL.EQ.13) THEN   
C...Z0 + jet:   
        MSUB(15)=1  
        MSUB(30)=1  
      ELSEIF(MSEL.EQ.14) THEN   
C...W+/- + jet: 
        MSUB(16)=1  
        MSUB(31)=1  
      ELSEIF(MSEL.EQ.15) THEN   
C...Z0 & W+/- pair production:  
        MSUB(19)=1  
        MSUB(20)=1  
        MSUB(22)=1  
        MSUB(23)=1  
        MSUB(25)=1  
      ELSEIF(MSEL.EQ.16) THEN   
C...H0 production:  
        MSUB(3)=1   
        MSUB(5)=1   
        MSUB(8)=1   
        MSUB(102)=1 
      ELSEIF(MSEL.EQ.17) THEN   
C...H0 & Z0 or W+/- pair production:    
        MSUB(24)=1  
        MSUB(26)=1  
      ELSEIF(MSEL.EQ.21) THEN   
C...Z'0 production: 
        MSUB(141)=1 
      ELSEIF(MSEL.EQ.22) THEN   
C...H+/- production:    
        MSUB(142)=1 
      ELSEIF(MSEL.EQ.23) THEN   
C...R production:   
        MSUB(143)=1 
      ENDIF 
    
C...Count number of subprocesses on.    
      MINT(44)=0    
      DO 120 ISUB=1,200 
      IF(MINT(43).LT.4.AND.ISUB.GE.91.AND.ISUB.LE.96.AND.   
     &MSUB(ISUB).EQ.1) THEN 
        WRITE(MSTU(11),1200) ISUB,CHLH(MINT(41)),CHLH(MINT(42)) 
        STOP    
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).EQ.-1) THEN 
        WRITE(MSTU(11),1300) ISUB   
        STOP    
      ELSEIF(MSUB(ISUB).EQ.1.AND.ISET(ISUB).LE.-2) THEN 
        WRITE(MSTU(11),1400) ISUB   
        STOP    
      ELSEIF(MSUB(ISUB).EQ.1) THEN  
        MINT(44)=MINT(44)+1 
      ENDIF 
  120 CONTINUE  
      IF(MINT(44).EQ.0) THEN    
        WRITE(MSTU(11),1500)    
        STOP    
      ENDIF 
      MINT(45)=MINT(44)-MSUB(91)-MSUB(92)-MSUB(93)-MSUB(94) 
    
C...Maximum 4 generations; set maximum number of allowed flavours.  
      MSTP(1)=MIN(4,MSTP(1))    
      MSTU(114)=MIN(MSTU(114),2*MSTP(1))    
      MSTP(54)=MIN(MSTP(54),2*MSTP(1))  
    
C...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton. 
      DO 140 I=-20,20   
      VINT(180+I)=0.    
      IA=IABS(I)    
      IF(IA.GE.1.AND.IA.LE.2*MSTP(1)) THEN  
        DO 130 J=1,MSTP(1)  
        IB=2*J-1+MOD(IA,2)  
        IPM=(5-ISIGN(1,I))/2    
        IDC=J+MDCY(IA,2)+2  
  130   IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.IPM) VINT(180+I)= 
     &  VINT(180+I)+VCKM((IA+1)/2,(IB+1)/2) 
      ELSEIF(IA.GE.11.AND.IA.LE.10+2*MSTP(1)) THEN  
        VINT(180+I)=1.  
      ENDIF 
  140 CONTINUE  
    
C...Choose Lambda value to use in alpha-strong. 
      MSTU(111)=MSTP(2) 
      IF(MSTP(3).GE.1) THEN 
        ALAM=PARP(1)    
        IF(MSTP(51).EQ.1) ALAM=0.2  
        IF(MSTP(51).EQ.2) ALAM=0.29 
        IF(MSTP(51).EQ.3) ALAM=0.2  
        IF(MSTP(51).EQ.4) ALAM=0.4  
        IF(MSTP(51).EQ.11) ALAM=0.16    
        IF(MSTP(51).EQ.12) ALAM=0.26    
        IF(MSTP(51).EQ.13) ALAM=0.36    
        PARP(1)=ALAM    
        PARP(61)=ALAM   
        PARU(112)=ALAM  
        PARJ(81)=ALAM   
      ENDIF 
    
C...Initialize widths and partial widths for resonances.    
      CALL PYHIINRE   
    
C...Reset variables for cross-section calculation.  
      DO 150 I=0,200    
      DO 150 J=1,3  
      NGEN(I,J)=0   
  150 XSEC(I,J)=0.  
      VINT(108)=0.  
    
C...Find parametrized total cross-sections. 
      IF(MINT(43).EQ.4) CALL PYHIXTOT 
    
C...Maxima of differential cross-sections.  
      IF(MSTP(121).LE.0) CALL PYHIMAXI    
    
C...Initialize possibility of overlayed events. 
      IF(MSTP(131).NE.0) CALL PYHIOVLY(1) 
    
C...Initialize multiple interactions with variable impact parameter.    
      IF(MINT(43).EQ.4.AND.(MINT(45).NE.0.OR.MSTP(131).NE.0).AND.   
     &MSTP(82).GE.2) CALL PYHIMULT(1) 
C      IF(MSTP(122).GE.1) WRITE(MSTU(11),1600)  
    
C...Formats for initialization information. 
 1000 FORMAT(///20X,'The Lund Monte Carlo - PYHITHIA version ',I1,
     &'.',I1/ 
     &20X,'**  Last date of change:  ',I2,1X,A3,1X,I4,'  **'/)  
 1100 FORMAT('1',18('*'),1X,'PYHIINIT: initialization of PYHITHIA ',    
     &'(hijing pythia) routines',1X,17('*'))    
 1200 FORMAT(1X,'Error: process number ',I3,' not meaningful for ',A6,  
     &'-',A6,' interactions.'/1X,'Execution stopped!')  
 1300 FORMAT(1X,'Error: requested subprocess',I4,' not implemented.'/   
     &1X,'Execution stopped!')  
 1400 FORMAT(1X,'Error: requested subprocess',I4,' not existing.'/  
     &1X,'Execution stopped!')  
 1500 FORMAT(1X,'Error: no subprocess switched on.'/    
     &1X,'Execution stopped.')  
 1600 FORMAT(/1X,22('*'),1X,'PYHIINIT: initialization completed',1X,  
     &22('*'))  
    
      RETURN    
      END   

    
C*********************************************************************  
    
      SUBROUTINE PYHIINRE 
    
C...Calculates full and effective widths of guage bosons, stores masses 
C...and widths, rescales coefficients to be used for resonance  
C...production generation.  
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYHISUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYHISUBS/ 
      COMMON/PYHIPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYHIPARS/ 
      COMMON/PYHIINT1/MINT(400),VINT(400) 
      SAVE /PYHIINT1/ 
      COMMON/PYHIINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYHIINT2/ 
      COMMON/PYHIINT4/WIDP(21:40,0:40),WIDE(21:40,0:40),WIDS(21:40,3) 
      SAVE /PYHIINT4/ 
      COMMON/PYHIINT6/PROC(0:200) 
      CHARACTER PROC*28 
      SAVE /PYHIINT6/ 
      DIMENSION WDTP(0:40),WDTE(0:40,0:5)   
    
C...Calculate full and effective widths of gauge bosons.    
      AEM=PARU(101) 
      XW=PARU(102)  
      DO 100 I=21,40    
      DO 100 J=0,40 
      WIDP(I,J)=0.  
  100 WIDE(I,J)=0.  
    
C...W+/-:   
      WMAS=PMAS(24,1)   
      WFAC=AEM/(24.*XW)*WMAS    
      CALL PYHIWIDT(24,WMAS,WDTP,WDTE)    
      WIDS(24,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(24,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(24,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 110 I=0,40 
      WIDP(24,I)=WFAC*WDTP(I)   
  110 WIDE(24,I)=WFAC*WDTE(I,0) 
    
C...H+/-:   
      HCMAS=PMAS(37,1)  
      HCFAC=AEM/(8.*XW)*(HCMAS/WMAS)**2*HCMAS   
      CALL PYHIWIDT(37,HCMAS,WDTP,WDTE)   
      WIDS(37,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(37,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(37,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 120 I=0,40 
      WIDP(37,I)=HCFAC*WDTP(I)  
  120 WIDE(37,I)=HCFAC*WDTE(I,0)    
    
C...Z0: 
      ZMAS=PMAS(23,1)   
      ZFAC=AEM/(48.*XW*(1.-XW))*ZMAS    
      CALL PYHIWIDT(23,ZMAS,WDTP,WDTE)    
      WIDS(23,1)=((WDTE(0,1)+WDTE(0,2))**2+ 
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(23,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(23,3)=0. 
      DO 130 I=0,40 
      WIDP(23,I)=ZFAC*WDTP(I)   
  130 WIDE(23,I)=ZFAC*WDTE(I,0) 
    
C...H0: 
      HMAS=PMAS(25,1)   
      HFAC=AEM/(8.*XW)*(HMAS/WMAS)**2*HMAS  
      CALL PYHIWIDT(25,HMAS,WDTP,WDTE)    
      WIDS(25,1)=((WDTE(0,1)+WDTE(0,2))**2+ 
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(25,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(25,3)=0. 
      DO 140 I=0,40 
      WIDP(25,I)=HFAC*WDTP(I)   
  140 WIDE(25,I)=HFAC*WDTE(I,0) 
    
C...Z'0:    
      ZPMAS=PMAS(32,1)  
      ZPFAC=AEM/(48.*XW*(1.-XW))*ZPMAS  
      CALL PYHIWIDT(32,ZPMAS,WDTP,WDTE)   
      WIDS(32,1)=((WDTE(0,1)+WDTE(0,2)+WDTE(0,3))**2+   
     &2.*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+   
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(32,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(32,3)=0. 
      DO 150 I=0,40 
      WIDP(32,I)=ZPFAC*WDTP(I)  
  150 WIDE(32,I)=ZPFAC*WDTE(I,0)    
    
C...R:  
      RMAS=PMAS(40,1)   
      RFAC=0.08*RMAS/((MSTP(1)-1)*(1.+6.*(1.+ULALPS(RMAS**2)/PARU(1)))) 
      CALL PYHIWIDT(40,RMAS,WDTP,WDTE)    
      WIDS(40,1)=((WDTE(0,1)+WDTE(0,2))*(WDTE(0,1)+WDTE(0,3))+  
     &(WDTE(0,1)+WDTE(0,2)+WDTE(0,1)+WDTE(0,3))*(WDTE(0,4)+WDTE(0,5))+  
     &2.*WDTE(0,4)*WDTE(0,5))/WDTP(0)**2    
      WIDS(40,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)    
      WIDS(40,3)=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)    
      DO 160 I=0,40 
      WIDP(40,I)=WFAC*WDTP(I)   
  160 WIDE(40,I)=WFAC*WDTE(I,0) 
    
C...Q:  
      KFLQM=1   
      DO 170 I=1,MIN(8,MDCY(21,3))  
      IDC=I+MDCY(21,2)-1    
      IF(MDME(IDC,1).LE.0) GOTO 170 
      KFLQM=I   
  170 CONTINUE  
      MINT(46)=KFLQM    
      KFPR(81,1)=KFLQM  
      KFPR(81,2)=KFLQM  
      KFPR(82,1)=KFLQM  
      KFPR(82,2)=KFLQM  
    
C...Set resonance widths and branching ratios in JETSET.    
      DO 180 I=1,6  
      IF(I.LE.3) KC=I+22    
      IF(I.EQ.4) KC=32  
      IF(I.EQ.5) KC=37  
      IF(I.EQ.6) KC=40  
      PMAS(KC,2)=WIDP(KC,0) 
      PMAS(KC,3)=MIN(0.9*PMAS(KC,1),10.*PMAS(KC,2)) 
      DO 180 J=1,MDCY(KC,3) 
      IDC=J+MDCY(KC,2)-1    
      BRAT(IDC)=WIDE(KC,J)/WIDE(KC,0)   
  180 CONTINUE  
    
C...Special cases in treatment of gamma*/Z0: redefine process name. 
      IF(MSTP(43).EQ.1) THEN    
        PROC(1)='f + fb -> gamma*'  
      ELSEIF(MSTP(43).EQ.2) THEN    
        PROC(1)='f + fb -> Z0'  
      ELSEIF(MSTP(43).EQ.3) THEN    
        PROC(1)='f + fb -> gamma*/Z0'   
      ENDIF 
    
C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name. 
      IF(MSTP(44).EQ.1) THEN    
        PROC(141)='f + fb -> gamma*'    
      ELSEIF(MSTP(44).EQ.2) THEN    
        PROC(141)='f + fb -> Z0'    
      ELSEIF(MSTP(44).EQ.3) THEN    
        PROC(141)='f + fb -> Z''0'  
      ELSEIF(MSTP(44).EQ.4) THEN    
        PROC(141)='f + fb -> gamma*/Z0' 
      ELSEIF(MSTP(44).EQ.5) THEN    
        PROC(141)='f + fb -> gamma*/Z''0'   
      ELSEIF(MSTP(44).EQ.6) THEN    
        PROC(141)='f + fb -> Z0/Z''0'   
      ELSEIF(MSTP(44).EQ.7) THEN    
        PROC(141)='f + fb -> gamma*/Z0/Z''0'    
      ENDIF 
    
      RETURN    
      END   

C
C
C
C     Modified for HIJING program
c
c    modification July 22, 1997  In pyremnn put an upper limit
c     on the total pt kick the parton can accumulate via multiple
C     scattering. Set the upper limit to be the sqrt(s)/2,
c     this is fix cronin bug for Pb+Pb events at SPS energy.
c
C
C Last modification Oct. 1993 to comply with non-vax
C machines' compiler 
C
C
      SUBROUTINE LU1ENT(IP,KF,PE,THE,PHI)   
    
C...Purpose: to store one parton/particle in commonblock LUJETS.    
      COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5) 
      SAVE /LUJETS/ 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
    
C...Standard checks.    
      MSTU(28)=0    
      IF(MSTU(12).GE.1) CALL LULIST(0)  
      IPA=MAX(1,IABS(IP))   
      IF(IPA.GT.MSTU(4)) CALL LUERRM(21,    
     &'(LU1ENT:) writing outside LUJETS memory')    
      KC=LUCOMP(KF) 
      IF(KC.EQ.0) CALL LUERRM(12,'(LU1ENT:) unknown flavour code')  
    
C...Find mass. Reset K, P and V vectors.    
      PM=0. 
      IF(MSTU(10).EQ.1) PM=P(IPA,5) 
      IF(MSTU(10).GE.2) PM=ULMASS(KF)   
      DO 100 J=1,5  
      K(IPA,J)=0    
      P(IPA,J)=0.   
  100 V(IPA,J)=0.   
    
C...Store parton/particle in K and P vectors.   
      K(IPA,1)=1    
      IF(IP.LT.0) K(IPA,1)=2    
      K(IPA,2)=KF   
      P(IPA,5)=PM   
      P(IPA,4)=MAX(PE,PM)   
      PA=SQRT(P(IPA,4)**2-P(IPA,5)**2)  
      P(IPA,1)=PA*SIN(THE)*COS(PHI) 
      P(IPA,2)=PA*SIN(THE)*SIN(PHI) 
      P(IPA,3)=PA*COS(THE)  
    
C...Set N. Optionally fragment/decay.   
      N=IPA 
      IF(IP.EQ.0) CALL LUEXEC   
    
      RETURN    
      END   

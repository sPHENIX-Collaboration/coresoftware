    
C*********************************************************************  
    
      SUBROUTINE PYHIFRAM(IFRAME) 
    
C...Performs transformations between different coordinate frames.   
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/PYHIPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYHIPARS/ 
      COMMON/PYHIINT1/MINT(400),VINT(400) 
      SAVE /PYHIINT1/ 
    
      IF(IFRAME.LT.1.OR.IFRAME.GT.2) THEN   
        WRITE(MSTU(11),1000) IFRAME,MINT(6) 
        RETURN  
      ENDIF 
      IF(IFRAME.EQ.MINT(6)) RETURN  
    
      IF(MINT(6).EQ.1) THEN 
C...Transform from fixed target or user specified frame to  
C...CM-frame of incoming particles. 
        CALL LUROBO(0.,0.,-VINT(8),-VINT(9),-VINT(10))  
        CALL LUROBO(0.,-VINT(7),0.,0.,0.)   
        CALL LUROBO(-VINT(6),0.,0.,0.,0.)   
        MINT(6)=2   
    
      ELSE  
C...Transform from particle CM-frame to fixed target or user specified  
C...frame.  
        CALL LUROBO(VINT(6),VINT(7),VINT(8),VINT(9),VINT(10))   
        MINT(6)=1   
      ENDIF 
      MSTI(6)=MINT(6)   
    
 1000 FORMAT(1X,'Error: illegal values in subroutine PYHIFRAM.',1X,   
     &'No transformation performed.'/1X,'IFRAME =',1X,I5,'; MINT(6) =', 
     &1X,I5)    
    
      RETURN    
      END   

      FVOIDP FUNCTION HIJCRDN_ADDRESS() 

      COMMON/HIJCRDN/YP(3,300),YT(3,300)
      SAVE  /HIJCRDN/
      EXTERNAL HIDATA ! All that is required to force loading of BLOCK DATA
      FVOIDP GETADDR

      HIJCRDN_ADDRESS = GETADDR(YP)
      
      RETURN
      END

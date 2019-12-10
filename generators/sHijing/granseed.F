      FVOIDP FUNCTION RANSEED_ADDRESS() 

      COMMON/RANSEED/NSEED
      SAVE  /RANSEED/
      EXTERNAL HIDATA ! All that is required to force loading of BLOCK DATA
      FVOIDP GETADDR

      RANSEED_ADDRESS = GETADDR(NSEED)
      
      RETURN
      END

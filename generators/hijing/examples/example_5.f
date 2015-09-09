 
C****************************************************************************
C	The program for checking conservation laws fulfillment in 
C	                  the HIJING model.
C	          Written by V.Uzhinsky, CERN, Oct. 2003
C***************************************************************************

      CHARACTER FRAME*8,PROJ*8,TARG*8

      COMMON/HIMAIN1/ NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/

      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
      SAVE  /HIMAIN2/

C         ********information of produced particles
C
      CHARACTER *1 KEY
      CHARACTER *12 FNAME

      PARAMETER (NWPAWC=50000)
      COMMON/PAWC/HMEMORY(NWPAWC)

      COMMON/RANSEED/NSEED                                    
      SAVE  /RANSEED/                                         
      NSEED=0                                                 

      write(6,*)'====================================================='
      write(6,*)'The program for checking conservation laws fulfillment'
      write(6,*)'              in the HIJING model                    '
      write(6,*)'  Only hadronic (hh) interactions are considered     '
      write(6,*)'====================================================='
      write(6,*)'        You can work in Lab. or CM systems           '
      write(6,*)' Would you like to use CM system? Y - Yes, N - No?   '
    
      read(5,1001) KEY
1001  format(A1)      
      
      if(KEY.eq.'Y'.or.KEY.eq.'y') then
        write(6,*)'CM system is used --------------------------------'
        FRAME='CMS'
      else
        write(6,*)'LAB system is used -------------------------------'      
        FRAME='LAB'
      endif
      
      write(6,*)'Enter the corresponding energy (GeV)'                              

      read(5,*) EFRM 
      
      write(6,*)
      write(6,*)'Enter a type of the projectile particle'
      write(6,*)
      write(6,*)' P proton,            PBAR anti-proton,' 
      write(6,*)' N neutron,           NBAR anti-neutron,'
      write(6,*)' PI+ - positive pion, PI- negative pion,'
      write(6,*)' K+ positive kaon,    K- negative kaon'      

      read(5,1002) PROJ
1002  format(A8)

      if(PROJ.ne.'A') then
        IAP=1
        if(PROJ.eq.'P'   ) IZP= 1
        if(PROJ.eq.'PBAR') IZP=-1
        if(PROJ.eq.'N'   ) IZP= 0
        if(PROJ.eq.'NBAR') IZP= 0
        if(PROJ.eq.'PI+' ) IZP= 1
        if(PROJ.eq.'PI-' ) IZP=-1
        if(PROJ.eq.'K+'  ) IZP= 1
        if(PROJ.eq.'K-'  ) IZP=-1
      else
        write(6,*)
        write(6,*)'Enter mass number and charge of the proj. nucleus'
        read(5,*)  IAP, IZP
      endif
        
      write(6,*)
      write(6,*)'Enter a type of the target particle (same notations)'      
      read(5,1002) TARG                                          

      if(TARG.ne.'A') then
        IAT=1
        if(TARG.eq.'P'   ) IZT= 1
        if(TARG.eq.'PBAR') IZT=-1
        if(TARG.eq.'N'   ) IZT= 0
        if(TARG.eq.'NBAR') IZT= 0
        if(TARG.eq.'PI+' ) IZT= 1
        if(TARG.eq.'PI-' ) IZT=-1
        if(TARG.eq.'K+'  ) IZT= 1
        if(TARG.eq.'K-'  ) IZT=-1
      else
        write(6,*)
        write(6,*)'Enter mass number and charge of the target nucleus'
        read(5,*)  IAT, IZT
      endif
      
      write(6,*)
      write(6,*)'Enter number of events'                                            
      read(5,*) N_events                                             

      write(6,*)'Enter FILENAME for HBOOK output'
      read(5,1003) FNAME
1003  format(A12)

c------------------------------------------------------------------
      CALL HLIMIT(NWPAWC)

      EFRMl=EFRM*(9./10.)
      EFRMh=EFRM*(11./10.)     
      CALL HBOOK1(10,' Energy conservation',100,EFRMl,EFRMh,0.)

      Pxl=-5.
      Pxh= 5.
      CALL HBOOK1(20,' Px momentum conservation',100,Pxl,Pxh,0.)

      Pyl=-5.         
      Pyh= 5.         
      CALL HBOOK1(30,' Py momentum conservation',100,Pyl,Pyh,0.)

      if(FRAME.eq.'CMS') then
        Pzl=-5.
        Pzh= 5.
      else
        if(PROJ.eq.'P'.or.PROJ.eq.'PBAR'.or.
     &     PROJ.eq.'N'.or.PROJ.eq.'NBAR'.or.PROJ.eq.'A') then
         Plab=sqrt(EFRM**2-0.88)
        endif

        if(PROJ.eq.'PI+'.or.PROJ.eq.'PI-') then
         Plab=sqrt(EFRM**2-0.0196)
        endif

        if(PROJ.eq.'K+'.or.PROJ.eq.'K-') then
         Plab=sqrt(EFRM-0.25)
        endif

        Pzl=Plab*(3./4.)
        Pzh=Plab*(5./4.)
      endif
 
      CALL HBOOK1(40,' Pz momentum conservation',100,Pzl,Pzh,0.)

      CALL HBOOK1(50,' Charge conservation'       ,100,-5.,5.,0.)

      CALL HBOOK1(60,' Baryon number conservation',100,-5.,5.,0.)

      CALL HBOOK1(70,' Lepton number conservation',100,-5.,5.,0.)

      CALL HBOOK1(80,' Azimuthal isotropy'      ,100,-3.2,3.2,0.)
C----------------------------------------------------------------

      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
C			********Initialize HIJING

      WRITE(6,*)' Simulation of interactions with'             
      WRITE(6,*)                                               
      WRITE(6,*)' Proj = ',PROJ,' and Targ = ',TARG            
      WRITE(6,*)' IAP  =',IAP  ,'            IAT  =',IAT       
      WRITE(6,*)' IZP  =',IZP  ,'            IZT  =',IZT       
      WRITE(6,*)                                               
      WRITE(6,*)' Reference frame -   ',FRAME                  
      WRITE(6,*)' ENERGY            ',EFRM,' GeV'                     
      WRITE(6,*)' Number of generated events -',N_events      
      WRITE(6,*)                                              

        

      BMIN=0.0
      BMAX=0.0

      DO 2000 I_event=1,N_events
C
        WRITE(6,*)' Event # ',I_event,' ------------------'
C
        CALL HIJING(FRAME,BMIN,BMAX)
C
        Sum_E =0.
        Sum_Px=0.
        Sum_Py=0. 
        Sum_Pz=0.
        Sum_Q =0.
        Sum_B =0.
        Sum_L =0. 
      
        DO 3000 I=1,NATT
                                                               
          if(KATT(I,2).eq. 0) go to 3000 ! reject non-interacting projectile
          if(KATT(I,2).eq. 1) go to 3000 ! reject elastic scattering
          if(KATT(I,2).eq.10) go to 3000 ! reject non-interacting target
          if(KATT(I,2).eq.11) go to 3000 ! reject elastic scattering

          ICH=LUCHGE(KATT(I,1))/3                              

          Amass=ULMASS(KATT(I,1))                             

          Sum_E =Sum_E + PATT(i,4)
          Sum_Px=Sum_Px+ PATT(i,1)
          Sum_Py=Sum_Py+ PATT(i,2)
          Sum_Pz=Sum_Pz+ PATT(i,3)
          Sum_Q =Sum_Q + ICH

          if(IABS(KATT(i,1)).gt.2000) then
            Sum_B=Sum_B+IABS(KATT(i,1))/KATT(i,1)
          endif

          if((IABS(KATT(i,1)).gt.10.).and.
     &       (IABS(KATT(i,1)).lt.20.))           then
            Sum_L=Sum_L+IABS(KATT(i,1))/KATT(i,1)
          endif

          Phi=ULANGL(PATT(i,1),PATT(i,2))
          CALL HF1(80,Phi,1.)

3000    CONTINUE

        CALL HF1(10,Sum_E ,1.) 
        CALL HF1(20,Sum_Px,1.)
        CALL HF1(30,Sum_Py,1.)
        CALL HF1(40,Sum_Pz,1.)
        CALL HF1(50,Sum_Q ,1.)
        CALL HF1(60,Sum_B ,1.)
        CALL HF1(70,Sum_L ,1.)
2000  CONTINUE

      CALL HRPUT(0,FNAME,'N')
      END

      FUNCTION RAN(NSEED)                                      
      RAN=RLU(NSEED)                                           
      RETURN                                                   
      END                                                     

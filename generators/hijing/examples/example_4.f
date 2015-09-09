C****************************************************************************
C          Program # 3 from Comp. Phys. Commun. 83 (1994) 307
C	           by M. Gyulassy and X-.N. Wang
C             Modified by V.Uzhinsky, CERN, Oct. 2003	
C***************************************************************************

      CHARACTER FRAME*8,PROJ*8,TARG*8

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/

      COMMON/RANSEED/NSEED                                     
      SAVE  /RANSEED/                                          

      CHARACTER *1 KEY

      NSEED=0                                                 

C....switch on triggered jet production:
      IHPR2(3)=1

      FRAME='CMS'
      
      write(6,*)'===================================================='
      write(6,*)' Simulation of events with triggered hard processes '
      write(6,*)'          in hh-, hA- and AA-collisions             '
      write(6,*)'    Calculation will be performed in CM System      '
      write(6,*)'===================================================='

      write(6,*)
      write(6,*)'Enter the energy per NN-collision (GeV)'
      read(5,*)  EFRM

      write(6,*)   
      write(6,*)'Enter a type of the "projectile" particle'
      write(6,*)
      write(6,*)' P proton,            PBAR anti-proton,'
      write(6,*)' N neutron,           NBAR anti-neutron,'
      write(6,*)' PI+ - positive pion, PI- negative pion,'
      write(6,*)' K+ positive kaon,    K- negative kaon'
      write(6,*)
      write(6,*)' A - nucleus --------------------------'

      read(5,1) PROJ
1     format(A8)
    
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
      write(6,*)'Enter a type of the "target" particle (same notations)'
      read(5,1) TARG

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
      
      if(PROJ.eq.'A'.or.TARG.eq.'A') then
        write(6,*)'Enter Min. and Max. values of impact parameter (fm)'
        read(5,*)  BMIN, BMAX
      else
        BMIN=0.0
        BMAX=0.0
      endif

      if(PROJ.eq.'A'.or.TARG.eq.'A') then
        write(6,*)
        write(6,*)' Would you like to take jet quenching into account?'
        write(6,*)'                 Y - Yes, N - No?   '

        read(5,2) KEY
2       format(A1)

        if(KEY.eq.'Y'.or.KEY.eq.'y') then
          IHPR2(4)=1 
        else
          IHPR2(4)=0 
        endif
      endif

C....set the pt range of the triggered jets:
      write(6,*)
      write(6,*)' Set the Pt of triggered jets (in GeV/c)'
      read(5,*)   Pt_trigger

      write(6,*)
      write(6,*)'Enter number of events'
      read(5,*)  N_events

C....initialize HIJING for requested interactions
      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)

      write(6,*)' Simulation of interactions with'             
      write(6,*)                                               
      write(6,*)' Proj = ',PROJ,' and  Targ = ',TARG            
      write(6,*)' IAP  =',IAP  ,'            IAT  =',IAT       
      write(6,*)' IZP  =',IZP  ,'             IZT  =',IZT       
      write(6,*)                                               
      write(6,*)' Reference frame -   ',FRAME                  
      write(6,*)' ENERGY            ',EFRM,' GeV'                     
      write(6,*)' Number of generated events',N_events       
      write(6,*)' Triggered jets Pt =>',Pt_trigger,' (GeV/c)'
      write(6,*)                                               

      DO 100 I_event=1,N_events

        write(6,*)'Event # ',I_event,' -------------------------'
        write(6,*)

        HIPR1(10)=-Pt_trigger

        CALL HIJING(FRAME,BMIN,BMAX)

C....print out flavor code of the first jet:
        write(6,*)'Flavor code of the first jet: ',IHNT2(9)

C....and its four momentum:
        write(6,*)'     Px           Py           Pz',
     ,'            E     (GeV/c, GeV)'
        write(6,3)HINT1(21),HINT1(22),HINT1(23),HINT1(24)
3       format(4(e11.4,2x))        
        write(6,*)

C....print out flavor code of the second jet:
        write(6,*)'Flavor code of the second jet:',IHNT2(10)

C....and its four momentum:
        write(6,*)'     Px           Py           Pz',
     ,'            E     (GeV/c, GeV)'
        write(6,3)HINT1(31),HINT1(32),HINT1(33),HINT1(34)
        write(6,*)

        PAUSE
100   CONTINUE
 
      STOP
      END    

      FUNCTION RAN(NSEED)                                      
      RAN=RLU(NSEED)                                           
      RETURN                                                   
      END                                                    

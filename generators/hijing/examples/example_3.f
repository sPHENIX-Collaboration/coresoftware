C****************************************************************************
C          Program # 2 from Comp. Phys. Commun. 83 (1994) 307
C	           by M. Gyulassy and X-.N. Wang
C             Modified by V.Uzhinsky, CERN, Oct. 2003	
C***************************************************************************

      CHARACTER FRAME*8,PROJ*8,TARG*8

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/

      COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/

      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
      SAVE  /HIMAIN2/
C
      DIMENSION GB(101), XB(101), DNDP(50)

      COMMON/RANSEED/NSEED                                     
      SAVE  /RANSEED/                                          

      NSEED=0                                                 

      do i=1,50
         DNDP(i)=0.
      enddo

      FRAME='CMS'
      
      write(6,*)'===================================================='
      write(6,*)'  Calculation of transverse momentum distribution   '
      write(6,*)'        of charged pions in minimum bias            '
      write(6,*)'            hh-, hA- and AA-collisions              '
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
      
      write(6,*)
      write(6,*)'Enter number of events per each value of impact',
     ,' parameter (e.a. 10)'

      read(5,*)  N_EVENT

C....initialize HIJING for requested interactions
      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)

      write(6,*)' Simulation of interactions with'             
      write(6,*)                                               
      write(6,*)' Proj = ',PROJ,' and  Targ = ',TARG            
      write(6,*)' IAP  =',IAP  ,'            IAT  =',IAT       
      write(6,*)' IZP  =',IZP  ,'            IZT  =',IZT       
      write(6,*)                                               
      write(6,*)' Reference frame -   ',FRAME                  
      write(6,*)' ENERGY            ',EFRM,' GeV'                     
      write(6,*)' Number of generated events per B interval -',N_event       
      write(6,*)                                               

C....set BMIN=0 and BMAX=R_A+R_B

      BMIN=0.0
      BMAX=HIPR1(34)+HIPR1(35)

C....calculate the Glauber probability and its integrated value:

      DIP=(BMAX-BMIN)/100.0
      GBTOT=0.0
      DO 100 I=1,101
         XB(I)=BMIN+(I-1)*DIP
         OV=PROFILE(XB(I))
         GB(I)=XB(I)*(1.0-exp(-HINT1(12)*OV))
         GBTOT=GBTOT+GB(I)
100   CONTINUE
      write(6,*)'Inelastic X-section (mb) ',GBTOT*DIP*10.*6.28
	
      PAUSE
C....generating N_EVENT for each of 100 impact parameter intervals:

      NONT=0
      GNORM=GBTOT
      
      DO 300 IB=1,100
         B1=XB(iB)
         B2=XB(IB+1)

C....normalized Glauber probability:

         W_GB=(GB(IB)+GB(IB+1))/2.0/GBTOT

         DO 200 IE=1,N_EVENT
            CALL HIJING(FRAME,B1,B2)

C....count number of events without any interaction
C....and renormalize the total Glauber prabability:

            IF(NP+NT .EQ. 0) THEN
              NONT=NONT+1
              GNORM=GNORM-GB(IB)/FLOAT(N_EVENT)
              GO TO 200
            ENDIF

C....calculate pt distribution of charged pions:

            DO 150 K=1,NATT
C....select charged pions only:
               IF(ABS(KATT(K,1)) .NE. 211) GO TO 150

C....calculate pt:
               PTR=SQRT(PATT(K,1)**2+PATT(K,2)**2)

C....calculate pt distribution and weight with normalized
C....Glauber probability to get minimum bias results:

               IF(PTR .GE. 10.0) GO TO 150
               IPT=1+PTR/0.2
               DNDP(IPT)=DNDP(IPT)+W_GB/FLOAT(N_EVENT)/0.2
150         CONTINUE
200      CONTINUE
300   CONTINUE   

C....renormalize the distribution by the renormalized Glauber
C....probability which excludes the events without any interaction:

      IF(NONT.NE.0) THEN
        DO 400 I=1,50
           DNDP(I)=DNDP(I)*GBTOT/GNORM
           if(DNDP(I).ne.0.) write(6,350)0.2*(i-1),DNDP(I)
350        format(1x,f5.1,2x,e11.4)
400     CONTINUE
      ENDIF

      STOP
      END    

      FUNCTION RAN(NSEED)                                      
      RAN=RLU(NSEED)                                           
      RETURN                                                   
      END                                                    

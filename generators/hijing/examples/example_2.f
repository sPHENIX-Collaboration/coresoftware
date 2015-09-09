C****************************************************************************
C          Program # 1 from Comp. Phys. Commun. 83 (1994) 307
C	           by M. Gyulassy and X-.N. Wang
C             Modified by V.Uzhinsky, CERN, Oct. 2003	
C***************************************************************************

      CHARACTER FRAME*8,PROJ*8,TARG*8

      DIMENSION DNDPT(50), DNDY(50)

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/

C....information of produced particles:

      COMMON/HIMAIN1/NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      SAVE  /HIMAIN1/

      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
      SAVE  /HIMAIN2/
C
C....information of produced partons:

      COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),PJPY(300,500)
     &		     ,PJPZ(300,500),PJPE(300,500),PJPM(300,500)
     &		     ,NTJ(300),KFTJ(300,500),PJTX(300,500),PJTY(300,500)
     &		     ,PJTZ(300,500),PJTE(300,500),PJTM(300,500)
      SAVE  /HIJJET1/

      COMMON/HIJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100)
     &         ,K2SG(900,100),PXSG(900,100),PYSG(900,100),PZSG(900,100)
     &         ,PESG(900,100),PMSG(900,100)
      SAVE  /HIJJET2/
C 
      COMMON/HISTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
      SAVE  /HISTRNG/

      COMMON/RANSEED/NSEED                                     
      SAVE  /RANSEED/                                          

      NSEED=0                                                 

      do k=1,50
        DNDPT(K)=0.
        DNDY(K) =0.
      enddo

      FRAME='CMS'
      
      write(6,*)'===================================================='
      write(6,*)'  Calculation of transverse momentum and rapidity   '
      write(6,*)'       distributions of charged particles in        '
      write(6,*)'hh-, hA- and AA-collisions at fixed impact parameter'
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
      write(6,*)'Enter number of events'

      read(5,*)  N_EVENT

      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)

      write(6,*)' Simulation of interactions with'             
      write(6,*)                                               
      write(6,*)' Proj = ',PROJ,' and  Targ = ',TARG            
      write(6,*)' IAP  =',IAP  ,'            IAT  =',IAT       
      write(6,*)' IZP  =',IZP  ,'            IZT  =',IZT       
      write(6,*)                                               
      write(6,*)' Reference frame -   ',FRAME                  
      write(6,*)' ENERGY            ',EFRM,' (GeV)'                     
      write(6,*)' Number of generated events -',N_event       
      write(6,*)                                               

      if(PROJ.eq.'A'.or.TARG.eq.'A') then
        write(6,*)'Enter Min. and Max. values of impact parameter (fm)'
        read(5,*)  BMIN, BMAX
      else
        BMIN=0.0
        BMAX=0.0
      endif

      DO 2000 J=1,N_EVENT

        write(6,*)' Event # ',J,' ------------------------------'

        CALL HIJING(FRAME,BMIN,BMAX)
C
C....calculate rapidity and transverse momentum distributions of
C....produced charged particles:  

        DO 1000 I=1,NATT

C....exclude beam nucleons as produced particles:

          if(KATT(I,2).EQ.0 .OR. KATT(I,2).EQ.10) GO TO 1000

C....select charged particles only:

          IF(LUCHGE(KATT(I,1)) .EQ. 0) GO TO 1000

          PTR=SQRT(PATT(I,1)**2+PATT(I,2)**2)
          IF (PTR .GE. 10.0) GO TO 100

          IPT=1+PTR/0.2
          DNDPT(IPT)=DNDPT(IPT)+1.0/FLOAT(N_EVENT)/0.2/2.0/PTR

100       Y=0.5*LOG((PATT(I,4)+PATT(I,3))/(PATT(I,4)-PATT(I,3)))
          IF(ABS(Y) .GE. 10.0) GO TO 1000

          IY=1+ABS(Y)/0.2
          DNDY(IY)=DNDY(IY)+1.0/FLOAT(N_EVENT)/0.2/2.0
1000    CONTINUE
2000  CONTINUE

C....print out the rapidity and transverse momentum distributions:

      do k=1,50
        WRITE(6,2)0.2*(K-1),DNDPT(K),DNDY(K)
2       format(1x,f5.1,2(2x,e11.4))
      enddo

      STOP
      END    

      FUNCTION RAN(NSEED)                                      
      RAN=RLU(NSEED)                                           
      RETURN                                                   
      END                                                    

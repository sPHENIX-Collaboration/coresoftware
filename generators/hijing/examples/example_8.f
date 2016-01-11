
C****************************************************************************
C     The program for calculation of the general properties of
C             partons in the interactions by the HIJING model.
C	          Written by V.Uzhinsky, CERN, Oct. 2003
C***************************************************************************

      CHARACTER FRAME*8,PROJ*8,TARG*8

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)     
      SAVE  /HIPARNT/
C
C         ********Event Options and parameters

C....information of produced partons:
      
      COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),PJPY(300,500)
     &               ,PJPZ(300,500),PJPE(300,500),PJPM(300,500)
     &               ,NTJ(300),KFTJ(300,500),PJTX(300,500),PJTY(300,500)
     &               ,PJTZ(300,500),PJTE(300,500),PJTM(300,500)
      SAVE  /HIJJET1/

      COMMON/HIJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100)
     &         ,K2SG(900,100),PXSG(900,100),PYSG(900,100),PZSG(900,100)
     &         ,PESG(900,100),PMSG(900,100)
      SAVE  /HIJJET2/
C
      CHARACTER *1 KEY
      CHARACTER *12 FNAME

      PARAMETER (NWPAWC=1000000)
      COMMON/PAWC/HMEMORY(NWPAWC)

      COMMON/RANSEED/NSEED
      SAVE  /RANSEED/
      NSEED=0

      write(6,*)'====================================================='
      write(6,*)'The program for calculation of general properties of '
      write(6,*)'           partons in the HIJING model.              '
      write(6,*)'====================================================='
      write(6,*)
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

      write(6,*)'Enter the corresponding energy per NN collisions (GeV)'

      read(5,*) EFRM

      write(6,*)
      write(6,*)'Enter a type of the projectile particle'
      write(6,*)
      write(6,*)' P proton,            PBAR anti-proton,'
      write(6,*)' N neutron,           NBAR anti-neutron,'
      write(6,*)' PI+ - positive pion, PI- negative pion,'
      write(6,*)' K+ positive kaon,    K- negative kaon'
      write(6,*)' A - nucleus --------------------------'

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
      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
C                       ********Initialize HIJING

      WRITE(6,*)' Simulation of interactions with'
      WRITE(6,*)
      WRITE(6,*)' Proj = ',PROJ,' and  Targ = ',TARG
      WRITE(6,*)' IAP  =',IAP  ,'            IAT  =',IAT
      WRITE(6,*)' IZP  =',IZP  ,'            IZT  =',IZT
      WRITE(6,*)
      WRITE(6,*)' Reference frame -   ',FRAME
      WRITE(6,*)' ENERGY            ',EFRM,' GeV'
      WRITE(6,*)' Number of generated events -',N_events
      WRITE(6,*)
      
      BMIN=0.   
      BMAX=HIPR1(34)+HIPR1(35)

      WRITE(6,*)
      
c------------------------------------------------------------------
      CALL HLIMIT(NWPAWC)
c------------------------------------------------------------------
c----------------- Projectile associated parton -------------------
c------------------------------------------------------------------

      if(PROJ.eq.'A'.or.TARG.eq.'A') then
        CH_max=10*MAX(IAP,IAT)
      else
        CH_max=99.5
      endif

      CALL HBOOK1(110,' Projectile parton multiplicity distribution',
     ,100,-0.5,CH_max,0.)

      if(FRAME.eq.'CMS') then
        Eta_l=-15.
        Eta_h=+15.
      else
        Eta_l=-2.0
        Eta_h=Alog(2.*EFRM)+2.
      endif
      NbinEta=(Eta_h-Eta_l)/0.2

      CALL HBOOK1(120,' Projectile parton pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      if(FRAME.eq.'CMS') then
        Y_l=-Alog(EFRM)-2.
        Y_h=+Alog(EFRM)+2.
      else
        Y_l=-2.0
        Y_h=Alog(2.*EFRM)+2.
      endif
      NbinY=(Y_h-Y_l)/0.2

      CALL HBOOK1(130,' Projectile parton rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(140,' Projectile parton Pt distribution',
     ,100,0.,10.,0.)

      if(FRAME.eq.'CMS') then
        E_l= 0.
        E_h=EFRM/2.
      else
        E_l= 0.
        E_h=EFRM
      endif
      NbinE=(E_h-E_l)/0.5

      CALL HBOOK1(150,' Projectile parton energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(160,' Projectile parton Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(170,' Projectile parton Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(180,' Phi - Eta correlation of Projectile parton',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

      CALL HBOOK1(190,' ID distribution of projectile partons',
     , 28,-6.5,21.5,0.)
C----------------------------------------------------------------
c----------------- Target associated partons  -------------------
c----------------------------------------------------------------
      CALL HBOOK1(210,' Target parton multiplicity distribution',
     ,100,-0.5,CH_max,0.)

      CALL HBOOK1(220,' Target parton pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(230,' Target parton rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(240,' Target parton Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(250,' Target parton energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(260,' Target parton Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(270,' Target parton Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(280,' Phi - Eta correlation of Target parton',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

      CALL HBOOK1(290,' ID distribution of target partons',
     , 28,-6.5,21.5,0.)
c----------------------------------------------------------------
c-------------- "Central" produced partons ----------------------
c----------------------------------------------------------------
      CALL HBOOK1(310,' Central parton multiplicity distribution',
     ,100,-0.5,CH_max,0.)

      CALL HBOOK1(320,' Central parton pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(330,' Central parton rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(340,' Central parton Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(350,' Central parton energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(360,' Central parton Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(370,' Central parton Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(380,' Phi - Eta correlation of Central parton',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

      CALL HBOOK1(390,' ID distribution of central partons',
     , 28,-6.5,21.5,0.)
c----------------------------------------------------------------
c----------------------------------------------------------------
      Pi=4.*ATAN(1.)                   ! 3.14159

      DO I_event=1,N_events

        WRITE(6,*)' Event # ',I_event,' ------------------'
C
        CALL HIJING(FRAME,BMIN,BMAX)
C
        N_prp=0
        N_trp=0
        N_crp=0

        DO I=1,IAP

          N_prp=N_prp+NPJ(I)
          
          do j=1,NPJ(I)

            ID=KFPJ(I,J)
            E=PJPE(I,J)
            Pz=PJPZ(I,J)
            Pt=sqrt(PJPX(I,J)**2+PJPY(I,J)**2)
            P=sqrt(Pz**2+Pt**2)

            if(P-Pz.ge.1.0e-4) then
              Eta= 0.5*Alog((P+Pz)/(P-Pz))
            elseif(P+Pz.ge.1.0e-4) then
              Eta=-0.5*Alog((P-Pz)/(P+pz))
            else
              Eta=15.
            endif

            if(E-Pz.ge.1.0e-4) then
              Y= 0.5*Alog((E+Pz)/(E-Pz))
            elseif(E+Pz.ge.1.0e-4) then
              Y=-0.5*Alog((E-Pz)/(E+pz))
            else
              Y=15.
            endif

            if(P.ge.1.0e-4) then
              Cos_Theta=Pz/P
            else
              Cos_Theta=1.
            endif

            Phi=ULANGL(PJPX(I,J),PJPY(I,J))*180./Pi
c==========================================================

            Call HF1(120,Eta,1.)
            Call HF1(130,  Y,1.)
            Call HF1(140, Pt,1.)
            Call HF1(150,  E,1.)
            Call HF1(160,Cos_Theta,1.)
            Call HF1(170,Phi,1.)
            Call HFILL(180,Eta,Phi,1.)
            Call Hf1(190,Float(Id),1.)
          enddo
        ENDDO

        DO I=1,IAT

          N_trp=N_trp+NTJ(I)
          
          do j=1,NTJ(I)

            ID=KFTJ(I,J)
            E=PJTE(I,J)
            Pz=PJTZ(I,J)
            Pt=sqrt(PJTX(I,J)**2+PJTY(I,J)**2)
            P=sqrt(Pz**2+Pt**2)

            if(P-Pz.ge.1.0e-4) then
              Eta= 0.5*Alog((P+Pz)/(P-Pz))
            elseif(P+Pz.ge.1.0e-4) then
              Eta=-0.5*Alog((P-Pz)/(P+pz))
            else
              Eta=15.
            endif

            if(E-Pz.ge.1.0e-4) then
              Y= 0.5*Alog((E+Pz)/(E-Pz))
            elseif(E+Pz.ge.1.0e-4) then
              Y=-0.5*Alog((E-Pz)/(E+pz))
            else
              Y=15.
            endif

            if(P.ge.1.0e-4) then
              Cos_Theta=Pz/P
            else
              Cos_Theta=1.
            endif

            Phi=ULANGL(PJTX(I,J),PJTY(I,J))*180./Pi
c==========================================================

            Call HF1(220,Eta,1.)
            Call HF1(230,  Y,1.)
            Call HF1(240, Pt,1.)
            Call HF1(250,  E,1.)
            Call HF1(260,Cos_Theta,1.)
            Call HF1(270,Phi,1.)
            Call HFILL(280,Eta,Phi,1.)
            Call Hf1(290,Float(Id),1.)
          enddo
        ENDDO  

        DO I=1,NSG

          N_crp=N_crp+NJSG(I)
          
          do j=1,NJSG(I)

            ID=K2SG(I,J)
            E=PESG(I,J)
            Pz=PZSG(I,J)
            Pt=sqrt(PXSG(I,J)**2+PYSG(I,J)**2)
            P=sqrt(Pz**2+Pt**2)

            if(P-Pz.ge.1.0e-4) then
              Eta= 0.5*Alog((P+Pz)/(P-Pz))
            elseif(P+Pz.ge.1.0e-4) then
              Eta=-0.5*Alog((P-Pz)/(P+pz))
            else
              Eta=15.
            endif

            if(E-Pz.ge.1.0e-4) then
              Y= 0.5*Alog((E+Pz)/(E-Pz))
            elseif(E+Pz.ge.1.0e-4) then
              Y=-0.5*Alog((E-Pz)/(E+pz))
            else
              Y=15.
            endif

            if(P.ge.1.0e-4) then
              Cos_Theta=Pz/P
            else
              Cos_Theta=1.
            endif

            Phi=ULANGL(PXSG(I,J),PYSG(I,J))*180./Pi
c==========================================================

            Call HF1(320,Eta,1.)
            Call HF1(330,  Y,1.)
            Call HF1(340, Pt,1.)
            Call HF1(350,  E,1.)
            Call HF1(360,Cos_Theta,1.)
            Call HF1(370,Phi,1.)
            Call HFILL(380,Eta,Phi,1.)
            Call Hf1(390,Float(Id),1.)
          enddo
        ENDDO
       
        CALL HF1(110,float(N_prp),1.)
        CALL HF1(210,float(N_trp),1.)
        CALL HF1(310,float(N_crp),1.)

      ENDDO

c================ Normalization ===========================

      C1=1./float(N_events)                        ! Multiplicity distr.
      C2=0.

      Call HOPERA(110,'+',110,110,C1,C2)
      Call HOPERA(210,'+',210,210,C1,C2)
      Call HOPERA(310,'+',310,310,C1,C2)

      C1=1./float(N_events)/((Eta_h-Eta_l)/NbinEta)! Eta distr.
      C2=0.

      Call HOPERA(120,'+',120,120,C1,C2)
      Call HOPERA(220,'+',220,220,C1,C2)
      Call HOPERA(320,'+',320,320,C1,C2)

      C1=1./float(N_events)/((Y_h-Y_l)/NbinY)      ! Y distr.
      C2=0.

      Call HOPERA(130,'+',130,130,C1,C2)
      Call HOPERA(230,'+',230,230,C1,C2)
      Call HOPERA(330,'+',330,330,C1,C2)

      C1=1./float(N_events)/0.1                    ! Pt distr.
      C2=0.

      Call HOPERA(140,'+',140,140,C1,C2)
      Call HOPERA(240,'+',240,240,C1,C2)
      Call HOPERA(340,'+',340,340,C1,C2)

      C1=1./float(N_events)/((E_h-E_l)/NbinE)      ! E distr.
      C2=0.

      Call HOPERA(150,'+',150,150,C1,C2)
      Call HOPERA(250,'+',250,250,C1,C2)
      Call HOPERA(350,'+',350,350,C1,C2)

      C1=1./float(N_events)/0.05                   ! Cos Theta distr.
      C2=0.

      Call HOPERA(160,'+',160,160,C1,C2)
      Call HOPERA(260,'+',260,260,C1,C2)
      Call HOPERA(360,'+',360,360,C1,C2)

      C1=1./float(N_events)                        ! Phi distr.
      C2=0.

      Call HOPERA(170,'+',170,170,C1,C2)
      Call HOPERA(270,'+',270,270,C1,C2)
      Call HOPERA(370,'+',370,370,C1,C2)

c================ Writing results =========================

      CALL HRPUT(0,FNAME,'N')
      END

      FUNCTION RAN(NSEED)                                            
      RAN=RLU(NSEED)                                           
      RETURN                                                   
      END

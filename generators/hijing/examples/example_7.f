
C****************************************************************************
C     The program for calculation of the general properties of
C    neutral particles in the interactions by the HIJING model.
C	          Written by V.Uzhinsky, CERN, Oct. 2003
C***************************************************************************

      CHARACTER FRAME*8,PROJ*8,TARG*8

      COMMON/HIMAIN1/ NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)

C         ********information of produced particles

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)     
      SAVE  /HIPARNT/
C
C         ********Event Options and parameters

      CHARACTER *1 KEY
      CHARACTER *12 FNAME

      PARAMETER (NWPAWC=1000000)
      COMMON/PAWC/HMEMORY(NWPAWC)

      COMMON/RANSEED/NSEED
      SAVE  /RANSEED/
      NSEED=0

      write(6,*)'====================================================='
      write(6,*)'The program for calculation of general properties of '
      write(6,*)'   neutral particles (Pi_0, Lambda0, K0_S, Phi       '
      write(6,*)'             in the HIJING model.                    '
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
c------------------------ Pi0 mesons ------------------------------
c------------------------------------------------------------------

      Nu_max=MIN(IAP*IAT+1,2000)

      if(PROJ.eq.'A'.or.TARG.eq.'A') then
        CH_max=5*Nu_max/2.
      else
        CH_max=99.5
      endif

      CALL HBOOK1(10,' Pi0 meson multiplicity distribution',
     ,100,-0.5,CH_max,0.)

      if(FRAME.eq.'CMS') then
        Eta_l=-15.
        Eta_h=+15.
      else
        Eta_l=-2.0
        Eta_h=Alog(2.*EFRM)+2.
      endif
      NbinEta=(Eta_h-Eta_l)/0.2

      CALL HBOOK1(20,' Pi0 meson pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      if(FRAME.eq.'CMS') then
        Y_l=-Alog(EFRM)-2.
        Y_h=+Alog(EFRM)+2.
      else
        Y_l=-2.0
        Y_h=Alog(2.*EFRM)+2.
      endif
      NbinY=(Y_h-Y_l)/0.2

      CALL HBOOK1(30,' Pi0 meson rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(40,' Pi0 meson Pt distribution',
     ,100,0.,10.,0.)

      if(FRAME.eq.'CMS') then
        E_l= 0.
        E_h=EFRM/2.
      else
        E_l= 0.
        E_h=EFRM
      endif
      NbinE=(E_h-E_l)/0.5

      CALL HBOOK1(50,' Pi0 meson energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(60,' Pi0 meson Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(70,' Pi0 meson Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(80,' Phi - Eta correlation of Pi0 meson',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

C----------------------------------------------------------------
c--------------------------- Lambda0 ---------------------------
c----------------------------------------------------------------
      CALL HBOOK1(110,' Lambda0 multiplicity distribution',
     ,100,-0.5,float(IAP+IAT),0.)

      CALL HBOOK1(120,' Lambda0 pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(130,' Lambda0 rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(140,' Lambda0 Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(150,' Lambda0 energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(160,' Lambda0 Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(170,' Lambda0 Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(180,' Phi - Eta correlation of Lambda0',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

c----------------------------------------------------------------
c------------------------------- K0_S ---------------------------
c----------------------------------------------------------------
      CALL HBOOK1(210,' K0S multiplicity distribution',
     ,100,-0.5,CH_max,0.)

      CALL HBOOK1(220,' K0S pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(230,' K0S rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(240,' K0S Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(250,' K0S energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(260,' K0S Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(270,' K0S Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(280,' Phi - Eta correlation of K0S',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

c----------------------------------------------------------------
c------------------------- Phi meson  ---------------------------
c----------------------------------------------------------------
      CALL HBOOK1(310,' Phi meson multiplicity distribution',
     ,100,-0.5,CH_max/3.,0.)

      CALL HBOOK1(320,' Phi meson pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(330,' Phi meson rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(340,' Phi meson Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(350,' Phi meson energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(360,' Phi meson Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(370,'  Phi meson Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(380,' Phi - Eta correlation of Phi meson',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

c----------------------------------------------------------------

      IHPR2(12)=1                   ! To suppress the particle decays
      CALL LUGIVE('MDCY(C333,1)=0') ! To suppress Phi decay
c     CALL LUGIVE('MDCY(C1114,1)=0;MDCY(C-1114,1)=0') ! To suppress Delta- Decay

c----------------------------------------------------------------
      Pi=4.*ATAN(1.)                   ! 3.14159

      DO 2000 I_event=1,N_events

        WRITE(6,*)' Event # ',I_event,' ------------------'
C
        CALL HIJING(FRAME,BMIN,BMAX)
C
        N_pi0=0
        N_lam=0
        N_k0s=0
        N_phi=0

c        write(6,*)' NATT --------- ',NATT
        DO 3000 I=1,NATT

c          write(6,*)i,KATT(I,1)

          if(KATT(I,2).eq. 0) go to 3000 ! reject non-interacting projectile
          if(KATT(I,2).eq. 1) go to 3000 ! reject elastic scattering
          if(KATT(I,2).eq.10) go to 3000 ! reject non-interacting target
          if(KATT(I,2).eq.11) go to 3000 ! reject elastic scattering

          ID=KATT(i,1)

          ICH=LUCHGE(KATT(I,1))/3

          Amass=ULMASS(KATT(I,1))

          E=PATT(i,4)
          Pz=PATT(i,3)
          Pt=sqrt(PATT(i,1)**2+PATT(i,2)**2)
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

          Phi=ULANGL(PATT(i,1),PATT(i,2))*180./Pi
c==========================================================

          if(ID.eq.111) then             ! Pi_0
            N_pi0=N_pi0+1

            Call HF1(20,Eta,1.)
            Call HF1(30,  Y,1.)
            Call HF1(40, Pt,1.)
            Call HF1(50,  E,1.)
            Call HF1(60,Cos_Theta,1.)
            Call HF1(70,Phi,1.)
            Call HFILL(80,Eta,Phi,1.)
            Call Hf1(90,Float(Id),1.)
          endif

          if(ID.eq.3122) then         ! Lambda0
            N_lam=N_lam+1

            Call HF1(120,Eta,1.)
            Call HF1(130,  Y,1.)
            Call HF1(140, Pt,1.)
            Call HF1(150,  E,1.)
            Call HF1(160,Cos_Theta,1.)
            Call HF1(170,Phi,1.)
            Call HFILL(180,Eta,Phi,1.)
          endif

          if(ID.eq.310) then          ! K0_S
            N_k0s=N_k0s+1

            Call HF1(220,Eta,1.)
            Call HF1(230,  Y,1.)
            Call HF1(240, Pt,1.)
            Call HF1(250,  E,1.)
            Call HF1(260,Cos_Theta,1.)
            Call HF1(270,Phi,1.)
            Call HFILL(280,Eta,Phi,1.)
          endif

          if(Id.eq.333) then          ! Phi
            N_phi=N_phi+1

            Call HF1(320,Eta,1.)
            Call HF1(330,  Y,1.)
            Call HF1(340, Pt,1.)
            Call HF1(350,  E,1.)
            Call HF1(360,Cos_Theta,1.)
            Call HF1(370,Phi,1.)
            Call HFILL(380,Eta,Phi,1.)
          endif

3000    CONTINUE

c       PAUSE       
        CALL HF1( 10,float(N_pi0),1.)
        CALL HF1(110,float(N_lam),1.)
        CALL HF1(210,float(N_k0s),1.)
        CALL HF1(310,float(N_phi),1.)

2000  CONTINUE

c================ Normalization ===========================

      C1=1./float(N_events)                        ! Multiplicity distr.
      C2=0.

      Call HOPERA( 10,'+', 10, 10,C1,C2)
      Call HOPERA(110,'+',110,110,C1,C2)
      Call HOPERA(210,'+',210,210,C1,C2)
      Call HOPERA(310,'+',310,310,C1,C2)

      C1=1./float(N_events)/((Eta_h-Eta_l)/NbinEta)! Eta distr.
      C2=0.

      Call HOPERA( 20,'+', 20, 20,C1,C2)
      Call HOPERA(120,'+',120,120,C1,C2)
      Call HOPERA(220,'+',220,220,C1,C2)
      Call HOPERA(320,'+',320,320,C1,C2)

      C1=1./float(N_events)/((Y_h-Y_l)/NbinY)      ! Y distr.
      C2=0.

      Call HOPERA( 30,'+', 30, 30,C1,C2)
      Call HOPERA(130,'+',130,130,C1,C2)
      Call HOPERA(230,'+',230,230,C1,C2)
      Call HOPERA(330,'+',330,330,C1,C2)

      C1=1./float(N_events)/0.1                    ! Pt distr.
      C2=0.

      Call HOPERA( 40,'+', 40, 40,C1,C2)
      Call HOPERA(140,'+',140,140,C1,C2)
      Call HOPERA(240,'+',240,240,C1,C2)
      Call HOPERA(340,'+',340,340,C1,C2)

      C1=1./float(N_events)/((E_h-E_l)/NbinE)      ! E distr.
      C2=0.

      Call HOPERA( 50,'+', 50, 50,C1,C2)
      Call HOPERA(150,'+',150,150,C1,C2)
      Call HOPERA(250,'+',250,250,C1,C2)
      Call HOPERA(350,'+',350,350,C1,C2)

      C1=1./float(N_events)/0.05                   ! Cos Theta distr.
      C2=0.

      Call HOPERA( 60,'+', 60, 60,C1,C2)
      Call HOPERA(160,'+',160,160,C1,C2)
      Call HOPERA(260,'+',260,260,C1,C2)
      Call HOPERA(360,'+',360,360,C1,C2)

      C1=1./float(N_events)                        ! Phi distr.
      C2=0.

      Call HOPERA( 70,'+', 70, 70,C1,C2)
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

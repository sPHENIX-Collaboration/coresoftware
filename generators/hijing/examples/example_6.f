
C****************************************************************************
C     The program for calculation of the general properties of
C                   interactions in the HIJING model.
C	          Written by V.Uzhinsky, CERN, Oct. 2003
C***************************************************************************

      CHARACTER FRAME*8,PROJ*8,TARG*8

      COMMON/HIMAIN1/ NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
      COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)

C         ********information of produced particles
C

      COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      SAVE  /HIPARNT/

      CHARACTER *1 KEY
      CHARACTER *12 FNAME

      PARAMETER (NWPAWC=1000000)
      COMMON/PAWC/HMEMORY(NWPAWC)

      COMMON/RANSEED/NSEED
      SAVE  /RANSEED/
      NSEED=0

      write(6,*)'====================================================='
      write(6,*)'The program for calculation of general properties of '
      write(6,*)'         interactions in the HIJING model.           '
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
C			********Initialize HIJING

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

c------------------------------------------------------------------
c----------------- Charged particles ------------------------------
c------------------------------------------------------------------
      CALL HLIMIT(NWPAWC)

      CALL HBOOK1(2,'Impact parameter distribution',
     ,100,BMIN,BMAX,0.)

      Nu_max=MIN(IAP*IAT+1,2000)
      CALL HBOOK1(4,'Nu-distribution',
     ,100,-0.5,float(Nu_max),0.)

      CALL HBOOK1(6,'Distribution on the number of wounded A nucleons',
     ,IAP+1,-0.5,float(IAP)+0.5,0.) 

      CALL HBOOK1(8,'Distribution on the number of wounded B nucleons',
     ,IAT+1,-0.5,float(IAT)+0.5,0.)

      if(PROJ.eq.'A'.or.TARG.eq.'A') then
        CH_max=5*Nu_max
        M_neg=CH_max/2.
      else
        CH_max=99.5
        M_neg=CH_max
      endif
      CALL HBOOK1(10,' Charged particle multiplicity distribution',
     ,100,-0.5,CH_max,0.)

      if(FRAME.eq.'CMS') then
        Eta_l=-15.
        Eta_h=+15.
      else
        Eta_l=-2.0
        Eta_h=Alog(2.*EFRM)+2.
      endif
      NbinEta=(Eta_h-Eta_l)/0.2

      CALL HBOOK1(20,' Charged particle pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      if(FRAME.eq.'CMS') then
        Y_l=-Alog(EFRM)-2.
        Y_h=+Alog(EFRM)+2.
      else
        Y_l=-2.0
        Y_h=Alog(2.*EFRM)+2.
      endif
      NbinY=(Y_h-Y_l)/0.2

      CALL HBOOK1(30,' Charged particle rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(40,' Charged particle Pt distribution',
     ,100,0.,10.,0.)

      if(FRAME.eq.'CMS') then
        E_l= 0.
        E_h=EFRM/2.
      else
        E_l= 0.
        E_h=EFRM
      endif
      NbinE=(E_h-E_l)/0.5

      CALL HBOOK1(50,' Charged particle energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(60,' Charged particle Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(70,' Charged particle Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(80,' Phi - Eta correlation of charged particles',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

      CALL HBOOK1(90,' Particle composition (IDs)',
     ,6600,-3300.,3300.,0.)

C----------------------------------------------------------------
c--------------- Negative charged particles ---------------------
c----------------------------------------------------------------
      CALL HBOOK1(110,' Negative particle multiplicity distribution',
     ,100,-0.5,float(M_neg),0.)

      CALL HBOOK1(120,' Negative particle pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(130,' Negative particle rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(140,' Negative particle Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(150,' Negative particle energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(160,' Negative particle Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(170,' Negative particle Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(180,' Phi - Eta correlation of negative particles',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

c----------------------------------------------------------------
c----------------Positive charged particles ---------------------
c----------------------------------------------------------------
      CALL HBOOK1(210,' Positive particle multiplicity distribution',
     ,100,-0.5,float(M_neg),0.)

      CALL HBOOK1(220,' Positive particle pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(230,' Positive particle rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(240,' Positive particle Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(250,' Positive particle energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(260,' Positive particle Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(270,' Positive particle Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(280,' Phi - Eta correlation of positive particles',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)

c----------------------------------------------------------------
c------------------------- Protons ------------------------------
c----------------------------------------------------------------
      CALL HBOOK1(310,' Proton multiplicity distribution',
     ,100,-0.5,float(IZP+IZT)+10.,0.)

      CALL HBOOK1(320,' Proton pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(330,' Proton rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(340,' Proton Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(350,' Proton energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(360,' Proton Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(370,' Proton Phi distribution',
     ,180,0.,180.,0.)
      CALL HBOOK2(380,' Phi - Eta correlation of protons',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)


c----------------------------------------------------------------
c----------------------- Gamma quanta ---------------------------
c----------------------------------------------------------------
      CALL HBOOK1(410,' Gamma multiplicity distribution',
     ,100,-0.5,float(M_neg),0.)

      CALL HBOOK1(420,' Gamma pseudo-rapidity distribution',
     ,NbinEta,Eta_l,Eta_h,0.)

      CALL HBOOK1(430,' Gamma rapidity distribution',
     ,NbinY,Y_l,Y_h,0.)

      CALL HBOOK1(440,' Gamma Pt distribution',
     ,100,0.,10.,0.)

      CALL HBOOK1(450,' Gamma energy distribution',
     ,NbinE,E_l,E_h,0.)

      CALL HBOOK1(460,' Gamma Cos(Theta) distribution',
     ,40,-1.,1.,0.)

      CALL HBOOK1(470,' Gamma Phi distribution',
     ,180,0.,180.,0.)

      CALL HBOOK2(480,' Phi - Eta correlation of Gammas',
     ,NbinEta,Eta_l,Eta_h,180,-180.,180.,0.)
c----------------------------------------------------------------

      Pi=4.*ATAN(1.)                   ! 3.14159

      DO 2000 I_event=1,N_events

        WRITE(6,*)' Event # ',I_event,' ------------------'
C
        CALL HIJING(FRAME,BMIN,BMAX)
C
        B=HINT1(19)
        Call HF1(2,B,1.)

        Nu=N0+N01+N10+N11
        Call HF1(4,float(Nu),1.)

        Call HF1(6,float(NP),1.)
        Call HF1(8,float(NT),1.)

        N_ch=0
        N_neg=0
        N_pos=0
        N_protons=0
        N_gammas=0

        DO 3000 I=1,NATT

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

          if(ICH.ne.0) then
            N_ch=N_ch+1

            Call HF1(20,Eta,1.)
            Call HF1(30,  Y,1.)
            Call HF1(40, Pt,1.)
            Call HF1(50,  E,1.)
            Call HF1(60,Cos_Theta,1.)
            Call HF1(70,Phi,1.)
            Call HFILL(80,Eta,Phi,1.)
            Call Hf1(90,Float(Id),1.)
          endif

          if(ICH.lt.0) then
            N_neg=N_neg+1

            Call HF1(120,Eta,1.)
            Call HF1(130,  Y,1.)
            Call HF1(140, Pt,1.)
            Call HF1(150,  E,1.)
            Call HF1(160,Cos_Theta,1.)
            Call HF1(170,Phi,1.)
            Call HFILL(180,Eta,Phi,1.)
          endif

          if(ICH.gt.0) then
            N_pos=N_pos+1

            Call HF1(220,Eta,1.)
            Call HF1(230,  Y,1.)
            Call HF1(240, Pt,1.)
            Call HF1(250,  E,1.)
            Call HF1(260,Cos_Theta,1.)
            Call HF1(270,Phi,1.)
            Call HFILL(280,Eta,Phi,1.)
          endif

          if(Id.eq.2212) then
            N_protons=N_protons+1

            Call HF1(320,Eta,1.)
            Call HF1(330,  Y,1.)
            Call HF1(340, Pt,1.)
            Call HF1(350,  E,1.)
            Call HF1(360,Cos_Theta,1.)
            Call HF1(370,Phi,1.)
            Call HFILL(380,Eta,Phi,1.)
          endif

          if(Id.eq.  22) then
            N_gammas=N_gammas+1

            Call HF1(420,Eta,1.)
            Call HF1(430,  Y,1.)
            Call HF1(440, Pt,1.)
            Call HF1(450,  E,1.)
            Call HF1(460,Cos_Theta,1.)
            Call HF1(470,Phi,1.)
            Call HFILL(480,Eta,Phi,1.)
          endif


3000    CONTINUE


        CALL HF1( 10,float(N_ch),1.)
        CALL HF1(110,float(N_neg),1.)
        CALL HF1(210,float(N_pos),1.)
        CALL HF1(310,float(N_protons),1.)
        CALL HF1(410,float(N_gammas ),1.)

2000  CONTINUE

c================ Normalization ===========================

      C1=1./float(N_events)                        ! Multiplicity distr.
      C2=0.

      Call HOPERA( 10,'+', 10, 10,C1,C2)
      Call HOPERA(110,'+',110,110,C1,C2)
      Call HOPERA(210,'+',210,210,C1,C2)
      Call HOPERA(310,'+',310,310,C1,C2)
      Call HOPERA(410,'+',410,410,C1,C2)

      C1=1./float(N_events)/((Eta_h-Eta_l)/NbinEta)! Eta distr.
      C2=0.

      Call HOPERA( 20,'+', 20, 20,C1,C2)
      Call HOPERA(120,'+',120,120,C1,C2)
      Call HOPERA(220,'+',220,220,C1,C2)
      Call HOPERA(320,'+',320,320,C1,C2)
      Call HOPERA(420,'+',420,420,C1,C2)

      C1=1./float(N_events)/((Y_h-Y_l)/NbinY)      ! Y distr.
      C2=0.

      Call HOPERA( 30,'+', 30, 30,C1,C2)
      Call HOPERA(130,'+',130,130,C1,C2)
      Call HOPERA(230,'+',230,230,C1,C2)
      Call HOPERA(330,'+',330,330,C1,C2)
      Call HOPERA(430,'+',430,430,C1,C2)

      C1=1./float(N_events)/0.1                    ! Pt distr.
      C2=0.

      Call HOPERA( 40,'+', 40, 40,C1,C2)
      Call HOPERA(140,'+',140,140,C1,C2)
      Call HOPERA(240,'+',240,240,C1,C2)
      Call HOPERA(340,'+',340,340,C1,C2)
      Call HOPERA(440,'+',440,440,C1,C2)

      C1=1./float(N_events)/((E_h-E_l)/NbinE)      ! E distr.
      C2=0.

      Call HOPERA( 50,'+', 50, 50,C1,C2)
      Call HOPERA(150,'+',150,150,C1,C2)
      Call HOPERA(250,'+',250,250,C1,C2)
      Call HOPERA(350,'+',350,350,C1,C2)
      Call HOPERA(450,'+',450,450,C1,C2)

      C1=1./float(N_events)/0.05                   ! Cos Theta distr.
      C2=0.

      Call HOPERA( 60,'+', 60, 60,C1,C2)
      Call HOPERA(160,'+',160,160,C1,C2)
      Call HOPERA(260,'+',260,260,C1,C2)
      Call HOPERA(360,'+',360,360,C1,C2)
      Call HOPERA(460,'+',460,460,C1,C2)

      C1=1./float(N_events)                        ! Phi distr.
      C2=0.

      Call HOPERA( 70,'+', 70, 70,C1,C2)
      Call HOPERA(170,'+',170,170,C1,C2)
      Call HOPERA(270,'+',270,270,C1,C2)
      Call HOPERA(370,'+',370,370,C1,C2)
      Call HOPERA(470,'+',470,470,C1,C2)

c================ Writing results =========================

      CALL HRPUT(0,FNAME,'N')
      END

      FUNCTION RAN(NSEED)                                            
      RAN=RLU(NSEED)                                           
      RETURN                                                   
      END

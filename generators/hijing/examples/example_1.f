C****************************************************************************
C
C
C
C	The following is an example program for calling HIJING. one should
C	include all the common blocks and the data values which are listed
C	below in his own program.
C***************************************************************************

        CHARACTER FRAME*8,PROJ*8,TARG*8
	COMMON/HIMAIN1/ NATT,EATT,JATT,NT,NP,N0,N01,N10,N11
        SAVE  /HIMAIN1/

	COMMON/HIMAIN2/KATT(130000,4),PATT(130000,4)
        SAVE  /HIMAIN2/

C         ********information of produced particles
C
	COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),PJPY(300,500)
     &		     ,PJPZ(300,500),PJPE(300,500),PJPM(300,500)
     &		     ,NTJ(300),KFTJ(300,500),PJTX(300,500),PJTY(300,500)
     &		     ,PJTZ(300,500),PJTE(300,500),PJTM(300,500)
        SAVE  /HIJJET1/

        COMMON/HIJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100)
     &         ,K2SG(900,100),PXSG(900,100),PYSG(900,100),PZSG(900,100)
     &         ,PESG(900,100),PMSG(900,100)
        SAVE  /HIJJET2/

C         ********information of produced partons

      COMMON/RANSEED/NSEED                                     ! Uzhi
      SAVE  /RANSEED/                                          ! Uzhi
      NSEED=0                                                  ! Uzhi

      EFRM =200.0
      FRAME='CMS'
      PROJ ='P'    ! 'A'                                       ! Uzhi
      TARG ='P'    ! 'A'                                       ! Uzhi
      IAP  =1      ! 197                                       ! Uzhi
      IZP  =1      !  79                                       ! Uzhi
      IAT  =1      ! 197                                       ! Uzhi
      IZT  =1      !  79                                       ! Uzhi
      N_events=100                                             ! Uzhi

C   Simulation of PP-interactions at \sqrt{s}=200 GeV          ! Uzhi


      CALL HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
C			********Initialize HIJING

      WRITE(6,*)' Simulation of interactions with'             ! Uzhi
      WRITE(6,*)                                               ! Uzhi
      WRITE(6,*)' Proj = ',PROJ,' and Targ = ',TARG            ! Uzhi
      WRITE(6,*)' IAP  =',IAP  ,'            IAT  =',IAT       ! Uzhi
      WRITE(6,*)' IZP  =',IZP  ,'            IZT  =',IZT       ! Uzhi
      WRITE(6,*)                                               ! Uzhi
      WRITE(6,*)' Reference frame -   ',FRAME                  ! Uzhi
      WRITE(6,*)' ENERGY            ',EFRM,' GeV'              ! Uzhi
      WRITE(6,*)' Number of generated events -',N_events       ! Uzhi
      WRITE(6,*)                                               ! Uzhi

      BMIN=0.0
      BMAX=0.0
      DO 2000 I_event=1,N_events

        WRITE(6,*)' Event # ',I_event,' ------------------'    ! Uzhi

        CALL HIJING(FRAME,BMIN,BMAX)
C
        WRITE(6,*)' Multiplicity of produced particles - ',NATT! Uzhi
        write(6,*)                                             ! Uzhi
        WRITE(6,*)'   ID Charge  Mass (GeV)  Px           Py', ! Uzhi
     &            '           Pz    (GeV/c)'                   ! Uzhi
        WRITE(6,*)' ----------------------------------------'  ! Uzhi


        DO 1000 I=1,NATT
                                                               ! Uzhi
C         IF(LUCHGE(KATT(I,1)).NE.0) THEN                      ! Uzhi
C           this select charged particles only
C	                           !information of produced particles
C			           !is stored in common blocks HIMAIN1 and
C					HIMAIN2
C         ENDIF                                                ! Uzhi

          ICH=LUCHGE(KATT(I,1))/3                              ! Uzhi

          Amass=ULMASS(KATT(I,1))                              ! Uzhi

          WRITE(6,900) KATT(I,1),ICH,Amass,                    ! Uzhi
     &                 PATT(I,1),PATT(I,2),PATT(I,3)           ! Uzhi
 900      FORMAT(1X,I6,I4,2x,F7.3,3(2X,E11.4))                 ! Uzhi

1000    CONTINUE

        PAUSE                                                  ! Uzhi

2000	CONTINUE
	STOP
	END

      FUNCTION RAN(NSEED)                                      ! Uzhi
      RAN=RLU(NSEED)                                           ! Uzhi
      RETURN                                                   ! Uzhi
      END                                                      ! Uzhi

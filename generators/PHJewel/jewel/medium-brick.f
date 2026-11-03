C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++ medium-brick.f -- static "brick" medium for JEWEL 2.6.0          ++
C++                                                                 ++
C++ Local addition (ProjectJEWEL260, 2026-09-30). Same interface as ++
C++ medium-simple.f, but the medium is a static, uniform plasma at  ++
C++ temperature TBRICK that exists for lab times 0 < t < LBRICK fm  ++
C++ everywhere in space and then vanishes. A parton produced at the ++
C++ origin at t=0 and moving at c therefore sees exactly LBRICK fm  ++
C++ of medium whatever its direction. No longitudinal flow, no      ++
C++ expansion, no geometry.                                         ++
C++                                                                 ++
C++ Parameter file labels: TBRICK [GeV], LBRICK [fm], NF,           ++
C++ MDFACTOR, MDSCALEFAC.                                           ++
C++ getcentrality() returns TBRICK and getb() returns LBRICK so     ++
C++ that both land in the EVINFOFILE 'centr' and 'b' columns.       ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE MEDINIT(FILE,id,etam,mass)
      IMPLICIT NONE
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      common/rapmax2/etamax2
      double precision etamax2
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      common/logfile/logfid
      integer logfid
      INTEGER I,LUN,POS,IOS,id,mass
      double precision etam
      CHARACTER*512 BUFFER,LABEL
      CHARACTER*512 FILE
      character firstchar
      logical fileexist

      etamax2 = etam
      logfid = id
      IOS=0
      LUN=77
C--defaults
      TBRICK=0.35d0
      LBRICK=3.d0
      D3=0.9d0
      ZETA3=1.2d0
      NF=3
      MDFACTOR=0.45d0
      MDSCALEFAC=1.d0

      write(logfid,*)
      inquire(file=FILE,exist=fileexist)
      if(fileexist)then
        write(logfid,*)'Reading brick medium parameters from ',FILE
        OPEN(unit=LUN,file=FILE,status='old',err=10)
        do 20 i=1,1000
          READ(LUN, '(A)', iostat=ios) BUFFER
          if (ios.ne.0) goto 30
          firstchar = buffer(1:1)
          if (firstchar.eq.'#') goto 20
          POS=SCAN(BUFFER,' ')
          LABEL=BUFFER(1:POS)
          BUFFER=BUFFER(POS+1:)
          IF (LABEL=="TBRICK")THEN
            READ(BUFFER,*,IOSTAT=IOS) TBRICK
          ELSE IF (LABEL=="LBRICK") THEN
            READ(BUFFER,*,IOSTAT=IOS) LBRICK
          ELSE IF (LABEL=="NF") THEN
            READ(BUFFER,*,IOSTAT=IOS) NF
          ELSE IF (LABEL=="MDFACTOR") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDFACTOR
          ELSE IF (LABEL=="MDSCALEFAC") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDSCALEFAC
          else
            write(logfid,*)'unknown label ',label
          endif
 20     continue
 30     close(LUN,status='keep')
        write(logfid,*)'...done'
        goto 40
 10     write(logfid,*)'Could not open medium parameter file, '//
     &  'will run with default settings.'
      else
        write(logfid,*)'No medium parameter file found, '//
     &  'will run with default settings.'
      endif

 40   write(logfid,*)'BRICK medium, using parameters:'
      write(logfid,*)'TBRICK     =',TBRICK
      write(logfid,*)'LBRICK     =',LBRICK
      write(logfid,*)'NF         =',NF
      write(logfid,*)'MDFACTOR   =',MDFACTOR
      write(logfid,*)'MDSCALEFAC =',MDSCALEFAC
      write(logfid,*)
      END



      SUBROUTINE MEDNEXTEVT
      IMPLICIT NONE
      END



      double precision function getcentrality()
      implicit none
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      getcentrality=TBRICK
      end



      double precision function getb()
      implicit none
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      getb=LBRICK
      end



      SUBROUTINE PICKVTX(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      X=0.d0
      Y=0.d0
      END



      SUBROUTINE SETB(BVAL)
      IMPLICIT NONE
      DOUBLE PRECISION BVAL
      END



      SUBROUTINE GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS)
      IMPLICIT NONE
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      common/logfile/logfid
      integer logfid
      DOUBLE PRECISION GETTEMP,GETMD,GETMS
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E,MD,TEMP
      INTEGER TYPE
      DOUBLE PRECISION R,PYR,pmax,wt,theta,phi,pi,p,E2
      DATA PI/3.141592653589793d0/

      R=PYR(0)
      IF(R.LT.(2.*12.*NF*D3/3.)/(2.*12.*NF*D3/3.+3.*16.*ZETA3/2.))THEN
         TYPE=2
         MS=GETMS(X,Y,Z,T)
      ELSE
         TYPE=21
         MS=GETMD(X,Y,Z,T)
      ENDIF
      TEMP=GETTEMP(X,Y,Z,T)
      pmax = 10.*temp

      IF(TEMP.LT.1.D-2)THEN
       write(logfid,*)'asking for a scattering centre without medium:'
       write(logfid,*)'at (x,y,z,t)=',X,Y,Z,T
       write(logfid,*)'making one up to continue but '//
     & 'something is wrong!'
       TYPE=21
       PX=0.d0
       PY=0.d0
       PZ=0.d0
       MS=GETMS(0.d0,0.d0,0.d0,0.d0)
       MD=GETMD(0.d0,0.d0,0.d0,0.d0)
       E=SQRT(PX**2+PY**2+PZ**2+MS**2)
       RETURN
      ENDIF

 10   p = pyr(0)**0.3333333*pmax
      E2 = sqrt(p**2+ms**2)
      if (type.eq.2) then
        wt = (exp(ms/temp)-1.)/(exp(E2/temp)-1.)
      else
        wt = (exp(ms/temp)+1.)/(exp(E2/temp)+1.)
      endif
      if (wt.gt.1.) write(logfid,*)'Error in getscatterer: weight = ',wt
      if (wt.lt.0.) write(logfid,*)'Error in getscatterer: weight = ',wt
      if (pyr(0).gt.wt) goto 10
      phi = pyr(0)*2.*pi
      theta = -acos(2.*pyr(0)-1.)+pi
      px  = p*sin(theta)*cos(phi)
      py  = p*sin(theta)*sin(phi)
      pz  = p*cos(theta)
      E   = E2
      END



      SUBROUTINE AVSCATCEN(X,Y,Z,T,PX,PY,PZ,E,m)
      IMPLICIT NONE
      double precision x,y,z,t,px,py,pz,e,getms,m
      m  = getms(x,y,z,t)
      e  = m
      px = 0.d0
      py = 0.d0
      pz = 0.d0
      end



      SUBROUTINE maxscatcen(PX,PY,PZ,E,m)
      IMPLICIT NONE
      double precision px,py,pz,e,getmsmax,m
      m  = getmsmax()
      e  = m
      px = 0.d0
      py = 0.d0
      pz = 0.d0
      end



      DOUBLE PRECISION FUNCTION GETMD(X1,Y1,Z1,T1)
      IMPLICIT NONE
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION X1,Y1,Z1,T1,GETTEMP
      GETMD=MDSCALEFAC*3.*GETTEMP(X1,Y1,Z1,T1)
      GETMD=MAX(GETMD,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMS(X2,Y2,Z2,T2)
      IMPLICIT NONE
      DOUBLE PRECISION X2,Y2,Z2,T2,GETMD
      GETMS=GETMD(X2,Y2,Z2,T2)/SQRT(2.)
      END



      DOUBLE PRECISION FUNCTION GETNEFF(X3,Y3,Z3,T3)
      IMPLICIT NONE
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      DOUBLE PRECISION X3,Y3,Z3,T3,PI,GETTEMP
      DATA PI/3.141592653589793d0/
      GETNEFF=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *GETTEMP(X3,Y3,Z3,T3)**3/PI**2
      END



      DOUBLE PRECISION FUNCTION GETTEMP(X4,Y4,Z4,T4)
      IMPLICIT NONE
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      DOUBLE PRECISION X4,Y4,Z4,T4
      GETTEMP=0.D0
      IF(T4.LT.0.d0) RETURN
      IF(T4.GT.LBRICK) RETURN
      GETTEMP=TBRICK
      END



      DOUBLE PRECISION FUNCTION GETTEMPMAX()
      IMPLICIT NONE
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      GETTEMPMAX=TBRICK
      END



      DOUBLE PRECISION FUNCTION GETMDMAX()
      IMPLICIT NONE
      COMMON/MDFAC/MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
      GETMDMAX=MDSCALEFAC*3.*GETTEMPMAX()
      GETMDMAX=MAX(GETMDMAX,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMDMIN()
      IMPLICIT NONE
      DOUBLE PRECISION GETMDMAX
      GETMDMIN=GETMDMAX()
      END



      DOUBLE PRECISION FUNCTION GETMSMAX()
      IMPLICIT NONE
      DOUBLE PRECISION GETMDMAX
      GETMSMAX=GETMDMAX()/SQRT(2.D0)
      END



      DOUBLE PRECISION FUNCTION GETNATMDMIN()
      IMPLICIT NONE
      DOUBLE PRECISION GETNEFFMAX
      GETNATMDMIN=GETNEFFMAX()
      END



      DOUBLE PRECISION FUNCTION GETLTIMEMAX()
      IMPLICIT NONE
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      GETLTIMEMAX=LBRICK
      END



      DOUBLE PRECISION FUNCTION GETNEFFMAX()
      IMPLICIT NONE
      COMMON/BRICKPAR/TBRICK,LBRICK,D3,ZETA3,NF
      DOUBLE PRECISION TBRICK,LBRICK,D3,ZETA3
      INTEGER NF
      DOUBLE PRECISION PI
      DATA PI/3.141592653589793d0/
      GETNEFFMAX=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *TBRICK**3/PI**2
      END



      DOUBLE PRECISION FUNCTION MEDDERIV(XVAL,W)
      IMPLICIT NONE
      DOUBLE PRECISION XVAL
      INTEGER W
      MEDDERIV=0.D0
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++ Copyright (C) 2017 Korinna C. Zapp [Korinna.Zapp@cern.ch]       ++
C++                                                                 ++
C++ This file is part of JEWEL 2.2.0                                ++
C++                                                                 ++
C++ The JEWEL homepage is jewel.hepforge.org                        ++
C++                                                                 ++
C++ The medium model was partly implemented by Jochen Klein.        ++
C++ Raghav Kunnawalkam Elayavalli helped with the implementation    ++
C++ of the V+jet processes.                                         ++
C++                                                                 ++
C++ Please follow the MCnet GUIDELINES and cite Eur.Phys.J. C74     ++
C++ (2014) no.2, 2762 [arXiv:1311.0048] for the code and            ++
C++ JHEP 1303 (2013) 080 [arXiv:1212.1599] and                      ++
C++ optionally EPJC 60 (2009) 617 [arXiv:0804.3568] for the         ++
C++ physics. The reference for V+jet processes is EPJC 76 (2016)    ++
C++ no.12 695 [arXiv:1608.03099] and for recoil effects it is       ++
C++ arXiv:1707.01539.
C++                                                                 ++
C++ JEWEL relies heavily on PYTHIA 6 for the event generation. The  ++
C++ modified version of PYTHIA 6.4.25 that is distributed with      ++
C++ JEWEL is, however, not an official PYTHIA release and must not  ++
C++ be used for anything else. Please refer to results as           ++
C++ "JEWEL+PYTHIA".                                                 ++
C++                                                                 ++
C++ JEWEL also uses code provided by S. Zhang and J. M. Jing        ++
C++ (Computation of Special Functions, John Wiley & Sons, New York, ++
C++ 1996 and http://jin.ece.illinois.edu) for computing the         ++
C++ exponential integral Ei(x).                                     ++
C++                                                                 ++
C++                                                                 ++
C++ JEWEL  is free software; you can redistribute it and/or         ++
C++ modify it under the terms of the GNU General Public License     ++
C++ as published by the Free Software Foundation; either version 2  ++
C++ of the License, or (at your option) any later version.          ++
C++                                                                 ++
C++ JEWEL is distributed in the hope that it will be useful,        ++
C++ but WITHOUT ANY WARRANTY; without even the implied warranty of  ++
C++ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    ++
C++ GNU General Public License for more details.                    ++
C++                                                                 ++
C++ You should have received a copy of the GNU General Public       ++  
C++ License along with this program; if not, write to the Free      ++
C++ Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, ++
C++ MA 02110-1301 USA                                               ++
C++                                                                 ++
C++ Linking JEWEL statically or dynamically with other modules is   ++
C++ making a combined work based on JEWEL. Thus, the terms and      ++
C++ conditions of the GNU General Public License cover the whole    ++
C++ combination.                                                    ++
C++                                                                 ++
C++ In addition, as a special exception, I give you permission to   ++
C++ combine JEWEL with the code for the computation of special      ++
C++ functions provided by S. Zhang and J. M. Jing. You may copy and ++
C++ distribute such a system following the terms of the GNU GPL for ++
C++ JEWEL and the licenses of the other code concerned, provided    ++
C++ that you include the source code of that other code when and as ++
C++ the GNU GPL requires distribution of source code.               ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE MEDINIT
      IMPLICIT NONE
      COMMON/MEDPARAM/BMIN,BMAX,CENTRMIN,CENTRMAX,BREAL,RAU
      DOUBLE PRECISION BMIN,BMAX,CENTRMIN,CENTRMAX,BREAL,RAU
      COMMON/MEDIUM/MEDIUM
      LOGICAL MEDIUM
      DATA MEDIUM/.FALSE./
      RAU=5.d0
      END



      SUBROUTINE MEDNEXTEVT
      IMPLICIT NONE
      END



      double precision function getcentrality()
      implicit none
      getcentrality=-1.d0
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



      SUBROUTINE GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS,MD,TEMP)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E,MD,TEMP
      INTEGER TYPE
      WRITE(*,*)'GETSCATTERER called although in vacuum'
      END



      SUBROUTINE AVSCATCEN(X,Y,Z,T,PX,PY,PZ,E,MS)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E
      WRITE(*,*)'AVSCATCEN called although in vacuum'
      END



      SUBROUTINE MAXSCATCEN(PX,PY,PZ,E,MS)
      IMPLICIT NONE
      DOUBLE PRECISION MS,PX,PY,PZ,E
      ms=0.5d0
	e=0.5d0
	px=0.d0
	py=0.d0
	pz=0.d0
      END



      DOUBLE PRECISION FUNCTION GETNEFF(X,Y,Z,T)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,T
      GETNEFF=0.d0
      END



      DOUBLE PRECISION FUNCTION GETTEMP(X,Y,Z,T)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,T
      GETTEMP=0.d0
      END


      
      DOUBLE PRECISION FUNCTION GETTEMPMAX()
      IMPLICIT NONE
      GETTEMPMAX=0.d0
      END



      DOUBLE PRECISION FUNCTION GETMDMAX()
      IMPLICIT NONE
      GETMDMAX=0.d0
      END



      DOUBLE PRECISION FUNCTION GETMDMIN()
      IMPLICIT NONE
      GETMDMIN=0.d0
      END



      DOUBLE PRECISION FUNCTION GETMSMAX()
      IMPLICIT NONE
      GETMSMAX=0.d0
      END



      DOUBLE PRECISION FUNCTION GETNEFFMAX()
      IMPLICIT NONE
      GETNEFFMAX=0.d0
      END



	DOUBLE PRECISION FUNCTION GETNATMDMIN()
	IMPLICIT NONE
      GETNATMDMIN=0.d0
	END



      DOUBLE PRECISION FUNCTION GETLTIMEMAX()
      IMPLICIT NONE
	GETLTIMEMAX=0.
	END



      DOUBLE PRECISION FUNCTION GETMD(X,Y,Z,T)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,T,GETTEMP
      GETMD=0.D0
      END



      DOUBLE PRECISION FUNCTION GETMS(X,Y,Z,T)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,T,GETMD
      GETMS=0.D0
      END



      DOUBLE PRECISION FUNCTION MEDDERIV(XVAL,W)
      MEDDERIV=0.D0
      END

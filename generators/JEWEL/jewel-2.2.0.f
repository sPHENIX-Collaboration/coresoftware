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

      PROGRAM JEWEL
	IMPLICIT NONE
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	INTEGER MSTU,MSTJ
	DOUBLE PRECISION PARU,PARJ
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	INTEGER MDCY,MDME,KFDP
	DOUBLE PRECISION BRAT
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
	INTEGER MSEL,MSELPD,MSUB,KFIN
	DOUBLE PRECISION CKIN 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
      COMMON/PYDATR/MRPY(6),RRPY(100)
	INTEGER MRPY
	DOUBLE PRECISION RRPY
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--use nuclear pdf?      
      COMMON/NPDF/MASS,NSET,EPS09,INITSTR
      INTEGER NSET
      DOUBLE PRECISION MASS
      LOGICAL EPS09
      CHARACTER*10 INITSTR
C--number of protons
	common/np/nproton
	integer nproton
C--organisation of event record
	common/evrecord/nsim,npart,offset,hadrotype,sqrts,collider,hadro,
     &shorthepmc,channel,isochannel
	integer nsim,npart,offset,hadrotype
	double precision sqrts
	character*4 collider,channel
	character*2 isochannel
	logical hadro,shorthepmc
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
C--number of scattering events
	COMMON/CHECK/NSCAT,NSCATEFF,NSPLIT
	DOUBLE PRECISION NSCAT,NSCATEFF,NSPLIT
C--number of extrapolations in tables
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
C--local variables
	integer j,i,kk,poissonian
      integer nsimpp,nsimpn,nsimnp,nsimnn,nsimsum,nsimchn
	double precision sumofweightstot,wdisctot,scalefac
	double precision gettemp,r,tau
	character*2 b1,b2

	call init()

	SUMOFWEIGHTSTOT=0.d0
      WDISCTOT=0.d0

C--e+ + e- event generation
	if (collider.eq.'EEJJ') then
	  b1 = 'e+'
	  b2 = 'e-'
	  write(logfid,*)
	  write(logfid,*)
     &'####################################################'
	  write(logfid,*)
	  write(logfid,*)'generating ',nsim,' events in ',b1,' + ',b2,
     &' channel'
	  write(logfid,*)
	  write(logfid,*)
     &'####################################################'
	  write(logfid,*)
	  SUMOFWEIGHTS=0.d0
        WDISC=0.d0
	  call initpythia(b1,b2)
	    write(logfid,*)
C--e+ + e- event loop
	  DO 100 J=1,NSIM
	    call genevent(j,b1,b2)
 100	  CONTINUE
	  sumofweightstot = sumofweightstot+sumofweights
	  wdisctot = wdisctot + wdisc
	  write(logfid,*)
	  write(logfid,*)'cross section in e+ + e- channel:',PARI(1),'mb'
	  write(logfid,*)'sum of event weights in e+ + e- channel:',
     &	sumofweights-wdisc
	  write(logfid,*)

	else
C--hadronic event generation
	  if (isochannel.eq.'PP') then
	    nsimpp = nsim
	    nsimpn = 0
	    nsimnp = 0
	    nsimnn = 0
	  elseif (isochannel.eq.'PN') then
	    nsimpp = 0
	    nsimpn = nsim
	    nsimnp = 0
	    nsimnn = 0
	  elseif (isochannel.eq.'NP') then
	    nsimpp = 0
	    nsimpn = 0
	    nsimnp = nsim
	    nsimnn = 0
	  elseif (isochannel.eq.'NN') then
	    nsimpp = 0
	    nsimpn = 0
	    nsimnp = 0
	    nsimnn = nsim
	  else
	    nsimpp = poissonian(nsim*nproton**2/mass**2)
	    nsimpn = poissonian(nsim*nproton*(mass-nproton*1.d0)/mass**2)
	    nsimnp = poissonian(nsim*nproton*(mass-nproton*1.d0)/mass**2)
	    nsimnn = poissonian(nsim*(mass-nproton*1.d0)**2/mass**2)
	    nsimsum = nsimpp + nsimpn + nsimnp + nsimnn
	    scalefac = nsim*1.d0/(nsimsum*1.d0)
	    nsimpp = int(nsimpp*scalefac)
	    nsimpn = int(nsimpn*scalefac)
	    nsimnp = int(nsimnp*scalefac)
	    nsimnn = int(nsimnn*scalefac)
	    nsimsum = nsimpp + nsimpn + nsimnp + nsimnn
	  endif
C--loop over channels
	  do 101 kk=1,4
	    if (kk.eq.1) then
	      b1 = 'p+'
	      b2 = 'p+'
	      nsimchn = nsimpp
	    elseif (kk.eq.2) then
	      b1 = 'p+'
	      b2 = 'n0'
	      nsimchn = nsimpn
	    elseif (kk.eq.3) then
	      b1 = 'n0'
	      b2 = 'p+'
	      nsimchn = nsimnp
	    else
	      b1 = 'n0'
	      b2 = 'n0'
	      nsimchn = nsimnn
	    endif
	    write(logfid,*)
	    write(logfid,*)
     &'####################################################'
	    write(logfid,*)
	    write(logfid,*)'generating ',nsimchn,' events in ',
     &b1,' + ',b2,' channel'
	    write(logfid,*)
	    write(logfid,*)
     &'####################################################'
	    write(logfid,*)
	    SUMOFWEIGHTS=0.d0
          WDISC=0.d0
	    call initpythia(b1,b2)
	    write(logfid,*)
C--event loop
	    DO 102 J=1,nsimchn
	      call genevent(j,b1,b2)
 102	    CONTINUE
	    sumofweightstot = sumofweightstot+sumofweights
	    wdisctot = wdisctot + wdisc
	    write(logfid,*)
	    write(logfid,*)'cross section in ',b1,' + ',b2,' channel:',
     &	PARI(1),'mb'
	    write(logfid,*)'sum of event weights in ',b1,' + ',b2,
     &	' channel:',sumofweights-wdisc
	    write(logfid,*)
 101	  continue
	endif
 
C--finish
	WRITE(HPMCFID,'(A)')'HepMC::IO_GenEvent-END_EVENT_LISTING'
	WRITE(HPMCFID,*)
	CLOSE(HPMCFID,status='keep')

	write(logfid,*)
	write(logfid,*)'mean number of scatterings:',
     &      NSCAT/(SUMOFWEIGHTSTOT-WDISCTOT)
	write(logfid,*)'mean number of effective scatterings:',
     &      NSCATEFF/(SUMOFWEIGHTSTOT-WDISCTOT)
	write(logfid,*)'mean number of splittings:',
     &      NSPLIT/(SUMOFWEIGHTSTOT-WDISCTOT)
	write(logfid,*)
	write(logfid,*)'number of extrapolations in splitting integral: ',
     &	noverspliti,' (',(noverspliti*1.d0)/(ntotspliti*1.d0),'%)'
	write(logfid,*)
     &	'number of extrapolations in splitting partonic PDFs: ',
     &	noverpdf,' (',(noverpdf*1.d0)/(ntotpdf*1.d0),'%)'
	write(logfid,*)
     &	'number of extrapolations in splitting cross sections: ',
     &	noverxsec,' (',(noverxsec*1.d0)/(ntotxsec*1.d0),'%)'
	write(logfid,*)
     &	'number of extrapolations in Sudakov form factor: ',
     &	noversuda,' (',(noversuda*1.d0)/(ntotsuda*1.d0),'%)'
	write(logfid,*)
	write(logfid,*)'number of good events: ',ngood
	write(logfid,*)'total number of discarded events: ',NDISC
	write(logfid,*)'number of events for which conversion '//
     &'to hepmc failed: ',NSTRANGE
	call printtime

	close(logfid,status='keep')

	END



***********************************************************************
***********************************************************************
***   END OF MAIN PROGRAM - NOW COME THE SUBROUTINES   ****************
***********************************************************************
***********************************************************************


***********************************************************************
***	  subroutine init
***********************************************************************
	subroutine init()
	implicit none
	INTEGER PYCOMP
	INTEGER NMXHEP
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	INTEGER MSTU,MSTJ
	DOUBLE PRECISION PARU,PARJ
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	INTEGER MDCY,MDME,KFDP
	DOUBLE PRECISION BRAT
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
	INTEGER MSEL,MSELPD,MSUB,KFIN
	DOUBLE PRECISION CKIN 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
      COMMON/PYDATR/MRPY(6),RRPY(100)
	INTEGER MRPY
	DOUBLE PRECISION RRPY
C--use nuclear pdf?      
      COMMON/NPDF/MASS,NSET,EPS09,INITSTR
      INTEGER NSET
      DOUBLE PRECISION MASS
      LOGICAL EPS09
      CHARACTER*10 INITSTR
C--pdfset
	common/pdf/pdfset
	integer pdfset
C--number of protons
	common/np/nproton
	integer nproton
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--splitting integral
      COMMON/SPLITINT/SPLITIGGV(1000,1000),SPLITIQQV(1000,1000),
     &SPLITIQGV(1000,1000),QVAL(1000),ZMVAL(1000),QMAX,ZMMIN,NPOINT
      INTEGER NPOINT
      DOUBLE PRECISION SPLITIGGV,SPLITIQQV,SPLITIQGV,
     &QVAL,ZMVAL,QMAX,ZMMIN
C--pdf common block
	COMMON/PDFS/QINQX(2,1000),GINQX(2,1000),QINGX(2,1000),
     &GINGX(2,1000)
	DOUBLE PRECISION QINQX,GINQX,QINGX,GINGX
C--cross secttion common block
	COMMON/XSECS/INTQ1(1001,101),INTQ2(1001,101),
     &INTG1(1001,101),INTG2(1001,101)
	DOUBLE PRECISION INTQ1,INTQ2,INTG1,INTG2
C--Sudakov common block
	COMMON/INSUDA/SUDAQQ(1000,2),SUDAQG(1000,2),SUDAGG(1000,2)
     &,SUDAGC(1000,2)
	DOUBLE PRECISION SUDAQQ,SUDAQG,SUDAGG,SUDAGC
C--exponential integral for negative arguments
      COMMON/EXPINT/EIX(3,1000),VALMAX,NVAL
      INTEGER NVAL
      DOUBLE PRECISION EIX,VALMAX
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--factor in front of formation times
	COMMON/FTIMEFAC/FTFAC
	DOUBLE PRECISION FTFAC
C--factor in front of alphas argument
	COMMON/ALPHASFAC/PTFAC
	DOUBLE PRECISION PTFAC
C--number of scattering events
	COMMON/CHECK/NSCAT,NSCATEFF,NSPLIT
	DOUBLE PRECISION NSCAT,NSCATEFF,NSPLIT
C--number of extrapolations in tables
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
C--event weight exponent
	COMMON/WEXPO/WEIGHTEX
	DOUBLE PRECISION WEIGHTEX
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--max rapidity
	common/rapmax/etamax
	double precision etamax
C--memory for error message from getdeltat
	common/errline/errl
	integer errl
C--organisation of event record
	common/evrecord/nsim,npart,offset,hadrotype,sqrts,collider,hadro,
     &shorthepmc,channel,isochannel
	integer nsim,npart,offset,hadrotype
	double precision sqrts
	character*4 collider,channel
	character*2 isochannel
	logical hadro,shorthepmc
C--extra storage for scattering centres before interactions
      common/storescatcen/nscatcen,maxnscatcen,scatflav(10000),
     &scatcen(10000,5),writescatcen,writedummies
	integer nscatcen,maxnscatcen,scatflav
	double precision scatcen
	logical writescatcen,writedummies
C--Pythia parameters
	common/pythiaparams/PTMIN,PTMAX,weighted
	double precision PTMIN,PTMAX
	LOGICAL WEIGHTED

C--Variables local to this program
	INTEGER NJOB,ios,pos,i,j,jj,intmass
	DOUBLE PRECISION GETLTIMEMAX,EOVEST,r,pyr
	character firstchar
	CHARACTER*2 SNSET
      CHARACTER*80 PDFFILE,XSECFILE,FILEMED,FILESPLIT,buffer,
     &label,value
      CHARACTER*100 HEPMCFILE,LOGFILE,FILENAME2
	CHARACTER(LEN=100) filename
	LOGICAL PDFEXIST,SPLITIEXIST,XSECEXIST

	data maxnscatcen/10000/

      HPMCFID = 4
	logfid = 3

C--default settings
	nsim = 10
	njob = 0
	logfile = 'out.log'
	hepmcfile = 'out.hepmc'
	filesplit = 'splitint.dat'
	pdffile = 'pdfs.dat'
	xsecfile = 'xsecs.dat'
	filemed = 'medium-params.dat'
	nf = 3
	lqcd = 0.4
	q0 = 1.5
	ptmin = 5.
	ptmax = 350.
	etamax = 3.1
	collider = 'PPJJ'
	isochannel = 'XX'
	channel = 'MUON'
	sqrts = 2760
	pdfset = 10042
	nset = 1
	mass = 208.
      nproton = 82
	weighted = .true.
	weightex = 5.
	angord = .true.
	allhad = .false.
	hadro = .true.
	hadrotype = 0
	shorthepmc = .true.
	compress = .true.
	writescatcen = .false.
	writedummies = .false.
	
	lps = lqcd
	scatrecoil = .false.
	if (.not.hadro) shorthepmc = .true.

	SCALEFACM=1.
	ptfac=1.
	ftfac=1.d0

	if (iargc().eq.0) then
	  write(*,*)'No parameter file given, '// 
     &'will run with default settings.'
	else
	  call getarg(1,filename)
	  write(*,*)'Reading parameters from ',filename
	  open(unit=1,file=filename,status='old',err=110)
	  do 120 i=1,1000
          read(1, '(A)', iostat=ios) buffer
	    if(ios.ne.0) goto 130
	    firstchar = buffer(1:1)
	    if (firstchar.eq.'#') goto 120
          pos=scan(buffer,' ')
          label=buffer(1:pos)
          value=buffer(pos+1:)
          if(label.eq."NEVENT")then
            read(value,*,iostat=ios) nsim
          elseif(label.eq."NJOB")then
            read(value,*,iostat=ios) njob
          elseif(label.eq."LOGFILE")then
            read(value,'(a)',iostat=ios) logfile
          elseif(label.eq."HEPMCFILE")then
            read(value,'(a)',iostat=ios) hepmcfile
          elseif(label.eq."SPLITINTFILE")then
            read(value,'(a)',iostat=ios) filesplit
          elseif(label.eq."PDFFILE")then
            read(value,'(a)',iostat=ios) pdffile
          elseif(label.eq."XSECFILE")then
            read(value,'(a)',iostat=ios) xsecfile
          elseif(label.eq."MEDIUMPARAMS")then
            read(value,'(a)',iostat=ios) filemed
          elseif(label.eq."NF")then
            read(value,*,iostat=ios) nf
          elseif(label.eq."LAMBDAQCD")then
            read(value,*,iostat=ios) lqcd
          elseif(label.eq."Q0")then
            read(value,*,iostat=ios) q0
          elseif(label.eq."PTMIN")then
            read(value,*,iostat=ios) ptmin
          elseif(label.eq."PTMAX")then
            read(value,*,iostat=ios) ptmax
          elseif(label.eq."ETAMAX")then
            read(value,*,iostat=ios) etamax
          elseif(label.eq."PROCESS")then
            read(value,*,iostat=ios) collider
          elseif(label.eq."ISOCHANNEL")then
            read(value,*,iostat=ios) isochannel
	    elseif(label.eq."CHANNEL")then
	    read(value,*,iostat=ios) channel
          elseif(label.eq."SQRTS")then
            read(value,*,iostat=ios) sqrts
          elseif(label.eq."PDFSET")then
            read(value,*,iostat=ios) pdfset
          elseif(label.eq."NSET")then
            read(value,*,iostat=ios) nset
          elseif(label.eq."MASS")then
            read(value,*,iostat=ios) mass
          elseif(label.eq."NPROTON")then
            read(value,*,iostat=ios) nproton
          elseif(label.eq."WEIGHTED")then
            read(value,*,iostat=ios) weighted
          elseif(label.eq."WEXPO")then
            read(value,*,iostat=ios) weightex
          elseif(label.eq."ANGORD")then
            read(value,*,iostat=ios) angord
          elseif(label.eq."KEEPRECOILS")then
            read(value,*,iostat=ios) allhad
          elseif(label.eq."HADRO")then
            read(value,*,iostat=ios) hadro
          elseif(label.eq."HADROTYPE")then
            read(value,*,iostat=ios) hadrotype
          elseif(label.eq."SHORTHEPMC")then
            read(value,*,iostat=ios) shorthepmc
          elseif(label.eq."COMPRESS")then
            read(value,*,iostat=ios) compress
          elseif(label.eq."WRITESCATCEN")then
            read(value,*,iostat=ios) writescatcen
          elseif(label.eq."WRITEDUMMIES")then
            read(value,*,iostat=ios) writedummies
	    else
	      write(*,*)'unknown label ',label
	    endif
 120	  continue


 110	  write(*,*)
     &		'Unable to open parameter file, will exit the run.'
	  call exit(1)

 130	  close(1,status='keep')
	  write(*,*)'...done'
	endif

	if (ptmin.lt.3.d0) ptmin = 3.d0
	if (.not.writescatcen) writedummies = .false.

	OPEN(unit=logfid,file=LOGFILE,status='unknown')
	MSTU(11)=logfid

	call printtime
	call printlogo(logfid)


	write(logfid,*)
	write(logfid,*)'parameters of the run:'
	write(logfid,*)'NEVENT       = ',nsim
	write(logfid,*)'NJOB         = ',njob
	write(logfid,*)'LOGFILE      = ',logfile
	write(logfid,*)'HEPMCFILE    = ',hepmcfile
	write(logfid,*)'SPLITINTFILE = ',filesplit
	write(logfid,*)'PDFFILE      = ',pdffile
	write(logfid,*)'XSECFILE     = ',xsecfile
	write(logfid,*)'MEDIUMPARAMS = ',filemed
	write(logfid,*)'NF           = ',nf
	write(logfid,*)'LAMBDAQCD    = ',lqcd
	write(logfid,*)'Q0           = ',q0
	write(logfid,*)'PTMIN        = ',ptmin
	write(logfid,*)'PTMAX        = ',ptmax
	write(logfid,*)'ETAMAX       = ',etamax
	write(logfid,*)'PROCESS      = ',collider
	write(logfid,*)'ISOCHANNEL   = ',isochannel
	write(logfid,*)'CHANNEL      = ',channel
	write(logfid,*)'SQRTS        = ',sqrts
	write(logfid,*)'PDFSET       = ',pdfset
	write(logfid,*)'NSET         = ',nset
	write(logfid,*)'MASS         = ',mass
	write(logfid,*)'NPROTON      = ',nproton
	write(logfid,*)'WEIGHTED     = ',weighted
	write(logfid,*)'WEXPO        = ',weightex
	write(logfid,*)'ANGORD       = ',angord
	write(logfid,*)'KEEPRECOILS  = ',allhad
	write(logfid,*)'HADRO        = ',hadro
	write(logfid,*)'HADROTYPE    = ',hadrotype
	write(logfid,*)'SHORTHEPMC   = ',shorthepmc
	write(logfid,*)'COMPRESS     = ',compress
	write(logfid,*)'WRITESCATCEN = ',writescatcen
	write(logfid,*)'WRITEDUMMIES = ',writedummies
	write(logfid,*)
	call flush(logfid)

	if ((collider.ne.'PPJJ').and.(collider.ne.'EEJJ')
     &	.and.(collider.ne.'PPYJ').and.(collider.ne.'PPYQ')
     &	.and.(collider.ne.'PPYG')
     &	.and.(collider.ne.'PPZJ').and.(collider.ne.'PPZQ')
     &	.and.(collider.ne.'PPZG').and.(collider.ne.'PPWJ')
     &	.and.(collider.ne.'PPWQ').and.(collider.ne.'PPWG')
     &      .and.(collider.ne.'PPDY')) then
	  write(logfid,*)'Fatal error: colliding system unknown, '//
     &	'will exit now'
	  call exit(1)
	endif

C--initialize medium
	intmass = int(mass)
      CALL MEDINIT(FILEMED,logfid,etamax,intmass)
      CALL MEDNEXTEVT

	OPEN(unit=HPMCFID,file=HEPMCFILE,status='unknown')
	WRITE(HPMCFID,*)
	WRITE(HPMCFID,'(A)')'HepMC::Version 2.06.05'
	WRITE(HPMCFID,'(A)')'HepMC::IO_GenEvent-START_EVENT_LISTING'

	NPART=2
	
	if(ptmax.gt.0.)then
	  EOVEST=MIN(1.5*(PTMAX+50.)*COSH(ETAMAX),sqrts/2.)
	else
	  EOVEST=sqrts/2.
	endif

  
	CALL EIXINT
	CALL INSUDAINT(EOVEST)

	write(logfid,*)
	 INQUIRE(file=FILESPLIT,exist=SPLITIEXIST)
	 IF(SPLITIEXIST)THEN
	  write(logfid,*)'read splitting integrals from ',FILESPLIT
	  OPEN(unit=10,file=FILESPLIT,status='old')
	  READ(10,*)QMAX,ZMMIN,NPOINT
	  DO 893 I=1,NPOINT+1
	   READ(10,*) QVAL(I),ZMVAL(I)
 893    CONTINUE	 
	  DO 891 I=1,NPOINT+1
	   DO 892 J=1,NPOINT+1
	    READ(10,*)SPLITIGGV(I,J),SPLITIQQV(I,J),SPLITIQGV(I,J)
 892	   CONTINUE
 891	  CONTINUE
	  CLOSE(10,status='keep')
	 ELSE
 	  write(logfid,*)'have to integrate splitting functions, '// 
     &'this may take some time'
	  CALL SPLITFNCINT(EOVEST)
	  INQUIRE(file=FILESPLIT,exist=SPLITIEXIST)
	  IF(.NOT.SPLITIEXIST)THEN
 	   write(logfid,*)'write splitting integrals to ',FILESPLIT
	   OPEN(unit=10,file=FILESPLIT,status='new')
	   WRITE(10,*)QMAX,ZMMIN,NPOINT
	   DO 896 I=1,NPOINT+1
	    WRITE(10,*) QVAL(I),ZMVAL(I)
 896     CONTINUE	 
	   DO 897 I=1,NPOINT+1
	    DO 898 J=1,NPOINT+1
	     WRITE(10,*)SPLITIGGV(I,J),SPLITIQQV(I,J),SPLITIQGV(I,J)
 898	    CONTINUE
 897	   CONTINUE
	   CLOSE(10,status='keep')
	  ENDIF 
	 ENDIF
	write(logfid,*)

	INQUIRE(file=PDFFILE,exist=PDFEXIST)
	IF(PDFEXIST)THEN
	write(logfid,*)'read pdfs from ',PDFFILE
	 OPEN(unit=10,file=PDFFILE,status='old')
	 DO 872 I=1,2
	  DO 873 J=1,1000
	   READ(10,*)QINQX(I,J),GINQX(I,J),QINGX(I,J),GINGX(I,J)
 873	  CONTINUE
 872	 CONTINUE
	 CLOSE(10,status='keep')
	ELSE
 	 write(logfid,*)'have to integrate pdfs, this may take some time'
	 CALL PDFINT(EOVEST)
	 INQUIRE(file=PDFFILE,exist=PDFEXIST)
	 IF(.NOT.PDFEXIST)THEN
 	  write(logfid,*)'write pdfs to ',PDFFILE
	  OPEN(unit=10,file=PDFFILE,status='new')
	  DO 876 I=1,2
	   DO 877 J=1,1000
	    WRITE(10,*)QINQX(I,J),GINQX(I,J),QINGX(I,J),GINGX(I,J)
 877	   CONTINUE
 876	  CONTINUE
	  CLOSE(10,status='keep')
	 ENDIF
	ENDIF 
	write(logfid,*)

	INQUIRE(file=XSECFILE,exist=XSECEXIST)
	IF(XSECEXIST)THEN
	write(logfid,*)'read cross sections from ',XSECFILE
	 OPEN(unit=10,file=XSECFILE,status='old')
	  DO 881 J=1,1001
         DO 885 JJ=1,101
	   READ(10,*)INTQ1(J,JJ),INTQ2(J,JJ),
     &INTG1(J,JJ),INTG2(J,JJ)
 885     CONTINUE
 881	  CONTINUE
	 CLOSE(10,status='keep')
	ELSE
	 write(logfid,*)'have to integrate cross sections, '//
     &'this may take some time'
	 CALL XSECINT(EOVEST)
	 INQUIRE(file=XSECFILE,exist=XSECEXIST)
	 IF(.NOT.XSECEXIST)THEN
	  write(logfid,*)'write cross sections to ',XSECFILE
	  OPEN(unit=10,file=XSECFILE,status='new')
	   DO 883 J=1,1001
          DO 884 JJ=1,101
	    WRITE(10,*)INTQ1(J,JJ),INTQ2(J,JJ),
     &INTG1(J,JJ),INTG2(J,JJ)
 884      CONTINUE
 883	   CONTINUE
	  CLOSE(10,status='keep')
	 ENDIF 
	ENDIF
	write(logfid,*)
	CALL FLUSH(3)



C--initialise random number generator status
      IF(NJOB.GT.0)THEN
       MRPY(1)=NJOB*1000
       MRPY(2)=0
      ENDIF

C--Call PYR once for initialization
	R=PYR(0)

	NDISC=0
      NGOOD=0
      NSTRANGE=0
      
	ERRCOUNT=0
	errl = 0

	NSCAT=0.d0
	NSCATEFF=0.d0
	NSPLIT=0.d0

	ntotspliti=0
	noverspliti=0
	ntotpdf=0
	noverpdf=0
	ntotxsec=0
	noverxsec=0
	ntotsuda=0
	noversuda=0

	IF(NSET.EQ.0)THEN
	 EPS09=.FALSE.
	ELSE
	 EPS09=.TRUE.
	 IF(NSET.LT.10)THEN
	  WRITE(SNSET,'(i1)') NSET
	 ELSE
	  WRITE(SNSET,'(i2)') NSET
	 ENDIF
	  INITSTR='EPS09LO,'//SNSET
	ENDIF 

	end



***********************************************************************
***	  subroutine initpythia
***********************************************************************
	subroutine initpythia(beam1,beam2)
	implicit none
	INTEGER PYCOMP
	INTEGER NMXHEP
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	INTEGER MSTU,MSTJ
	DOUBLE PRECISION PARU,PARJ
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	INTEGER MDCY,MDME,KFDP
	DOUBLE PRECISION BRAT
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
	INTEGER MSEL,MSELPD,MSUB,KFIN
	DOUBLE PRECISION CKIN 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
      COMMON/PYDATR/MRPY(6),RRPY(100)
	INTEGER MRPY
	DOUBLE PRECISION RRPY
C--use nuclear pdf?      
      COMMON/NPDF/MASS,NSET,EPS09,INITSTR
      INTEGER NSET
      DOUBLE PRECISION MASS
      LOGICAL EPS09
      CHARACTER*10 INITSTR
C--pdfset
	common/pdf/pdfset
	integer pdfset
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
C--event weight exponent
	COMMON/WEXPO/WEIGHTEX
	DOUBLE PRECISION WEIGHTEX
C--memory for error message from getdeltat
	common/errline/errl
	integer errl
C--organisation of event record
	common/evrecord/nsim,npart,offset,hadrotype,sqrts,collider,hadro,
     &shorthepmc,channel,isochannel
	integer nsim,npart,offset,hadrotype
	double precision sqrts
	character*4 collider,channel
	character*2 isochannel
	logical hadro,shorthepmc
C--Pythia parameters
	common/pythiaparams/PTMIN,PTMAX,weighted
	double precision PTMIN,PTMAX
	LOGICAL WEIGHTED

C--Variables local to this program
	character*2 beam1,beam2


C--initialise PYTHIA
C--no multiple interactions
	 MSTP(81) = 0
C--initial state radiation
	 MSTP(61)=1
C--switch off final state radiation
	 MSTP(71)=0
C--No hadronisation (yet)
       MSTP(111)=0
C--parameter affecting treatment of string corners
       PARU(14)=1.
C--Min shat in simulation
       CKIN(1)=2.      
C--pT-cut
       CKIN(3)=PTMIN
       CKIN(4)=PTMAX
C--use LHAPDF
	 MSTP(52)=2
C--choose pdf: CTEQ6ll (LO fit/LO alphas) - 10042
C	         MSTW2008 (LO central) - 21000
	 MSTP(51)=PDFSET
	 IF(COLLIDER.EQ.'PPYQ')THEN
	  MSEL=0
	  MSUB(29)=1
	 ELSEIF(COLLIDER.EQ.'PPYG')THEN
	  MSEL=0
	  MSUB(14)=1
	  MSUB(115)=1
	 ELSEIF(COLLIDER.EQ.'PPYJ')THEN
	  MSEL=0
	  MSUB(14)=1
	  MSUB(29)=1
	  MSUB(115)=1
	 ELSEIF((COLLIDER.EQ.'PPZJ').or.(COLLIDER.EQ.'PPZQ')
     &	.or.(COLLIDER.EQ.'PPZG')
     &      .or.(collider.eq.'PPDY'))THEN
	  MSEL=0
	  IF((COLLIDER.EQ.'PPZJ').or.(COLLIDER.EQ.'PPZQ')) MSUB(30)=1
	  IF((COLLIDER.EQ.'PPZJ').or.(COLLIDER.EQ.'PPZG')) MSUB(15)=1
	  IF(COLLIDER.EQ.'PPDY') MSUB(1)=1
	  MDME(174,1)=0          !Z decay into d dbar', 
	  MDME(175,1)=0          !Z decay into u ubar', 
	  MDME(176,1)=0          !Z decay into s sbar', 
	  MDME(177,1)=0          !Z decay into c cbar', 
	  MDME(178,1)=0          !Z decay into b bbar', 
	  MDME(179,1)=0          !Z decay into t tbar', 
	  MDME(182,1)=0          !Z decay into e- e+', 
	  MDME(183,1)=0          !Z decay into nu_e nu_ebar', 
	  MDME(184,1)=0          !Z decay into mu- mu+', 
	  MDME(185,1)=0          !Z decay into nu_mu nu_mubar', 
	  MDME(186,1)=0          !Z decay into tau- tau+', 
	  MDME(187,1)=0          !Z decay into nu_tau nu_taubar',
	  if (channel.EQ.'ELEC')THEN
	    MDME(182,1)=1
	  ELSEIF(channel.EQ.'MUON')THEN
	    MDME(184,1)=1
	  ENDIF
	 ELSEIF((COLLIDER.EQ.'PPWJ').or.(COLLIDER.EQ.'PPWQ')
     &	.or.(COLLIDER.EQ.'PPWG'))THEN
	  MSEL=0
	  IF((COLLIDER.EQ.'PPWJ').or.(COLLIDER.EQ.'PPWQ')) MSUB(31)=1
	  IF((COLLIDER.EQ.'PPWJ').or.(COLLIDER.EQ.'PPWG')) MSUB(16)=1
	  MDME(190,1)=0          ! W+ decay into dbar u,
	  MDME(191,1)=0          ! W+ decay into dbar c,
	  MDME(192,1)=0          ! W+ decay into dbar t,
	  MDME(194,1)=0          ! W+ decay into sbar u,
	  MDME(195,1)=0          ! W+ decay into sbar c,
	  MDME(196,1)=0          ! W+ decay into sbar t,
	  MDME(198,1)=0          ! W+ decay into bbar u,
	  MDME(199,1)=0          ! W+ decay into bbar c,
	  MDME(200,1)=0          ! W+ decay into bbar t,
	  MDME(202,1)=0          ! W+ decay into b'bar u,
	  MDME(203,1)=0          ! W+ decay into b'bar c,
	  MDME(204,1)=0          ! W+ decay into b'bar t,
	  MDME(206,1)=0          ! W+ decay into e+ nu_e,
	  MDME(207,1)=0          ! W+ decay into mu+ nu_mu,
	  MDME(208,1)=0          ! W+ decay into tau+ nu_tau,
	  MDME(209,1)=0      ! W+ decay into tau'+ nu'_tau,
	  if (channel.EQ.'ELEC')THEN
	   MDME(206,1)=1
	  ELSEIF(channel.EQ.'MUON')THEN
	   MDME(207,1)=1
	  ENDIF
	 ELSE
C--All QCD processes are active
        MSEL=1
	 ENDIF
!	 MSEL=0
!	 MSUB(11)=1
!	 MSUB(12)=1
!	 MSUB(53)=1
!	 MSUB(13)=1
!	 MSUB(68)=1
!	 MSUB(28)=1

C--weighted events
       IF(WEIGHTED) MSTP(142)=1

C--number of errors to be printed
	 MSTU(22)=MAX(10,INT(5.*NSIM/100.))

C--number of lines in event record
	MSTU(4)=23000
	MSTU(5)=23000

C--switch off pi0 decay
      MDCY(PYCOMP(111),1)=0
C--initialisation call
	 IF(COLLIDER.EQ.'EEJJ')THEN
	  OFFSET=9
	  CALL PYINIT('CMS',beam1,beam2,sqrts)
	 ELSEIF((COLLIDER.EQ.'PPJJ').OR.(COLLIDER.EQ.'PPYJ').OR.
     & 		(COLLIDER.EQ.'PPYG').OR.(COLLIDER.EQ.'PPYQ'))THEN
	  OFFSET=8
	  CALL PYINIT('CMS',beam1,beam2,sqrts)
	 ELSEIF((COLLIDER.EQ.'PPWJ').OR.(COLLIDER.EQ.'PPZJ').or.
     &	(COLLIDER.EQ.'PPWQ').OR.(COLLIDER.EQ.'PPZQ').or.
     &	(COLLIDER.EQ.'PPWG').OR.(COLLIDER.EQ.'PPZG'))THEN
	  OFFSET=10
	  CALL PYINIT('CMS',beam1,beam2,sqrts)
	 elseif (collider.eq.'PPDY') then
	  CALL PYINIT('CMS',beam1,beam2,sqrts)
	 ENDIF

	end



***********************************************************************
***	  subroutine genevent
***********************************************************************
	subroutine genevent(j,b1,b2)
	implicit none
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
	INTEGER PYCOMP
	INTEGER NMXHEP
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	INTEGER MSTU,MSTJ
	DOUBLE PRECISION PARU,PARJ
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	INTEGER MDCY,MDME,KFDP
	DOUBLE PRECISION BRAT
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
	INTEGER MSEL,MSELPD,MSUB,KFIN
	DOUBLE PRECISION CKIN 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
      COMMON/PYDATR/MRPY(6),RRPY(100)
	INTEGER MRPY
	DOUBLE PRECISION RRPY
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--factor in front of formation times
	COMMON/FTIMEFAC/FTFAC
	DOUBLE PRECISION FTFAC
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
C--number of scattering events
	COMMON/CHECK/NSCAT,NSCATEFF,NSPLIT
	DOUBLE PRECISION NSCAT,NSCATEFF,NSPLIT
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
C--event weight exponent
	COMMON/WEXPO/WEIGHTEX
	DOUBLE PRECISION WEIGHTEX
C--max rapidity
	common/rapmax/etamax
	double precision etamax
C--production point
	common/jetpoint/x0,y0
	double precision x0,y0
C--organisation of event record
	common/evrecord/nsim,npart,offset,hadrotype,sqrts,collider,hadro,
     &shorthepmc,channel,isochannel
	integer nsim,npart,offset,hadrotype
	double precision sqrts
	character*4 collider,channel
	character*2 isochannel
	logical hadro,shorthepmc
C--extra storage for scattering centres before interactions
      common/storescatcen/nscatcen,maxnscatcen,scatflav(10000),
     &scatcen(10000,5),writescatcen,writedummies
	integer nscatcen,maxnscatcen,scatflav
	double precision scatcen
	logical writescatcen,writedummies

C--Variables local to this program
	INTEGER NOLD,PID,IPART,LME1,LME2,j,i,LME1ORIG,LME2ORIG,llep1,
     &llep2,lv
	DOUBLE PRECISION PYR,ENI,QMAX1,R,GETMASS,PYP,Q1,Q2,P21,P22,ETOT,
     &QMAX2,POLD,EN1,EN2,BETA(3),ENEW1,ENEW2,emax,lambda,x1,x2,x3,
     &MEWEIGHT,PSWEIGHT,WEIGHT,EPS1,EPS2,THETA1,THETA2,Z1,Z2,
     &getltimemax,pi,m1,m2
	character*2 b1,b2
	CHARACTER*2 TYPE1,TYPE2
	LOGICAL FIRSTTRIP,WHICH1,WHICH2,ISDIQUARK
	DATA PI/3.141592653589793d0/

	 N=0
	 COLMAX=600
	 DISCARD=.FALSE.
       DO 91 I=1,23000
        MV(I,1)=0.d0
        MV(I,2)=0.d0
        MV(I,3)=0.d0
        MV(I,4)=0.d0
        MV(I,5)=0.d0
 91    CONTINUE
	 nscatcen = 0

       CALL MEDNEXTEVT

C--initialisation with matrix element	 
C--production vertex
        CALL PICKVTX(X0,Y0)
        LTIME=GETLTIMEMAX()

 99	  CALL PYEVNT
        NPART=N-OFFSET
        EVWEIGHT=PARI(10)
	  SUMOFWEIGHTS=SUMOFWEIGHTS+EVWEIGHT
	  IF((COLLIDER.EQ.'EEJJ').AND.(ABS(K(8,2)).GT.6))THEN
	   WDISC=WDISC+EVWEIGHT
	   NDISC=NDISC+1
	   GOTO 102
	  ELSE
	   NGOOD=NGOOD+1
	  ENDIF 

C--DY: don't have to do anything
	  if (collider.eq.'PPDY') then
	    CALL PYEXEC
	    call CONVERTTOHEPMC(HPMCFID,NGOOD,PID,b1,b2)
	    goto 102
	  endif


C--   prepare event record
	  if((COLLIDER.EQ.'PPZJ').OR.(COLLIDER.EQ.'PPZQ').or.
     &	(COLLIDER.EQ.'PPZG').or.(COLLIDER.EQ.'PPWJ').or.
     &	(COLLIDER.EQ.'PPWQ').or.(COLLIDER.EQ.'PPWG'))THEN 
             LME1ORIG=7
             LME2ORIG=8
	       if(abs(k(7,2)).gt.21) then
	         lv=7
		 else
	         lv=8
	       endif
          ELSE
             LME1ORIG=OFFSET-1
             LME2ORIG=OFFSET
          ENDIF
        DO 180 IPART=OFFSET+1, OFFSET+NPART
C--find decay leptons in V+jet events
	  if((COLLIDER.EQ.'PPZJ').OR.(COLLIDER.EQ.'PPZQ').or.
     &	(COLLIDER.EQ.'PPZG').or.(COLLIDER.EQ.'PPWJ').or.
     &	(COLLIDER.EQ.'PPWQ').or.(COLLIDER.EQ.'PPWG'))THEN 
	     if(k(ipart,3).eq.offset-1) llep1=ipart
	     if(k(ipart,3).eq.offset) llep2=ipart
	   endif
         IF(K(IPART,3).EQ.(LME1ORIG))THEN
          LME1=IPART
	    IF(K(IPART,2).EQ.21)THEN
	     TYPE1='GC'
	    ELSE
	     TYPE1='QQ'
	    ENDIF
         ELSEIF(K(IPART,3).EQ.LME2ORIG)THEN
          LME2=IPART        
	    IF(K(IPART,2).EQ.21)THEN
	     TYPE2='GC'
	    ELSE
	     TYPE2='QQ'
	    ENDIF
	   ELSE
	    TRIP(IPART)=0
	    ANTI(IPART)=0
	    ZD(IPART)=0.d0
	    THETAA(IPART)=0.d0
	   ENDIF 
C--assign colour indices
         IF(K(IPART,1).EQ.2)THEN
	    IF(K(IPART-1,1).EQ.2)THEN
C--in middle of colour singlet
	     IF(FIRSTTRIP)THEN
	      TRIP(IPART)=COLMAX+1
	      ANTI(IPART)=TRIP(IPART-1)
	     ELSE
	      TRIP(IPART)=ANTI(IPART-1)
	      ANTI(IPART)=COLMAX+1
	     ENDIF
	     COLMAX=COLMAX+1
	    ELSE
C--beginning of colour singlet
	     IF(((ABS(K(IPART,2)).LT.10).AND.(K(IPART,2).GT.0))
     &	    .OR.(ISDIQUARK(K(IPART,2)).AND.(K(IPART,2).LT.0)))THEN
	      TRIP(IPART)=COLMAX+1
	      ANTI(IPART)=0
	      FIRSTTRIP=.TRUE.
	     ELSE
	      TRIP(IPART)=0
	      ANTI(IPART)=COLMAX+1
	      FIRSTTRIP=.FALSE.
	     ENDIF
	     COLMAX=COLMAX+1
	    ENDIF
	   ENDIF 
         IF(K(IPART,1).EQ.1)THEN
C--end of colour singlet
	    IF(FIRSTTRIP)THEN
	     TRIP(IPART)=0
	     ANTI(IPART)=TRIP(IPART-1)
	    ELSE
	     TRIP(IPART)=ANTI(IPART-1)
	     ANTI(IPART)=0
	    ENDIF
	   ENDIF
 180    CONTINUE
	  if (k(lme1,1).lt.11) K(LME1,1)=1
	  if (k(lme2,1).lt.11) K(LME2,1)=1
	  PID=K(LME1,2)
	  ENI=MAX(P(LME1,4),P(LME2,4))
	  DO 183 IPART=OFFSET+1, OFFSET+NPART
	   IF((IPART.NE.LME1).AND.(IPART.NE.LME2).AND.(K(IPART,1).LT.11))
     &	   K(IPART,1)=4
	   if (k(ipart,2).eq.22) k(ipart,1)=4
 183    CONTINUE	  

C--find virtualities and adapt four-vectors
	  if((COLLIDER.EQ.'PPZJ').OR.(COLLIDER.EQ.'PPZQ').or.
     &	(COLLIDER.EQ.'PPZG').or.(COLLIDER.EQ.'PPWJ').or.
     &	(COLLIDER.EQ.'PPWQ').or.(COLLIDER.EQ.'PPWG'))THEN 
	    if (abs(k(lme1,2)).gt.21) then
           QMAX1=0.d0
           QMAX2=sqrt(pari(18)+p(lme1,5)**2)
	    else
           QMAX1=sqrt(pari(18)+p(lme2,5)**2)
           QMAX2=0.d0
	    endif
           EMAX=P(LME1,4)+P(LME2,4)
           THETA1=-1.d0
           THETA2=-1.d0
        ELSEIF(COLLIDER.EQ.'PPJJ'.OR.COLLIDER.EQ.'PPYJ'
     &          .OR.COLLIDER.EQ.'PPYQ'.OR.COLLIDER.EQ.'PPYG')THEN
	     if (k(lme1,1).eq.4) then
	       qmax1 = 0.d0
	     else
             QMAX1=pari(17)
	     endif
	     if (k(lme2,1).eq.4) then
	       qmax2 = 0.d0
	     else
             QMAX2=pari(17)
	     endif
!        QMAX1=PYP(LME1,10)*exp(0.3*abs(pyp(lme1,17)-pyp(lme2,17))/2.)/2.
!        QMAX2=PYP(LME2,10)*exp(0.3*abs(pyp(lme1,17)-pyp(lme2,17))/2.)/2.
         EMAX=P(LME1,4)+P(LME2,4)
         THETA1=-1.d0
         THETA2=-1.d0
        ENDIF 
        EN1=P(LME1,4)
        EN2=P(LME2,4)
        BETA(1)=(P(LME1,1)+P(LME2,1))/(P(LME1,4)+P(LME2,4))
        BETA(2)=(P(LME1,2)+P(LME2,2))/(P(LME1,4)+P(LME2,4))
        BETA(3)=(P(LME1,3)+P(LME2,3))/(P(LME1,4)+P(LME2,4))
        CALL PYROBO(LME1,LME1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
        CALL PYROBO(LME2,LME2,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
	  ETOT=P(LME1,4)+P(LME2,4)
	  IF(COLLIDER.EQ.'EEJJ')THEN
         QMAX1=ETOT
         QMAX2=ETOT
	   EMAX=P(LME1,4)+P(LME2,4)
	   THETA1=-1.d0
	   THETA2=-1.d0
        ENDIF
C--   find virtuality
        Q1=GETMASS(0.d0,QMAX1,THETA1,EMAX,TYPE1,EMAX,.FALSE.,
     &       Z1,WHICH1)
        Q2=GETMASS(0.d0,QMAX2,THETA2,EMAX,TYPE2,EMAX,.FALSE.,
     &       Z2,WHICH2)
 182	  if (abs(k(lme1,2)).gt.21) then
	    m1=p(lme1,5)
	  else
	    m1=q1
	  endif
 	  if (abs(k(lme2,2)).gt.21) then
	    m2=p(lme2,5)
	  else
	    m2=q2
	  endif
        ENEW1=ETOT/2.d0 + (m1**2-m2**2)/(2.*ETOT)
        ENEW2=ETOT/2.d0 - (m1**2-m2**2)/(2.*ETOT)
	  P21 = (ETOT/2.d0 + (m1**2-m2**2)/(2.*ETOT))**2 - m1**2
	  P22 = (ETOT/2.d0 - (m1**2-m2**2)/(2.*ETOT))**2 - m2**2
	  WEIGHT=1.d0
	  IF((PYR(0).GT.WEIGHT).OR.(P21.LT.0.d0).OR.(P22.LT.0.d0)
     &	.OR.(ENEW1.LT.0.d0).OR.(ENEW2.LT.0.d0)
     &	)THEN
	   IF(Q1.GT.Q2)THEN
          Q1=GETMASS(0.d0,Q1,THETA1,EMAX,TYPE1,EMAX,.FALSE.,
     &	Z1,WHICH1)
	   ELSE
          Q2=GETMASS(0.d0,Q2,THETA2,EMAX,TYPE2,EMAX,.FALSE.,
     &	Z2,WHICH2)
	   ENDIF
	   GOTO 182
	  ENDIF
        POLD=PYP(LME1,8)
	  P(LME1,1)=P(LME1,1)*SQRT(P21)/POLD
	  P(LME1,2)=P(LME1,2)*SQRT(P21)/POLD
	  P(LME1,3)=P(LME1,3)*SQRT(P21)/POLD
	  P(LME1,4)=ENEW1
	  P(LME1,5)=m1
        POLD=PYP(LME2,8)
	  P(LME2,1)=P(LME2,1)*SQRT(P22)/POLD
	  P(LME2,2)=P(LME2,2)*SQRT(P22)/POLD
	  P(LME2,3)=P(LME2,3)*SQRT(P22)/POLD
	  P(LME2,4)=ENEW2
	  P(LME2,5)=m2
        CALL PYROBO(LME1,LME1,0d0,0d0,BETA(1),BETA(2),BETA(3))
        CALL PYROBO(LME2,LME2,0d0,0d0,BETA(1),BETA(2),BETA(3))
C--correct for overestimated energy
	  IF(Q1.GT.0.d0)THEN
	   EPS1=0.5-0.5*SQRT(1.-Q0**2/Q1**2)
     &	   *SQRT(1.-Q1**2/P(LME1,4)**2)
	   IF((Z1.LT.EPS1).OR.(Z1.GT.(1.-EPS1)))THEN
          Q1=GETMASS(0.d0,Q1,THETA1,EMAX,TYPE1,EMAX,.FALSE.,
     &	Z1,WHICH1)
          CALL PYROBO(LME1,LME1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
          CALL PYROBO(LME2,LME2,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
	    GOTO 182
	   ENDIF
	  ENDIF 
	  IF(Q2.GT.0.d0)THEN
	   EPS2=0.5-0.5*SQRT(1.-Q0**2/Q2**2)
     &	   *SQRT(1.-Q2**2/P(LME2,4)**2)
         IF((Z2.LT.EPS2).OR.(Z2.GT.(1.-EPS2)))THEN
          Q2=GETMASS(0.d0,Q2,THETA2,EMAX,TYPE2,EMAX,.FALSE.,
     &	Z2,WHICH2)
          CALL PYROBO(LME1,LME1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
          CALL PYROBO(LME2,LME2,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
	    GOTO 182
         ENDIF
        ENDIF
        
C--correct to ME for first parton
	  IF(COLLIDER.EQ.'EEJJ')THEN
         BETA(1)=(P(LME1,1)+P(LME2,1))/(P(LME1,4)+P(LME2,4))
         BETA(2)=(P(LME1,2)+P(LME2,2))/(P(LME1,4)+P(LME2,4))
         BETA(3)=(P(LME1,3)+P(LME2,3))/(P(LME1,4)+P(LME2,4))
         CALL PYROBO(LME1,LME1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
         CALL PYROBO(LME2,LME2,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
         IF(Q1.GT.0.d0)THEN
C--generate z value      
	    X1=Z1*(ETOT**2+Q1**2)/ETOT**2
	    X2=(ETOT**2-Q1**2)/ETOT**2
	    X3=(1.-Z1)*(ETOT**2+Q1**2)/ETOT**2
	    PSWEIGHT=(1.-X1)*(1.+(X1/(2.-X2))**2)/X3
     &	+ (1.-X2)*(1.+(X2/(2.-X1))**2)/X3 
	    MEWEIGHT=X1**2+X2**2
	    WEIGHT=MEWEIGHT/PSWEIGHT
	    IF(PYR(0).GT.WEIGHT)THEN
 184	     Q1=GETMASS(0.d0,Q1,THETA1,EMAX,TYPE1,EMAX,.FALSE.,
     &	Z1,WHICH1)
 	    ENDIF
 	   ENDIF 
C--correct to ME for second parton
	   IF(Q2.GT.0.d0)THEN
C--generate z value      
	    X1=(ETOT**2-Q2**2)/ETOT**2
	    X2=Z2*(ETOT**2+Q2**2)/ETOT**2
	    X3=(1.-Z2)*(ETOT**2+Q2**2)/ETOT**2
	    PSWEIGHT=(1.-X1)*(1.+(X1/(2.-X2))**2)/X3
     &	+ (1.-X2)*(1.+(X2/(2.-X1))**2)/X3 
	    MEWEIGHT=X1**2+X2**2
	    WEIGHT=MEWEIGHT/PSWEIGHT
	    IF(PYR(0).GT.WEIGHT)THEN
 185	     Q2=GETMASS(0.d0,Q2,THETA2,EMAX,TYPE2,EMAX,.FALSE.,
     &	Z2,WHICH2)
	    ENDIF
	   ENDIF
 186     ENEW1=ETOT/2.d0 + (Q1**2-Q2**2)/(2.*ETOT)
         ENEW2=ETOT/2.d0 - (Q1**2-Q2**2)/(2.*ETOT)
	   P21 = (ETOT/2.d0 + (Q1**2-Q2**2)/(2.*ETOT))**2 - Q1**2
	   P22 = (ETOT/2.d0 - (Q1**2-Q2**2)/(2.*ETOT))**2 - Q2**2
         POLD=PYP(LME1,8)
	   P(LME1,1)=P(LME1,1)*SQRT(P21)/POLD
	   P(LME1,2)=P(LME1,2)*SQRT(P21)/POLD
	   P(LME1,3)=P(LME1,3)*SQRT(P21)/POLD
	   P(LME1,4)=ENEW1
	   P(LME1,5)=Q1
         POLD=PYP(LME2,8)
	   P(LME2,1)=P(LME2,1)*SQRT(P22)/POLD
	   P(LME2,2)=P(LME2,2)*SQRT(P22)/POLD
	   P(LME2,3)=P(LME2,3)*SQRT(P22)/POLD
	   P(LME2,4)=ENEW2
	   P(LME2,5)=Q2
         CALL PYROBO(LME1,LME1,0d0,0d0,BETA(1),BETA(2),BETA(3))
         CALL PYROBO(LME2,LME2,0d0,0d0,BETA(1),BETA(2),BETA(3))
C--correct for overestimated energy
	   IF(Q1.GT.0.d0)THEN
	   EPS1=0.5-0.5*SQRT(1.-Q0**2/Q1**2)
     &	   *SQRT(1.-Q1**2/P(LME1,4)**2)
	    IF((Z1.LT.EPS1).OR.(Z1.GT.(1.-EPS1)))THEN
           Q1=GETMASS(0.d0,Q1,THETA1,EMAX,TYPE1,EMAX,.FALSE.,
     &	Z1,WHICH1)
           CALL PYROBO(LME1,LME1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
           CALL PYROBO(LME2,LME2,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
	     GOTO 186
	    ENDIF
	   ENDIF 
	   IF(Q2.GT.0.d0)THEN
	   EPS2=0.5-0.5*SQRT(1.-Q0**2/Q2**2)
     &	   *SQRT(1.-Q2**2/P(LME2,4)**2)
          IF((Z2.LT.EPS2).OR.(Z2.GT.(1.-EPS2)))THEN
           Q2=GETMASS(0.d0,Q2,THETA2,EMAX,TYPE2,EMAX,.FALSE.,
     &	Z2,WHICH2)
           CALL PYROBO(LME1,LME1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
           CALL PYROBO(LME2,LME2,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
	     GOTO 186
          ENDIF
         ENDIF 
	  ENDIF

C--transfer recoil to decay leptons in V+jet
	  if((COLLIDER.EQ.'PPZJ').OR.(COLLIDER.EQ.'PPZQ').or.
     &	(COLLIDER.EQ.'PPZG').or.(COLLIDER.EQ.'PPWJ').or.
     &	(COLLIDER.EQ.'PPWQ').or.(COLLIDER.EQ.'PPWG'))THEN 
	    beta(1)=p(lv,1)/p(lv,4)
	    beta(2)=p(lv,2)/p(lv,4)
	    beta(3)=p(lv,3)/p(lv,4)
          CALL PYROBO(llep1,llep1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
          CALL PYROBO(llep2,llep2,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
	    if (abs(k(lme1,2)).gt.21) then
	      beta(1)=p(lme1,1)/p(lme1,4)
	      beta(2)=p(lme1,2)/p(lme1,4)
	      beta(3)=p(lme1,3)/p(lme1,4)
	    else
	      beta(1)=p(lme2,1)/p(lme2,4)
	      beta(2)=p(lme2,2)/p(lme2,4)
	      beta(3)=p(lme2,3)/p(lme2,4)
	    endif
          CALL PYROBO(llep1,llep1,0d0,0d0,BETA(1),BETA(2),BETA(3))
          CALL PYROBO(llep2,llep2,0d0,0d0,BETA(1),BETA(2),BETA(3))
	  endif

  
        ZA(LME1)=1.d0
        ZA(LME2)=1.d0
	  THETAA(LME1)=P(LME1,5)/(SQRT(Z1*(1.-Z1))*P(LME1,4))
	  THETAA(LME2)=P(LME2,5)/(SQRT(Z2*(1.-Z2))*P(LME2,4))
	  ZD(LME1)=Z1
	  ZD(LME2)=Z2
	  QQBARD(LME1)=WHICH1
	  QQBARD(LME2)=WHICH2

        MV(LME1,1)=X0
        MV(LME1,2)=Y0
        MV(LME1,3)=0.d0
        MV(LME1,4)=0.d0
        IF(P(LME1,5).GT.0.d0)THEN
         LAMBDA=1.d0/(FTFAC*P(LME1,4)*0.2/Q1**2)
          MV(LME1,5)=-LOG(1.d0-PYR(0))/LAMBDA
        ELSE
         MV(LME1,5)=LTIME
        ENDIF
         
        MV(LME2,1)=X0
        MV(LME2,2)=Y0
        MV(LME2,3)=0.d0
        MV(LME2,4)=0.d0
        IF(P(LME2,5).GT.0.d0)THEN
         LAMBDA=1.d0/(FTFAC*P(LME2,4)*0.2/Q2**2)
          MV(LME2,5)=-LOG(1.d0-PYR(0))/LAMBDA
        ELSE
         MV(LME2,5)=LTIME
        ENDIF

C--develop parton shower
	 CALL MAKECASCADE
	 IF(DISCARD) THEN
	  NGOOD=NGOOD-1
 	  WDISC=WDISC+EVWEIGHT
	  NDISC=NDISC+1
        write(logfid,*)'discard event',J
	  GOTO 102
	 ENDIF

       IF(.NOT.ALLHAD)THEN
        DO 86 I=1,N
         IF(K(I,1).EQ.3) K(I,1)=22
 86     CONTINUE
       ENDIF
       IF(HADRO)THEN
        CALL MAKESTRINGS(HADROTYPE)
	  IF(DISCARD) THEN
         write(logfid,*)'discard event',J
	   WDISC=WDISC+EVWEIGHT
	   NDISC=NDISC+1
	   NGOOD=NGOOD-1
	   GOTO 102
	  ENDIF
        CALL PYEXEC
	  IF(MSTU(30).NE.ERRCOUNT)THEN
         write(logfid,*)'PYTHIA discards event',J,
     &	'  (error number',MSTU(30),')'
	   ERRCOUNT=MSTU(30)
	   WDISC=WDISC+EVWEIGHT
	   NDISC=NDISC+1
	   NGOOD=NGOOD-1
	   GOTO 102
	  ENDIF
       ENDIF

	 IF(MSTU(30).NE.ERRCOUNT)THEN
	  ERRCOUNT=MSTU(30)
	 ELSE 
	  CALL CONVERTTOHEPMC(HPMCFID,NGOOD,PID,b1,b2)
	 ENDIF

C--write message to log-file
 102  IF(NSIM.GT.100)THEN
       IF(MOD(J,NSIM/100).EQ.0)THEN
 	  write(logfid,*) 'done with event number ',J
 	 ENDIF
	else
 	  write(logfid,*) 'done with event number ',J
      ENDIF
	call flush(logfid)
	end



***********************************************************************
***	  subroutine makestrings
***********************************************************************
	SUBROUTINE MAKESTRINGS(WHICH)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
	INTEGER WHICH
	IF(WHICH.EQ.0)THEN
	 CALL MAKESTRINGS_VAC
	ELSEIF(WHICH.EQ.1)THEN
	 CALL MAKESTRINGS_MINL
	ELSE
	WRITE(logfid,*)'error: unknown hadronisation type in MAKESTRINGS'
	ENDIF
	END


***********************************************************************
***	  subroutine makestrings_vac
***********************************************************************
      SUBROUTINE MAKESTRINGS_VAC
      IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--local variables
      INTEGER NOLD,I,J,LQUARK,LMATCH,LLOOSE,NOLD1
      DOUBLE PRECISION EADDEND,PYR,DIR
      LOGICAL ISDIQUARK,compressevent,roomleft
      DATA EADDEND/10.d0/
	
	i = 0
	if (compress) roomleft = compressevent(i)
      NOLD1=N
C--remove all active lines that are leptons, gammas, hadrons etc.
	DO 52 I=1,NOLD1
	 IF((K(I,1).EQ.4).AND.(TRIP(I).EQ.0).AND.(ANTI(I).EQ.0))THEN
C--copy line to end of event record
        N=N+1
        IF(N.GT.22990) THEN
         write(logfid,*)'event too long for event record'
         DISCARD=.TRUE.
         RETURN
        ENDIF
        K(N,1)=11
        K(N,2)=K(I,2)
        K(N,3)=I
        K(N,4)=0
        K(N,5)=0
        P(N,1)=P(I,1)
        P(N,2)=P(I,2)
        P(N,3)=P(I,3)
        P(N,4)=P(I,4)
        P(N,5)=P(I,5)
        K(I,1)=17
        K(I,4)=N
        K(I,5)=N
	  TRIP(N)=TRIP(I)
	  ANTI(N)=ANTI(I)
	 ENDIF
 52	CONTINUE
      NOLD=N
C--first do strings with existing (anti)triplets
C--find string end (=quark or antiquark)
 43   LQUARK=0
      DO 40 I=1,NOLD
       IF((K(I,1).EQ.11).OR.(K(I,1).EQ.12).OR.(K(I,1).EQ.13)
     &            .OR.(K(I,1).EQ.14)) K(I,1)=17
       IF(((K(I,1).EQ.1).OR.(K(I,1).EQ.3).OR.(K(I,1).EQ.4).OR.
     &   (K(I,1).EQ.5)).AND.((K(I,2).LT.6).OR.ISDIQUARK(K(I,2))))THEN
        LQUARK=I
	  GOTO 41
       ENDIF
 40   CONTINUE
	GOTO 50
 41	CONTINUE
C--copy string end to end of event record
      N=N+1
      IF(N.GT.22990) THEN
       write(logfid,*)'event too long for event record'
       DISCARD=.TRUE.
       RETURN
      ENDIF
      K(N,1)=2
      K(N,2)=K(LQUARK,2)
      K(N,3)=LQUARK
      K(N,4)=0
      K(N,5)=0
      P(N,1)=P(LQUARK,1)
      P(N,2)=P(LQUARK,2)
      P(N,3)=P(LQUARK,3)
      P(N,4)=P(LQUARK,4)
      P(N,5)=P(LQUARK,5)
      K(LQUARK,1)=16
      K(LQUARK,4)=N
      K(LQUARK,5)=N
	TRIP(N)=TRIP(LQUARK)
	ANTI(N)=ANTI(LQUARK)
C--append matching colour partner
	LMATCH=0
	DO 44 J=1,10000000
	 DO 42 I=1,NOLD
	  IF(((K(I,1).EQ.1).OR.(K(I,1).EQ.3).OR.(K(I,1).EQ.4)
     &						.OR.(K(I,1).EQ.5))
     &      .AND.(((TRIP(I).EQ.ANTI(N)).AND.(TRIP(I).NE.0))
     &		.OR.((ANTI(I).EQ.TRIP(N)).AND.(ANTI(I).NE.0))))THEN
         N=N+1
         IF(N.GT.22990) THEN
          write(logfid,*)'event too long for event record'
          DISCARD=.TRUE.
          RETURN
         ENDIF
         K(N,2)=K(I,2)
         K(N,3)=I
         K(N,4)=0
         K(N,5)=0
         P(N,1)=P(I,1)
         P(N,2)=P(I,2)
         P(N,3)=P(I,3)
         P(N,4)=P(I,4)
         P(N,5)=P(I,5)
	   TRIP(N)=TRIP(I)
	   ANTI(N)=ANTI(I)
         K(I,1)=16
         K(I,4)=N
         K(I,5)=N
         IF(K(I,2).EQ.21)THEN
          K(N,1)=2
          GOTO 44
         ELSE
          K(N,1)=1
          GOTO 43
         ENDIF
	  ENDIF
 42	 CONTINUE
C--no matching colour partner found
	 write(logfid,*)'Error in MAKESTRINGS_VAC: failed to reconstruct '//
     &'colour singlet system, will discard event'
	 discard = .true.
	 return
 44	CONTINUE
C--now take care of purely gluonic remainder system
C-----------------------------------------
C--find gluon where anti-triplet is not matched
 50   LLOOSE=0
      DO 45 I=1,NOLD
       IF(((K(I,1).EQ.1).OR.(K(I,1).EQ.3).OR.(K(I,1).EQ.4)
     &					.OR.(K(I,1).EQ.5)))THEN
	  DO 46 J=1,NOLD
	   IF(((K(I,1).EQ.1).OR.(K(I,1).EQ.3).OR.(K(I,1).EQ.4)
     &					.OR.(K(I,1).EQ.5)))THEN
	    IF(ANTI(I).EQ.TRIP(J)) GOTO 45
	   ENDIF
 46	  CONTINUE
        LLOOSE=I
	  GOTO 47
       ENDIF
 45   CONTINUE
	GOTO 51
 47	CONTINUE
C--generate artificial triplet end
	 write(logfid,*)'Error in MAKESTRINGS_VAC: failed to reconstruct '//
     &'colour singlet system, will discard event'
	 discard = .true.
	 return
C--copy loose gluon to end of event record
      N=N+1
      IF(N.GT.22990) THEN
       write(logfid,*)'event too long for event record'
       DISCARD=.TRUE.
       RETURN
      ENDIF
      K(N,1)=2
      K(N,2)=K(LLOOSE,2)
      K(N,3)=LLOOSE
      K(N,4)=0
      K(N,5)=0
      P(N,1)=P(LLOOSE,1)
      P(N,2)=P(LLOOSE,2)
      P(N,3)=P(LLOOSE,3)
      P(N,4)=P(LLOOSE,4)
      P(N,5)=P(LLOOSE,5)
      K(LLOOSE,1)=16
      K(LLOOSE,4)=N
      K(LLOOSE,5)=N
	TRIP(N)=TRIP(LLOOSE)
	ANTI(N)=ANTI(LLOOSE)
C--append matching colour partner
	LMATCH=0
	DO 48 J=1,10000000
	 DO 49 I=1,NOLD
	  IF(((K(I,1).EQ.1).OR.(K(I,1).EQ.3).OR.(K(I,1).EQ.4)
     &				.OR.(K(I,1).EQ.5))
     &		.AND.(ANTI(I).EQ.TRIP(N)))THEN
         N=N+1
         IF(N.GT.22990) THEN
          write(logfid,*)'event too long for event record'
          DISCARD=.TRUE.
          RETURN
         ENDIF
         K(N,2)=K(I,2)
         K(N,3)=I
         K(N,4)=0
         K(N,5)=0
         P(N,1)=P(I,1)
         P(N,2)=P(I,2)
         P(N,3)=P(I,3)
         P(N,4)=P(I,4)
         P(N,5)=P(I,5)
	   TRIP(N)=TRIP(I)
	   ANTI(N)=ANTI(I)
         K(I,1)=16
         K(I,4)=N
         K(I,5)=N
         K(N,1)=2
         GOTO 48
	  ENDIF
 49	 CONTINUE
C--no matching colour partner found, add artificial end point
	 write(logfid,*)'Error in MAKESTRINGS_VAC: failed to reconstruct '//
     &'colour singlet system, will discard event'
	 discard = .true.
	 return
 48	CONTINUE
 51	CONTINUE
	CALL CLEANUP(NOLD1)
	END


***********************************************************************
***	  subroutine makestrings_minl
***********************************************************************
      SUBROUTINE MAKESTRINGS_MINL
      IMPLICIT NONE
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
C--local variables
      INTEGER NOLD,I,J,LMAX,LMIN,LEND,nold1
      DOUBLE PRECISION EMAX,MINV,MMIN,Z,GENERATEZ,MCUT,EADDEND,PYR,DIR,
     &pyp
      DATA MCUT/1.d8/
      DATA EADDEND/10.d0/
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
	logical compressevent,roomleft

	 i = 0
	 if (compress) roomleft = compressevent(i)
      NOLD1=N
C--remove all active lines that are leptons, gammas, hadrons etc.
	DO 52 I=1,NOLD1
	 IF((K(I,1).EQ.4).AND.(TRIP(I).EQ.0).AND.(ANTI(I).EQ.0))THEN
C--copy line to end of event record
        N=N+1
        IF(N.GT.22990) THEN
         write(logfid,*)'event too long for event record'
         DISCARD=.TRUE.
         RETURN
        ENDIF
        K(N,1)=11
        K(N,2)=K(I,2)
        K(N,3)=I
        K(N,4)=0
        K(N,5)=0
        P(N,1)=P(I,1)
        P(N,2)=P(I,2)
        P(N,3)=P(I,3)
        P(N,4)=P(I,4)
        P(N,5)=P(I,5)
        K(I,1)=17
        K(I,4)=N
        K(I,5)=N
	  TRIP(N)=TRIP(I)
	  ANTI(N)=ANTI(I)
	 ENDIF
 52	CONTINUE
       NOLD=N
C--find most energetic unfragmented parton in event
 43    EMAX=0
       LMAX=0
       DO 40 I=1,NOLD
        IF((K(I,1).EQ.11).OR.(K(I,1).EQ.12).OR.(K(I,1).EQ.13)
     &            .OR.(K(I,1).EQ.14)) K(I,1)=17
        if (abs(pyp(I,17)).gt.4.d0) k(i,1)=17
        IF(((K(I,1).EQ.1).OR.(K(I,1).EQ.3).OR.(K(I,1).EQ.4)
     &	.OR.(K(I,1).EQ.5)).AND.(P(I,4).GT.EMAX))THEN
         EMAX=P(I,4)
         LMAX=I
        ENDIF
 40    CONTINUE
C--if there is non, we are done
       IF(LMAX.EQ.0) GOTO 50
C--check if highest energy parton is (anti)quark or gluon
       IF(K(LMAX,2).EQ.21)THEN
C--split gluon in qqbar pair and store one temporarily in line 1
C--make new line in event record for string end
        N=N+2
        IF(N.GT.22990) THEN
         write(logfid,*)'event too long for event record'
         DISCARD=.TRUE.
         RETURN
        ENDIF
	  IF((N-2).GT.NOLD)THEN
         DO 47 J=NOLD,N-3
          K(N+NOLD-J,1)=K(N+NOLD-J-2,1)
          K(N+NOLD-J,2)=K(N+NOLD-J-2,2)
          IF(K(N+NOLD-J-2,3).GT.NOLD) THEN
           K(N+NOLD-J,3)=K(N+NOLD-J-2,3)+2
          ELSE
           K(N+NOLD-J,3)=K(N+NOLD-J-2,3)
          ENDIF
          K(N+NOLD-J,4)=0
          K(N+NOLD-J,5)=0
          P(N+NOLD-J,1)=P(N+NOLD-J-2,1)
          P(N+NOLD-J,2)=P(N+NOLD-J-2,2)
          P(N+NOLD-J,3)=P(N+NOLD-J-2,3)
          P(N+NOLD-J,4)=P(N+NOLD-J-2,4)
          P(N+NOLD-J,5)=P(N+NOLD-J-2,5)
          K(K(N+NOLD-J-2,3),4)=K(K(N+NOLD-J-2,3),4)+2
          K(K(N+NOLD-J-2,3),5)=K(K(N+NOLD-J-2,3),5)+2
 47      CONTINUE
	  ENDIF
        NOLD=NOLD+2
        K(LMAX,1)=18
        Z=GENERATEZ(0.d0,0.d0,1.d-3,'QG')
        IF(Z.GT.0.5)THEN
         K(NOLD-1,2)=1
         K(NOLD,2)=-1
        ELSE
         Z=1.-Z
         K(NOLD-1,2)=-1
         K(NOLD,2)=1
        ENDIF
        K(NOLD-1,1)=1
        K(NOLD-1,3)=LMAX
        K(NOLD-1,4)=0
        K(NOLD-1,5)=0
        P(NOLD-1,1)=(1.-Z)*P(LMAX,1)
        P(NOLD-1,2)=(1.-Z)*P(LMAX,2)
        P(NOLD-1,3)=(1.-Z)*P(LMAX,3)
        P(NOLD-1,4)=(1.-Z)*P(LMAX,4)
        P(NOLD-1,5)=P(LMAX,5)
        K(NOLD,1)=1
        K(NOLD,3)=LMAX
        K(NOLD,4)=0
        K(NOLD,5)=0
        P(NOLD,1)=Z*P(LMAX,1)
        P(NOLD,2)=Z*P(LMAX,2)
        P(NOLD,3)=Z*P(LMAX,3)
        P(NOLD,4)=Z*P(LMAX,4)
        P(NOLD,5)=P(LMAX,5)
        K(LMAX,1)=18
        K(LMAX,4)=NOLD-1
        K(LMAX,5)=NOLD
        LMAX=NOLD
       ENDIF
       N=N+1
       IF(N.GT.22990) THEN
        write(logfid,*)'event too long for event record'
        DISCARD=.TRUE.
        RETURN
       ENDIF
       K(N,1)=2
       K(N,2)=K(LMAX,2)
       K(N,3)=LMAX
       K(N,4)=0
       K(N,5)=0
       P(N,1)=P(LMAX,1)
       P(N,2)=P(LMAX,2)
       P(N,3)=P(LMAX,3)
       P(N,4)=P(LMAX,4)
       P(N,5)=P(LMAX,5)
       K(LMAX,1)=16
       K(LMAX,4)=N
       K(LMAX,5)=N
       LEND=LMAX
C--find closest partner
 42    MMIN=1.d10
       LMIN=0
       DO 41 I=1,NOLD
        IF(((K(I,1).EQ.1).OR.(K(I,1).EQ.3).OR.(K(I,1)
     &			.EQ.4).OR.(K(I,1).EQ.5))
     &      .AND.((K(I,2).EQ.21).OR.((K(I,2)*K(LEND,2).LT.0.d0).AND.
     &		(K(I,3).NE.K(LEND,3))))
     &      .AND.(P(I,1)*P(LEND,1).GT.0.d0))THEN
         MINV=P(I,4)*P(LMAX,4)-P(I,1)*P(LMAX,1)-P(I,2)*P(LMAX,2)
     &            -P(I,3)*P(LMAX,3)
         IF((MINV.LT.MMIN).AND.(MINV.GT.0.d0).AND.(MINV.LT.MCUT))THEN
          MMIN=MINV
          LMIN=I
         ENDIF
        ENDIF
 41    CONTINUE
C--if no closest partner can be found, generate artificial end point for string
       IF(LMIN.EQ.0)THEN
        N=N+1
        IF(N.GT.22990) THEN
         write(logfid,*)'event too long for event record'
         DISCARD=.TRUE.
         RETURN
        ENDIF
        K(N,1)=1
        K(N,2)=-K(LEND,2)
        K(N,3)=0
        K(N,4)=0
        K(N,5)=0
        P(N,1)=0.d0
        P(N,2)=0.d0
        IF(PYR(0).LT.0.5)THEN
         DIR=1.d0
        ELSE
         DIR=-1.d0
        ENDIF
        P(N,3)=DIR*EADDEND
        P(N,4)=EADDEND
        P(N,5)=0.d0
        GOTO 43
       ELSE
C--else build closest partner in string
        N=N+1
        IF(N.GT.22990) THEN
         write(logfid,*)'event too long for event record'
         DISCARD=.TRUE.
         RETURN
        ENDIF
        K(N,2)=K(LMIN,2)
        K(N,3)=LMIN
        K(N,4)=0
        K(N,5)=0
        P(N,1)=P(LMIN,1)
        P(N,2)=P(LMIN,2)
        P(N,3)=P(LMIN,3)
        P(N,4)=P(LMIN,4)
        P(N,5)=P(LMIN,5)
        K(LMIN,1)=16
        K(LMIN,4)=N
        K(LMIN,5)=N
        IF(K(LMIN,2).EQ.21)THEN
         K(N,1)=2
         LMAX=LMIN
         GOTO 42
        ELSE
         K(N,1)=1
         GOTO 43
        ENDIF
       ENDIF
 50    CONTINUE
       CALL CLEANUP(NOLD)
      END


***********************************************************************
***	  subroutine cleanup
***********************************************************************
	SUBROUTINE CLEANUP(NFIRST)
	IMPLICIT NONE
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--local variables
	INTEGER NFIRST,NLAST,I,J
	
	NLAST=N
	DO 21 I=1,NLAST-NFIRST
	 DO 22 J=1,5
	  K(I,J)=K(NFIRST+I,J)
	  P(I,J)=P(NFIRST+I,J)
	  V(I,J)=V(NFIRST+I,J)
 22	 CONTINUE
	 K(I,3)=0	 
 21	CONTINUE
      N=NLAST-NFIRST
	END


***********************************************************************
***	  subroutine makecascade
***********************************************************************
	SUBROUTINE MAKECASCADE
      IMPLICIT NONE
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc

C--local variables
	INTEGER NOLD,I
	LOGICAL CONT

 10	NOLD=N
	CONT=.FALSE.
 	DO 11 I=2,NOLD
	 if (i.gt.n) goto 10
C--check if parton may evolve, i.e. do splitting or scattering
	 IF((K(I,1).EQ.1).OR.(K(I,1).EQ.2))THEN
	  CONT=.TRUE.
	  CALL MAKEBRANCH(I)
	  IF(DISCARD) GOTO 12
	 ENDIF
 11	CONTINUE
 	IF(CONT) GOTO 10
 12	END


***********************************************************************
***	  subroutine makebranch
***********************************************************************
      SUBROUTINE MAKEBRANCH(L)
      IMPLICIT NONE
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--number of scattering events
	COMMON/CHECK/NSCAT,NSCATEFF,NSPLIT
	DOUBLE PRECISION NSCAT,NSCATEFF,NSPLIT
C--variables for coherent scattering
	COMMON/COHERENT/NSTART,NEND,ALLQS(10000,6),SCATCENTRES(10000,10),
     &QSUMVEC(4),QSUM2
	INTEGER NSTART,NEND
	DOUBLE PRECISION ALLQS,SCATCENTRES,QSUMVEC,QSUM2
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--extra storage for scattering centres before interactions
       common/storescatcen/nscatcen,maxnscatcen,scatflav(10000),
     & scatcen(10000,5),writescatcen,writedummies
	 integer nscatcen,maxnscatcen,scatflav
	 double precision scatcen
	 logical writescatcen,writedummies
C--local variables
      INTEGER L,LINE,NOLD,TYPI,LINEOLD,LKINE,nendold,nscatcenold
      DOUBLE PRECISION THETA,PHI,PYP,FORMTIME,STARTTIME,TLEFT,
     &TSUM,DELTAT,NEWMASS,GETMASS,Q,GETMS,ZDEC,X,DTCORR
	LOGICAL OVERQ0,QQBARDEC
	CHARACTER TYP
	LOGICAL RADIATION,RETRYSPLIT,MEDIND,roomleft,compressevent

	LINE=L
	NSTART=0
	NEND=0
	STARTTIME=MV(LINE,4)
	TSUM=0.d0
	QSUM2=0.d0
	QSUMVEC(1)=0.d0
	QSUMVEC(2)=0.d0
	QSUMVEC(3)=0.d0
	QSUMVEC(4)=0.d0
	RETRYSPLIT=.FALSE.
      MEDIND=.FALSE.
	X=0.d0
	Q=0.d0
	TYPI=0

      IF ((N.GT.20000).and.compress) roomleft = compressevent(line)

20	IF(DISCARD) RETURN
	IF(((K(LINE,1).EQ.1).AND.(P(LINE,5).GT.0.d0))
     &	.OR.((K(LINE,1).EQ.2).AND.(zd(line).gt.0.d0)))THEN
       IF(MEDIND)THEN
        FORMTIME=starttime
       ELSE 
	  FORMTIME=MIN(MV(LINE,5),LTIME)
	 ENDIF
	 RADIATION=.TRUE.
	ELSE
	 FORMTIME=LTIME
	 RADIATION=.FALSE.
	ENDIF
	TLEFT=FORMTIME-STARTTIME
      IF(K(LINE,2).EQ.21)THEN
       TYP='G'
      ELSE
       TYP='Q'
      ENDIF
      MEDIND=.FALSE.

      IF(TLEFT.LE.1.d-10)THEN
C--no scattering
	 IF(RADIATION)THEN
C--if there is radiation associated with the parton then form it now
C--rotate such that momentum points in z-direction
        NOLD=N
        nscatcenold=nscatcen
        THETA=PYP(LINE,13)
        PHI=PYP(LINE,15)
        CALL PYROBO(LINE,LINE,0d0,-PHI,0d0,0d0,0d0)
        CALL PYROBO(LINE,LINE,-THETA,0d0,0d0,0d0,0d0)
        CALL MAKESPLITTING(LINE)
C--rotate back
        CALL PYROBO(LINE,LINE,THETA,0d0,0d0,0d0,0d0)
        CALL PYROBO(LINE,LINE,0d0,PHI,0d0,0d0,0d0)
        IF(DISCARD) RETURN
        CALL PYROBO(N-1,N,THETA,0d0,0d0,0d0,0d0)
        CALL PYROBO(N-1,N,0d0,PHI,0d0,0d0,0d0)
C--set the production vertices: x_mother + (tprod - tprod_mother) * beta_mother
        MV(N-1,1)=MV(LINE,1)
     &	+(MV(N-1,4)-MV(LINE,4))*P(LINE,1)/max(pyp(line,8),P(LINE,4))
        MV(N-1,2)=MV(LINE,2)
     &	+(MV(N-1,4)-MV(LINE,4))*P(LINE,2)/max(pyp(line,8),P(LINE,4))
        MV(N-1,3)=MV(LINE,3)
     &	+(MV(N-1,4)-MV(LINE,4))*P(LINE,3)/max(pyp(line,8),P(LINE,4))
        MV(N,  1)=MV(LINE,1)
     &	+(MV(N,  4)-MV(LINE,4))*P(LINE,1)/max(pyp(line,8),P(LINE,4))
        MV(N,  2)=MV(LINE,2)
     &	+(MV(N,  4)-MV(LINE,4))*P(LINE,2)/max(pyp(line,8),P(LINE,4))
        MV(N,  3)=MV(LINE,3)
     &	+(MV(N,  4)-MV(LINE,4))*P(LINE,3)/max(pyp(line,8),P(LINE,4))

	  LINE=N
	  NSTART=0
	  NEND=0
	  STARTTIME=MV(N,4)
	  QSUMVEC(1)=0.d0
	  QSUMVEC(2)=0.d0
	  QSUMVEC(3)=0.d0
	  QSUMVEC(4)=0.d0
	  QSUM2=0.d0
	  TSUM=0.d0
	  GOTO 21
	 ELSE
	  NSTART=0
	  NEND=0
	  STARTTIME=FORMTIME
	  QSUMVEC(1)=0.d0
	  QSUMVEC(2)=0.d0
	  QSUMVEC(3)=0.d0
	  QSUMVEC(4)=0.d0
	  QSUM2=0.d0
	  TSUM=0.d0
	  GOTO 21
	 ENDIF
	ELSE
C--do scattering
C--find delta t for the scattering
	 DELTAT=TLEFT
	 OVERQ0=.FALSE.
	 CALL DOINSTATESCAT(LINE,X,TYPI,Q,STARTTIME+TSUM,DELTAT,
     &		OVERQ0,.FALSE.)
	 TSUM=TSUM+DELTAT
	 TLEFT=TLEFT-DELTAT
C--do initial state splitting if there is one
	 NOLD=N
	 LINEOLD=LINE
	 ZDEC=ZD(LINE)
	 QQBARDEC=QQBARD(LINE)
        nscatcenold=nscatcen
 25	 IF(X.LT.1.d0) THEN
	  CALL MAKEINSPLIT(LINE,X,QSUM2,Q,TYPI,STARTTIME+TSUM,DELTAT)
        IF(DISCARD) RETURN
	  IF(X.LT.1.d0)THEN
	   LINE=N
	   LKINE=N
	   IF(K(LINE,2).EQ.21)THEN
	    NEWMASS=GETMASS(0.d0,SCALEFACM*SQRT(-QSUM2),-1.d0,P(LINE,4),
     &			'GC',SQRT(-QSUM2),.FALSE.,ZDEC,QQBARDEC)
          IF(ZDEC.GT.0.d0)THEN
           THETAA(LINE)=NEWMASS/(SQRT(ZDEC*(1.-ZDEC))*P(LINE,4))
          ELSE
           THETAA(LINE)=0.d0
          ENDIF 
	    ZD(LINE)=ZDEC
	    QQBARD(LINE)=QQBARDEC
	   ELSE	
	    NEWMASS=GETMASS(0.d0,SCALEFACM*SQRT(-QSUM2),-1.d0,P(LINE,4),
     &			'QQ',SQRT(-QSUM2),.FALSE.,ZDEC,QQBARDEC)
	    IF(ZDEC.GT.0.d0)THEN
           THETAA(LINE)=NEWMASS/(SQRT(ZDEC*(1.-ZDEC))*P(LINE,4))
          ELSE
           THETAA(LINE)=0.d0
          ENDIF 
	    ZD(LINE)=ZDEC
	    QQBARD(LINE)=QQBARDEC
	   ENDIF
	   ZDEC=ZD(LINE)
	   QQBARDEC=QQBARD(LINE)
	  ELSE
	   LKINE=LINE
	   NEND=NSTART
	   QSUM2=ALLQS(NEND,1)
	   QSUMVEC(1)=ALLQS(NEND,2)
	   QSUMVEC(2)=ALLQS(NEND,3)
	   QSUMVEC(3)=ALLQS(NEND,4)
	   QSUMVEC(4)=ALLQS(NEND,5)
	   IF(-ALLQS(NEND,1).GT.Q0**2/SCALEFACM**2)THEN
	    OVERQ0=.TRUE.
	   ELSE
	    OVERQ0=.FALSE.
	   ENDIF
	   tleft = starttime+tsum+tleft-allqs(1,6)
	   tsum = allqs(1,6)-starttime
	  ENDIF 
	 ENDIF
	 IF(X.EQ.1.d0)THEN
	  NEWMASS=0.d0
	  IF(NEND.GT.0)THEN
	   CALL DOFISTATESCAT(LINE,STARTTIME+TSUM,TLEFT,DELTAT,
     &		NEWMASS,OVERQ0,ZDEC,QQBARDEC)
	   IF(NEWMASS.GT.(P(LINE,5)*(1.d0+1.d-6)))THEN
	    MEDIND=.TRUE.
	   ELSE
	    MEDIND=.FALSE.
	    ZDEC=ZD(LINE)
	    QQBARDEC=QQBARD(LINE)
	   ENDIF 
	   TSUM=TSUM+DELTAT
	   TLEFT=TLEFT-DELTAT
	   LKINE=LINE
	  ENDIF
	 ENDIF
C--do kinematics
	 RETRYSPLIT=.FALSE.
	 IF(NEND.GT.0) THEN
	  nendold=nend
	  CALL DOKINEMATICS(LKINE,lineold,NSTART,NEND,NEWMASS,RETRYSPLIT,
     &		STARTTIME+TSUM,X,ZDEC,QQBARDEC)
	  IF(RETRYSPLIT) THEN
	   tleft = starttime+tsum+tleft-allqs(1,6)
	   tsum = allqs(1,6)-starttime
	   if (x.lt.1.d0) then
	     NEND=NSTART
	     QSUM2=ALLQS(NEND,1)
	     QSUMVEC(1)=ALLQS(NEND,2)
	     QSUMVEC(2)=ALLQS(NEND,3)
	     QSUMVEC(3)=ALLQS(NEND,4)
	     QSUMVEC(4)=ALLQS(NEND,5)
	     TYPI=K(L,2)
	     IF(-ALLQS(NEND,1).GT.Q0**2/SCALEFACM**2)THEN
	       OVERQ0=.TRUE.
	     ELSE
	       OVERQ0=.FALSE.
	     ENDIF
	     N=NOLD
	     LINE=LINEOLD
	     X=1.d0
	     K(LINE,1)=1
	     nscatcen=nscatcenold
	     NSPLIT=NSPLIT-EVWEIGHT
	     GOTO 25
	   else
	     LINE=N
	     STARTTIME=STARTTIME+TSUM
	     TSUM=0.d0
	   endif
	  ELSE
	   LINE=N
	   STARTTIME=STARTTIME+TSUM
	   TSUM=0.d0
	  ENDIF
	 ELSE
	  STARTTIME=STARTTIME+TSUM
	  TSUM=0.d0
	 ENDIF
	 IF(P(LINE,5).GT.0.d0) RADIATION=.TRUE.
	ENDIF

 21   IF(((K(LINE,1).EQ.1).AND.(P(LINE,5).GT.0.d0))
     &	.OR.((K(LINE,1).EQ.2).AND.(zd(line).gt.0.d0))
     &	.OR.(STARTTIME.LT.LTIME))THEN
	 GOTO 20
	ENDIF
	IF((K(LINE,1).EQ.1).AND.(P(LINE,5).EQ.0.d0)) K(LINE,1)=4
	IF((K(LINE,1).EQ.2).AND.(zd(line).lt.0.d0)) K(LINE,1)=5
      END



***********************************************************************
***	  subroutine makesplitting
***********************************************************************
	SUBROUTINE MAKESPLITTING(L)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--factor in front of formation times
	COMMON/FTIMEFAC/FTFAC
	DOUBLE PRECISION FTFAC
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--number of scattering events
	COMMON/CHECK/NSCAT,NSCATEFF,NSPLIT
	DOUBLE PRECISION NSCAT,NSCATEFF,NSPLIT
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights

C--local variables
	INTEGER L,DIR
	DOUBLE PRECISION PHIQ,PYR,PI,GENERATEZ,BMAX1,CMAX1,PTS,MB,MC,
     &GETMASS,PZ,EPS,QH,Z,R,LAMBDA,WEIGHT,ZDECB,ZDECC,XDEC(3),THETA,
     &GETTEMP
      LOGICAL QUARK,QQBAR,QQBARDECB,QQBARDECC
	integer bin
	DATA PI/3.141592653589793d0/

      IF((N+2).GT.22990) THEN
       write(logfid,*)'event too long for event record'
       DISCARD=.TRUE.
       RETURN
      ENDIF

      XDEC(1)=MV(L,1)+(MV(L,5)-MV(L,4))*P(L,1)/P(L,4)
      XDEC(2)=MV(L,2)+(MV(L,5)-MV(L,4))*P(L,2)/P(L,4)
      XDEC(3)=MV(L,3)+(MV(L,5)-MV(L,4))*P(L,3)/P(L,4)
	IF(GETTEMP(XDEC(1),XDEC(2),XDEC(3),MV(L,5)).GT.0.d0)THEN
	 THETA=-1.d0
	ELSE
	 THETA=THETAA(L)
	ENDIF 

C--on-shell partons cannot split
	IF((P(L,5).EQ.0d0).OR.(K(L,1).EQ.11).OR.(K(L,1).EQ.12)
     &  .OR.(K(L,1).EQ.13).OR.(K(L,1).EQ.14).OR.(K(L,1).EQ.3)
     &  .or.(zd(l).lt.0.d0)) GOTO 31
C--quark or gluon?
	IF(K(L,2).EQ.21)THEN
	 QUARK=.FALSE.
	ELSE
	 QUARK=.TRUE.
	 QQBAR=.FALSE.
	ENDIF
C--if gluon decide on kind of splitting
	QQBAR=QQBARD(L)
C--if g->gg splitting decide on colour order
	IF(QUARK.OR.QQBAR)THEN
	 DIR=0
	ELSE
	 IF(PYR(0).LT.0.5)THEN
	  DIR=1
	 ELSE
	  DIR=-1
	 ENDIF
	ENDIF
	Z=ZD(L)
	IF(Z.EQ.0.d0)THEN
	 write(logfid,*)'makesplitting: z=0',L
	 goto 36
	ENDIF  
	GOTO 35
C--generate z value
 36	IF(ANGORD.AND.(ZA(L).NE.1.d0))THEN
C--additional z constraint due to angular ordering
	 QH=4.*P(L,5)**2*(1.-ZA(L))/(ZA(L)*P(K(L,3),5)**2)
	 IF(QH.GT.1)THEN
	  write(logfid,*)L,': reject event: angular ordering
     &      conflict in medium'
	  CALL PYLIST(3)
	  DISCARD=.TRUE.
	  GOTO 31
	 ENDIF
	 EPS=0.5-0.5*SQRT(1.-QH)
	ELSE
	 EPS=0d0
	ENDIF
 	IF(QUARK)THEN
	 Z=GENERATEZ(P(L,5)**2,P(L,4),EPS,'QQ')
	ELSE
	 IF(QQBAR)THEN
	  Z=GENERATEZ(P(L,5)**2,P(L,4),EPS,'QG')
	 ELSE
	  Z=GENERATEZ(P(L,5)**2,P(L,4),EPS,'GG')
	 ENDIF
 	ENDIF
 35	CONTINUE
C--maximum virtualities for daughters
	BMAX1=MIN(P(L,5),Z*P(L,4))
      CMAX1=MIN(P(L,5),(1.-Z)*P(L,4))
C--generate mass of quark or gluon (particle b) from Sudakov FF
 30	IF(QUARK.OR.QQBAR)THEN
 	 MB=GETMASS(0.d0,BMAX1,THETA,Z*P(L,4),'QQ',
     &      BMAX1,.FALSE.,ZDECB,QQBARDECB)
	ELSE
 	 MB=GETMASS(0.d0,BMAX1,THETA,Z*P(L,4),'GC',
     &      BMAX1,.FALSE.,ZDECB,QQBARDECB)
 	ENDIF
C--generate mass gluon (particle c) from Sudakov FF
 	IF(QUARK.OR.(.NOT.QQBAR))THEN
       MC=GETMASS(0.d0,CMAX1,THETA,(1.-Z)*P(L,4),'GC',
     &	CMAX1,.FALSE.,ZDECC,QQBARDECC)
	ELSE
       MC=GETMASS(0.d0,CMAX1,THETA,(1.-Z)*P(L,4),'QQ',
     &	CMAX1,.FALSE.,ZDECC,QQBARDECC)
	ENDIF
C--quark (parton b) momentum
 182	PZ=(2.*Z*P(L,4)**2-P(L,5)**2-MB**2+MC**2)/(2.*P(L,3))
	PTS=Z**2*(P(L,4)**2)-PZ**2-MB**2
C--if kinematics doesn't work out, generate new virtualities
C     for daughters
C--massive phase space weight	
      IF((MB.EQ.0.d0).AND.(MC.EQ.0.d0).AND.(PTS.LT.0.d0)) GOTO 36
 	WEIGHT=1.d0
	IF((PYR(0).GT.WEIGHT).OR.(PTS.LT.0.d0)
     &	.OR.((MB+MC).GT.P(L,5)))THEN
	 IF(MB.GT.MC)THEN
 	  IF(QUARK.OR.QQBAR)THEN
 	   MB=GETMASS(0.d0,MB,THETA,Z*P(L,4),'QQ',
     &      BMAX1,.FALSE.,ZDECB,QQBARDECB)
	  ELSE
 	   MB=GETMASS(0.d0,MB,THETA,Z*P(L,4),'GC',
     &      BMAX1,.FALSE.,ZDECB,QQBARDECB)
 	  ENDIF
	 ELSE
 	  IF(QUARK.OR.(.NOT.QQBAR))THEN
         MC=GETMASS(0.d0,MC,THETA,(1.-Z)*P(L,4),'GC',
     &	CMAX1,.FALSE.,ZDECC,QQBARDECC)
	  ELSE
         MC=GETMASS(0.d0,MC,THETA,(1.-Z)*P(L,4),'QQ',
     &	CMAX1,.FALSE.,ZDECC,QQBARDECC)
	  ENDIF
	 ENDIF
	 GOTO 182
	ENDIF
	N=N+2
C--take care of first daughter (radiated gluon or antiquark)
	K(N-1,1)=K(L,1)
	IF(QQBAR)THEN
	 K(N-1,2)=-1
	 TRIP(N-1)=0
	 ANTI(N-1)=ANTI(L)
	ELSE
	 K(N-1,2)=21
	 IF((K(L,2).GT.0).AND.(DIR.GE.0))THEN
	  TRIP(N-1)=TRIP(L)
	  ANTI(N-1)=COLMAX+1
	 ELSE
	  TRIP(N-1)=COLMAX+1
	  ANTI(N-1)=ANTI(L)
	 ENDIF
	 COLMAX=COLMAX+1
	ENDIF
	K(N-1,3)=L
	K(N-1,4)=0
	K(N-1,5)=0
	P(N-1,4)=(1-Z)*P(L,4)
	P(N-1,5)=MC
	ZA(N-1)=1.-Z
	IF(ZDECC.GT.0.d0)THEN
	 THETAA(N-1)=P(N-1,5)/(SQRT(ZDECC*(1.-ZDECC))*P(N-1,4))
	ELSE
	 THETAA(N-1)=0.d0
	ENDIF 
	ZD(N-1)=ZDECC
	QQBARD(N-1)=QQBARDECC
C--take care of second daughter (final quark or gluon or quark from 
C	 gluon splitting)
	K(N,1)=K(L,1)
	IF(QUARK)THEN
	 K(N,2)=K(L,2)
	 IF(K(N,2).GT.0)THEN
	  TRIP(N)=ANTI(N-1)
	  ANTI(N)=0
	 ELSE
	  TRIP(N)=0
	  ANTI(N)=TRIP(N-1)
	 ENDIF
	ELSEIF(QQBAR)THEN
	 K(N,2)=1
	 TRIP(N)=TRIP(L)
	 ANTI(N)=0
	ELSE
	 K(N,2)=21
	 IF(DIR.EQ.1)THEN
	  TRIP(N)=ANTI(N-1)
	  ANTI(N)=ANTI(L)
	 ELSE
	  TRIP(N)=TRIP(L)
	  ANTI(N)=TRIP(N-1)
	 ENDIF
	ENDIF
	K(N,3)=L
	K(N,4)=0
	K(N,5)=0
	P(N,3)=PZ
	P(N,4)=Z*P(L,4)
	P(N,5)=MB
	ZA(N)=Z
	IF(ZDECB.GT.0.d0)THEN
	 THETAA(N)=P(N,5)/(SQRT(ZDECB*(1.-ZDECB))*P(N,4))
	ELSE 
	 THETAA(N)=0.d0
	ENDIF 
	ZD(N)=ZDECB
	QQBARD(N)=QQBARDECB
C--azimuthal angle
	PHIQ=2*PI*PYR(0)
	P(N,1)=SQRT(PTS)*COS(PHIQ)
	P(N,2)=SQRT(PTS)*SIN(PHIQ)
C--gluon momentum
	P(N-1,1)=P(L,1)-P(N,1)
	P(N-1,2)=P(L,2)-P(N,2)
	P(N-1,3)=P(L,3)-P(N,3)
      MV(N-1,4)=MV(L,5)
      IF(P(N-1,5).GT.0.d0)THEN
       LAMBDA=1.d0/(FTFAC*P(N-1,4)*0.2/P(N-1,5)**2)
	 MV(N-1,5)=MV(L,5)-LOG(1.d0-PYR(0))/LAMBDA
      ELSE
      MV(N-1,5)=0.d0
      ENDIF
      MV(N,4)=MV(L,5)
      IF(P(N,5).GT.0.d0)THEN
       LAMBDA=1.d0/(FTFAC*P(N,4)*0.2/P(N,5)**2)
	 MV(N,5)=MV(L,5)-LOG(1.d0-PYR(0))/LAMBDA
      ELSE
       MV(N,5)=0.d0
      ENDIF
C--take care of initial quark (or gluon)
      IF(K(L,1).EQ.2)THEN
       K(L,1)=13
      ELSE
	 K(L,1)=11
      ENDIF
	K(L,4)=N-1
	K(L,5)=N
	NSPLIT=NSPLIT+EVWEIGHT
 31	CONTINUE
 	END


***********************************************************************
***	  subroutine makeinsplit
***********************************************************************
	SUBROUTINE MAKEINSPLIT(L,X,TSUM,VIRT,TYPI,TIME,TAURAD)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--factor in front of formation times
	COMMON/FTIMEFAC/FTFAC
	DOUBLE PRECISION FTFAC
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--number of scattering events
	COMMON/CHECK/NSCAT,NSCATEFF,NSPLIT
	DOUBLE PRECISION NSCAT,NSCATEFF,NSPLIT
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights

C--local variables
	INTEGER L,TYPI,NOLD,DIR
	DOUBLE PRECISION X,VIRT,MB2,MC2,GETMASS,PZ,KT2,THETA,PHI,PI,
     &PHIQ,PYP,PYR,R,TIME,TSUM,TAURAD,LAMBDA,ZDEC
      LOGICAL QQBARDEC
	CHARACTER*2 TYP2,TYPC
	integer bin
	DATA PI/3.141592653589793d0/

      IF((N+2).GT.22990) THEN
       write(logfid,*)'event too long for event record'
       DISCARD=.TRUE.
       RETURN
      ENDIF

	IF(K(L,2).EQ.21)THEN
	 IF(TYPI.EQ.21)THEN
	  TYP2='GG'
	  TYPC='GC'
	 ELSE
	  TYP2='QG'
	  TYPC='QQ'
	 ENDIF
	ELSE
	 IF(TYPI.EQ.21)THEN
	  TYP2='GQ'
	  TYPC='QQ'
	 ELSE
	  TYP2='QQ'
	  TYPC='GC'
	 ENDIF
	ENDIF

C--if g->gg decide on colour configuration
	IF(TYP2.EQ.'GG')THEN
	 IF(PYR(0).LT.0.5)THEN
	  DIR=1
	 ELSE
	  DIR=-1
	 ENDIF
	ELSE
	 DIR=0
	ENDIF

	MB2=VIRT**2
	MB2=P(L,5)**2-MB2
	MC2=GETMASS(0.d0,SCALEFACM*SQRT(-TSUM),-1.d0,
     &	(1.-X)*P(L,4),TYPC,(1.-X)*P(L,4),
     &      .FALSE.,ZDEC,QQBARDEC)**2

C--rotate such that momentum points in z-direction
      NOLD=N
      THETA=PYP(L,13)
      PHI=PYP(L,15)
      CALL PYROBO(L,L,0d0,-PHI,0d0,0d0,0d0)
      CALL PYROBO(L,L,-THETA,0d0,0d0,0d0,0d0)
	PZ=(2*X*P(L,4)**2-P(L,5)**2-MB2+MC2)/(2*P(L,3))
	KT2=X**2*(P(L,4)**2)-PZ**2-MB2
	IF(KT2.LT.0.d0)THEN
	 MC2=0.d0
	 PZ=(2*X*P(L,4)**2-P(L,5)**2-MB2+MC2)/(2*P(L,3))
	 KT2=X**2*(P(L,4)**2)-PZ**2-MB2
	 IF(KT2.LT.0.d0)THEN
        CALL PYROBO(L,L,THETA,0d0,0d0,0d0,0d0)
        CALL PYROBO(L,L,0d0,PHI,0d0,0d0,0d0)
        X=1.d0
	  RETURN
	 ENDIF
	ENDIF	
	N=N+2
C--take care of first daughter (radiated gluon or antiquark)
	K(N-1,1)=K(L,1)
	IF(TYP2.EQ.'QG')THEN
	 K(N-1,2)=-TYPI
	 IF(K(N-1,2).GT.0)THEN
	  TRIP(N-1)=TRIP(L)
	  ANTI(N-1)=0
	 ELSE
	  TRIP(N-1)=0
	  ANTI(N-1)=ANTI(L)
	 ENDIF
	ELSEIF(TYP2.EQ.'GQ')THEN
	 K(N-1,2)=K(L,2)
       IF(K(N-1,2).GT.0)THEN
	  TRIP(N-1)=COLMAX+1
	  ANTI(N-1)=0
	 ELSE
	  TRIP(N-1)=0
	  ANTI(N-1)=COLMAX+1
	 ENDIF
	 COLMAX=COLMAX+1
	ELSE
	 K(N-1,2)=21
	 IF((K(L,2).GT.0).AND.(DIR.GE.0))THEN
	  TRIP(N-1)=TRIP(L)
	  ANTI(N-1)=COLMAX+1
	 ELSE
	  TRIP(N-1)=COLMAX+1
	  ANTI(N-1)=ANTI(L)
	 ENDIF
	 COLMAX=COLMAX+1
	ENDIF
	K(N-1,3)=L
	K(N-1,4)=0
	K(N-1,5)=0
	P(N-1,4)=(1.-X)*P(L,4)
	P(N-1,5)=SQRT(MC2)
C--take care of second daughter (final quark or gluon or quark from 
C	 gluon splitting)
	K(N,1)=K(L,1)
	IF(TYP2.EQ.'QG')THEN
	 K(N,2)=TYPI
	 IF(K(N,2).GT.0)THEN
	  TRIP(N)=TRIP(L)
	  ANTI(N)=0
	 ELSE
	  TRIP(N)=0
	  ANTI(N)=ANTI(L)
	 ENDIF
	ELSEIF(TYPI.NE.21)THEN
	 K(N,2)=K(L,2)
       IF(K(N,2).GT.0)THEN
	  TRIP(N)=ANTI(N-1)
	  ANTI(N)=0
	 ELSE
	  TRIP(N)=0
	  ANTI(N)=TRIP(N-1)
	 ENDIF
	ELSE
	 K(N,2)=21
	 IF(K(N-1,2).EQ.21)THEN
	  IF(DIR.EQ.1)THEN
	   TRIP(N)=ANTI(N-1)
	   ANTI(N)=ANTI(L)
	  ELSE
	   TRIP(N)=TRIP(L)
	   ANTI(N)=TRIP(N-1)
	  ENDIF
	 ELSEIF(K(N-1,2).GT.0)THEN
	  TRIP(N)=TRIP(L)
	  ANTI(N)=TRIP(N-1)
	 ELSE
	  TRIP(N)=ANTI(N-1)
	  ANTI(N)=ANTI(L)
	 ENDIF
	ENDIF
	K(N,3)=L
	K(N,4)=0
	K(N,5)=0
	P(N,3)=PZ
	P(N,4)=X*P(L,4)
	IF(MB2.LT.0.d0)THEN
	 P(N,5)=-SQRT(-MB2)
	ELSE
	 P(N,5)=SQRT(MB2)
	ENDIF
C--azimuthal angle
	PHIQ=2*PI*PYR(0)
	P(N,1)=SQRT(KT2)*COS(PHIQ)
	P(N,2)=SQRT(KT2)*SIN(PHIQ)
C--gluon momentum
	P(N-1,1)=P(L,1)-P(N,1)
	P(N-1,2)=P(L,2)-P(N,2)
	P(N-1,3)=P(L,3)-P(N,3)
	MV(L,5)=TIME-TAURAD
      MV(N-1,4)=MV(L,5)
      IF(P(N-1,5).GT.0.d0)THEN
       LAMBDA=1.d0/(FTFAC*P(N-1,4)*0.2/P(N-1,5)**2)
	 MV(N-1,5)=MV(L,5)-LOG(1.d0-PYR(0))/LAMBDA
      ELSE
       MV(N-1,5)=0.d0
      ENDIF
      MV(N,4)=MV(L,5)
      IF(P(N,5).GT.0.d0)THEN
	 MV(N,5)=TIME
      ELSE
       MV(N,5)=0.d0
      ENDIF
	ZA(N-1)=1.d0
      THETAA(N-1)=-1.d0
	ZD(N-1)=ZDEC
	QQBARD(N-1)=QQBARDEC
	ZA(N)=1.d0
	THETAA(N)=-1.d0
	ZD(N)=0.d0
	QQBARD(N)=.FALSE.
C--take care of initial quark (or gluon)
      IF(K(L,1).EQ.2)THEN
       K(L,1)=13
      ELSE
	 K(L,1)=11
      ENDIF
	K(L,4)=N-1
	K(L,5)=N
	NSPLIT=NSPLIT+EVWEIGHT
      CALL PYROBO(L,L,THETA,0d0,0d0,0d0,0d0)
      CALL PYROBO(N-1,N,THETA,0d0,0d0,0d0,0d0)
      CALL PYROBO(L,L,0d0,PHI,0d0,0d0,0d0)
      CALL PYROBO(N-1,N,0d0,PHI,0d0,0d0,0d0)

C--set the production vertices: x_mother + (tprod - tprod_mother) * beta_mother
      MV(N-1,1)=MV(L,1)+(MV(N-1,4)-MV(L,4))*P(L,1)/max(pyp(l,8),P(L,4))
      MV(N-1,2)=MV(L,2)+(MV(N-1,4)-MV(L,4))*P(L,2)/max(pyp(l,8),P(L,4))
      MV(N-1,3)=MV(L,3)+(MV(N-1,4)-MV(L,4))*P(L,3)/max(pyp(l,8),P(L,4))
      MV(N,  1)=MV(L,1)+(MV(N,  4)-MV(L,4))*P(L,1)/max(pyp(l,8),P(L,4))
      MV(N,  2)=MV(L,2)+(MV(N,  4)-MV(L,4))*P(L,2)/max(pyp(l,8),P(L,4))
      MV(N,  3)=MV(L,3)+(MV(N,  4)-MV(L,4))*P(L,3)/max(pyp(l,8),P(L,4))

	END


***********************************************************************
***	  subroutine doinstatescat
***********************************************************************
	SUBROUTINE DOINSTATESCAT(L,X,TYPI,Q,TSTART,DELTAT,OVERQ0,
     &				RETRYSPLIT)
	IMPLICIT NONE
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--factor in front of formation times
	COMMON/FTIMEFAC/FTFAC
	DOUBLE PRECISION FTFAC
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--variables for coherent scattering
	COMMON/COHERENT/NSTART,NEND,ALLQS(10000,6),SCATCENTRES(10000,10),
     &QSUMVEC(4),QSUM2
	INTEGER NSTART,NEND
	DOUBLE PRECISION ALLQS,SCATCENTRES,QSUMVEC,QSUM2
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--local variables
	INTEGER L,TYPI,COUNTER,COUNTMAX,COUNT2
	DOUBLE PRECISION X,DELTAT,DELTAL,PYR,R,PNORAD,GETPNORAD1,GETNOSCAT,
     &WEIGHT,LOW,FMAX,GETPDF,SIGMATOT,GETSSCAT,PFCHANGE,PI,TNOW,TLEFT,
     &XMAX,PQQ,PQG,PGQ,PGG,ALPHAS,TSTART,TSUM,Q,QOLD,Q2OLD,GETNEWMASS,
     &GENERATEZ,TMAX,TMAXNEW,DT,XSC,YSC,ZSC,TSC,MS1,MD1,GETMS,GETMD,
     &GETTEMP,GETNEFF,LAMBDA,RTAU,PHI,TAUEST,QSUMVECOLD(4),ZDUM,WEIGHT,
     &pyp
	LOGICAL FCHANGE,NORAD,OVERQ0,NOSCAT,GETDELTAT,RETRYSPLIT,
     &QQBARDUM	
	CHARACTER TYP
	CHARACTER*2 TYP2
	DATA PI/3.141592653589793d0/
	DATA COUNTMAX/10000/

	COUNTER=0
	
      XSC=MV(L,1)+(TSTART-MV(L,4))*P(L,1)/P(L,4)
      YSC=MV(L,2)+(TSTART-MV(L,4))*P(L,2)/P(L,4)
      ZSC=MV(L,3)+(TSTART-MV(L,4))*P(L,3)/P(L,4)
      TSC=TSTART
      MD1=GETMD(XSC,YSC,ZSC,TSC)
      MS1=GETMS(XSC,YSC,ZSC,TSC)

      IF(MD1.LE.1.D-4.OR.MS1.LE.1.D-4)THEN
       write(logfid,*)'problem!',GETTEMP(XSC,YSC,ZSC,TSC),
     &GETNEFF(XSC,YSC,ZSC,TSC)
      ENDIF

C--check for scattering
      NOSCAT=.NOT.GETDELTAT(L,TSTART,DELTAT,DT)
	IF(NOSCAT.AND.(.NOT.RETRYSPLIT)) GOTO 116

C--decide whether there will be radiation
	PNORAD=GETPNORAD1(L,xsc,ysc,zsc,tsc)
	IF((PYR(0).LT.PNORAD).OR.(P(L,4).LT.1.001*Q0))THEN
	 NORAD=.TRUE.
	ELSE
	 NORAD=.FALSE.
	ENDIF

C--decide whether q or g is to be scattered
      IF(K(L,2).EQ.21)THEN
       TYP='G'
       TYP2='GC'
	 SIGMATOT=GETSSCAT(P(L,4),p(l,1),p(l,2),p(l,3),P(L,5),
     &	Q0,'G','C',xsc,ysc,zsc,tsc,0)
	 IF((SIGMATOT.EQ.0.d0).OR.(PNORAD.EQ.1.d0))THEN
	  PFCHANGE=0.d0
	 ELSE
	  PFCHANGE=GETSSCAT(P(L,4),p(l,1),p(l,2),p(l,3),P(L,5),
     &	Q0,'G','Q',xsc,ysc,zsc,tsc,0)
     &	/SIGMATOT
	 ENDIF
	 SIGMATOT=GETSSCAT(P(L,4),p(l,1),p(l,2),p(l,3),P(L,5),
     &	0.d0,'G','C',xsc,ysc,zsc,tsc,0)
      ELSE
       TYP='Q'
       TYP2='QQ'
	 SIGMATOT=GETSSCAT(P(L,4),p(l,1),p(l,2),p(l,3),P(L,5),
     &	Q0,'Q','C',xsc,ysc,zsc,tsc,0)
	 IF((SIGMATOT.EQ.0.d0).OR.(PNORAD.EQ.1.d0))THEN
	  PFCHANGE=0.d0
	 ELSE
	  PFCHANGE=GETSSCAT(P(L,4),p(l,1),p(l,2),p(l,3),P(L,5),
     &	Q0,'Q','G',xsc,ysc,zsc,tsc,0)
     &	/SIGMATOT
	 ENDIF
	 SIGMATOT=GETSSCAT(P(L,4),p(l,1),p(l,2),p(l,3),P(L,5),
     &	0.d0,'Q','C',xsc,ysc,zsc,tsc,0)
      ENDIF
	IF((PFCHANGE.LT.-1.d-4).OR.(PFCHANGE.GT.1.d0+1.d-4)) THEN
      write(logfid,*)'error: flavour change probability=',
     &	PFCHANGE,'for ',TYP
	ENDIF
	IF(PYR(0).LT.PFCHANGE)THEN
	 FCHANGE=.TRUE.
	ELSE
	 FCHANGE=.FALSE.
	ENDIF
      IF (NORAD) FCHANGE=.FALSE.
C--set TYPI
	IF(TYP.EQ.'G')THEN
	 IF(FCHANGE)THEN
	  TYPI=INT(SIGN(2.d0,PYR(0)-0.5))
	 ELSE
	  TYPI=K(L,2)
	 ENDIF
	ELSE
	 IF(FCHANGE)THEN
	  TYPI=21
	 ELSE
	  TYPI=K(L,2)
	 ENDIF
	ENDIF
	LOW=Q0**2/SCALEFACM**2
	TMAX=4.*(P(L,4)**2-P(L,5)**2)
	XMAX=1.-Q0**2/(SCALEFACM**2*4.*TMAX)

	IF(SIGMATOT.EQ.0.d0) GOTO 116

	RTAU=PYR(0)

C--generate a trial emission
C--pick a x value from splitting function
 112	COUNTER=COUNTER+1
	IF(TYP.EQ.'G')THEN
	 IF(FCHANGE)THEN
	  X=GENERATEZ(0.d0,0.d0,1.-XMAX,'QG')
	 ELSE
	  X=GENERATEZ(0.d0,0.d0,1.-XMAX,'GG')
	 ENDIF
	ELSE
	 IF(FCHANGE)THEN
	  X=1.-GENERATEZ(0.d0,0.d0,1.-XMAX,'QQ')
	 ELSE
	  X=GENERATEZ(0.d0,0.d0,1.-XMAX,'QQ')
	 ENDIF
	ENDIF
      IF(NORAD) X=1.d0
C--initialisation
      TMAXNEW=(X*P(L,4))**2
	PHI=0.d0
	TLEFT=DELTAT
	TNOW=TSTART
	QSUMVEC(1)=0.d0
	QSUMVEC(2)=0.d0
	QSUMVEC(3)=0.d0
	QSUMVEC(4)=0.d0
	QSUM2=-1.d-10
	OVERQ0=.FALSE.
	Q=P(L,5)
	QOLD=P(L,5)
      TAUEST=DELTAT
C--generate first momentum transfer
	DELTAL=DT
	NSTART=1
	NEND=1
	TNOW=TNOW+DELTAL
	TSUM=DELTAL
	TLEFT=TLEFT-DELTAL
	ALLQS(NEND,6)=TNOW
	Q2OLD=QSUM2
C--get new momentum transfer
	COUNT2=0
 118	CALL GETQVEC(L,NEND,TNOW-MV(L,4),X)
	IF(-QSUM2.GT.P(L,4)**2)THEN
	 QSUMVEC(1)=0.d0
	 QSUMVEC(2)=0.d0
	 QSUMVEC(3)=0.d0
	 QSUMVEC(4)=0.d0
	 QSUM2=Q2OLD
	 IF(COUNT2.LT.100)THEN
	  COUNT2=COUNT2+1
	  GOTO 118
	 ELSE
	  ALLQS(NEND,1)=0.d0
	  ALLQS(NEND,2)=0.d0
	  ALLQS(NEND,3)=0.d0
	  ALLQS(NEND,4)=0.d0
	  ALLQS(NEND,5)=0.d0
	 ENDIF
	ENDIF
C--update OVERQ0
	IF(-ALLQS(NEND,1).GT.LOW) OVERQ0=.TRUE.
C--get new virtuality
	 IF(OVERQ0.AND.(.NOT.NORAD))THEN
	  Q=GETNEWMASS(L,SCALEFACM**2*QSUM2,SCALEFACM**2*Q2OLD,0.d0,
     &	  .TRUE.,X,ZDUM,QQBARDUM)
	 ELSE
	  Q=0.d0
	 ENDIF

C--estimate formation time
 111	IF((Q.EQ.0.d0).OR.(Q.EQ.P(L,5)))THEN
 	 TAUEST=DELTAT
	ELSE
 	 TAUEST=FTFAC*(1.-PHI)*0.2*X*P(L,4)/Q**2
	ENDIF
	LAMBDA=1.d0/TAUEST
	TAUEST=-LOG(1.d0-RTAU)/LAMBDA

C--find number, position and momentum transfers of further scatterings
	NOSCAT=.NOT.GETDELTAT(L,TNOW,MIN(TLEFT,TAUEST),DELTAL)
	IF((.NOT.NOSCAT).AND.(.NOT.NORAD))THEN
C--add a momentum transfer
	 NEND=NEND+1
	 IF(NEND.GE.100)THEN
	  nend=nend-1
	  goto 114
	 ENDIF
	 TNOW=TNOW+DELTAL
	 TSUM=TSUM+DELTAL
	 TLEFT=TLEFT-DELTAL
C--update phase
	 IF((Q.NE.0.d0).AND.(Q.NE.P(L,5)))THEN
	  PHI=PHI+5.*DELTAL*Q**2/(1.*X*P(L,4))
	 ENDIF
C--get new momentum transfer
	 ALLQS(NEND,6)=TNOW
	 Q2OLD=QSUM2
	 QSUMVECOLD(1)=QSUMVEC(1)
	 QSUMVECOLD(2)=QSUMVEC(2)
	 QSUMVECOLD(3)=QSUMVEC(3)
	 QSUMVECOLD(4)=QSUMVEC(4)
	 COUNT2=0
 119	 CALL GETQVEC(L,NEND,TNOW-MV(L,4),X)
	 IF(-QSUM2.GT.P(L,4)**2)THEN
	  QSUMVEC(1)=QSUMVECOLD(1)
	  QSUMVEC(2)=QSUMVECOLD(2)
	  QSUMVEC(3)=QSUMVECOLD(3)
	  QSUMVEC(4)=QSUMVECOLD(4)
	  QSUM2=Q2OLD
	  IF(COUNT2.LT.100)THEN
	   COUNT2=COUNT2+1
	   GOTO 119
	  ELSE
	   ALLQS(NEND,1)=0.d0
	   ALLQS(NEND,2)=0.d0
	   ALLQS(NEND,3)=0.d0
	   ALLQS(NEND,4)=0.d0
	   ALLQS(NEND,5)=0.d0
	  ENDIF
	 ENDIF
C--update OVERQ0
	 IF((-QSUM2.GT.LOW)
     &	.OR.(-ALLQS(NEND,1).GT.LOW)) OVERQ0=.TRUE.
C--get new virtuality
	 QOLD=Q
	 IF(OVERQ0.AND.(.NOT.NORAD))THEN
	  Q=GETNEWMASS(L,SCALEFACM**2*QSUM2,SCALEFACM**2*Q2OLD,0.d0,
     &	  .TRUE.,X,ZDUM,QQBARDUM)
	 ELSE
	  Q=0.d0
	 ENDIF
	 GOTO 111
	ENDIF

C--do reweighting
 114	TMAXNEW=X**2*P(L,4)**2
	IF(NORAD)THEN
	 WEIGHT=1.d0
	 Q=0.d0
	 X=1.d0
	ELSEIF((-QSUM2.LT.LOW).OR.(Q.EQ.0.d0))THEN
	 WEIGHT=0.d0
	ELSEIF(-QSUM2.GT.P(L,4)**2)THEN
	 WEIGHT=0.d0
	ELSE	 
	 IF(TYP.EQ.'G')THEN
 	  FMAX=2.*LOG(-SCALEFACM**2*QSUM2/Q0**2)
     & 	  *ALPHAS(Q0**2/4.,LPS)/(2.*PI)
	  IF(QSUM2.EQ.0.d0)THEN
	   WEIGHT=0.d0
	   NORAD=.TRUE.
	  ELSE
	   IF(FCHANGE)THEN
	    WEIGHT=2.*GETPDF(X,SCALEFACM*SQRT(-QSUM2),'QG')/(PQG(X)*FMAX)
	    IF((WEIGHT.GT.1.d0+1.d-4).OR.(WEIGHT.LT.-1.d-4))THEN
	      write(logfid,*)'x,sqrt(qsum^2),getpdf,fmax:',X,
     &	SQRT(-QSUM2),GETPDF(X,SCALEFACM*SQRT(-QSUM2),'QG'),'qg',
     &	FMAX
          ENDIF
	   ELSE
	    WEIGHT=GETPDF(X,SCALEFACM*SQRT(-QSUM2),'GG')/(PGG(X)*FMAX)
	    IF((WEIGHT.GT.1.d0+1.d-4).OR.(WEIGHT.LT.-1.d-4))THEN
	      write(logfid,*)'x,sqrt(qsum^2),getpdf,fmax:',X,
     &	SQRT(-QSUM2),GETPDF(X,SCALEFACM*SQRT(-QSUM2),'GG'),'gg',
     &	FMAX
          ENDIF
	   ENDIF
	  ENDIF
	 ELSE
 	  FMAX=LOG(-SCALEFACM**2*QSUM2/Q0**2)
     & 	  *ALPHAS(Q0**2/4.,LPS)/(2.*PI)
	  IF(QSUM2.EQ.0.d0)THEN
	   WEIGHT=0.d0
	   NORAD=.TRUE.
	  ELSE
	   IF(FCHANGE)THEN
	    WEIGHT=GETPDF(X,SCALEFACM*SQRT(-QSUM2),'GQ')/(PGQ(X)*FMAX)
	    IF((WEIGHT.GT.1.d0+1.d-4).OR.(WEIGHT.LT.-1.d-4))THEN
	     write(logfid,*)'x,sqrt(qsum^2),getpdf:,fmax',X,
     &	SQRT(-QSUM2),GETPDF(X,SCALEFACM*SQRT(-QSUM2),'GQ'),'gq',
     &	FMAX
          ENDIF
	   ELSE
	    WEIGHT=GETPDF(X,SCALEFACM*SQRT(-QSUM2),'QQ')/(PQQ(X)*FMAX)
	    IF((WEIGHT.GT.1.d0+1.d-4).OR.(WEIGHT.LT.-1.d-4))THEN
	     write(logfid,*)'x,sqrt(qsum^2),getpdf,fmax:',X,
     &	SQRT(-QSUM2),GETPDF(X,SCALEFACM*SQRT(-QSUM2),'QQ'),'qq',
     &	FMAX
          ENDIF
	   ENDIF
	  ENDIF
	 ENDIF
	ENDIF
	IF((WEIGHT.GT.1.d0+1.d-4).OR.(WEIGHT.LT.-1.d-4))
     &	write(logfid,*)'error: weight=',WEIGHT
 115	IF(PYR(0).GT.WEIGHT)THEN
	 IF(COUNTER.LT.COUNTMAX)THEN
	  GOTO 112
	 ELSE
	  Q=0.d0
	  X=1.d0
	  NEND=NSTART
	  QSUM2=ALLQS(NEND,1)
	  QSUMVEC(1)=ALLQS(NEND,2)
	  QSUMVEC(2)=ALLQS(NEND,3)
	  QSUMVEC(3)=ALLQS(NEND,4)
	  QSUMVEC(4)=ALLQS(NEND,5)
	  TYPI=K(L,2)
	  IF(-ALLQS(NEND,1).GT.LOW)THEN
	   OVERQ0=.TRUE.
	  ELSE
	   OVERQ0=.FALSE.
	  ENDIF
        DELTAT=ALLQS(NEND,6)-TSTART
	  TNOW=ALLQS(1,6)
	  RETURN
	 ENDIF
	ENDIF
C--found meaningful configuration, now do final checks
C--check if phase is unity and weight with 1/Nscat
      IF(((TLEFT.LT.TAUEST).OR.(PYR(0).GT.1.d0/(NEND*1.d0)))
     &			.AND.(.NOT.NORAD))THEN
	 Q=0.d0
	 X=1.d0
	 NEND=NSTART
	 QSUM2=ALLQS(NEND,1)
	 QSUMVEC(1)=ALLQS(NEND,2)
	 QSUMVEC(2)=ALLQS(NEND,3)
	 QSUMVEC(3)=ALLQS(NEND,4)
	 QSUMVEC(4)=ALLQS(NEND,5)
	 TYPI=K(L,2)
	 IF(-ALLQS(NEND,1).GT.LOW)THEN
	  OVERQ0=.TRUE.
	 ELSE
	  OVERQ0=.FALSE.
	 ENDIF
       DELTAT=ALLQS(NEND,6)-TSTART
	 TNOW=ALLQS(1,6)
	ELSE
       IF(.NOT.NORAD)THEN
	  TLEFT=TLEFT-TAUEST
	  TNOW=TNOW+TAUEST
	  TSUM=TSUM+TAUEST
	 ENDIF
       DELTAT=TSUM
	ENDIF
	RETURN
C--exit in case of failure
 116	Q=0.d0
	X=1.d0
	NSTART=0
	NEND=0
	QSUMVEC(1)=0.d0
	QSUMVEC(2)=0.d0
	QSUMVEC(3)=0.d0
	QSUMVEC(4)=0.d0
	QSUM2=0.d0
	OVERQ0=.FALSE.
	TYPI=K(L,2)
	RETURN
	END


***********************************************************************
***	  subroutine dofistatescat
***********************************************************************
	SUBROUTINE DOFISTATESCAT(L,TNOW,DTLEFT,DELTAT,NEWMASS,
     &		OVERQ0,Z,QQBAR)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--factor in front of formation times
	COMMON/FTIMEFAC/FTFAC
	DOUBLE PRECISION FTFAC
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--variables for coherent scattering
	COMMON/COHERENT/NSTART,NEND,ALLQS(10000,6),SCATCENTRES(10000,10),
     &QSUMVEC(4),QSUM2
	INTEGER NSTART,NEND
	DOUBLE PRECISION ALLQS,SCATCENTRES,QSUMVEC,QSUM2
C--local variables
	INTEGER L,COUNTER,COUNTMAX,COUNT2
	DOUBLE PRECISION TNOW,DELTAT,NEWMASS,TLEFT,DELTAL,Q2OLD,
     &GETNEWMASS,PYR,TSUM,QSUMVECOLD(4),RTAU,LAMBDA,DTLEFT,PHI,
     &TAUEST,LOW,Z,pyp
	LOGICAL OVERQ0,NOSCAT,GETDELTAT,QQBAR
	CHARACTER TYP
	DATA COUNTMAX/100/
	DELTAL=0.d0

	IF(-QSUM2.GT.P(L,4)**2)
     & write(logfid,*) 'DOFISTATESCAT has a problem:',-QSUM2,P(L,4)**2

      IF(K(L,2).EQ.21)THEN
       TYP='G'
	ELSE
	 TYP='Q'
	ENDIF
	LOW=Q0**2/SCALEFACM**2

	TSUM=0.d0
	PHI=0.d0
	DELTAT=0.d0

C--check for radiation with first (given) momentum transfer
	Q2OLD=0.d0
	IF(OVERQ0.OR.(-QSUM2.GT.LOW))THEN
	 NEWMASS=GETNEWMASS(L,SCALEFACM**2*QSUM2,SCALEFACM**2*Q2OLD,
     &	NEWMASS,.FALSE.,1.d0,Z,QQBAR)
	 OVERQ0=.TRUE.
	ELSE
	 NEWMASS=P(L,5)
	ENDIF

	RTAU=PYR(0)

	TLEFT=DTLEFT
 222	IF((NEWMASS.EQ.0.d0).OR.(NEWMASS.EQ.P(L,5)))THEN
 	 TAUEST=TLEFT
	ELSE
 	 TAUEST=FTFAC*(1.-PHI)*0.2*P(L,4)/NEWMASS**2
	ENDIF
	LAMBDA=1.d0/TAUEST
	TAUEST=-LOG(1.d0-RTAU)/LAMBDA
      NOSCAT=.NOT.GETDELTAT(L,TNOW+TSUM,MIN(TAUEST,TLEFT),DELTAL)
	IF(.NOT.NOSCAT)THEN
C--do scattering
	 NEND=NEND+1
	 IF(NEND.gt.countmax)THEN
	  nend=nend-1
	  goto 218
	 ENDIF
	 IF(NSTART.EQ.0) NSTART=1
	 TSUM=TSUM+DELTAL
	 TLEFT=TLEFT-DELTAL
	 IF((NEWMASS.NE.0.d0).AND.(NEWMASS.NE.P(L,5)))THEN
	  PHI=PHI+5.*DELTAL*NEWMASS**2/(1.*P(L,4))
	 ENDIF
	 ALLQS(NEND,6)=TNOW+TSUM
	 QSUMVECOLD(1)=QSUMVEC(1)
	 QSUMVECOLD(2)=QSUMVEC(2)
	 QSUMVECOLD(3)=QSUMVEC(3)
	 QSUMVECOLD(4)=QSUMVEC(4)
	 Q2OLD=QSUM2
C--get new momentum transfer
	 COUNT2=0
 219	 CALL GETQVEC(L,NEND,TNOW+TSUM-MV(L,4),1.d0)
	 IF(-QSUM2.GT.P(L,4)**2)THEN
	  QSUMVEC(1)=QSUMVECOLD(1)
	  QSUMVEC(2)=QSUMVECOLD(2)
	  QSUMVEC(3)=QSUMVECOLD(3)
	  QSUMVEC(4)=QSUMVECOLD(4)
	  QSUM2=Q2OLD
	  IF(COUNT2.LT.100)THEN
	   COUNT2=COUNT2+1
	   GOTO 219
	  ELSE
	   ALLQS(NEND,1)=0.d0
	   ALLQS(NEND,2)=0.d0
	   ALLQS(NEND,3)=0.d0
	   ALLQS(NEND,4)=0.d0
	   ALLQS(NEND,5)=0.d0
	  ENDIF
	 ENDIF
C--figure out new virtuality
	 IF(OVERQ0.OR.(-QSUM2.GT.LOW))THEN
	  NEWMASS=GETNEWMASS(L,SCALEFACM**2*QSUM2,SCALEFACM**2*Q2OLD,
     &	  NEWMASS,.FALSE.,1.d0,Z,QQBAR)
	  OVERQ0=.TRUE.
	 ENDIF
	 GOTO 222
	ENDIF
C--no more scattering
 218	if ((newmass**2.gt.low).and.(newmass.ne.p(l,5))) then
	  if ((TLEFT.LT.TAUEST).OR.(PYR(0).GT.1.d0/(NEND*1.d0))) then
	    if (nend.eq.countmax) then
	      deltat=tsum
	    else if (TLEFT.LT.TAUEST) then
	      DELTAT=TSUM+tleft
	    else
	      DELTAT=TSUM+tauest
	    endif
	    NEWMASS=P(L,5)
	  ELSE
	    DELTAT=TSUM+TAUEST
	  ENDIF
	else  
	  DELTAT=0.d0
	  NSTART=1
	  NEND=1
	  QSUM2=ALLQS(NEND,1)
	  QSUMVEC(1)=ALLQS(NEND,2)
	  QSUMVEC(2)=ALLQS(NEND,3)
	  QSUMVEC(3)=ALLQS(NEND,4)
	  QSUMVEC(4)=ALLQS(NEND,5)
	  IF(-ALLQS(NEND,1).GT.LOW)THEN
	    OVERQ0=.TRUE.
	  ELSE
	    OVERQ0=.FALSE.
	  ENDIF
	  NEWMASS=P(L,5)
	endif
	return
	END


***********************************************************************
***	  function getnewmass
***********************************************************************
	DOUBLE PRECISION FUNCTION GETNEWMASS(L,Q2,QOLD2,MASS,IN,X,
     &	ZDEC,QQBARDEC)
	IMPLICIT NONE
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	INTEGER L
	DOUBLE PRECISION Q2,QOLD2,R,PYR,PNOSPLIT1,PNOSPLIT2,Z,QA,
     &GETSUDAKOV,GETMASS,PKEEP,X,MASS,ZDEC,QTMP,ZOLD
	LOGICAL IN,QQBARDEC,QQBAROLD
	CHARACTER*2 TYP	

	IF(x*P(L,4).LT.Q0)THEN
	 GETNEWMASS=0.d0
	 ZDEC=0.d0
	 QQBARDEC=.FALSE.
	 RETURN
	ENDIF
	IF (-Q2.LT.Q0**2)THEN
	 GETNEWMASS=0.d0
	 RETURN
	ENDIF
      IF(K(L,2).EQ.21)THEN
       TYP='GC'
      ELSE
       TYP='QQ'
      ENDIF
	IF(SQRT(-QOLD2).LE.Q0)THEN
	 IF(IN)THEN
	  GETNEWMASS=GETMASS(0.d0,SQRT(-Q2),-1.d0,
     &	X*P(L,4),TYP,X*P(L,4),IN,ZDEC,QQBARDEC)
	 ELSE
	  GETNEWMASS=GETMASS(0.d0,SQRT(-Q2),-1.d0,P(L,4),TYP,
     &	SQRT(-Q2),IN,ZDEC,QQBARDEC)
	 ENDIF
	 GETNEWMASS=MIN(GETNEWMASS,X*P(L,4))
	 RETURN
	ENDIF
	Z=1.d0
	QA=1.d0	
	IF(MAX(P(L,5),MASS).GT.0.d0)THEN
	   IF(-Q2.GT.-QOLD2)THEN
	      ZOLD=ZDEC
	      QQBAROLD=QQBARDEC
	      QTMP=GETMASS(0.d0,SQRT(-Q2),-1.d0,X*P(L,4),TYP,
     &		SQRT(-Q2),IN,ZDEC,QQBARDEC)
	      IF(QTMP.LT.SQRT(-QOLD2))THEN
	        GETNEWMASS=MASS
	        ZDEC=ZOLD
              QQBARDEC=QQBAROLD
	      ELSE
	         GETNEWMASS=QTMP
	      ENDIF
	   ELSE
	     PNOSPLIT1=GETSUDAKOV(SQRT(-QOLD2),QA,Q0,Z,X*P(L,4),
     &      TYP,MV(L,4),IN)
	     PNOSPLIT2=GETSUDAKOV(SQRT(-Q2),QA,Q0,Z,X*P(L,4),
     &      TYP,MV(L,4),IN)
	     PKEEP=(1.-PNOSPLIT2)/(1.-PNOSPLIT1)
	     IF(PYR(0).LT.PKEEP)THEN
	       IF(P(L,5).LT.SQRT(-Q2))THEN
		   GETNEWMASS=MASS
		 ELSE
 55		   GETNEWMASS=GETMASS(Q0,SQRT(-Q2),-1.d0,X*P(L,4),TYP,
     &		SQRT(-Q2),IN,ZDEC,QQBARDEC)
		   IF((GETNEWMASS.EQ.0.d0).AND.(X*P(L,4).GT.Q0)) GOTO 55
		 ENDIF
	     ELSE
	       GETNEWMASS=0.d0
	       ZDEC=0.d0
	       QQBARDEC=.FALSE.
	     ENDIF
	   ENDIF
	 ELSE
	   IF(-Q2.GT.-QOLD2)THEN
	     GETNEWMASS=GETMASS(0.d0,SQRT(-Q2),-1.d0,
     &        X*P(L,4),TYP,X*P(L,4),IN,ZDEC,QQBARDEC)
           if(getnewmass.lt.SQRT(-QOLD2))then
	       GETNEWMASS=0.d0
	       ZDEC=0.d0
	       QQBARDEC=.FALSE.
           endif
	   ELSE
	     GETNEWMASS=0.d0
	     ZDEC=0.d0
	     QQBARDEC=.FALSE.
	   ENDIF
	 ENDIF
	 GETNEWMASS=MIN(GETNEWMASS,x*P(L,4))
	END	


***********************************************************************
***	  function getpnorad1
***********************************************************************
	DOUBLE PRECISION FUNCTION GETPNORAD1(LINE,x,y,z,t)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	INTEGER LINE
	DOUBLE PRECISION UP,LOW,CCOL,SIGMATOT,GETSSCAT,GETXSECINT,
     &SCATPRIMFUNC,MS1,MD1,shat,pcms2,avmom(5),x,y,z,t,getmd
	
	md1 = getmd(x,y,z,t)
	call avscatcen(x,y,z,t,
     &avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))
	ms1 = avmom(5)
	shat = avmom(5)**2 + p(line,5)**2 + 2.*(avmom(4)*p(line,4)
     &       -avmom(1)*p(line,1)-avmom(2)*p(line,2)-avmom(3)*p(line,3))
	pcms2 = (shat+p(line,5)**2-ms1**2)**2/(4.*shat)-p(line,5)**2
	up = 4.*pcms2
	 LOW=Q0**2/SCALEFACM**2
	 IF((UP.LE.LOW).OR.(P(LINE,4).LT.Q0/SCALEFACM))THEN
	  GETPNORAD1=1.d0
	  RETURN
	 ENDIF
	 IF(K(LINE,2).EQ.21)THEN
	  CCOL=3./2.
C--probability for no initial state radiation
	  SIGMATOT=GETSSCAT(P(LINE,4),p(line,1),p(line,2),p(line,3),
     &		P(LINE,5),0.d0,'G','C',x,y,z,t,0)
	  IF(SIGMATOT.EQ.0.d0)THEN
	   GETPNORAD1=-1.d0
	   RETURN
	  ENDIF
	   GETPNORAD1=(CCOL*(SCATPRIMFUNC(LOW,MD1)-
     &SCATPRIMFUNC(0.d0,MD1))
     &		+ GETXSECINT(UP,MD1,'GB'))/SIGMATOT
	 ELSE
	  CCOL=2./3.
C--probability for no initial state radiation
	  SIGMATOT=GETSSCAT(P(LINE,4),p(line,1),p(line,2),p(line,3),
     &		P(LINE,5),0.d0,'Q','C',x,y,z,t,0)
	  IF(SIGMATOT.EQ.0.d0)THEN
	   GETPNORAD1=1.d0
	   RETURN
	  ENDIF
	   GETPNORAD1=(CCOL*(SCATPRIMFUNC(LOW,MD1)-
     &SCATPRIMFUNC(0.d0,MD1))
     &		+ GETXSECINT(UP,MD1,'QB'))/SIGMATOT
	 ENDIF
	IF((GETPNORAD1.LT.-1.d-4).OR.(GETPNORAD1.GT.1.d0+1.d-4))THEN
       write(logfid,*)'error: P_norad=',GETPNORAD1,
     &	P(LINE,4),P(LINE,5),LOW,UP,K(LINE,2),MD1
	ENDIF
	END


***********************************************************************
***	  subroutine getqvec
***********************************************************************
	SUBROUTINE GETQVEC(L,J,DT,X)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--variables for coherent scattering
	COMMON/COHERENT/NSTART,NEND,ALLQS(10000,6),SCATCENTRES(10000,10),
     &QSUMVEC(4),QSUM2
	INTEGER NSTART,NEND
	DOUBLE PRECISION ALLQS,SCATCENTRES,QSUMVEC,QSUM2
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	INTEGER L,J,COUNTER,COUNTMAX,COUNT2,i
      DOUBLE PRECISION XSC,YSC,ZSC,TSC,GETMD,GETTEMP,DT,X,PYR,NEWMOM(4),
     &T,PT,MAXT,PHI2,BETA(3),PHI,THETA,GETT,PYP,PI,PT2,GETMS,
     &savemom(5),theta2,mb2,pz,kt2,phiq,maxt2,xi,md,shat,pcms2,
     &avmom(5)
	CHARACTER TYPS
	DATA PI/3.141592653589793d0/
	DATA COUNTMAX/1000/

      IF (J.GT.10000)THEN
       discard = .true.
	 return
      ENDIF

	COUNTER=0
	COUNT2=0

      XSC=MV(L,1)+DT*P(L,1)/P(L,4)
      YSC=MV(L,2)+DT*P(L,2)/P(L,4)
      ZSC=MV(L,3)+DT*P(L,3)/P(L,4)
      TSC=MV(L,4)+DT
	md = GETMD(XSC,YSC,ZSC,TSC)

	call AVSCATCEN(xsc,ysc,zsc,tsc,
     &avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))

	do 210 i=1,5
	  savemom(i) = p(l,i)
 210	continue

	xi = sqrt(max(x**2*p(l,4)**2,p(l,5)**2) - p(l,5)**2)/pyp(l,8)
	p(l,1) = xi*p(l,1)
	p(l,2) = xi*p(l,2)
	p(l,3) = xi*p(l,3)
	p(l,4) = max(x*p(l,4),p(l,5))


 444  CALL GETSCATTERER(XSC,YSC,ZSC,TSC,
     &K(1,2),P(1,1),P(1,2),P(1,3),P(1,4),P(1,5))
      MV(1,1)=XSC
      MV(1,2)=YSC
      MV(1,3)=ZSC
      MV(1,4)=TSC
      TYPS='Q'
      IF(K(1,2).EQ.21)TYPS='G'

	shat = avmom(5)**2 + savemom(5)**2 + 2.*(avmom(4)*savemom(4)
     &    -avmom(1)*savemom(1)-avmom(2)*savemom(2)-avmom(3)*savemom(3))
	pcms2 = (shat+savemom(5)**2-avmom(5)**2)**2/(4.*shat)
     &	-savemom(5)**2
	maxt = 4.*pcms2

      K(1,1)=13
	SCATCENTRES(J,1)=K(1,2)
	SCATCENTRES(J,2)=P(1,1)
	SCATCENTRES(J,3)=P(1,2)
	SCATCENTRES(J,4)=P(1,3)
	SCATCENTRES(J,5)=P(1,4)
	SCATCENTRES(J,6)=P(1,5)
	SCATCENTRES(J,7)=MV(1,1)
	SCATCENTRES(J,8)=MV(1,2)
	SCATCENTRES(J,9)=MV(1,3)
	SCATCENTRES(J,10)=MV(1,4)
C--transform to scattering centre's rest frame and rotate such that parton momentum is in z-direction
      BETA(1)=P(1,1)/P(1,4)
      BETA(2)=P(1,2)/P(1,4)
      BETA(3)=P(1,3)/P(1,4)
      CALL PYROBO(L,L,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
      CALL PYROBO(1,1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
      THETA=PYP(L,13)
      PHI=PYP(L,15)
      CALL PYROBO(L,L,0d0,-PHI,0d0,0d0,0d0)
      CALL PYROBO(1,1,0d0,-PHI,0d0,0d0,0d0)
      CALL PYROBO(L,L,-THETA,0d0,0d0,0d0,0d0)
      CALL PYROBO(1,1,-THETA,0d0,0d0,0d0,0d0)
C--pick a t from differential scattering cross section
 204  T=-GETT(0.d0,MAXT,md)
 202	NEWMOM(4)=P(L,4)+T/(2.*p(1,5))
	NEWMOM(3)=(T-2.*P(L,5)**2+2.*p(l,4)*NEWMOM(4))/(2.*P(L,3))
	PT2=NEWMOM(4)**2-NEWMOM(3)**2-P(L,5)**2
	IF(DABS(PT2).LT.1.d-10) PT2=0.d0	
	IF(T.EQ.0.d0) PT2=0.d0
	IF(PT2.LT.0.d0)THEN
	 T=0.d0
	 GOTO 202
	ENDIF
	PT=SQRT(PT2)
      PHI2=PYR(0)*2*PI
	NEWMOM(1)=PT*COS(PHI2)
	NEWMOM(2)=PT*SIN(PHI2)
	P(1,1)=NEWMOM(1)-P(L,1)
	P(1,2)=NEWMOM(2)-P(L,2)
	P(1,3)=NEWMOM(3)-P(L,3)
	P(1,4)=NEWMOM(4)-P(L,4)
	P(1,5)=0.d0
C--transformation to lab
      CALL PYROBO(L,L,THETA,0d0,0d0,0d0,0d0)
      CALL PYROBO(1,1,THETA,0d0,0d0,0d0,0d0)
      CALL PYROBO(L,L,0d0,PHI,0d0,0d0,0d0)
      CALL PYROBO(1,1,0d0,PHI,0d0,0d0,0d0)
      CALL PYROBO(L,L,0d0,0d0,BETA(1),BETA(2),BETA(3))
      CALL PYROBO(1,1,0d0,0d0,BETA(1),BETA(2),BETA(3))
	ALLQS(J,1)=T
	ALLQS(J,2)=P(1,1)
	ALLQS(J,3)=P(1,2)
	ALLQS(J,4)=P(1,3)
	ALLQS(J,5)=P(1,4)
	QSUMVEC(1)=QSUMVEC(1)+ALLQS(NEND,2)
	QSUMVEC(2)=QSUMVEC(2)+ALLQS(NEND,3)
	QSUMVEC(3)=QSUMVEC(3)+ALLQS(NEND,4)
	QSUMVEC(4)=QSUMVEC(4)+ALLQS(NEND,5)
	QSUM2=QSUMVEC(4)**2-QSUMVEC(1)**2-QSUMVEC(2)**2-QSUMVEC(3)**2
	IF(QSUM2.GT.0.d0)THEN
	 QSUMVEC(1)=QSUMVEC(1)-ALLQS(NEND,2)
	 QSUMVEC(2)=QSUMVEC(2)-ALLQS(NEND,3)
	 QSUMVEC(3)=QSUMVEC(3)-ALLQS(NEND,4)
	 QSUMVEC(4)=QSUMVEC(4)-ALLQS(NEND,5)
	 QSUM2=QSUMVEC(4)**2-QSUMVEC(1)**2-QSUMVEC(2)**2-QSUMVEC(3)**2
	 IF(COUNTER.GT.COUNTMAX)THEN
	  write(logfid,*)'GETQVEC unable to find q vector'
	  ALLQS(J,1)=0.d0
	  ALLQS(J,2)=0.d0
	  ALLQS(J,3)=0.d0
	  ALLQS(J,4)=0.d0
	  ALLQS(J,5)=0.d0
	 ELSE
	  COUNTER=COUNTER+1
	  GOTO 444
	 ENDIF
	ENDIF
	do 211 i=1,5
	  p(l,i) = savemom(i)
 211	continue
	END

***********************************************************************
***	  subroutine dokinematics
***********************************************************************
      SUBROUTINE DOKINEMATICS(L,lold,N1,N2,NEWM,RETRYSPLIT,
     &	TIME,X,Z,QQBAR)
      IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--factor in front of formation times
	COMMON/FTIMEFAC/FTFAC
	DOUBLE PRECISION FTFAC
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--discard event flag
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--variables for coherent scattering
	COMMON/COHERENT/NSTART,NEND,ALLQS(10000,6),SCATCENTRES(10000,10),
     &QSUMVEC(4),QSUM2
	INTEGER NSTART,NEND
	DOUBLE PRECISION ALLQS,SCATCENTRES,QSUMVEC,QSUM2
C--number of scattering events
	COMMON/CHECK/NSCAT,NSCATEFF,NSPLIT
	DOUBLE PRECISION NSCAT,NSCATEFF,NSPLIT
C--event weight
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
C--extra storage for scattering centres before interactions
      common/storescatcen/nscatcen,maxnscatcen,scatflav(10000),
     &scatcen(10000,5),writescatcen,writedummies
	integer nscatcen,maxnscatcen,scatflav
	double precision scatcen
	logical writescatcen,writedummies
C--local variables
      INTEGER L,LINE,N1,N2,J,DIR,lold,nold,colmaxold,statold,nscatcenold
      DOUBLE PRECISION PYR,PI,BETA(3),THETA,PHI,PYP,PHI2,MAXT,T,
     &NEWMASS,DELTAM,DM,TTOT,DMLEFT,LAMBDA,TIME,ENDTIME,X,tmp,
     &m32,newm2,shat,theta2,z,gettemp,E3new,E4new,p32,p42,p3old,
     &newm,mass2,enew,pt2,pt,pl,m12,firsttime,pcms2
      CHARACTER*2 TYP
	LOGICAL RETRYSPLIT,QQBAR,QQBARDEC,rejectt,redokin,reshuffle
	DATA PI/3.141592653589793d0/

      IF((N+2*(n2-n1+1)).GT.22990)THEN
        write(logfid,*)'event too long for event record'
        DISCARD=.TRUE.
        RETURN
      ENDIF

	firsttime = mv(l,5)

	redokin = .false.

	newm2=newm
	nold=n
	colmaxold=colmax
	statold=k(l,1)
 204	DELTAM=NEWM2-P(L,5)
 	DMLEFT=DELTAM

	TTOT=0.d0
	DO 220 J=N1,N2
	 TTOT=TTOT+ALLQS(J,1)
 220  CONTINUE

	LINE=L

	DO 222 J=N1,N2
	
C--projectile type
	 IF(K(LINE,2).EQ.21)THEN
	  TYP='GC'
	  IF(PYR(0).LT.0.5)THEN
	   DIR=1
	  ELSE
	   DIR=-1
	  ENDIF
	 ELSE
	  TYP='QQ'
	  DIR=0
	 ENDIF
       K(1,1)=6
	 K(1,2)=SCATCENTRES(J,1)
	 P(1,1)=SCATCENTRES(J,2)
	 P(1,2)=SCATCENTRES(J,3)
	 P(1,3)=SCATCENTRES(J,4)
	 P(1,4)=SCATCENTRES(J,5)
	 P(1,5)=SCATCENTRES(J,6)
       MV(1,1)=SCATCENTRES(J,7)
       MV(1,2)=SCATCENTRES(J,8)
       MV(1,3)=SCATCENTRES(J,9)
       MV(1,4)=SCATCENTRES(J,10)
	 T=ALLQS(J,1)
	 if (t.eq.0.d0) then
	   rejectt = .true.
	 else 
	   rejectt = .false.
	 endif

C--transform to c.m.s. and rotate such that parton momentum is in z-direction
       BETA(1)=(P(1,1)+p(line,1))/(P(1,4)+p(line,4))
       BETA(2)=(P(1,2)+p(line,2))/(P(1,4)+p(line,4))
       BETA(3)=(P(1,3)+p(line,3))/(P(1,4)+p(line,4))
       IF ((BETA(1).GT.1.d0).OR.(BETA(2).GT.1.d0).OR.(BETA(3).GT.1.d0)
     &	.or.(sqrt(beta(1)**2+beta(2)**2+beta(3)**2).gt.1.d0))THEN
	   reshuffle = .false.
	 else 
	   reshuffle = .true.
	 endif
 205	 if (.not.reshuffle) then
         BETA(1)=P(1,1)/P(1,4)
         BETA(2)=P(1,2)/P(1,4)
         BETA(3)=P(1,3)/P(1,4)
         CALL PYROBO(LINE,LINE,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
         CALL PYROBO(1,1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
         THETA=PYP(LINE,13)
         PHI=PYP(LINE,15)
         CALL PYROBO(LINE,LINE,0d0,-PHI,0d0,0d0,0d0)
         CALL PYROBO(1,1,0d0,-PHI,0d0,0d0,0d0)
         CALL PYROBO(LINE,LINE,-THETA,0d0,0d0,0d0,0d0)
         CALL PYROBO(1,1,-THETA,0d0,0d0,0d0,0d0)

	   maxt = -2.*p(1,5)*p(line,4)
	   if (t.lt.maxt) then
	     t=0.d0
	     rejectt = .true.
	   endif
	   m12 = -p(line,5)**2
 203	   enew = p(line,4)+t/(2.*p(1,5))
	   pl = (t+2.*p(line,4)*enew-2.*m12)/(2.*p(line,3))
	   pt2 = enew**2-pl**2-m12
	   if (t.eq.0.d0) pt2 = 0.d0
	   if (dabs(pt2).lt.1.d-8) pt2 = 0.d0
	   if (pt2.lt.0.d0) then
	     write(logfid,*)' This should not have happened: pt^2<0!'
	     write(logfid,*)t,enew,pl,pt2
	     t = 0.d0
	     rejectt = .true.
	     goto 203
	   endif
	   pt = sqrt(pt2)
	   phi2 = pyr(0)*2.*pi
	   n=n+2
	   p(n,1)=pt*cos(phi2)
	   p(n,2)=pt*sin(phi2)
	   p(n,3)=pl
	   p(n,4)=enew
	   p(n,5)=p(line,5)
!---------------------------------       
         P(N-1,1)=P(1,1)+P(LINE,1)-P(N,1)
         P(N-1,2)=P(1,2)+P(LINE,2)-P(N,2)
         P(N-1,3)=P(1,3)+P(LINE,3)-P(N,3)
         P(N-1,4)=P(1,4)+P(LINE,4)-P(N,4)
	   mass2 = P(N-1,4)**2-P(N-1,1)**2-P(N-1,2)**2-P(N-1,3)**2
	   if ((mass2.lt.0.d0).and.(mass2.gt.-1.-6))  mass2=0.d0
         if (mass2.lt.0.d0)  
     &	write(logfid,*)'messed up scattering centres mass^2: ',
     &	mass2,p(1,5)**2
         P(N-1,5)=SQRT(mass2)
	   if (abs(p(n-1,5)-p(1,5)).gt.1.d-6)
     &	write(logfid,*)'messed up scattering centres mass: ',
     &	p(n-1,5),p(1,5),p(l,5)
	   call flush(logfid)
!---------------------------------       
!        P(N-1,1)=P(1,1)
!        P(N-1,2)=P(1,2)
!        P(N-1,3)=P(1,3)
!        P(N-1,4)=P(1,4)
!        P(N-1,5)=P(1,5)
!---------------------------------       
	 else 
         CALL PYROBO(LINE,LINE,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
         CALL PYROBO(1,1,0d0,0d0,-BETA(1),-BETA(2),-BETA(3))
	   if ((p(1,4).lt.0.d0).or.(p(line,4).lt.0.d0)) then
           CALL PYROBO(1,1,0d0,0d0,BETA(1),BETA(2),BETA(3))
           CALL PYROBO(LINE,LINE,0d0,0d0,BETA(1),BETA(2),BETA(3))
	     reshuffle = .false.
	     goto 205
	   endif
         THETA=PYP(LINE,13)
         PHI=PYP(LINE,15)
         CALL PYROBO(LINE,LINE,0d0,-PHI,0d0,0d0,0d0)
         CALL PYROBO(1,1,0d0,-PHI,0d0,0d0,0d0)
         CALL PYROBO(LINE,LINE,-THETA,0d0,0d0,0d0,0d0)
         CALL PYROBO(1,1,-THETA,0d0,0d0,0d0,0d0)
	   shat = (p(1,4)+p(line,4))**2
	   p3old = p(line,3)

	   maxt = -4.*p(line,3)**2
	   if (t.lt.maxt) then
	     t=0.d0
	     rejectt = .true.
	   endif
	   theta2 = acos(1.d0+t/(2.*p(line,3)**2))
	   phi2 = pyr(0)*2.*pi
	   n=n+2
	   p(n,1)=p(line,3)*sin(theta2)*cos(phi2)
	   p(n,2)=p(line,3)*sin(theta2)*sin(phi2)
	   p(n,3)=p(line,3)*cos(theta2)
	   p(n,4)=p(line,4)
	   p(n,5)=p(line,5)
!---------------------------------       
         P(N-1,1)=P(1,1)+P(LINE,1)-P(N,1)
         P(N-1,2)=P(1,2)+P(LINE,2)-P(N,2)
         P(N-1,3)=P(1,3)+P(LINE,3)-P(N,3)
         P(N-1,4)=P(1,4)+P(LINE,4)-P(N,4)
	   mass2 = P(N-1,4)**2-P(N-1,1)**2-P(N-1,2)**2-P(N-1,3)**2
	   if ((mass2.lt.0.d0).and.(mass2.gt.-1.-6))  mass2=0.d0
         if (mass2.lt.0.d0)  
     &	write(logfid,*)'messed up scattering centres mass^2: ',
     &	mass2,p(1,5)**2
         P(N-1,5)=SQRT(mass2)
	   if (abs(p(n-1,5)-p(1,5)).gt.1.d-6)
     &	write(logfid,*)'messed up scattering centres mass: ',
     &	p(n-1,5),p(1,5),p(l,5)
	   call flush(logfid)
!---------------------------------       
!        P(N-1,1)=P(1,1)
!        P(N-1,2)=P(1,2)
!        P(N-1,3)=P(1,3)
!        P(N-1,4)=P(1,4)
!        P(N-1,5)=P(1,5)
!---------------------------------       
	 endif
C--outgoing projectile
       ZA(N)=1.d0
	 THETAA(N)=-1.d0
       ZD(N)=Z
       QQBARD(N)=QQBAR
       K(N,1)=K(LINE,1)
       K(N,2)=K(LINE,2)
	 K(N,3)=L
	 K(N,4)=0
	 K(N,5)=0
	 IF(ALLHAD.and.(.not.rejectt))THEN
	  IF(K(N,2).EQ.21)THEN
	   IF(DIR.EQ.1)THEN
	    TRIP(N)=COLMAX+1
	    ANTI(N)=ANTI(LINE)
	   ELSE
	    TRIP(N)=TRIP(LINE)
	    ANTI(N)=COLMAX+1
	   ENDIF
	  ELSEIF(K(N,2).GT.0)THEN
	   TRIP(N)=COLMAX+1	
	   ANTI(N)=0
	  ELSE
	   TRIP(N)=0
	   ANTI(N)=COLMAX+1
	  ENDIF
	  COLMAX=COLMAX+1
	 ELSE
	  TRIP(N)=TRIP(LINE)
	  ANTI(N)=ANTI(LINE)
	 ENDIF
C--take care of incoming projectile
       IF(K(LINE,1).EQ.1)THEN
	  K(LINE,1)=12
       ELSE
        K(LINE,1)=14
       ENDIF
	 K(LINE,4)=N-1
	 K(LINE,5)=N
C--outgoing scattering centre
       ZA(N-1)=1.d0
	 THETAA(N-1)=-1.d0
       ZD(N-1)=-1.d0
       QQBARD(N-1)=.false.
C--temporary status code, will be overwritten later
       K(N-1,1)=3
	 K(N-1,2)=21
	 K(N-1,3)=0
	 K(N-1,4)=0
	 K(N-1,5)=0
	 IF(ALLHAD.and.(.not.rejectt))THEN
	  IF((K(N,2).GT.0).AND.(DIR.GE.0))THEN
	   TRIP(N-1)=TRIP(LINE)
	   ANTI(N-1)=TRIP(N)
	  ELSE
	   TRIP(N-1)=ANTI(N)
	   ANTI(N-1)=ANTI(LINE)
	  ENDIF
	 ELSE
	  TRIP(N-1)=0
	  ANTI(N-1)=0
	 ENDIF

	 if (reshuffle.and.(dm.gt.0.d0)) then
C--adjust mass and re-shuffle momenta

	   IF(TTOT.EQ.0.d0)THEN
	    DM=0.d0
	   ELSE
	    if (dmleft.lt.0.d0) then
	      DM=max(DMLEFT*T/TTOT*1.5d0,dmleft)
	    else
	      DM=min(DMLEFT*T/TTOT*1.5d0,dmleft)
	    endif
 	   ENDIF
	   TTOT=TTOT-ALLQS(J,1)

	   newmass = p(n,5)+dm
	   if (newmass.lt.0.d0) then
	     m32 = -NEWMASS**2
	   else
	     m32 = NEWMASS**2
	   endif
	   E3new = (shat + m32 - p(1,5)**2)/(2.d0*sqrt(shat))
	   E4new = (shat - m32 + p(1,5)**2)/(2.d0*sqrt(shat))
	   p32 = E3new**2 - m32
	   p42 = E4new**2 - p(1,5)**2
	   if ((p32.lt.0.d0).or.(p42.lt.0.d0).or.
     &       (E3new.lt.0.d0).or.(E4new.lt.0.d0)) then
	     p32 = 0.d0
	     p42 = 0.d0
	     E4new = p(n-1,5)
	     E3new = sqrt(shat) - E4new
	     m32 = E3new**2
	     if ((E3new.lt.0.d0).or.(E4new.lt.0.d0)) then
	       E3new = p(n,4)
	       E4new = p(n-1,4)
	       p32 = p3old**2
	       p42 = p3old**2
	   	 if (p(n,5).lt.0.d0) then
	     	   m32 = -p(n,5)**2
	   	 else
	     	   m32 = p(n,5)**2
	   	 endif 
	     endif
	   endif
	   p(n,1) = sqrt(p32)*p(n,1)/p3old
	   p(n,2) = sqrt(p32)*p(n,2)/p3old
	   p(n,3) = sqrt(p32)*p(n,3)/p3old
	   p(n,4) = E3new
	   p(n,5) = sign(sqrt(abs(m32)),newmass)
	   tmp = p(n,4)**2-p(n,1)**2-p(n,2)**2-p(n,3)**2
	   if (abs(tmp-m32).gt.1.d-6) 
     &	write(logfid,*) 'Oups, messed up projectiles mass:',
     &	tmp,m32,p(n,5)
!---------------------------------       
	   p(n-1,1) = sqrt(p42)*p(n-1,1)/p3old
	   p(n-1,2) = sqrt(p42)*p(n-1,2)/p3old
	   p(n-1,3) = sqrt(p42)*p(n-1,3)/p3old
	   p(n-1,4) = E4new
	   tmp = p(n-1,4)**2-p(n-1,1)**2-p(n-1,2)**2-p(n-1,3)**2
     &	-p(n-1,5)**2
	   if (abs(tmp).gt.1.d-6) 
     &	write(logfid,*) 'Oups, messed up scattering centres mass:',
     &	tmp,p3old,p(n-1,1),p(n-1,2),p(n-1,3),p(n-1,4),p(n-1,5)
	   if ((abs(p(n,1)+p(n-1,1)).gt.1.d-6).or.
     &     (abs(p(n,2)+p(n-1,2)).gt.1.d-6).or.
     &     (abs(p(n,3)+p(n-1,3)).gt.1.d-6)) 
     &	write(logfid,*) 'Oups, momentum not conserved', 
     &	p(n,1)+p(n-1,1),p(n,2)+p(n-1,2),p(n,3)+p(n-1,3)
!---------------------------------       
!        P(N-1,1)=P(1,1)
!        P(N-1,2)=P(1,2)
!        P(N-1,3)=P(1,3)
!        P(N-1,4)=P(1,4)
!        P(N-1,5)=P(1,5)
!---------------------------------       
	 endif

C--transformation to lab
       CALL PYROBO(N-1,N,THETA,0d0,0d0,0d0,0d0)
       CALL PYROBO(LINE,LINE,THETA,0d0,0d0,0d0,0d0)
       CALL PYROBO(N-1,N,0d0,PHI,0d0,0d0,0d0)
       CALL PYROBO(LINE,LINE,0d0,PHI,0d0,0d0,0d0)
       CALL PYROBO(N-1,N,0d0,0d0,BETA(1),BETA(2),BETA(3))
       CALL PYROBO(LINE,LINE,0d0,0d0,BETA(1),BETA(2),BETA(3))
       CALL PYROBO(1,1,THETA,0d0,0d0,0d0,0d0)
       CALL PYROBO(1,1,0d0,PHI,0d0,0d0,0d0)
       CALL PYROBO(1,1,0d0,0d0,BETA(1),BETA(2),BETA(3))
      if (.not.allhad) then
	  k(n-1,1)=13
	 else
        IF(SCATRECOIL.AND.(P(N-1,4).GT.(10.*3.*
     &GETTEMP(MV(1,1),MV(1,2),MV(1,3),MV(1,4)))))THEN
         K(N-1,1)=2
        ELSE
         K(N-1,1)=3
        ENDIF
	 endif
	 if (rejectt) k(n-1,1)=11
       MV(N,4)=MV(1,4)
       MV(N-1,4)=MV(1,4)
C--set the production vertices: x_mother + (tprod - tprod_mother) * beta_mother
       MV(N-1,1)=MV(line,1)
     &	+(MV(N-1,4)-MV(line,4))*P(line,1)/max(pyp(line,8),P(line,4))
       MV(N-1,2)=MV(line,2)
     &	+(MV(N-1,4)-MV(line,4))*P(line,2)/max(pyp(line,8),P(line,4))
       MV(N-1,3)=MV(line,3)
     &	+(MV(N-1,4)-MV(line,4))*P(line,3)/max(pyp(line,8),P(line,4))
       MV(N,  1)=MV(line,1)
     &	+(MV(N,  4)-MV(line,4))*P(line,1)/max(pyp(line,8),P(line,4))
       MV(N,  2)=MV(line,2)
     &	+(MV(N,  4)-MV(line,4))*P(line,2)/max(pyp(line,8),P(line,4))
       MV(N,  3)=MV(line,3)
     &	+(MV(N,  4)-MV(line,4))*P(line,3)/max(pyp(line,8),P(line,4))
	 IF(P(N-1,5).GT.P(1,5))THEN
	   LAMBDA=1.d0/(FTFAC*0.2*P(N-1,4)/P(N-1,5)**2)
	   MV(N-1,5)=MV(N-1,4)-LOG(1.d0-PYR(0))/LAMBDA
	 ELSE
        MV(N-1,5)=0.d0
	 ENDIF
	 IF(J.LT.N2)THEN
        MV(N,5)=SCATCENTRES(J+1,10)
	 ELSE
	  IF(P(N,5).GT.0.d0)THEN
	   IF(DELTAM.EQ.0.d0)THEN
	    ENDTIME=firsttime
	   ELSE
	    IF(X.LT.1.d0)THEN
           LAMBDA=1.d0/(FTFAC*P(N,4)*0.2/P(N,5)**2)
	     ENDTIME=SCATCENTRES(J,10)-LOG(1.d0-PYR(0))/LAMBDA
	    ELSE
	     ENDTIME=TIME
	    ENDIF
	   ENDIF
	   MV(N,5)=ENDTIME
	  ELSE
         MV(N,5)=0.d0
	  ENDIF
	 ENDIF
	 MV(LINE,5)=ALLQS(J,6)


C--store scattering centre before interaction in separate common block
	 if (writescatcen.and.(.not.rejectt).and.
     &		(nscatcen.lt.maxnscatcen)) then
	  nscatcen = nscatcen+1
	  if (nscatcen.le.maxnscatcen) then
	   scatflav(nscatcen) = k(1,2)
	   scatcen(nscatcen,1) = p(1,1)
	   scatcen(nscatcen,2) = p(1,2)
	   scatcen(nscatcen,3) = p(1,3)
	   scatcen(nscatcen,4) = p(1,4)
	   scatcen(nscatcen,5) = p(1,5)
	  else
	   write(logfid,*) 
     &'WARNING: no room left to store further scattering centres'
	  endif
	 endif

!	if ((p(line,4).gt.100.d0).and.(p(n,4)-p(line,4).gt.1.d0)) then
!	  write(*,*)p(line,1),p(line,2),p(line,3),p(line,4),p(line,5)
!	  write(*,*)p(n,1),p(n,2),p(n,3),p(n,4),p(n,5)
!	  write(*,*)p(1,1),p(1,2),p(1,3),p(1,4),p(1,5)
!	  write(*,*)p(n-1,1),p(n-1,2),p(n-1,3),p(n-1,4),p(n-1,5)
!	  write(*,*)t
!	  write(*,*)GETTEMP(MV(1,1),MV(1,2),MV(1,3),MV(1,4))
!	  write(*,*)
!	endif

	 DMLEFT=DMLEFT-(p(n,5)-P(LINE,5))
	 LINE=N
	 tmp = abs(p(n,4)**2-p(n,1)**2-p(n,2)**2-p(n,3)**2)-p(n,5)**2
	 if (abs(tmp).ge.1.d-6) 
     &	write(logfid,*)tmp,j,p(l,5),p(line,5),p(n,5)
 222	CONTINUE
	if (p(n,5).lt.0.d0) then
	  RETRYSPLIT=.TRUE.
	  return
	endif
	if (p(n,5).ne.newm2) then
	  RETRYSPLIT=.TRUE.
	  redokin = .true.
	  n=nold
	  colmax=colmaxold
	  k(l,1)=statold
	  if (p(l,5).le.0.d0) then
	    newm2 = 0.d0
	  else
          if (p(l,5).lt.q0) then
            if ((newm2.eq.newm).and.(newm.ne.q0+1.d-6)) then
              newm2=q0+1.d-6
            else    
              RETRYSPLIT=.TRUE.
              return
            endif
          else
            newm2=p(l,5)
          endif
          n2=n1
        endif
	  goto 204
	endif
	if ((k(n,1).eq.1).and.
     &	((p(n,5).lt.0.d0).or.((p(n,5).gt.0.d0).and.(p(n,5).lt.q0))))
     &write(logfid,*)'dokinematics did not reach sensible mass: ',
     &p(n,5),newm,p(l,5),newm2
	NSCATEFF=NSCATEFF+EVWEIGHT
      END



***********************************************************************
***	  function getproba
***********************************************************************
	DOUBLE PRECISION FUNCTION GETPROBA(QI,QF,QAA,ZAA,EBB,TYPE,
     &	T1,INS2)
	IMPLICIT NONE
C--variables for Sudakov integration
	COMMON/SUDAINT/QA,ZA2,EB,T,INSTATE,TYP
	DOUBLE PRECISION QA,ZA2,EB,T
	CHARACTER*2 TYP
	LOGICAL INSTATE
C--local variables
	DOUBLE PRECISION QI,QF,QAA,ZAA,EBB,GETSUDAKOV,DERIV,T1
	CHARACTER*2 TYPE
	LOGICAL INS2

	QA=QAA
	ZA2=ZAA
	EB=EBB
	TYP=TYPE
	T=T1
	INSTATE=INS2
	GETPROBA=GETSUDAKOV(QI,QAA,QF,ZAA,EBB,TYPE,T1,INS2)
     &      *DERIV(QF,1)
	END


***********************************************************************
***	  function getsudakov
***********************************************************************
	DOUBLE PRECISION FUNCTION GETSUDAKOV(QMAX1,QA1,QB1,ZA1,EB1,
     &                                                TYPE3,T2,INS)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--variables for Sudakov integration
	COMMON/SUDAINT/QA,ZA2,EB,T,INSTATE,TYP
	DOUBLE PRECISION QA,ZA2,EB,T
	CHARACTER*2 TYP
	LOGICAL INSTATE
C--local variables
	DOUBLE PRECISION QMAX1,QA1,QB1,ZA1,EB1,TMAX,TB,YSTART,EPSI,
     &HFIRST,T2,GETINSUDAFAST,QB2
	CHARACTER*2 TYPE3
	LOGICAL INS
      DATA EPSI/1.d-4/

	QB2=QB1
	IF(INS)THEN
       IF(QB2.LT.Q0) write(logfid,*) 'error: Q < Q0',QB2,QMAX1
       IF(QB2.LT.(Q0+1.d-10)) QB2=QB2+1.d-10
      ELSE 
       IF(QB2.LT.Q0) write(logfid,*) 'error: Q < min',QB2,QMAX1
       IF(QB2.LT.(Q0+1.d-10)) QB2=QB2+1.d-10
      ENDIF 
      IF(QB2.GE.(QMAX1-1.d-10)) THEN
       GETSUDAKOV=1.d0
      ELSE
	 IF(INS)THEN
	  GETSUDAKOV=GETINSUDAFAST(QB1,QMAX1,TYPE3)
	 ELSE
	  QA=QA1
	  ZA2=ZA1
	  EB=EB1
	  TYP=TYPE3
	  T=T2
	  INSTATE=.FALSE.
        HFIRST=0.01*(QMAX1-QB1)
        YSTART=0.d0
        CALL ODEINT(YSTART,QB2,QMAX1,EPSI,HFIRST,0.d0,1)
        GETSUDAKOV=EXP(-YSTART)
	 ENDIF
      ENDIF
	END


***********************************************************************
***	  function getinsudakov
***********************************************************************
	DOUBLE PRECISION FUNCTION GETINSUDAKOV(QB,QMAX1,TYPE3)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--variables for Sudakov integration
	COMMON/SUDAINT/QA,ZA2,EB,T,INSTATE,TYP
	DOUBLE PRECISION QA,ZA2,EB,T
	CHARACTER*2 TYP
	LOGICAL INSTATE
C--local variables
	DOUBLE PRECISION QMAX1,QB,QB1,ZA1,EA1,YSTART,EPSI,
     &HFIRST
	CHARACTER*2 TYPE3
      DATA EPSI/1.d-4/

      QB1=QB
      IF(QB1.LT.Q0) write(logfid,*) 'error: Q < Q0',QB1,QMAX1
      IF(QB1.LT.(Q0+1.d-12)) QB1=QB1+1.d-12
      IF(QB1.GE.(QMAX1-1.d-12)) THEN
       GETINSUDAKOV=1.d0
      ELSE
	 TYP=TYPE3
       HFIRST=0.01*(QMAX1-QB1)
       YSTART=0.d0
       CALL ODEINT(YSTART,QB1,QMAX1,EPSI,HFIRST,0.d0,6)
       GETINSUDAKOV=EXP(-YSTART)
      ENDIF
	END


***********************************************************************
***	  function deriv
***********************************************************************
      DOUBLE PRECISION FUNCTION DERIV(XVAL,W4)
      IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--variables for splitting function integration
	COMMON/INTSPLITF/QQUAD,FM
	DOUBLE PRECISION QQUAD,FM
C--variables for Sudakov integration
	COMMON/SUDAINT/QA,ZA2,EB,T,INSTATE,TYP
	DOUBLE PRECISION QA,ZA2,EB,T
	CHARACTER*2 TYP
	LOGICAL INSTATE
C--variables for pdf integration
	COMMON/PDFINTV/XMAX,Z
	DOUBLE PRECISION XMAX,Z
C--variables for cross section integration 
	COMMON/XSECV/QLOW,MDX
	DOUBLE PRECISION QLOW,MDX
C--local variables
	INTEGER W4
      DOUBLE PRECISION XVAL,GETSPLITI,PI,ALPHAS,GETINSPLITI,
     &GETINSUDAFAST,SCATPRIMFUNC,PQQ,PQG,PGG,PGQ,
     &MEDDERIV
	DATA PI/3.141592653589793d0/

	IF(W4.EQ.1)THEN
C--Sudakov integration
	 IF(INSTATE)THEN
        DERIV=2.*GETINSPLITI(XVAL,TYP)/XVAL
	 ELSE
        DERIV=2.*GETSPLITI(QA,XVAL,ZA2,EB,TYP)/XVAL
	 ENDIF
	ELSEIF(W4.EQ.2)THEN
C--P(q->qg) integration
	 DERIV=(1.+FM)*ALPHAS(XVAL*(1.-XVAL)*QQUAD/1.,LPS)*
     &		PQQ(XVAL)/(2.*PI)
	ELSEIF(W4.EQ.3)THEN
C--P(g->gg) integration
       DERIV=(1.+FM)*ALPHAS(XVAL*(1.-XVAL)*QQUAD/1.,LPS)
     &           *PGG(XVAL)/(2.*PI)
	ELSEIF(W4.EQ.4)THEN
C--P(g->qq) integration
	 DERIV=(1.+FM)*ALPHAS(XVAL*(1-XVAL)*QQUAD/1.,LPS)*
     &	PQG(XVAL)/(2.*PI)	
	ELSEIF(W4.EQ.5)THEN
	 DERIV=EXP(-XVAL)/XVAL
	ELSEIF(W4.EQ.6)THEN
       DERIV=2.*GETINSPLITI(XVAL,TYP)/XVAL
	ELSEIF(W4.EQ.7)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'QQ')
     &	*ALPHAS((1.-Z)*XVAL**2/1.,LPS)
     &	*PQQ(Z)/(2.*PI*XVAL)
	ELSEIF(W4.EQ.8)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'GC')
     &	*ALPHAS((1.-Z)*XVAL**2/1.,LPS)
     &	*PGQ(Z)/(2.*PI*XVAL)
	ELSEIF(W4.EQ.9)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'QQ')
     &	*ALPHAS((1.-Z)*XVAL**2/1.,LPS)
     &	*PQG(Z)/(2.*PI*XVAL)	
	ELSEIF(W4.EQ.10)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'GC')
     &	*ALPHAS((1.-Z)*XVAL**2/1.,LPS)*
     &      *2.*PGG(Z)/(2.*PI*XVAL)
	ELSEIF(W4.EQ.11)THEN
	 DERIV=3.*GETINSPLITI(SCALEFACM*SQRT(XVAL),'GQ')
     &	*SCATPRIMFUNC(XVAL,MDX)/(2.*XVAL)
	ELSEIF(W4.EQ.12)THEN
	 DERIV=2.*GETINSPLITI(SCALEFACM*SQRT(XVAL),'QG')
     &	*SCATPRIMFUNC(XVAL,MDX)/(3.*XVAL)
	ELSEIF(W4.EQ.13)THEN
	 DERIV=GETINSUDAFAST(QLOW,SCALEFACM*SQRT(XVAL),'GC')
     &	*3.*2.*PI*ALPHAS(XVAL+MDX**2,LQCD)**2/(2.*(XVAL+MDX**2)**2)
	ELSEIF(W4.EQ.14)THEN
	 DERIV=GETINSUDAFAST(QLOW,SCALEFACM*SQRT(XVAL),'QQ')
     &	*2.*2.*PI*ALPHAS(XVAL+MDX**2,LQCD)**2/(3.*(XVAL+MDX**2)**2)
	ELSEIF(W4.EQ.21)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'QQ')*GETINSPLITI(XVAL,'QQ')
     &	/XVAL
	ELSEIF(W4.EQ.22)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'GC')*GETINSPLITI(XVAL,'GQ')
     &	/XVAL
	ELSEIF(W4.EQ.23)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'QQ')*GETINSPLITI(XVAL,'QG')
     &	/XVAL
	ELSEIF(W4.EQ.24)THEN
	 DERIV=2.*GETINSUDAFAST(XVAL,XMAX,'GC')*2.
     &	*GETINSPLITI(XVAL,'GG')/XVAL
      ELSE
       DERIV=MEDDERIV(XVAL,W4-100)
      ENDIF
      END


***********************************************************************
***	  function getspliti
***********************************************************************
	DOUBLE PRECISION FUNCTION GETSPLITI(QA,QB,ZETA,EB,TYPE1)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--splitting integral
      COMMON/SPLITINT/SPLITIGGV(1000,1000),SPLITIQQV(1000,1000),
     &SPLITIQGV(1000,1000),QVAL(1000),ZMVAL(1000),QMAX,ZMMIN,NPOINT
      INTEGER NPOINT
      DOUBLE PRECISION SPLITIGGV,SPLITIQQV,SPLITIQGV,
     &QVAL,ZMVAL,QMAX,ZMMIN
C--variables for splitting function integration
	COMMON/INTSPLITF/QQUAD,FM
	DOUBLE PRECISION QQUAD,FM
C--number of extrapolations in tables
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
C--local variables
	INTEGER I,J,LT,QLMAX,ZLMAX,QLINE,ZLINE
	DOUBLE PRECISION QA,QB,ZETA,EB,LOW,X1A(2),X2A(2),YA(2,2),Y,
     &SPLITINTGG,SPLITINTQG,A,B,YB(2)
	CHARACTER*2 TYPE1	

	ntotspliti=ntotspliti+1
	if (qb.gt.qmax) then
	  noverspliti=noverspliti+1
	  if (noverspliti.le.25) 
     &	write(logfid,*)'WARNING in getspliti: need to extrapolate: ',
     &	qb,qmax
	endif

C--find boundaries for z integration
      IF(ANGORD.AND.(ZETA.NE.1.d0))THEN
       LOW=MAX(0.5-0.5*SQRT(1.-Q0**2/QB**2)
     &	*SQRT(1.-QB**2/EB**2),
     &     0.5-0.5*SQRT(1.-4.*QB**2*(1.-ZETA)/(ZETA*QA**2)))
      ELSE
       LOW=0.5-0.5*SQRT(1.-Q0**2/QB**2)
     &	*SQRT(1.-QB**2/EB**2)
      ENDIF
C--find values in array
        QLMAX=INT((QB-QVAL(1))*NPOINT/(QVAL(1000)-QVAL(1))+1)
        QLINE=MAX(QLMAX,1)
        QLINE=MIN(QLINE,NPOINT)
        ZLMAX=INT((LOG(LOW)-LOG(ZMVAL(1)))*NPOINT/
     &        (LOG(ZMVAL(1000))-LOG(ZMVAL(1)))+1)
        ZLINE=MAX(ZLMAX,1)
        ZLINE=MIN(ZLINE,NPOINT)
	  IF((QLINE.GT.999).OR.(ZLINE.GT.999).OR.
     &	(QLINE.LT.1).OR.(ZLINE.LT.1))THEN 
         write(logfid,*)'ERROR in GETSPLITI: line number out of bound',
     &	QLINE,ZLINE
	  ENDIF
        IF((TYPE1.EQ.'GG').OR.(TYPE1.EQ.'GC'))THEN
         DO 17 I=1,2
          X1A(I)=QVAL(QLINE-1+I)
          X2A(I)=ZMVAL(ZLINE-1+I)
          DO 16 J=1,2
           YA(I,J)=SPLITIGGV(QLINE-1+I,ZLINE-1+J)
 16       CONTINUE
 17      CONTINUE
 	   DO 30 I=1,2
	    A=(YA(I,2)-YA(I,1))/(X2A(2)-X2A(1))
	    B=YA(I,1)-A*X2A(1)
	    YB(I)=A*LOW+B
 30	   CONTINUE
	   IF(X1A(1).EQ.X1A(2))THEN
	    Y=(YB(1)+YB(2))/2.
	   ELSE
	    A=(YB(2)-YB(1))/(X1A(2)-X1A(1))
	    B=YB(1)-A*X1A(1)
	    Y=A*QB+B
	   ENDIF
         IF(TYPE1.EQ.'GG')THEN
          GETSPLITI=MIN(Y,10.d0)
         ELSE
          SPLITINTGG=MIN(Y,10.d0)
         ENDIF
        ENDIF
        IF((TYPE1.EQ.'QG').OR.(TYPE1.EQ.'GC'))THEN
         DO 19 I=1,2
          X1A(I)=QVAL(QLINE-1+I)
          X2A(I)=ZMVAL(ZLINE-1+I)
          DO 18 J=1,2
           YA(I,J)=SPLITIQGV(QLINE-1+I,ZLINE-1+J)
 18       CONTINUE
 19      CONTINUE
 	   DO 31 I=1,2
	    A=(YA(I,2)-YA(I,1))/(X2A(2)-X2A(1))
	    B=YA(I,1)-A*X2A(1)
	    YB(I)=A*LOW+B
 31	   CONTINUE
	   IF(X1A(1).EQ.X1A(2))THEN
	    Y=(YB(1)+YB(2))/2.
	   ELSE
	    A=(YB(2)-YB(1))/(X1A(2)-X1A(1))
	    B=YB(1)-A*X1A(1)
	    Y=A*QB+B
	   ENDIF
         IF(TYPE1.EQ.'QG')THEN
          GETSPLITI=NF*MIN(Y,10.d0)
         ELSE
          SPLITINTQG=NF*MIN(Y,10.d0)
         ENDIF
        ENDIF
        IF(TYPE1.EQ.'QQ')THEN
         DO 21 I=1,2
          X1A(I)=QVAL(QLINE-1+I)
          X2A(I)=ZMVAL(ZLINE-1+I)
          DO 20 J=1,2
           YA(I,J)=SPLITIQQV(QLINE-1+I,ZLINE-1+J)
 20       CONTINUE
 21      CONTINUE
 	   DO 32 I=1,2
	    A=(YA(I,2)-YA(I,1))/(X2A(2)-X2A(1))
	    B=YA(I,1)-A*X2A(1)
	    YB(I)=A*LOW+B
 32	   CONTINUE
	   IF(X1A(1).EQ.X1A(2))THEN
	    Y=(YB(1)+YB(2))/2.
	   ELSE
	    A=(YB(2)-YB(1))/(X1A(2)-X1A(1))
	    B=YB(1)-A*X1A(1)
	    Y=A*QB+B
	   ENDIF
         GETSPLITI=MIN(Y,10.d0)
        ENDIF
        IF(TYPE1.EQ.'GC') GETSPLITI=SPLITINTGG+SPLITINTQG
      END


***********************************************************************
***	  function getinspliti
***********************************************************************
	DOUBLE PRECISION FUNCTION GETINSPLITI(QB,TYPE1)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION QB,LOW,PI,Y,SPLITINTGG,SPLITINTQG,UP,EI
	CHARACTER*2 TYPE1	
	DATA PI/3.141592653589793d0/

C--find boundaries for z integration
	 UP = 1. - Q0**2/(4.*QB**2)
       IF((TYPE1.EQ.'GG').OR.(TYPE1.EQ.'GC'))THEN
	  LOW=1.d0-UP
	  IF (UP.LE.LOW) THEN
	   GETINSPLITI=0.d0
	   RETURN
	  ENDIF
	  Y = 2.* ( LOG(LOG((1.-LOW)*QB**2/LPS**2))
     &	- LPS**2*EI(LOG((1.-LOW)*QB**2/LPS**2))/QB**2
     &	+ LPS**4*EI(2.*LOG((1.-LOW)*QB**2/LPS**2))/QB**4
     &	- LPS**6*EI(3.*LOG((1.-LOW)*QB**2/LPS**2))/QB**6
     &      - LOG(LOG((1.-UP)*QB**2/LPS**2))
     &	+ LPS**2*EI(LOG((1.-UP)*QB**2/LPS**2))/QB**2
     &	- LPS**4*EI(2.*LOG((1.-UP)*QB**2/LPS**2))/QB**4
     &	+ LPS**6*EI(3.*LOG((1.-UP)*QB**2/LPS**2))/QB**6
     &	+ LOW - LOG(LOW) - UP + LOG(UP) )
     &	*3.*12.*PI/(2.*PI*(33.-2.*NF))
        IF(TYPE1.EQ.'GG')THEN
         GETINSPLITI=Y
        ELSE
         SPLITINTGG=Y
        ENDIF
       ENDIF
       IF((TYPE1.EQ.'QG').OR.(TYPE1.EQ.'GC'))THEN
	  LOW=0.d0
	  IF (UP.LE.LOW) THEN
	   GETINSPLITI=0.d0
	   RETURN
	  ENDIF
	  Y = ( 2.*LPS**6*EI(3.*LOG((1.-LOW)*QB**2/LPS**2))/QB**6
     &	- 2.*LPS**4*EI(2.*LOG((1.-LOW)*QB**2/LPS**2))/QB**4
     &	+ 2.*LPS**2*EI(LOG((1.-LOW)*QB**2/LPS**2))/QB**2
     &	- 2.*LPS**6*EI(3.*LOG((1.-UP)*QB**2/LPS**2))/QB**6
     &	+ 2.*LPS**4*EI(2.*LOG((1.-UP)*QB**2/LPS**2))/QB**4
     &	- 2.*LPS**2*EI(LOG((1.-UP)*QB**2/LPS**2))/QB**2 )
     &	*12.*PI/(2.*2.*PI*(33.-2.*NF))
        IF(TYPE1.EQ.'QG')THEN
         GETINSPLITI=NF*Y
        ELSE
         SPLITINTQG=NF*Y
        ENDIF
       ENDIF
       IF(TYPE1.EQ.'QQ')THEN
	  LOW=0.d0
	  IF (UP.LE.LOW) THEN
	   GETINSPLITI=0.d0
	   RETURN
	  ENDIF
	  Y = ( 2.*LOG(LOG((1.-LOW)*QB**2/LPS**2))
     &	- 2.*LPS**2*EI(LOG((1.-LOW)*QB**2/LPS**2))/QB**2
     &	+ LPS**4*EI(2.*LOG((1.-LOW)*QB**2/LPS**2))/QB**4
     &	- 2.*LOG(LOG((1.-UP)*QB**2/LPS**2))
     &	+ 2.*LPS**2*EI(LOG((1.-UP)*QB**2/LPS**2))/QB**2
     &	- LPS**4*EI(2.*LOG((1.-UP)*QB**2/LPS**2))/QB**4 ) 
     &	*4.*12.*PI/(3.*2.*PI*(33.-2.*NF))
        GETINSPLITI=Y
       ENDIF
       IF(TYPE1.EQ.'GQ')THEN
	  LOW=1.d0-UP
	  IF (UP.LE.LOW) THEN
	   GETINSPLITI=0.d0
	   RETURN
	  ENDIF
	  Y = (UP**2/2.-2.*UP+2.*LOG(UP)-LOW**2/2.+2.*LOW- 2.*LOG(LOW)) 
     &	*4.*12.*PI/(3.*2.*PI*(33.-2.*NF)*LOG(QB**2/LPS**2))
        GETINSPLITI=Y
       ENDIF
       IF(TYPE1.EQ.'GC') GETINSPLITI=SPLITINTGG+SPLITINTQG
      END


***********************************************************************
***	  function getpdf
***********************************************************************
	DOUBLE PRECISION FUNCTION GETPDF(X,Q,TYP)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--pdf common block
	COMMON/PDFS/QINQX(2,1000),GINQX(2,1000),QINGX(2,1000),
     &GINGX(2,1000)
	DOUBLE PRECISION QINQX,GINQX,QINGX,GINGX
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--variables for pdf integration
	COMMON/PDFINTV/XMAX,Z
	DOUBLE PRECISION XMAX,Z
C--local variables
	DOUBLE PRECISION X,Q,QLOW,QHIGH,YSTART,EPSI,HFIRST
	CHARACTER*2 TYP
	DATA EPSI/1.d-4/	

	IF((X.LT.0.d0).OR.(X.GT.1.d0).OR.(Q.LT.Q0))THEN
	 write(logfid,*)'error in GETPDF: parameter out of bound',X,Q
	 GETPDF=0.d0
	 RETURN
	ENDIF

	IF(TYP.EQ.'QQ')THEN
	  Z=X
	  XMAX=Q
C--f_q^q
	  QLOW=MAX(Q0,Q0/(2.*SQRT(1.-X)))
	  QHIGH=Q
	  IF((QLOW.GE.QHIGH*(1.d0-1.d-10)).OR.(X.GT.1.d0-1.d-10))THEN
	   YSTART=0.d0
	  ELSE
         HFIRST=0.01*(QHIGH-QLOW)
         YSTART=0.d0
         CALL ODEINT(YSTART,QLOW,QHIGH,EPSI,HFIRST,0.d0,7)
	  ENDIF
	  GETPDF=YSTART
	ELSEIF(TYP.EQ.'GQ')THEN
	  Z=X
	  XMAX=Q
C--f_q^g
	  QLOW=MAX(Q0,MAX(Q0/(2.*SQRT(X)),Q0/(2.*SQRT(1.-X))))
	  QHIGH=Q
	  IF((QLOW.GE.QHIGH*(1.d0-1.d-10)).OR.(X.LT.0.d0+1.d-10)
     &	.OR.(X.GT.1.d0-1.d-10))THEN
	   YSTART=0.d0
	  ELSE
         HFIRST=0.01*(QHIGH-QLOW)
         YSTART=0.d0
         CALL ODEINT(YSTART,QLOW,QHIGH,EPSI,HFIRST,0.d0,8)
	  ENDIF
	  GETPDF=YSTART
	ELSEIF(TYP.EQ.'QG')THEN
	  Z=X
	  XMAX=Q
C--f_q^g
	  QLOW=MAX(Q0,Q0/(2.*SQRT(1.-X)))
	  QHIGH=Q
	  IF((QLOW.GE.QHIGH*(1.d0-1.d-10)).OR.(X.GT.1.d0-1.d-10))THEN
	   YSTART=0.d0
	  ELSE
         HFIRST=0.01*(QHIGH-QLOW)
         YSTART=0.d0
         CALL ODEINT(YSTART,QLOW,QHIGH,EPSI,HFIRST,0.d0,9)
	  ENDIF
	  GETPDF=YSTART
	ELSEIF(TYP.EQ.'GG')THEN
	  Z=X
	  XMAX=Q
C--f_q^q
	QLOW=MAX(Q0,MAX(Q0/(2.*SQRT(X)),Q0/(2.*SQRT(1.-X))))
	  QHIGH=Q
	  IF((QLOW.GE.QHIGH*(1.d0-1.d-10)).OR.(X.LT.0.d0+1.d-10)
     &	.OR.(X.GT.1.d0-1d-10))THEN
	   YSTART=0.d0
	  ELSE
         HFIRST=0.01*(QHIGH-QLOW)
         YSTART=0.d0
         CALL ODEINT(YSTART,QLOW,QHIGH,EPSI,HFIRST,0.d0,10)
	  ENDIF
	  GETPDF=YSTART
	ELSE
	 write(logfid,*)'error: pdf-type ',TYP,' does not exist'
	 GETPDF=0.d0
	ENDIF
	END

***********************************************************************
***	  function getpdfxint
***********************************************************************
	DOUBLE PRECISION FUNCTION GETPDFXINT(Q,TYP)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--pdf common block
	COMMON/PDFS/QINQX(2,1000),GINQX(2,1000),QINGX(2,1000),
     &GINGX(2,1000)
	DOUBLE PRECISION QINQX,GINQX,QINGX,GINGX
C--number of extrapolations in tables
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
C--local variables
	INTEGER J,Q2CLOSE,Q2LINE
	DOUBLE PRECISION Q,XA(2),YA(2),Y,A,B
	CHARACTER*2 TYP

	ntotpdf=ntotpdf+1
	if (q**2.gt.QINQX(1,1000)) then
	  noverpdf=noverpdf+1
	  if (noverpdf.le.25) 
     &	write(logfid,*)'WARNING in getpdfxint: need to extrapolate: ',
     &	q**2,QINQX(1,1000)
	endif

      Q2CLOSE=INT((LOG(Q**2)-LOG(QINQX(1,1)))*999.d0/
     &	(LOG(QINQX(1,1000))-LOG(QINQX(1,1)))+1)
      Q2LINE=MAX(Q2CLOSE,1)
      Q2LINE=MIN(Q2LINE,999)
	IF((Q2LINE.GT.999).OR.(Q2LINE.LT.1))THEN
       write(logfid,*)'ERROR in GETPDFXINT: line number out of bound',
     &	Q2LINE
	ENDIF

      IF(TYP.EQ.'QQ')THEN
       DO 11 J=1,2
        XA(J)=QINQX(1,Q2LINE-1+J)
        YA(J)=QINQX(2,Q2LINE-1+J)
 11    CONTINUE
      ELSEIF(TYP.EQ.'GQ')THEN
       DO 13 J=1,2
        XA(J)=GINQX(1,Q2LINE-1+J)
        YA(J)=GINQX(2,Q2LINE-1+J)
 13    CONTINUE
      ELSEIF(TYP.EQ.'QG')THEN
       DO 15 J=1,2
        XA(J)=QINGX(1,Q2LINE-1+J)
        YA(J)=QINGX(2,Q2LINE-1+J)
 15    CONTINUE
      ELSEIF(TYP.EQ.'GG')THEN
       DO 17 J=1,2
        XA(J)=GINGX(1,Q2LINE-1+J)
        YA(J)=GINGX(2,Q2LINE-1+J)
 17    CONTINUE
	ELSE
	 write(logfid,*)'error in GETPDFXINT: unknown integral type ',TYP
	ENDIF
	A=(YA(2)-YA(1))/(XA(2)-XA(1))
	B=YA(1)-A*XA(1)
	Y=A*Q**2+B
	GETPDFXINT=Y
	END


***********************************************************************
***	  subroutine getpdfxintexact
***********************************************************************
	DOUBLE PRECISION FUNCTION GETPDFXINTEXACT(Q,TYP)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--variables for pdf integration
	COMMON/PDFINTV/XMAX,Z
	DOUBLE PRECISION XMAX,Z
C--local variables
	DOUBLE PRECISION Q,EPSI,YSTART,HFIRST
	CHARACTER*2 TYP
	DATA EPSI/1.d-4/
	
      HFIRST=0.01d0
      YSTART=0.d0
	XMAX=Q
	Z=0.d0
	IF(TYP.EQ.'QQ')THEN
       CALL ODEINT(YSTART,Q0,Q,EPSI,HFIRST,0.d0,21)
	ELSEIF(TYP.EQ.'QG')THEN
       CALL ODEINT(YSTART,Q0,Q,EPSI,HFIRST,0.d0,23)
	ELSEIF(TYP.EQ.'GQ')THEN
       CALL ODEINT(YSTART,Q0,Q,EPSI,HFIRST,0.d0,22)
	ELSEIF(TYP.EQ.'GG')THEN
       CALL ODEINT(YSTART,Q0,Q,EPSI,HFIRST,0.d0,24)
	ENDIF
	GETPDFXINTEXACT=YSTART 
	END


***********************************************************************
***	  function getxsecint
***********************************************************************
	DOUBLE PRECISION FUNCTION GETXSECINT(TM,MD,TYP2)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--cross secttion common block
	COMMON/XSECS/INTQ1(1001,101),INTQ2(1001,101),
     &INTG1(1001,101),INTG2(1001,101)
	DOUBLE PRECISION INTQ1,INTQ2,INTG1,INTG2
C--variables for cross section integration 
	COMMON/XSECV/QLOW,MDX
	DOUBLE PRECISION QLOW,MDX
C--number of extrapolations in tables
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
C--local variables
	INTEGER TLINE,TCLOSE,MDCLOSE,MDLINE,I,J
	DOUBLE PRECISION TM,X1A(2),X2A(2),YA(2,2),Y,MD,YB(2),A,B
	CHARACTER*2 TYP2

	ntotxsec=ntotxsec+1
	if (tm.gt.intq1(1000,101)) then
	  noverxsec=noverxsec+1
	  if (noverpdf.le.25) 
     &	write(logfid,*)'WARNING in getxsecint: need to extrapolate: ',
     &	tm,intq1(1000,101)
	endif

       TCLOSE=INT((LOG(TM)-LOG(INTQ1(1,101)))*999.d0/
     &	(LOG(INTQ1(1000,101))-LOG(INTQ1(1,101)))+1)
       TLINE=MAX(TCLOSE,1)
       TLINE=MIN(TLINE,999)
       MDCLOSE=INT((MD-INTQ1(1001,1))*99.d0/
     &(INTQ1(1001,100)-INTQ1(1001,1))+1)
       MDLINE=MAX(MDCLOSE,1)
       MDLINE=MIN(MDLINE,99)
	 IF((TLINE.GT.999).OR.(MDLINE.GT.99)
     &  .OR.(TLINE.LT.1).OR.(MDLINE.LT.1)) THEN
      write(logfid,*)'ERROR in GETXSECINT: line number out of bound',
     &	TLINE,MDLINE
	 ENDIF

       IF(TYP2.EQ.'QA')THEN
C--first quark integral
        DO 12 I=1,2
         X1A(I)=INTQ1(1001,MDLINE-1+I)
         X2A(I)=INTQ1(TLINE-1+I,101)
         DO 11 J=1,2
          YA(I,J)=INTQ1(TLINE-1+J,MDLINE-1+I)
 11      CONTINUE
 12     CONTINUE
	 ELSEIF(TYP2.EQ.'QB')THEN
C--second quark integral
        DO 18 I=1,2
         X1A(I)=INTQ2(1001,MDLINE-1+I)
         X2A(I)=INTQ2(TLINE-1+I,101)
         DO 17 J=1,2
          YA(I,J)=INTQ2(TLINE-1+J,MDLINE-1+I)
 17      CONTINUE
 18     CONTINUE
	 ELSEIF(TYP2.EQ.'GA')THEN
C--first gluon integral
        DO 14 I=1,2
         X1A(I)=INTG1(1001,MDLINE-1+I)
         X2A(I)=INTG1(TLINE-1+I,101)
         DO 13 J=1,2
          YA(I,J)=INTG1(TLINE-1+J,MDLINE-1+I)
 13      CONTINUE
 14     CONTINUE
	 ELSEIF(TYP2.EQ.'GB')THEN
C--second gluon integral
        DO 16 I=1,2
         X1A(I)=INTG2(1001,MDLINE-1+I)
         X2A(I)=INTG2(TLINE-1+I,101)
         DO 15 J=1,2
          YA(I,J)=INTG2(TLINE-1+J,MDLINE-1+I)
 15      CONTINUE
 16     CONTINUE
	 ELSE
	  write(logfid,*)'error in GETXSECINT: unknown integral type ',
     &										TYP2
	 ENDIF
	 DO 19 I=1,2
	  A=(YA(I,2)-YA(I,1))/(X2A(2)-X2A(1))
	  B=YA(I,1)-A*X2A(1)
	  YB(I)=A*TM+B
 19	 CONTINUE
	 IF(X1A(1).EQ.X1A(2))THEN
	  Y=YB(1)
	 ELSE
	  A=(YB(2)-YB(1))/(X1A(2)-X1A(1))
	  B=YB(1)-A*X1A(1)
	  Y=A*MD+B
	 ENDIF
	 GETXSECINT=Y
	END


***********************************************************************
***	  function getinsudafast
***********************************************************************
	DOUBLE PRECISION FUNCTION GETINSUDAFAST(Q1,Q2,TYP)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION Q1,Q2,GETINSUDARED
	CHARACTER*2 TYP
	
	IF(Q2.LE.Q1)THEN
	 GETINSUDAFAST=1.d0
	ELSEIF(Q1.LE.Q0)THEN
	 GETINSUDAFAST=GETINSUDARED(Q2,TYP)
	ELSE
	 GETINSUDAFAST=GETINSUDARED(Q2,TYP)/GETINSUDARED(Q1,TYP)
	ENDIF
      IF(GETINSUDAFAST.GT.1.d0) GETINSUDAFAST=1.d0
	IF(GETINSUDAFAST.LT.(-1.d-10))THEN
	 write(logfid,*)'ERROR: GETINSUDAFAST < 0:',
     &	GETINSUDAFAST,' for',Q1,' ',Q2,' ',TYP
	ENDIF
	if (getinsudafast.lt.0.d0) getinsudafast = 0.d0
	END


***********************************************************************
***	  function getinsudared
***********************************************************************
	DOUBLE PRECISION FUNCTION GETINSUDARED(Q,TYP2)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--Sudakov common block
	COMMON/INSUDA/SUDAQQ(1000,2),SUDAQG(1000,2),SUDAGG(1000,2),
     &SUDAGC(1000,2)
	DOUBLE PRECISION SUDAQQ,SUDAQG,SUDAGG,SUDAGC
C--number of extrapolations in tables
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
C--local variables
	INTEGER QCLOSE,QBIN,I
	DOUBLE PRECISION Q,XA(2),YA(2),Y,A,B
	CHARACTER*2 TYP2

	ntotsuda=ntotsuda+1
	if (q.gt.sudaqq(1000,1)) then
	  noversuda=noversuda+1
	  if (noversuda.le.25) 
     &	write(logfid,*)'WARNING in getinsudared: need to extrapolate: ',
     &	q,sudaqq(1000,1)
	endif

      QCLOSE=INT((LOG(Q)-LOG(SUDAQQ(1,1)))*999.d0
     &	/(LOG(SUDAQQ(1000,1))-LOG(SUDAQQ(1,1)))+1)
      QBIN=MAX(QCLOSE,1)
      QBIN=MIN(QBIN,999)
	IF((QBIN.GT.999).OR.(QBIN.LT.1)) THEN
       write(logfid,*)
     &	'ERROR in GETINSUDARED: line number out of bound',QBIN
	ENDIF
	IF(TYP2.EQ.'QQ')THEN
       DO 16 I=1,2
        XA(I)=SUDAQQ(QBIN-1+I,1)
        YA(I)=SUDAQQ(QBIN-1+I,2)
 16    CONTINUE
	ELSEIF(TYP2.EQ.'QG')THEN
       DO 17 I=1,2
        XA(I)=SUDAQG(QBIN-1+I,1)
        YA(I)=SUDAQG(QBIN-1+I,2)
 17    CONTINUE
	ELSEIF(TYP2.EQ.'GG')THEN
       DO 18 I=1,2
        XA(I)=SUDAGG(QBIN-1+I,1)
        YA(I)=SUDAGG(QBIN-1+I,2)
 18    CONTINUE
	ELSEIF(TYP2.EQ.'GC')THEN
       DO 19 I=1,2
        XA(I)=SUDAGC(QBIN-1+I,1)
        YA(I)=SUDAGC(QBIN-1+I,2)
 19    CONTINUE
	ELSE
	 write(logfid,*)'error in GETINSUDARED: unknown type ',TYP2
	ENDIF
	A=(YA(2)-YA(1))/(XA(2)-XA(1))
	B=YA(1)-A*XA(1)
	Y=A*Q+B
	GETINSUDARED=Y
	IF(GETINSUDARED.LT.(-1.d-10))THEN
	 write(logfid,*) 'ERROR: GETINSUDARED < 0:',GETINSUDARED,Q,TYP2
	ENDIF
	if (getinsudared.lt.0.d0) getinsudared = 0.d0
	END


***********************************************************************
***	  function getsscat
***********************************************************************
      DOUBLE PRECISION FUNCTION GETSSCAT(EN,px,py,PZ,MP,LW,TYPE1,TYPE2,
     &	x,y,z,t,mode)
      IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--variables for cross section integration 
	COMMON/XSECV/QLOW,MDX
	DOUBLE PRECISION QLOW,MDX
C--local variables
	integer mode
      DOUBLE PRECISION UP,EN,LW,SCATPRIMFUNC,CCOL,MP,
     &LOW,GETPDFXINT,GETXSECINT,MDEB,pz,pcms2,shat,gettemp,
     &x,y,z,t,getmd,avmom(5),px,py,getmdmin,getmdmax,pproj,psct
      CHARACTER TYPE1,TYPE2

       IF(TYPE1.EQ.'Q')THEN
        CCOL=2./3.
       ELSE
        CCOL=3./2.
       ENDIF 
	 if (mode.eq.0) then
	   mdeb = getmd(x,y,z,t)
	   call avscatcen(x,y,z,t,
     &	avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))
	   shat = avmom(5)**2 + mp**2 + 
     &	2.*(avmom(4)*en - avmom(1)*px - avmom(2)*py - avmom(3)*pz)
	   pcms2 = (shat+mp**2-avmom(5)**2)**2/(4.*shat)-mp**2
	   up = 4.*pcms2
	 else
	   if (mode.eq.1) then
	     mdeb = getmdmin()
	   else 
	     mdeb = getmdmax()
	   endif 
	   call maxscatcen(avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))
	   psct = sqrt(avmom(1)**2+avmom(2)**2+avmom(3)**2)
	   pproj = sqrt(px**2+py**2+pz**2)
	   shat = avmom(5)**2 + mp**2 + 2.*(en*avmom(4) + pproj*psct)
	   pcms2 = (shat+mp**2-avmom(5)**2)**2/(4.*shat)-mp**2
	   up = 4.*pcms2
	 endif
	 LOW=LW**2
	 IF(LOW.GT.UP)THEN
	  GETSSCAT=0.d0
	  RETURN
	 ENDIF
	 IF((TYPE2.EQ.'C').OR.
     &	((TYPE1.EQ.'Q').AND.(TYPE2.EQ.'Q')).OR.
     &		((TYPE1.EQ.'G').AND.(TYPE2.EQ.'G')))THEN
        GETSSCAT=CCOL*(SCATPRIMFUNC(UP,MDEB)-SCATPRIMFUNC(LOW,MDEB))
	 ELSE
	  GETSSCAT=0.d0
	 ENDIF
	 LOW=Q0**2/SCALEFACM**2
	 IF(UP.GT.LOW)THEN
        IF(TYPE1.EQ.'Q')THEN
	   IF((TYPE2.EQ.'C').OR.(TYPE2.EQ.'G'))THEN
	    GETSSCAT=GETSSCAT+GETPDFXINT(SCALEFACM*SQRT(UP),'GQ')
     &	*3.*SCATPRIMFUNC(UP,MDEB)/2.
	    GETSSCAT=GETSSCAT-GETXSECINT(UP,MDEB,'QA')
	   ENDIF
	  ELSE
	   IF((TYPE2.EQ.'C').OR.(TYPE2.EQ.'G'))THEN
	    GETSSCAT=GETSSCAT+CCOL*(SCATPRIMFUNC(UP,MDEB)-
     &			SCATPRIMFUNC(LOW,MDEB))
     &		- GETXSECINT(UP,MDEB,'GB')
	   ENDIF
	   IF((TYPE2.EQ.'C').OR.(TYPE2.EQ.'Q'))THEN
	    GETSSCAT=GETSSCAT+2.*GETPDFXINT(SCALEFACM*SQRT(UP),'QG')
     &	*2.*SCATPRIMFUNC(UP,MDEB)/3.
	    GETSSCAT=GETSSCAT-2.*GETXSECINT(UP,MDEB,'GA')
	   ENDIF
	  ENDIF
	 ENDIF
	IF(GETSSCAT.LT.-1.d-4)
     &    write(logfid,*) 'error: cross section < 0',GETSSCAT,'for',
     &	EN,MP,LW,TYPE1,TYPE2,LW**2,UP
	GETSSCAT=MAX(GETSSCAT,0.d0)
      END



***********************************************************************
***	  function getmass
***********************************************************************
	DOUBLE PRECISION FUNCTION GETMASS(QBMIN,QBMAX,THETA,EP,TYPE,
     &                                   MAX2,INS,ZDEC,QQBARDEC)
	IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--Common block of Pythia
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	INTEGER MSTU,MSTJ
	DOUBLE PRECISION PARU,PARJ
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
	INTEGER MDCY,MDME,KFDP
	DOUBLE PRECISION BRAT
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--factor in front of alphas argument
	COMMON/ALPHASFAC/PTFAC
	DOUBLE PRECISION PTFAC
C--local variables
	DOUBLE PRECISION qbmin,qbmax,theta,ep,max2,zdec,
     &q2min,alphmax,alphas,log14,pref,q2max,sudaover,gmin,
     &gmax,arg,cand,eps,trueeps,trueval,oest,weight,getinspliti,
     &r,pyr,z,rz,thetanew,r2,pi,pqq,pgg,pqg,rmin
      CHARACTER*2 TYPE
	LOGICAL INS,QQBARDEC
      DATA PI/3.141592653589793d0/
	
	q2min = q0**2

	alphmax = alphas(3.*ptfac*q2min/16.,lps)
	log14 = log(0.25)

      IF(TYPE.EQ.'QQ')THEN
	 pref=4.*alphmax/(3.*2.*PI)
      ELSE
	 pref=29.*alphmax/(8.*2.*PI)
      ENDIF

C--check if phase space available, return 0.d0 otherwise
	IF((qbmax.LE.QBMIN).OR.(EP.LT.QBMIN)) THEN
	 getmass=0.d0
	 ZDEC=0.d0
	 QQBARDEC=.FALSE.
	 RETURN
	ENDIF

      q2max = qbmax**2
! 21	sudaover = exp(-pref*(log(q2min/(4.*q2max))**2 - log14**2))
!	IF(pyr(0).LE.sudaover)THEN
 21   if (q2max-qbmin**2.lt.1e-4)then
	    getmass=qbmin
	    zdec=0.5
	    IF(TYPE.EQ.'QQ')THEN
	      QQBARDEC=.FALSE.
	    ELSE
	      IF(PYR(0).LT.PQG(0.5d0)/(PQG(0.5d0)+PGG(0.5d0)))THEN
	        QQBARDEC=.TRUE.
	      ELSE 
	        QQBARDEC=.FALSE.
	      ENDIF
	    endif
	    return
        endif
        gmax = pref*log(q2min/(4.*q2max))**2
        if (qbmin.gt.0.d0) then
          rmin = exp(pref*log(q2min/(4.*qbmin**2))**2-gmax)
        else
	    rmin = 0.d0
	  endif  
	  
       r=pyr(0)*(1.d0-rmin)+rmin
       arg=gmax+log(r)
       if(arg.lt.0.d0)then
	 getmass=0.d0
	 ZDEC=0.d0
	 QQBARDEC=.FALSE.
	 RETURN
	endif
!	r=pyr(0)
!	gmin = pref*log14**2
!	gmax = pref*log(q2min/(4.*q2max))**2
!	arg = log(r*exp(gmax)+(1.-r)*exp(gmin))
	cand = q2min*exp(sqrt(arg/pref))/4.
	eps = q2min/(4.*cand)

	if ((cand.lt.q2min).or.(cand.lt.qbmin**2)) then
	 getmass=0.d0
	 ZDEC=0.d0
	 QQBARDEC=.FALSE.
	 RETURN
	endif

	IF((CAND.GT.MAX2**2).OR.(CAND.GT.EP**2))THEN
	 q2max=cand
	 goto 21
	ENDIF

	if (ins) then
	  trueval=getinspliti(sqrt(cand),type)
	  oest = -2.*pref*log(eps)
        weight = trueval/oest
	else
C--find true z interval
        TRUEEPS=0.5-0.5*SQRT(1.-q2min/cand)
     &	*SQRT(1.-cand/EP**2)
        IF(TRUEEPS.LT.EPS)
     &	WRITE(logfid,*)'error in getmass: true eps < eps',TRUEEPS,EPS
	  RZ=PYR(0)
	  z = 1.-eps**rz
	  if ((z.lt.trueeps).or.(z.gt.(1.-trueeps))) then
	    weight = 0.
	  else
	    if (type.eq.'QQ')then
!	      if (ins) then
!                trueval = alphas(ptfac*(1.-z)*cand,lps)*pqq(z)/(2.*pi)
!              else
	        trueval = alphas(ptfac*z*(1.-z)*cand,lps)*pqq(z)/(2.*pi)
!              endif
	      oest = 2.*pref/(1.-z)
	      weight = trueval/oest
	    else
	      if (pyr(0).lt.(17./29.)) z = 1.-z
!	      if (ins)then
!	        trueval = alphas(ptfac*(1.-z)*cand,lps)
!     &			*(pgg(z)+pqg(z))/(2.*pi)
!              else
	        trueval = alphas(ptfac*z*(1.-z)*cand,lps)
     &			*(pgg(z)+pqg(z))/(2.*pi)
!              endif
	      oest = alphmax*(17./(4.*z)+3./(1.-z))/(2.*pi)
	      weight = trueval/oest
	    endif
	    thetanew = sqrt(cand/(z*(1.-z)))/ep
	    if (angord.and.(theta.gt.0.).and.(thetanew.gt.theta)) 
     &								weight = 0.d0
	  endif
	endif
	IF (WEIGHT.GT.1.d0) WRITE(logfid,*) 
     &	'problem in getmass: weight> 1',
     &		WEIGHT,TYPE,EPS,TRUEEPS,Z,CAND
	R2=PYR(0)
	IF(R2.GT.WEIGHT)THEN
	 q2max=cand
	 GOTO 21
	ELSE
	 getmass=sqrt(cand)
	 if (.not.ins) then
	   ZDEC=Z
	   IF(TYPE.EQ.'QQ')THEN
	     QQBARDEC=.FALSE.
	   ELSE
	     IF(PYR(0).LT.PQG(Z)/(PQG(Z)+PGG(Z)))THEN
	       QQBARDEC=.TRUE.
	     ELSE 
	       QQBARDEC=.FALSE.
	     ENDIF
	   ENDIF
	  endif
	ENDIF
 	END



***********************************************************************
***	  function generatez
***********************************************************************
	DOUBLE PRECISION FUNCTION GENERATEZ(TI,EA,EPSI,TYPE)
      IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
      DOUBLE PRECISION TI,EA,EPS,PYR,X,R,HELP,R1,EPSI
	CHARACTER*2 TYPE

      IF(TI.EQ.0.d0)THEN
       EPS=EPSI
      ELSE
       EPS=MAX(0.5-0.5*SQRT(1.-Q0**2/TI)
     &      *SQRT(1.-TI/EA**2),EPSI)
      ENDIF
      IF(EPS.GT.0.5)THEN
       GENERATEZ=0.5
       GOTO 61
      ENDIF
 60   R=PYR(0)
 	IF(TYPE.EQ.'QQ')THEN
       X=1.-(1.-EPS)*(EPS/(1.-EPS))**R
       R=PYR(0)
       IF(R.LT.((1.+X**2)/2.))THEN
        GENERATEZ=X
       ELSE
        GOTO 60
       ENDIF
	ELSEIF(TYPE.EQ.'GG')THEN
       X=1./(1.+((1.-EPS)/EPS)**(1.-2.*R))
       R=PYR(0)
	 HELP=((1.-X)/X+X/(1.-X)+X*(1.-X))/(1./(1.-X)+1./X)
       IF(R.LT.HELP)THEN
        GENERATEZ=X
       ELSE
        GOTO 60
       ENDIF
	ELSE
	 R=PYR(0)*(1.-2.*EPS)+EPS
	 R1=PYR(0)/2.
	 HELP=0.5*(R**2+(1.-R)**2)
	 IF(R1.LT.HELP)THEN
	  GENERATEZ=R
	 ELSE
	  GOTO 60
	 ENDIF
	ENDIF
 61	END



***********************************************************************
***	  function scatprimfunc
***********************************************************************
      DOUBLE PRECISION FUNCTION SCATPRIMFUNC(T,MDEB)
      IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
      DOUBLE PRECISION T,PI,S,EI,ALPHAS,T1,MDEB
      DATA PI/3.141592653589793d0/

	 SCATPRIMFUNC = 2.*PI*(12.*PI)**2*(
     &	- EI(-LOG((T+MDEB**2)/LQCD**2))/LQCD**2
     &	- 1./((T+MDEB**2)*LOG((T+MDEB**2)/LQCD**2)))/(33.-2.*NF)**2
      END



***********************************************************************
***	  function intpqq
***********************************************************************
	DOUBLE PRECISION FUNCTION INTPQQ(Z,Q)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION Z,Q

	INTPQQ=6.*4.*(-2.*LOG(LOG(Q**2/LPS**2)
     &	+LOG(1.-Z)))/((33.-2.*NF)*3.)
	END



***********************************************************************
***	  function intpgglow
***********************************************************************
	DOUBLE PRECISION FUNCTION INTPGGLOW(Z,Q)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION Z,Q

	INTPGGLOW=6.*3.*(LOG(LOG(Q**2/LPS**2)+LOG(Z)))/(33.-2.*NF)
	END
	


***********************************************************************
***	  function intpgghigh
***********************************************************************
	DOUBLE PRECISION FUNCTION INTPGGHIGH(Z,Q)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION Z,Q

	INTPGGHIGH=-6.*3.*(LOG(LOG(Q**2/LPS**2)+LOG(1.-Z)))/(33.-2.*NF)
	END
	


***********************************************************************
***	  function intpqglow
***********************************************************************
	DOUBLE PRECISION FUNCTION INTPQGLOW(Z,Q)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION Z,Q,EI

	INTPQGLOW=6.*(LPS**2*EI(LOG(Q**2/LPS**2)+LOG(Z))/Q**2 
     & - 2.*LPS**4*EI(2.*(LOG(Q**2/LPS**2)+LOG(Z)))/Q**4
     & + 2.*LPS**6*EI(3.*(LOG(Q**2/LPS**2)+LOG(Z)))/Q**6)/
     &((33.-2.*NF)*2.)
	END
	


***********************************************************************
***	  function intpqghigh
***********************************************************************
	DOUBLE PRECISION FUNCTION INTPQGHIGH(Z,Q)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION Z,Q,EI

	INTPQGHIGH=-6.*(LPS**2*EI(LOG(Q**2/LPS**2)+LOG(1.-Z))/Q**2 
     & - 2.*LPS**4*EI(2.*(LOG(Q**2/LPS**2)+LOG(1.-Z)))/Q**4
     & + 2.*LPS**6*EI(3.*(LOG(Q**2/LPS**2)+LOG(1.-Z)))/Q**6)/
     &((33.-2.*NF)*2.)
	END



***********************************************************************
***	  function gett
***********************************************************************
 	DOUBLE PRECISION FUNCTION GETT(MINT,MAXT,MDEB)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION TMIN,TMAX,MAXI,PYR,R1,R2,ALPHAS,PI,Y,MAXT,
     &MDEB,MINT,T
	DATA PI/3.141592653589793d0/

	TMAX=MAXT+MDEB**2
	TMIN=MINT+MDEB**2
	IF(TMIN.GT.TMAX) THEN
	 GETT=0.d0
	 RETURN
	ENDIF
 20	R1=PYR(0)
	T=TMAX*TMIN/(TMAX+R1*(TMIN-TMAX))
	R2=PYR(0)
	IF(R2.LT.ALPHAS(T,LQCD)**2/ALPHAS(TMIN,LQCD)**2)THEN
	 GETT=T-MDEB**2
	ELSE
	 GOTO 20
	ENDIF

	END



***********************************************************************
***	  function ei
***********************************************************************
      DOUBLE PRECISION FUNCTION EI(X)
      IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--exponential integral for negative arguments
      COMMON/EXPINT/EIXS(3,1000),VALMAX,NVAL
      INTEGER NVAL
      DOUBLE PRECISION EIXS,VALMAX
C--local variables
      INTEGER K,LINE,LMAX
      DOUBLE PRECISION X,R,GA,XA(2),YA(2),Y,DY,A,B
	DOUBLE PRECISION YSTART,EPSI,HFIRST
	DATA EPSI/1.e-5/
	
	IF(DABS(X).GT.VALMAX)
     &	write(logfid,*)'warning: value out of array in Ei(x)',X,VALMAX

      IF(X.GE.0.d0)THEN
       LMAX=INT(X*NVAL/VALMAX)
       LINE=MAX(LMAX,1)
       LINE=MIN(LINE,999)
	 IF((LINE.GT.999).OR.(LINE.LT.1)) THEN
        write(logfid,*)'ERROR in EI: line number out of bound',LINE
	 ENDIF
       DO 26 K=1,2
        XA(K)=EIXS(1,LINE-1+K)
        YA(K)=EIXS(3,LINE-1+K)
 26    CONTINUE
	 A=(YA(2)-YA(1))/(XA(2)-XA(1))
	 B=YA(1)-A*XA(1)
	 Y=A*X+B
      ELSE
       LMAX=INT(-X*NVAL/VALMAX)
       LINE=MAX(LMAX,1)
       LINE=MIN(LINE,999)
	 IF((LINE.GT.999).OR.(LINE.LT.1)) THEN
        write(logfid,*)'ERROR in EI: line number out of bound',LINE
	 ENDIF
       DO 27 K=1,2
        XA(K)=EIXS(1,LINE-1+K)
        YA(K)=EIXS(2,LINE-1+K)
 27    CONTINUE
	 A=(YA(2)-YA(1))/(XA(2)-XA(1))
	 B=YA(1)-A*XA(1)
	 Y=-A*X+B
      ENDIF
      EI=Y
      END



***********************************************************************
***	  function pqq
***********************************************************************
	DOUBLE PRECISION FUNCTION PQQ(Z)
	IMPLICIT NONE
	DOUBLE PRECISION Z
	PQQ=4.*(1.+Z**2)/(3.*(1.-Z))
	END



***********************************************************************
***	  function pgq
***********************************************************************
	DOUBLE PRECISION FUNCTION PGQ(Z)
	IMPLICIT NONE
	DOUBLE PRECISION Z
	PGQ=4.*(1.+(1.-Z)**2)/(3.*Z)
	END



***********************************************************************
***	  function pgg
***********************************************************************
	DOUBLE PRECISION FUNCTION PGG(Z)
	IMPLICIT NONE
	DOUBLE PRECISION Z
	PGG=3.*((1.-Z)/Z + Z/(1.-Z) + Z*(1.-Z))
	END



***********************************************************************
***	  function pqg
***********************************************************************
	DOUBLE PRECISION FUNCTION PQG(Z)
	IMPLICIT NONE
	DOUBLE PRECISION Z
	PQG=0.5*(Z**2 + (1.-Z)**2)
	END



***********************************************************************
***	  function alphas
***********************************************************************
	DOUBLE PRECISION FUNCTION ALPHAS(T,LAMBDA)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--local variables
	DOUBLE PRECISION T,L0,PI,LAMBDA
	DATA PI/3.141592653589793d0/

	 ALPHAS=4.*PI/((11.-2.*NF/3.)*LOG(T/LAMBDA**2))
	END



***********************************************************************
***	  subroutine splitfncint
***********************************************************************
	SUBROUTINE SPLITFNCINT(EMAX)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--splitting integral
      COMMON/SPLITINT/SPLITIGGV(1000,1000),SPLITIQQV(1000,1000),
     &SPLITIQGV(1000,1000),QVAL(1000),ZMVAL(1000),QMAX,ZMMIN,NPOINT
      INTEGER NPOINT
      DOUBLE PRECISION SPLITIGGV,SPLITIQQV,SPLITIQGV,
     &QVAL,ZMVAL,QMAX,ZMMIN
C--variables for splitting function integration
	COMMON/INTSPLITF/QQUAD,FM
	DOUBLE PRECISION QQUAD,FM
C--max rapidity
	common/rapmax/etamax
	double precision etamax
C--local variables
	INTEGER NSTEP,I,J
	DOUBLE PRECISION EMAX,ZMMAX,EPSI,HFIRST,YSTART,LNZMMIN,
     &LNZMMAX,ZM,ZM2,Q,GETMSMAX,avmom(5),shat,pcms2
      DATA ZMMAX/0.5/
      DATA NSTEP/999/
	DATA EPSI/1.d-5/

	call maxscatcen(avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))
	shat = avmom(5)**2 +
     &    2.*emax*(avmom(4)+sqrt(avmom(1)**2+avmom(2)**2+avmom(3)**2))
	pcms2 = (shat-avmom(5)**2)**2/(4.*shat)
	qmax = sqrt(scalefacm*4.*pcms2)

	ZMMIN=Q0/EMAX

      LNZMMIN=LOG(ZMMIN)
      LNZMMAX=LOG(ZMMAX)

	NPOINT=NSTEP

	DO 100 I=1,NSTEP+1
	 Q=(I-1)*(QMAX-Q0)/NSTEP+Q0
       QVAL(I)=Q
	 QQUAD=Q**2
       DO 110 J=1,NSTEP+1
        ZM=EXP((J-1)*(LNZMMAX-LNZMMIN)/NSTEP+LNZMMIN)
        ZMVAL(J)=ZM
	  IF(Q**2.LT.Q0**2)THEN
	   ZM2=0.5
	  ELSE 
	   ZM2=0.5-0.5*SQRT(1.-Q0**2/Q**2)
	  ENDIF 
	  ZM=MAX(ZM,ZM2)
	  IF(ZM.EQ.0.5)THEN	
	   SPLITIQQV(I,J)=0.d0
	   SPLITIGGV(I,J)=0.d0
	   SPLITIQGV(I,J)=0.d0
	  ELSE
	   YSTART=0d0
	   HFIRST=0.01
	   FM=0.d0
	   CALL ODEINT(YSTART,ZM,1.-ZM,EPSI,HFIRST,0d0,2)
	   SPLITIQQV(I,J)=YSTART
	   YSTART=0d0
	   HFIRST=0.01
	   FM=0.d0
	   CALL ODEINT(YSTART,ZM,1.-ZM,EPSI,HFIRST,0d0,3)
	   SPLITIGGV(I,J)=YSTART
	   YSTART=0d0
	   HFIRST=0.01
	   FM=0.d0
	   CALL ODEINT(YSTART,ZM,1.-ZM,EPSI,HFIRST,0d0,4)
	   SPLITIQGV(I,J)=YSTART
	  ENDIF
 110   CONTINUE
 100	CONTINUE

	END



***********************************************************************
***	  subroutine pdfint
***********************************************************************
	SUBROUTINE PDFINT(EMAX)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--pdf common block
	COMMON/PDFS/QINQX(2,1000),GINQX(2,1000),QINGX(2,1000),
     &GINGX(2,1000)
	DOUBLE PRECISION QINQX,GINQX,QINGX,GINGX
C--variables for pdf integration
	COMMON/PDFINTV/XMAX,Z
	DOUBLE PRECISION XMAX,Z
C--max rapidity
	common/rapmax/etamax
	double precision etamax
C--local variables
	INTEGER I,J
	DOUBLE PRECISION EMAX,Q2,GETPDFXINTEXACT,YSTART,HFIRST,EPSI,
     &Q2MAX,DELTAQ2,avmom(5),shat,pcms2
	DATA EPSI/1.d-4/

	call maxscatcen(avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))
	shat = avmom(5)**2 +
     &    2.*emax*(avmom(4)+sqrt(avmom(1)**2+avmom(2)**2+avmom(3)**2))
	pcms2 = (shat-avmom(5)**2)**2/(4.*shat)
	q2max = scalefacm*4.*pcms2

	DELTAQ2=LOG(Q2MAX)-LOG(Q0**2)
	QINQX(1,1)=Q0**2
	GINQX(1,1)=Q0**2
	QINGX(1,1)=Q0**2
	GINGX(1,1)=Q0**2
	QINQX(2,1)=0.d0
	GINQX(2,1)=0.d0
	QINGX(2,1)=0.d0
	GINGX(2,1)=0.d0
	 DO 12 J=2,1000
	  Q2 = EXP((J-1)*DELTAQ2/999.d0 + LOG(Q0**2))
	  QINQX(1,J)=Q2
	  GINQX(1,J)=Q2
	  QINGX(1,J)=Q2
	  GINGX(1,J)=Q2
	  QINQX(2,J)=GETPDFXINTEXACT(SQRT(Q2),'QQ')
	  GINQX(2,J)=GETPDFXINTEXACT(SQRT(Q2),'GQ')
	  QINGX(2,J)=GETPDFXINTEXACT(SQRT(Q2),'QG')
	  GINGX(2,J)=GETPDFXINTEXACT(SQRT(Q2),'GG')
 12	 CONTINUE
	END



***********************************************************************
***	  subroutine xsecint
***********************************************************************
	SUBROUTINE XSECINT(EMAX)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--cross secttion common block
	COMMON/XSECS/INTQ1(1001,101),INTQ2(1001,101),
     &INTG1(1001,101),INTG2(1001,101)
	DOUBLE PRECISION INTQ1,INTQ2,INTG1,INTG2
C--variables for cross section integration 
	COMMON/XSECV/QLOW,MDX
	DOUBLE PRECISION QLOW,MDX
C--max rapidity
	common/rapmax/etamax
	double precision etamax
C--local variables
	INTEGER J,K
	DOUBLE PRECISION EMAX,TMAX,TMAXMAX,DELTATMAX,YSTART,HFIRST,EPSI,
     &GETMSMAX,GETMDMAX,MDMIN,MDMAX,DELTAMD,GETMDMIN,avmom(5),shat,pcms2
	DATA EPSI/1.d-4/

	call maxscatcen(avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))
	shat = avmom(5)**2 +
     &    2.*emax*(avmom(4)+sqrt(avmom(1)**2+avmom(2)**2+avmom(3)**2))
	pcms2 = (shat-avmom(5)**2)**2/(4.*shat)
	tmaxmax = scalefacm*4.*pcms2
	DELTATMAX=(LOG(TMAXMAX)-
     &	LOG(Q0**2*(1.d0+1.d-6)/SCALEFACM**2))/999.d0
      MDMIN=GETMDMIN()
      MDMAX=MAX(MDMIN,GETMDMAX())
      DELTAMD=(MDMAX-MDMIN)/99.d0

	 DO 12 J=1,1000
	  TMAX = EXP((J-1)*DELTATMAX
     &	  + LOG(Q0**2*(1.d0+1.d-6)/SCALEFACM**2))
	  INTQ1(J,101)=TMAX
	  INTQ2(J,101)=TMAX
	  INTG1(J,101)=TMAX
	  INTG2(J,101)=TMAX
        DO 13 K=1,100
         MDX=MDMIN+(K-1)*DELTAMD
         INTQ1(1001,K)=MDX
         INTQ2(1001,K)=MDX
         INTG1(1001,K)=MDX
         INTG2(1001,K)=MDX
	  IF(TMAX.LT.Q0**2/SCALEFACM**2)THEN
	   INTQ1(J,K)=0.d0
	   INTQ2(J,K)=0.d0
	   INTG1(J,K)=0.d0
	   INTG2(J,K)=0.d0
	  ELSE
C--first quark integral
	   QLOW=Q0
  	   HFIRST=0.01*(TMAX-Q0**2/SCALEFACM**2)
         YSTART=0.d0
        CALL ODEINT(YSTART,Q0**2/SCALEFACM**2,TMAX,EPSI,HFIRST
     &        ,0.d0,11)
	   INTQ1(J,K)=YSTART
C--second quark integral
	   QLOW=Q0
  	   HFIRST=0.01*(TMAX-Q0**2/SCALEFACM**2)
         YSTART=0.d0
        CALL ODEINT(YSTART,Q0**2/SCALEFACM**2,TMAX,EPSI,HFIRST
     &        ,0.d0,14)
	   INTQ2(J,K)=YSTART
C--first gluon integral
	   QLOW=Q0
         YSTART=0.d0
        CALL ODEINT(YSTART,Q0**2/SCALEFACM**2,TMAX,EPSI,HFIRST
     &        ,0.d0,12)
	   INTG1(J,K)=YSTART
C--second gluon integral
	   QLOW=Q0
         YSTART=0.d0
        CALL ODEINT(YSTART,Q0**2/SCALEFACM**2,TMAX,EPSI,HFIRST
     &        ,0.d0,13)
	   INTG2(J,K)=YSTART
	  ENDIF
 13     CONTINUE
 12	 CONTINUE
	END



***********************************************************************
***	  function insudaint
***********************************************************************
	SUBROUTINE INSUDAINT(EMAX)
	IMPLICIT NONE
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--Sudakov common block
	COMMON/INSUDA/SUDAQQ(1000,2),SUDAQG(1000,2),SUDAGG(1000,2),
     &SUDAGC(1000,2)
	DOUBLE PRECISION SUDAQQ,SUDAQG,SUDAGG,SUDAGC
C--max rapidity
	common/rapmax/etamax
	double precision etamax
C--local variables
	INTEGER I
	DOUBLE PRECISION QMAX,Q,GETINSUDAKOV,DELTA,EMAX,avmom(5),
     &shat,pcms2
	
	call maxscatcen(avmom(1),avmom(2),avmom(3),avmom(4),avmom(5))
	shat = avmom(5)**2 +
     &    2.*emax*(avmom(4)+sqrt(avmom(1)**2+avmom(2)**2+avmom(3)**2))
	pcms2 = (shat-avmom(5)**2)**2/(4.*shat)
	qmax = sqrt(scalefacm*4.*pcms2)
	DELTA=(LOG(3.*QMAX)-LOG(Q0**2*(1.d0+1.d-6)))/999.d0
	DO 22 I=1,1000
	 Q = EXP((I-1)*DELTA + LOG(Q0**2*(1.d0+1.d-6)))
	 SUDAQQ(I,1)=Q
	 SUDAQG(I,1)=Q
	 SUDAGG(I,1)=Q
	 SUDAGC(I,1)=Q
	 SUDAQQ(I,2)=GETINSUDAKOV(Q0,Q,'QQ')
	 SUDAQG(I,2)=GETINSUDAKOV(Q0,Q,'QG')
	 SUDAGG(I,2)=GETINSUDAKOV(Q0,Q,'GG')
	 SUDAGC(I,2)=GETINSUDAKOV(Q0,Q,'GC')
 22	CONTINUE
	END



***********************************************************************
***	  function eixint
***********************************************************************
	SUBROUTINE EIXINT
	IMPLICIT NONE
C--exponential integral for negative arguments
      COMMON/EXPINT/EIXS(3,1000),VALMAX,NVAL
      INTEGER NVAL
      DOUBLE PRECISION EIXS,VALMAX
C-local variables
	INTEGER I,K
	DOUBLE PRECISION X,EPSI,HFIRST,YSTART,EI,GA,R 
	DATA	EPSI/1.d-5/

	NVAL=1000
	VALMAX=55.

      DO 10 I=1,NVAL
       X=I*VALMAX/(NVAL*1.d0)
       EIXS(1,I)=X
C--do negative arguments first
	 YSTART=0d0
	 HFIRST=0.01
	 CALL ODEINT(YSTART,X,1000.d0,EPSI,HFIRST,0.d0,5)
       EIXS(2,I)=-YSTART
C--now do the positive arguments
       IF (X.EQ.0.0) THEN
        EI=-1.0D+300
       ELSE IF (X.LE.40.0) THEN
        EI=1.0D0
        R=1.0D0
        DO 15 K=1,100
         R=R*K*X/(K+1.0D0)**2
         EI=EI+R
         IF (DABS(R/EI).LE.1.0D-15) GO TO 20
15      CONTINUE
20      GA=0.5772156649015328D0
        EI=GA+DLOG(X)+X*EI
       ELSE
        EI=1.0D0
        R=1.0D0
        DO 25 K=1,20
         R=R*K/X
25       EI=EI+R
         EI=DEXP(X)/X*EI
       ENDIF
	 EIXS(3,I)=EI
 10   CONTINUE
	END



***********************************************************************
***	  function odeint
***********************************************************************
	subroutine odeint(ystart,a,b,eps,h1,hmin,w1)
	implicit none
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--local variables
	integer nmax,nstep,w1
	double precision ystart,a,b,eps,h1,hmin,x,h,y,dydx,
     &deriv,yscale,hdid,hnew
	data nmax/100000/

	x = a
	y = ystart
	h = sign(h1,b-a)
	do 20 nstep=1,nmax
	  dydx = deriv(x,w1)
	  yscale = abs(y) + abs(h*dydx) + 1.e-25
	  if (((x + h - b)*h).gt.0.) h = b-x
	  call rkstepper(x,y,dydx,h,hdid,hnew,yscale,eps,w1)
	  if ((x - b)*h.ge.0) then
	    ystart = y
	    return
	  endif
	  h = hnew
	  if (abs(h).lt.abs(hmin)) then
	    write(logfid,*)'Error in odeint: stepsize too small',w1
     &	,ystart,a,b,h1
	    return
	  endif	  
 20	continue
	write(logfid,*)'Error in odeint: too many steps',w1
     &	,ystart,a,b,h1
	end



***********************************************************************
***	  function rkstepper
***********************************************************************
	subroutine rkstepper(x,y,dydx,htest,hdid,hnew,yscale,eps,w1)
	implicit none
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--local variables
	integer w1
	double precision x,y,dydx,htest,hdid,hnew,yscale,eps,
     &yhalf,y1,y2,rk4step,dydxhalf,xnew,delta,err,h,safety, powerdown,
     &powerup,maxup,maxdown,deriv,fac
	logical reject
	data powerdown/0.25/
	data powerup/0.2/
	data safety/0.9/
	data maxdown/10./
	data maxup/5./

	reject = .false.
	h = htest
 10	xnew = x + h
	if (x.eq.xnew) then
	  write(logfid,*)'Error in rkstepper: step size not significant'
	  return
	endif
	yhalf = rk4step(x,y,dydx,h/2.,w1)
	dydxhalf = deriv(x+h/2.,w1)
	y2 = rk4step(x+h/2.,yhalf,dydxhalf,h/2.,w1)
	y1 = rk4step(x,y,dydx,h,w1)
	delta = y2-y1
	err = abs(delta)/(yscale*eps)
	if (err.gt.1.) then
	  reject = .true.
	  fac = max(1./maxdown,safety/err**powerdown)
	  h = h*fac
	  goto 10 
	else
	  if (reject) then
	    hnew = h
	  else
	    fac = min(maxup,safety/err**powerup)
	    hnew = fac*h
	  endif
	  x = xnew
	  y = y2 + delta/15.
	  hdid = h
	endif
	end



***********************************************************************
***	  function rk4step
***********************************************************************
	double precision function rk4step(x,y,dydx,h,w1)
	implicit none
	integer w1
	double precision x,y,dydx,h,k1,k2,k4,yout,deriv
	k1 = h*dydx
	k2 = h*deriv(x+h/2.,w1)
	k4 = h*deriv(x+h,w1)
	yout = y+k1/6.+2.*k2/3.+k4/6.
	rk4step = yout
	end



***********************************************************************
***	  function getdeltat
***********************************************************************
      LOGICAL FUNCTION GETDELTAT(LINE,TSTART,DTMAX1,DELTAT)
      IMPLICIT NONE
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--pythia common block
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--max rapidity
	common/rapmax/etamax
	double precision etamax
C--memory for error message from getdeltat
	common/errline/errl
	integer errl
C--local variables
      INTEGER LINE,I,NNULL
      DOUBLE PRECISION DTMAX,SIGMAMAX,NEFFMAX,LINVMAX,PYR,
     &R,TOFF,XS,YS,ZS,TS,GETSSCAT,GETMSMAX,GETMDMIN,MSMAX,MDMIN,
     &XSTART,YSTART,ZSTART,WEIGHT,MS,MD,NEFF,SIGMA,GETNEFF,
     &GETNEFFMAX,GETMS,GETMD,TAU,MDMAX,GETMDMAX,GETNATMDMIN,
     &SIGMAMIN,NEFFMIN,TSTART,DTMAX1,DELTAT
	CHARACTER PTYPE
	LOGICAL STOPNOW

C--initialization
	GETDELTAT=.FALSE.
      DELTAT=0.D0
	DTMAX=DTMAX1
	IF(K(LINE,2).EQ.21)THEN
	 PTYPE='G'
	ELSE
	 PTYPE='Q'
	ENDIF

	NNULL=0
	STOPNOW=.FALSE.

C--check for upper bound from plasma lifetime
      IF((TSTART+DTMAX).GE.LTIME)DTMAX=LTIME-TSTART
      IF(DTMAX.LT.0.D0) RETURN
	
C--calculate time relative to production of the considered parton
      TOFF=TSTART-MV(LINE,4)
	XSTART=MV(LINE,1)+TOFF*P(LINE,1)/P(LINE,4)
	YSTART=MV(LINE,2)+TOFF*P(LINE,2)/P(LINE,4)
	ZSTART=MV(LINE,3)+TOFF*P(LINE,3)/P(LINE,4)

C--calculate upper limit for density*cross section
	SIGMAMAX=GETSSCAT(P(LINE,4),p(line,1),p(line,2),p(line,3),
!     &	xstart,ystart,-sign(abs(zstart),p(line,3)),zstart+1.d-6)
     &	P(LINE,5),0.d0,PTYPE,'C',xstart,ystart,zstart,tstart,1)
	SIGMAMIN=GETSSCAT(P(LINE,4),p(line,1),p(line,2),p(line,3),
!     &	xstart,ystart,-sign(abs(zstart),p(line,3)),zstart+1.d-6)
     &	P(LINE,5),0.d0,PTYPE,'C',xstart,ystart,zstart,tstart,2)
	NEFFMAX=GETNEFFMAX()
	NEFFMIN=GETNATMDMIN()
	LINVMAX=5.d0*MAX(NEFFMIN*SIGMAMAX,NEFFMAX*SIGMAMIN)
	if(linvmax.eq.0.d0) return

	DO 333 I=1,1000000
	 DELTAT=DELTAT-LOG(PYR(0))/LINVMAX
	 XS=XSTART+DELTAT*P(LINE,1)/P(LINE,4)
	 YS=YSTART+DELTAT*P(LINE,2)/P(LINE,4)
	 ZS=ZSTART+DELTAT*P(LINE,3)/P(LINE,4)
	 TS=TSTART+DELTAT
	 IF(TS.LT.ZS)THEN
	  TAU=-1.d0
	 ELSE
	  TAU=SQRT(TS**2-ZS**2)
	 ENDIF
	 NEFF=GETNEFF(XS,YS,ZS,TS)
	 IF((TAU.GT.1.d0).AND.(NEFF.EQ.0.d0))THEN
	  IF(NNULL.GT.4)THEN
	   STOPNOW=.TRUE.
	  ELSE 
	   NNULL=NNULL+1
	  ENDIF
	 ELSE
	  NNULL=0
	 ENDIF
	 IF((DELTAT.GT.DTMAX).OR.STOPNOW) THEN
	  DELTAT=DTMAX
	  RETURN
	 ENDIF
	 IF(NEFF.GT.0.d0)THEN
	  SIGMA=GETSSCAT(P(LINE,4),p(line,1),p(line,2),p(line,3),
     &	P(LINE,5),0.d0,PTYPE,'C',xs,ys,zs,ts,0)
	 ELSE
	  SIGMA=0.d0
	 ENDIF
	 WEIGHT=5.d0*NEFF*SIGMA/LINVMAX
	 IF(WEIGHT.GT.1.d0+1d-6) then
	   if (line.ne.errl) then
     	     write(logfid,*)'error in GETDELTAT: weight > 1',WEIGHT,
     &	 NEFF*SIGMA/(NEFFMAX*SIGMAMIN),NEFF*SIGMA/(NEFFMIN*SIGMAMAX),
     &       p(line,4)
	     errl=line
	   endif
	 endif
       R=PYR(0)
	 IF(R.LT.WEIGHT)THEN
	  GETDELTAT=.TRUE.
	  RETURN
	 ENDIF
 333	CONTINUE
	END


	integer function poissonian(lambda)
	implicit none
	integer n
	double precision lambda,disc,p,pyr,u,v,pi
	data pi/3.141592653589793d0/
	
	if (lambda.gt.745.d0) then
	  u = pyr(0);
	  v = pyr(0);
	  poissonian = 
     &	int(sqrt(lambda)*sqrt(-2.*log(u))*cos(2.*pi*v)+lambda)
	else
	 disc=exp(-lambda)
	 p=1.d0
	 n=0	
 800   p = p*pyr(0)
	 if (p.gt.disc) then
	   n = n+1
	   goto 800
	 endif
	 poissonian=n
	endif
	end


***********************************************************************
***	  function ishadron
***********************************************************************
	LOGICAL FUNCTION ISHADRON(ID)
	IMPLICIT NONE
C--local variables
	INTEGER ID	
	IF(ABS(ID).LT.100) THEN
	 ISHADRON=.FALSE.
	ELSE
	 IF(MOD(INT(ABS(ID)/10.),10).EQ.0) THEN
	  ISHADRON = .FALSE.
	 ELSE
	  ISHADRON = .TRUE.
       ENDIF
      ENDIF
      END



***********************************************************************
***	  function isdiquark
***********************************************************************
	LOGICAL FUNCTION ISDIQUARK(ID)
	IMPLICIT NONE
C--local variables
	INTEGER ID	
	IF(ABS(ID).LT.1000) THEN
	 ISDIQUARK=.FALSE.
	ELSE 
	 IF(MOD(INT(ID/10),10).EQ.0) THEN
	  ISDIQUARK = .TRUE.
	 ELSE
	  ISDIQUARK = .FALSE.
       ENDIF
      ENDIF 
      END

***********************************************************************
***	  function islepton
***********************************************************************
      LOGICAL FUNCTION ISLEPTON(ID)
      IMPLICIT NONE
C--   local variables
      INTEGER ID
      IF((ABS(ID).EQ.11).OR.(ABS(ID).EQ.13).OR.(ABS(ID).EQ.15)) THEN
         ISLEPTON=.TRUE.
      ELSE
         ISLEPTON=.FALSE.
      ENDIF
      END
      
***********************************************************************
***	  function isparton
***********************************************************************
	LOGICAL FUNCTION ISPARTON(ID)
	IMPLICIT NONE
C--local variables
	INTEGER ID	
	LOGICAL ISDIQUARK
	IF((ABS(ID).LT.6).OR.(ID.EQ.21).OR.ISDIQUARK(ID)) THEN
	 ISPARTON=.TRUE.
	ELSE 
	 ISPARTON=.FALSE.
      ENDIF 
      END      



***********************************************************************
***	  function isprimstring
***********************************************************************
      logical function isprimstring(l)
      implicit none
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--local variables
	integer l
	logical isparton
	if ((K(l,2).ne.91).and.(K(l,2).ne.92)) then
	  isprimstring=.false.
	  return
	endif
	if ((K(K(l,3),3).eq.0).or.(isparton(K(K(K(l,3),3),2)))) then
        isprimstring=.true.
	else 
        isprimstring=.false.
	endif
	end



***********************************************************************
***	  function issecstring
***********************************************************************
      logical function issecstring(l)
      implicit none
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--local variables
	integer l
	logical isparton,isprimstring
	if ((K(l,2).ne.91).and.(K(l,2).ne.92)) then
	  issecstring = .false.
	  return
	endif
	if (isprimstring(l)) then
	  issecstring = .false.
	  return
	endif
	if (isparton(K(K(K(l,3),3),2))) then 
	  issecstring = .false.
	else
	  issecstring = .true.
	endif
	end



***********************************************************************
***	  function isprimhadron
***********************************************************************
      logical function isprimhadron(l)
      implicit none
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--local variables
	integer l
	logical isprimstring,isparton
	if (((K(K(l,3),2).EQ.91).OR.(K(K(l,3),2).EQ.92))
     &	.and.isprimstring(K(l,3))
     &	.and.(.not.isparton(K(l,2)))) then
	  isprimhadron=.true.
	else 
        isprimhadron=.false.
	endif
	if (k(l,1).eq.17) isprimhadron=.true.
	end



***********************************************************************
***	  function compressevent
***********************************************************************
	logical function compressevent(l1)
	implicit none
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
C--local variables
	integer l1,i,j,nold,nnew,nstart
	
	nold = n

	do 777 i=2,nold
	  if (((k(i,1).eq.11).or.(k(i,1).eq.12).or.(k(i,1).eq.13)).and.
     &	(i.ne.l1)) then
	    nnew = i
	    goto 778
	  endif
 777	continue
	compressevent = .false.
	return
 778	continue
	nstart = nnew
	do 779 i=nstart,nold
	  if (((k(i,1).ne.11).and.(k(i,1).ne.12).and.(k(i,1).ne.13)).or.
     &	(i.eq.l1)) then
	    do 780 j=1,5
	      p(nnew,j)=p(i,j)
	      v(nnew,j)=v(i,j)
	      mv(nnew,j)=mv(i,j)
 780	    continue
	    trip(nnew)=trip(i)
	    anti(nnew)=anti(i)
	    za(nnew)=za(i)
	    zd(nnew)=zd(i)
	    thetaa(nnew)=thetaa(i)
	    qqbard(nnew)=qqbard(i)
	    k(nnew,1)=k(i,1)
	    k(nnew,2)=k(i,2)
	    k(nnew,3)=0
	    k(nnew,4)=0
	    k(nnew,5)=0
	    if (l1.eq.i) l1=nnew
	    nnew=nnew+1
	  endif
 779	continue
	n=nnew-1
	if ((nold-n).le.10) then
	  compressevent = .false.
	else
	  compressevent = .true.
	endif
	do 781 i=nnew,nold
	  do 782 j=1,5
	    k(i,j)=0
	    p(i,j)=0.d0
	    v(i,j)=0.d0
	    mv(i,j)=0.d0
 782	  continue
	  trip(i)=0
	  anti(i)=0
	  za(i)=0.d0
	  zd(i)=0.d0
	  thetaa(i)=0.d0
	  qqbard(i)=.false.
 781	continue
	if (n.gt.23000) write(logfid,*)'Error in compressevent: n = ',n 
	if (l1.gt.n) write(logfid,*)'Error in compressevent: l1 = ',l1  
	call flush(logfid)
	return
	end



***********************************************************************
***	  subroutine pevrec
***********************************************************************
      SUBROUTINE PEVREC(NUM,COL)
C--identifier of file for hepmc output and logfile
	implicit none
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
C--variables for angular ordering
      COMMON/ANGOR/ZA(23000),ZD(23000),THETAA(23000),QQBARD(23000)
	DOUBLE PRECISION ZA,ZD,THETAA
      LOGICAL QQBARD
C--time common block
      COMMON/TIME/MV(23000,5)
      DOUBLE PRECISION MV
C--colour index common block
	COMMON/COLOUR/TRIP(23000),ANTI(23000),COLMAX
	INTEGER TRIP,ANTI,COLMAX
	INTEGER NUM,i
	LOGICAL COL

      DO 202 I=1,N
       V(I,1)=MV(I,1)
       V(I,2)=MV(I,2)
       V(I,3)=MV(I,3)
       V(I,4)=MV(I,4)
       V(I,5)=MV(I,5)
	 IF(COL) write(logfid,*)I,' (',TRIP(I),',',ANTI(I),')    [',
     &K(I,3),K(I,4),K(I,5),' ]  {',K(I,2),K(I,1),' } ',	 
     &ZD(I),THETAA(I)
 202  CONTINUE
      CALL PYLIST(NUM)

      END



***********************************************************************
***	  subroutine converttohepmc
***********************************************************************
	SUBROUTINE CONVERTTOHEPMC(J,EVNUM,PID,beam1,beam2)
	IMPLICIT NONE
      COMMON/PYJETS/N,NPAD,K(23000,5),P(23000,5),V(23000,5)
	INTEGER N,NPAD,K
	DOUBLE PRECISION P,V
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
C--Parameter common block
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,ANGORD,SCATRECOIL,
     &ALLHAD,compress,NF
      INTEGER NF
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM
      LOGICAL ANGORD,SCATRECOIL,ALLHAD,compress
C--organisation of event record
	common/evrecord/nsim,npart,offset,hadrotype,sqrts,collider,hadro,
     &shorthepmc,channel,isochannel
	integer nsim,npart,offset,hadrotype
	double precision sqrts
	character*4 collider,channel
	character*2 isochannel
	logical hadro,shorthepmc
C--extra storage for scattering centres before interactions
      common/storescatcen/nscatcen,maxnscatcen,scatflav(10000),
     &scatcen(10000,5),writescatcen,writedummies
	integer nscatcen,maxnscatcen,scatflav
	double precision scatcen
	logical writescatcen,writedummies
C--local variables
	INTEGER EVNUM,PBARCODE,VBARCODE,CODELIST(25000),I,PID,NSTART,
     &NFIRST,NVERTEX,NTOT,J,CODEFIRST
      DOUBLE PRECISION mproton,mneutron,pdummy,pscatcen
      LOGICAL ISHADRON,ISDIQUARK,ISPARTON,isprimhadron,isprimstring,
     &issecstring
	character*2 beam1,beam2
	data mproton/0.9383/
	data mneutron/0.9396/
	data pdummy/1.d-6/  
	
 5000 FORMAT(A2,I10,I3,3E14.6,2I2,I6,4I2,E14.6)
 5100 FORMAT(A2,2E14.6)
 5200 FORMAT(A2,6I7,2I2,1I7,4E14.6)
 5300 FORMAT(A2,2I2,5E14.6,2I2)
 5400 FORMAT(A2,I6,6I2,I6,I2)
 5500 FORMAT(A2,I6,I6,5E14.6,3I2,I6,I2)

	PBARCODE=0
	VBARCODE=0

	if (shorthepmc) then
C--short output
        IF(COLLIDER.EQ.'EEJJ')THEN
          NVERTEX=3
	    PBARCODE=5
        ELSE
          NVERTEX=1
	    PBARCODE=2
        ENDIF
	  nfirst = 0
	  do 131 i=1,N
	    if (((k(i,1).lt.6).or.(k(i,1).eq.17)))
     &	nfirst = nfirst+1
 131	  continue
	  if(writescatcen) NFIRST=NFIRST+nscatcen
	  if(writedummies) NFIRST=NFIRST+nscatcen

	  WRITE(J,5000)'E ',EVNUM,-1,0.d0,0.d0,0.d0,0,0,NVERTEX,1,2,0,1,
     &PARI(10)
	  WRITE(J,'(A2,I2,A5)')'N ',1,'"0"' 
	  WRITE(J,'(A)')'U GEV MM'
	  WRITE(J,5100)'C ',PARI(1)*1.d9,0.d0
  	  WRITE(J,5200)'H ',0,0,0,0,0,0,0,0,0,0.d0,0.d0,0.d0,0.d0
	  WRITE(J,5300)'F ',0,0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,0,0
C--write out vertex line	  
	  IF(COLLIDER.EQ.'EEJJ')THEN
	    WRITE(J,5400)'V ',-1,0,0,0,0,0,2,1,0
	    WRITE(J,5500)'P ',1,-11,0.d0,0.d0,sqrts/2.,sqrts/2.,
     &	0.00051,2,0,0,-1,0
	    WRITE(J,5500)'P ',2,11,0.d0,0.d0,-sqrts/2.,sqrts/2.,
     &	0.00051,2,0,0,-1,0
	    WRITE(J,5500)'P ',3,23,0.d0,0.d0,0.d0,sqrts,
     &	91.2,2,0,0,-2,0
	    WRITE(J,5400)'V ',-2,0,0,0,0,0,0,2,0
	    WRITE(J,5500)'P ',4,PID,sqrts/2.,0.d0,0.d0,sqrts/2.,
     &	0.000,2,0,0,-3,0
	    WRITE(J,5500)'P ',5,-PID,-sqrts/2.,0.d0,0.d0,sqrts/2.,
     &	0.000,2,0,0,-3,0
	    WRITE(J,5400)'V ',-3,0,0,0,0,0,0,NFIRST,0
        ELSE
	    WRITE(J,5400)'V ',-1,0,0,0,0,0,2,NFIRST,0
	    if (beam1.eq.'p+') then
	  	WRITE(J,5500)'P ',1,2212,0.d0,0.d0,
     &	sqrt(sqrts**2/4.-mproton**2),sqrts/2.,mproton,2,0,0,-1,0
	    else
	  	WRITE(J,5500)'P ',1,2112,0.d0,0.d0,
     &	sqrt(sqrts**2/4.-mneutron**2),sqrts/2.,mneutron,2,0,0,-1,0
	    endif
	    if (beam2.eq.'p+') then
	      WRITE(J,5500)'P ',2,2212,0.d0,0.d0,
     &	-sqrt(sqrts**2/4.-mproton**2),sqrts/2.,mproton,2,0,0,-1,0
	    else
	      WRITE(J,5500)'P ',2,2112,0.d0,0.d0,
     &	-sqrt(sqrts**2/4.-mneutron**2),sqrts/2.,mneutron,2,0,0,-1,0
	    endif
	  ENDIF
C--write out scattering centres
	if(writescatcen) then
	    do 133 i=1,nscatcen
	      pbarcode=pbarcode+1
	      WRITE(J,5500)'P ',pbarcode,scatflav(i),scatcen(I,1),
     &	  scatcen(I,2),scatcen(I,3),scatcen(I,4),scatcen(I,5),
     &	  3,0,0,0,0
 133	    continue
	  endif	  
C--write out dummy particles
	  if(writedummies) then
	    do 135 i=1,nscatcen
	      pbarcode=pbarcode+1
	      pscatcen=sqrt(scatcen(I,1)**2+scatcen(I,2)**2+
     &		scatcen(I,3)**2)
	      WRITE(J,5500)'P ',pbarcode,111,pdummy*scatcen(I,1)/pscatcen,
     &	  pdummy*scatcen(I,2)/pscatcen,pdummy*scatcen(I,3)/pscatcen,
     &	  pdummy,0.d0,1,0,0,0,0
 135	    continue
	  endif	  
C--write out particle lines
	  do 132 i=1,N
	    if(((k(i,1).lt.6).or.(k(i,1).eq.17))) then
	      pbarcode=pbarcode+1
		if((k(i,1).eq.3).or.(k(i,1).eq.5)) then
	        WRITE(J,5500)'P ',PBARCODE,K(I,2),P(I,1),P(I,2),P(I,3),
     &		P(I,4),P(I,5),4,0,0,0,0
	      else
	        WRITE(J,5500)'P ',PBARCODE,K(I,2),P(I,1),P(I,2),P(I,3),
     &		P(I,4),P(I,5),1,0,0,0,0
		endif
	    endif
 132	  continue

	else
C--long output
	  if (hadro) then
C--hadronised events
	    NFIRST=0
          IF(COLLIDER.EQ.'EEJJ')THEN
            NVERTEX=3
          ELSE
            NVERTEX=1
          ENDIF
	    DO 123 I=1,N
	      IF(K(i,3).ne.0)THEN
	        NSTART=I
	        GOTO 124
	      ENDIF
 123	    CONTINUE	 
 124	    CONTINUE	 
	    nstart=0

          DO 126 I=NSTART+1,N
	      IF(isprimhadron(i)) NFIRST=NFIRST+1
	      IF((ISHADRON(K(I,2)).OR.(ABS(K(I,2)).EQ.15))
     &	  .AND.(K(I,4).NE.0)) NVERTEX=NVERTEX+1
 126	    CONTINUE	 
 127	    CONTINUE	 

	    if(writescatcen) NFIRST=NFIRST+nscatcen
	    if(writedummies) NFIRST=NFIRST+nscatcen

	    WRITE(J,5000)'E ',EVNUM,-1,0.d0,0.d0,0.d0,0,0,NVERTEX,
     &1,2,0,1,PARI(10)
	    WRITE(J,'(A2,I2,A5)')'N ',1,'"0"' 
	    WRITE(J,'(A)')'U GEV MM'
	    WRITE(J,5100)'C ',PARI(1)*1.d9,0.d0
	    WRITE(J,5200)'H ',0,0,0,0,0,0,0,0,0,0.d0,0.d0,0.d0,0.d0
	    WRITE(J,5300)'F ',0,0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,0,0

C--write out vertex line	  
          IF(COLLIDER.EQ.'EEJJ')THEN
	      VBARCODE=-3
	      PBARCODE=5
	    ELSE
	      VBARCODE=-1
	      PBARCODE=2
	    ENDIF
	    IF(COLLIDER.EQ.'EEJJ')THEN
	      WRITE(J,5400)'V ',-1,0,0,0,0,0,2,1,0
	      WRITE(J,5500)'P ',1,-11,0.d0,0.d0,sqrts/2.,sqrts/2.,
     &	0.00051,2,0,0,-1,0
	      WRITE(J,5500)'P ',2,11,0.d0,0.d0,-sqrts/2.,sqrts/2.,
     &	0.00051,2,0,0,-1,0
	      WRITE(J,5500)'P ',3,23,0.d0,0.d0,0.d0,sqrts,
     &	91.2,2,0,0,-2,0
	      WRITE(J,5400)'V ',-2,0,0,0,0,0,0,2,0
	      WRITE(J,5500)'P ',4,PID,sqrts/2.,0.d0,0.d0,sqrts/2.,
     &	0.000,2,0,0,-3,0
		WRITE(J,5500)'P ',5,-PID,-sqrts/2.,0.d0,0.d0,sqrts/2.,
     &	0.000,2,0,0,-3,0
		WRITE(J,5400)'V ',VBARCODE,0,0,0,0,0,0,NFIRST,0
          ELSE
	      WRITE(J,5400)'V ',-1,0,0,0,0,0,2,NFIRST,0
	    if (beam1.eq.'p+') then
	  	WRITE(J,5500)'P ',1,2212,0.d0,0.d0,
     &	sqrt(sqrts**2/4.-mproton**2),sqrts/2.,mproton,2,0,0,-1,0
	    else
	  	WRITE(J,5500)'P ',1,2112,0.d0,0.d0,
     &	sqrt(sqrts**2/4.-mneutron**2),sqrts/2.,mneutron,2,0,0,-1,0
	    endif
	    if (beam2.eq.'p+') then
	      WRITE(J,5500)'P ',2,2212,0.d0,0.d0,
     &	-sqrt(sqrts**2/4.-mproton**2),sqrts/2.,mproton,2,0,0,-1,0
	    else
	      WRITE(J,5500)'P ',2,2112,0.d0,0.d0,
     &	-sqrt(sqrts**2/4.-mneutron**2),sqrts/2.,mneutron,2,0,0,-1,0
	    endif
	    ENDIF
       
	    CODEFIRST=NFIRST+PBARCODE

C--write out scattering centres
	  if(writescatcen) then
	    do 134 i=1,nscatcen
	      pbarcode=pbarcode+1
	      WRITE(J,5500)'P ',PBARCODE,scatflav(I),scatcen(I,1),
     &	  scatcen(I,2),scatcen(I,3),scatcen(I,4),scatcen(I,5),
     &	  3,0,0,0,0
 134	    continue
	  endif	  
C--write out dummy particles
	  if(writedummies) then
	    do 136 i=1,nscatcen
	      pbarcode=pbarcode+1
	      pscatcen=sqrt(scatcen(I,1)**2+scatcen(I,2)**2+
     &		scatcen(I,3)**2)
	      WRITE(J,5500)'P ',pbarcode,111,pdummy*scatcen(I,1)/pscatcen,
     &	  pdummy*scatcen(I,2)/pscatcen,pdummy*scatcen(I,3)/pscatcen,
     &	  pdummy,0.d0,1,0,0,0,0
 136	    continue
	  endif	  

C--first write out all particles coming directly from string or cluster decays
	     DO 125 I=NSTART+1,N
	       IF(.not.isprimhadron(i))THEN
	         GOTO 125
	       ELSE
	         IF (PBARCODE.EQ.CODEFIRST) GOTO 130
	         PBARCODE=PBARCODE+1
C--write out particle line	  
	         IF(K(I,4).GT.0)THEN
	           VBARCODE=VBARCODE-1
	           CODELIST(I)=VBARCODE
	          WRITE(J,5500)'P ',PBARCODE,K(I,2),P(I,1),P(I,2),P(I,3),
     &	     P(I,4),P(I,5),2,0,0,VBARCODE,0
	         ELSE 
	          WRITE(J,5500)'P ',PBARCODE,K(I,2),P(I,1),P(I,2),P(I,3),
     &	     P(I,4),P(I,5),1,0,0,0,0
	         ENDIF	    
	       ENDIF   
 125	     CONTINUE	   
 130	     CONTINUE	
C--now write out all other particles and vertices	
	     DO 129 I=NSTART+1,N
	       if (isprimhadron(i).or.isprimstring(i)) goto 129
	       if (isparton(K(i,2))) then
	         if (ishadron(K(K(i,3),2))) codelist(i)=codelist(K(i,3))
	         goto 129
	       endif
	       if (issecstring(i)) then
	         codelist(i)=codelist(K(i,3))
	         goto 129
	       endif
	       PBARCODE=PBARCODE+1
	       IF((K(I,3).NE.K(I-1,3)))THEN
C--write out vertex line	  
	         WRITE(J,5400)'V ',CODELIST(K(I,3)),0,0,0,0,0,0,
     &    		K(K(I,3),5)-K(K(I,3),4)+1,0
	       ENDIF 
C--write out particle line	  
	       IF(K(I,4).GT.0)THEN
	         VBARCODE=VBARCODE-1
	         CODELIST(I)=VBARCODE
	         WRITE(J,5500)'P ',PBARCODE,K(I,2),P(I,1),P(I,2),P(I,3),
     &		P(I,4),P(I,5),2,0,0,VBARCODE,0
	       ELSE 
	         WRITE(J,5500)'P ',PBARCODE,K(I,2),P(I,1),P(I,2),P(I,3),
     &		P(I,4),P(I,5),1,0,0,0,0
	       ENDIF	    
 129	     CONTINUE

	  else
C--partonic events
	  endif
	endif
	call flush(j)
	END
	


***********************************************************************
***	  subroutine printlogo
***********************************************************************
	subroutine printlogo(fid)
	implicit none
	integer fid

	write(fid,*)
	write(fid,*)'                   _______________'//
     &'__________________________                  '
	write(fid,*)'                  |               '//
     &'                          |                 '
	write(fid,*)'                  |  JJJJJ  EEEEE '//
     &' W       W  EEEEE  L      |                  '
	write(fid,*)'                  |      J  E     '//
     &' W       W  E      L      |                  '
	write(fid,*)' _________________|      J  EEE   '//
     &'  W  W  W   EEE    L      |_________________ '
	write(fid,*)'|                 |  J   J  E     '//
     &'  W W W W   E      L      |                 |'
	write(fid,*)'|                 |   JJJ   EEEEE '//
     &'   W   W    EEEEE  LLLLL  |                 |'
	write(fid,*)'|                 |_______________'//
     &'__________________________|                 |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'|                            '//
     &'this is JEWEL 2.1.0                              |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'| Copyright Korinna C. Zapp (2016)'//
     &'  [Korinna.Zapp@cern.ch]                    |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'| The JEWEL homepage is jewel.hepforge.org '//
     &'                                   |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'| The medium model was partly '//
     &'implemented by Jochen Klein                     |'
	write(fid,*)'| [Jochen.Klein@cern.ch]. Raghav '//
     &'Kunnawalkam Elayavalli helped with the       |'
	write(fid,*)'| implementation of the V+jet processes '//
     &'[raghav.k.e@cern.ch].                 |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'| Please cite JHEP 1303 (2013) '//
     &'080 [arXiv:1212.1599] and optionally           |'
	write(fid,*)'| EPJC C60 (2009) 617 [arXiv:0804.3568] '//
     &'for the physics and arXiv:1311.0048   |'
	write(fid,*)'| for the code. The reference for '//
     &'V+jet processes is EPJC 76 (2016) no.12 695 |'
       write(fid,*)'| [arXiv:1608.03099] and for recoil effects'//
     &' it is arXiv:1707.01539.           |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'| JEWEL contains code provided by '//
     &'S. Zhang and J. M. Jing                     |'
	write(fid,*)'| (Computation of Special Functions, '//
     &'John Wiley & Sons, New York, 1996 and    |'
	write(fid,*)'| http://jin.ece.illinois.edu) for '//
     &'computing the exponential integral Ei(x).  |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'| JEWEL relies heavily on PYTHIA 6'//
     &' for the event generation. The modified     |'
	write(fid,*)'| version of PYTHIA 6.4.25 that is'//
     &' shipped with JEWEL is, however, not an     |'
	write(fid,*)'| official PYTHIA release and must'//
     &' not be used for anything else. Please      |'
	write(fid,*)'| refer to results as "JEWEL+PYTHIA".'//
     &'                                         |'
	write(fid,*)'|                                 '//
     &'                                            |'
	write(fid,*)'|_________________________________'//
     &'____________________________________________|'
	write(fid,*)
	write(fid,*)
	end


***********************************************************************
***	  subroutine printtime
***********************************************************************
	subroutine printtime
	implicit none
C--identifier of file for hepmc output and logfile
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
C--local variables
	integer*4 date(3),time(3)

 1000 format (i2.2, '.', i2.2, '.', i4.4, ', ',
     &         i2.2, ':', i2.2, ':', i2.2 )
	call idate(date)
	call itime(time)
	write(logfid,1000)date,time
	end


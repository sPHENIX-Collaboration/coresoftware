C  jewel-api.f -- the PHJewel entry points into JEWEL 2.6.0.
C
C  PROGRAM JEWEL is   init(); per isospin channel: initpythia(b1,b2),
C  genevent(j,b1,b2) N times; finish.  These routines expose exactly
C  those pieces, so the caller (JewelCore.h) owns the event loop.
C  Strings never cross the C/Fortran boundary: the caller fills the
C  /JWAPI/ and /JWAPIC/ commons directly (jewel_capi.cc).
C
C  Called only from jewel_capi.cc.  Not part of the exported interface.

***********************************************************************
***	  JWFINIT: JEWEL init() in API mode
***********************************************************************
	SUBROUTINE JWFINIT()
	IMPLICIT NONE
	common/jwapi/apimode,apiwrite,apiseed
	logical apimode,apiwrite
	integer apiseed
	common/jwapich/cb1,cb2,ichan
	character*2 cb1,cb2
	integer ichan
	apimode = .true.
	ichan = 0
	cb1 = '  '
	cb2 = '  '
	call init()
	END

***********************************************************************
***	  JWFCHAN: switch isospin channel (1 pp, 2 pn, 3 np, 4 nn)
***	  as the channel loop of PROGRAM JEWEL does
***********************************************************************
	SUBROUTINE JWFCHAN(ICH)
	IMPLICIT NONE
	INTEGER ICH
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
	common/jwapich/cb1,cb2,ichan
	character*2 cb1,cb2
	integer ichan
	call jwfchanlog()
	if (ich.eq.1) then
	  cb1 = 'p+'
	  cb2 = 'p+'
	elseif (ich.eq.2) then
	  cb1 = 'p+'
	  cb2 = 'n0'
	elseif (ich.eq.3) then
	  cb1 = 'n0'
	  cb2 = 'p+'
	else
	  cb1 = 'n0'
	  cb2 = 'n0'
	endif
	ichan = ich
	write(logfid,*)
	write(logfid,*)
     &'####################################################'
	write(logfid,*)
	write(logfid,*)'PHJewel: generating events in ',
     &cb1,' + ',cb2,' channel'
	write(logfid,*)
	write(logfid,*)
     &'####################################################'
	write(logfid,*)
	SUMOFWEIGHTS=0.d0
	WDISC=0.d0
	call initpythia(cb1,cb2)
	write(logfid,*)
	call flush(logfid)
	END

***********************************************************************
***	  JWFNEXT: one pass of the event loop (may produce no event)
***********************************************************************
	SUBROUTINE JWFNEXT(IEV)
	IMPLICIT NONE
	INTEGER IEV
	common/jwapich/cb1,cb2,ichan
	character*2 cb1,cb2
	integer ichan
	call genevent(IEV,cb1,cb2)
	END

***********************************************************************
***	  JWFINFO: per-event quantities of the event just generated
***********************************************************************
	SUBROUTINE JWFINFO(D)
	IMPLICIT NONE
	DOUBLE PRECISION D(16)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
      common/storescatcen/nscatcen,maxnscatcen,scatflav(23000),
     &scatcen(23000,5),writescatcen,writedummies,dosubtraction,
     &putonshell,Rsub
	integer nscatcen,maxnscatcen,scatflav,putonshell
	double precision scatcen,Rsub
	logical writescatcen,writedummies,dosubtraction
	common/jwapich/cb1,cb2,ichan
	character*2 cb1,cb2
	integer ichan
	double precision getb,getcentrality
	D(1) = PARI(10)
	D(2) = PARI(1)
	D(3) = getb()
	D(4) = getcentrality()
	D(5) = MSTI(1)
	D(6) = MSTI(15)
	D(7) = MSTI(16)
	D(8) = PARI(33)
	D(9) = PARI(34)
	D(10) = PARI(17)
	D(11) = PARI(22)
	D(12) = NGOOD
	D(13) = nscatcen
	D(14) = ichan
	D(15) = 0.d0
	D(16) = 0.d0
	END

***********************************************************************
***	  JWFCOUNTS: run bookkeeping
***********************************************************************
	SUBROUTINE JWFCOUNTS(D)
	IMPLICIT NONE
	DOUBLE PRECISION D(16)
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
	common/totnumevents/nevdone
	integer nevdone
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	D(1) = NGOOD
	D(2) = NDISC
	D(3) = NSTRANGE
	D(4) = nevdone
	D(5) = sumofweights
	D(6) = wdisc
	D(7) = ntotspliti
	D(8) = noverspliti
	D(9) = ntotpdf
	D(10) = noverpdf
	D(11) = ntotxsec
	D(12) = noverxsec
	D(13) = ntotsuda
	D(14) = noversuda
	D(15) = 0.d0
	D(16) = 0.d0
	END

***********************************************************************
***	  JWFCONFIG: what JEWEL itself parsed from the card, so the caller
***	  can check it against its own reading of the same file
***********************************************************************
	SUBROUTINE JWFCONFIG(IC,DC)
	IMPLICIT NONE
	INTEGER IC(16)
	DOUBLE PRECISION DC(8)
	common/np/nproton,mass
	integer nproton,mass
	common/evrecord/nsim,npart,offset,hadrotype,sqrts,collider,hadro,
     &shorthepmc,isrscat,channel,isochannel
	integer nsim,npart,offset,hadrotype
	double precision sqrts
	character*4 collider,channel
	character*2 isochannel
	logical hadro,shorthepmc,isrscat
	common/pythiaparams/PTMIN,PTMAX,weighted,mpi,pi0dec
	double precision PTMIN,PTMAX
	LOGICAL WEIGHTED,mpi,pi0dec
	common/rapmax/etamax
	double precision etamax
	COMMON/PARAM/Q0,LPS,LQCD,LTIME,SCALEFACM,
     &RECSOFTCUT,RECHARDCUT,
     &ANGORD,aofirstsplit,SCATRECOIL,ALLHAD,compress,mpifsr,
     &NF,KINMODE,recmode
      INTEGER ANGORD,NF,KINMODE,recmode
	DOUBLE PRECISION Q0,LQCD,LTIME,LPS,SCALEFACM,
     &RECSOFTCUT,RECHARDCUT
      LOGICAL SCATRECOIL,ALLHAD,compress,mpifsr,aofirstsplit
      common/storescatcen/nscatcen,maxnscatcen,scatflav(23000),
     &scatcen(23000,5),writescatcen,writedummies,dosubtraction,
     &putonshell,Rsub
	integer nscatcen,maxnscatcen,scatflav,putonshell
	double precision scatcen,Rsub
	logical writescatcen,writedummies,dosubtraction
	integer i
	do 10 i=1,16
	  IC(i) = 0
 10	continue
	IC(1) = mass
	IC(2) = nproton
	IC(3) = -1
	if (isochannel.eq.'XX') IC(3) = 0
	if (isochannel.eq.'PP') IC(3) = 1
	if (isochannel.eq.'PN') IC(3) = 2
	if (isochannel.eq.'NP') IC(3) = 3
	if (isochannel.eq.'NN') IC(3) = 4
C--PROCESS as four character codes packed into one integer
	IC(4) = ichar(collider(1:1))*16777216
     &	+ ichar(collider(2:2))*65536
     &	+ ichar(collider(3:3))*256 + ichar(collider(4:4))
	if (weighted) IC(5) = 1
	if (allhad) IC(6) = 1
	if (writescatcen) IC(7) = 1
	if (writedummies) IC(8) = 1
	if (dosubtraction) IC(9) = 1
	if (hadro) IC(10) = 1
	if (shorthepmc) IC(11) = 1
	if (isrscat) IC(12) = 1
	if (mpi) IC(13) = 1
	IC(14) = NF
	DC(1) = sqrts
	DC(2) = PTMIN
	DC(3) = PTMAX
	DC(4) = etamax
	DC(5) = Q0
	DC(6) = LQCD
	DC(7) = LPS
	DC(8) = 0.d0
	END

***********************************************************************
***	  JWFFINISH: the end of PROGRAM JEWEL (without out.ran)
***********************************************************************
	SUBROUTINE JWFFINISH()
	IMPLICIT NONE
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
	common/evinfoid/evinfid,doevinfo
	integer evinfid
	logical doevinfo
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
	common/extrapolations/ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	integer ntotspliti,noverspliti,ntotpdf,noverpdf,
     &ntotxsec,noverxsec,ntotsuda,noversuda
	common/jwapi/apimode,apiwrite,apiseed
	logical apimode,apiwrite
	integer apiseed
C--close the last channel's bookkeeping in the log
	call jwfchanlog()
	if (apiwrite) then
	  WRITE(HPMCFID,'(A)')'HepMC::IO_GenEvent-END_EVENT_LISTING'
	  WRITE(HPMCFID,*)
	  CLOSE(HPMCFID,status='keep')
	endif
	if (doevinfo) CLOSE(EVINFID,status='keep')
	write(logfid,*)
	write(logfid,*)'number of extrapolations in splitting integral: ',
     &	noverspliti,' (',(noverspliti*1.d0)/(ntotspliti*1.d0),'%)'
	write(logfid,*)'number of extrapolations in pdfs: ',
     &	noverpdf,' (',(noverpdf*1.d0)/(ntotpdf*1.d0),'%)'
	write(logfid,*)
     &	'number of extrapolations in splitting cross sections: ',
     &	noverxsec,' (',(noverxsec*1.d0)/(ntotxsec*1.d0),'%)'
	write(logfid,*)
     &	'number of extrapolations in Sudakov form factor: ',
     &	noversuda,' (',(noversuda*1.d0)/(ntotsuda*1.d0),'%)'
	write(logfid,*)
	write(logfid,*)'total number of discarded events: ',NDISC
	write(logfid,*)'number of events for which conversion '//
     &'to hepmc failed: ',NSTRANGE
	call printtime
	close(logfid,status='keep')
	END

	SUBROUTINE JWFCHANLOG()
	IMPLICIT NONE
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
	INTEGER MSTP,MSTI
	DOUBLE PRECISION PARP,PARI
	common/hepmcid/hpmcfid,logfid
	integer hpmcfid,logfid
	COMMON/DISC/NDISC,NSTRANGE,NGOOD,errcount,wdisc,DISCARD
	LOGICAL DISCARD
	INTEGER NDISC,NSTRANGE,NGOOD,errcount
	double precision wdisc
	COMMON/WEIGHT/EVWEIGHT,sumofweights
	double precision EVWEIGHT,sumofweights
	common/jwapich/cb1,cb2,ichan
	character*2 cb1,cb2
	integer ichan
	if (ichan.le.0) return
	write(logfid,*)
	write(logfid,*)'cross section in ',cb1,' + ',cb2,' channel:',
     &	PARI(1),'mb'
	write(logfid,*)'sum of event weights in ',cb1,' + ',cb2,
     &	' channel:',sumofweights-wdisc
	write(logfid,*)
	END

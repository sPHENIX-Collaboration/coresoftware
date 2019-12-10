
C
C
C
C
	SUBROUTINE HIJFRG(JTP,NTP,IERROR)
C	NTP=1, fragment proj string, NTP=2, targ string, 
C       NTP=3, independent 
C	strings from jets.  JTP is the line number of the string
C*******Fragment all leading strings of proj and targ**************
C	IHNT2(1)=atomic #, IHNT2(2)=proton #(=-1 if anti-proton)  *
C******************************************************************
	COMMON/HIPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
	SAVE  /HIPARNT/
        COMMON/HIJDAT/HIDAT0(10,10),HIDAT(10)
	SAVE  /HIJDAT/
	COMMON/HISTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
	SAVE  /HISTRNG/
	COMMON/HIJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
	SAVE  /HIJJET1/
	COMMON/HIJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
     &		K2SG(900,100),PXSG(900,100),PYSG(900,100),
     &		PZSG(900,100),PESG(900,100),PMSG(900,100)
	SAVE  /HIJJET2/
C
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
	SAVE  /LUJETS/
	COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
	SAVE  /LUDAT1/
	COMMON/RANSEED/NSEED
	SAVE  /RANSEED/
	
	IERROR=0
	CALL LUEDIT(0)
	N=0
C			********initialize the document lines
	IF(NTP.EQ.3) THEN
		ISG=JTP
		N=NJSG(ISG)
		DO 100 I=1,NJSG(ISG)
			K(I,1)=K1SG(ISG,I)
			K(I,2)=K2SG(ISG,I)
			P(I,1)=PXSG(ISG,I)
			P(I,2)=PYSG(ISG,I)
			P(I,3)=PZSG(ISG,I)
			P(I,4)=PESG(ISG,I)
			P(I,5)=PMSG(ISG,I)
C BAC+++
C
C       Clear the starting point information in the Pythia arrays
C
			V(I,1)=0
			V(I,2)=0
			V(I,3)=0
			V(I,4)=0
			V(I,5)=0
C BAC---

100		CONTINUE
C		IF(IHPR2(1).GT.0) CALL ATTRAD(IERROR)
c		IF(IERROR.NE.0) RETURN
C		CALL LULIST(1)
		CALL LUEXEC
		RETURN
	ENDIF
C
	IF(NTP.EQ.2) GO TO 200
	IF(JTP.GT.IHNT2(1))   RETURN
	IF(NFP(JTP,5).NE.3.AND.NFP(JTP,3).NE.0
     &	    .AND.NPJ(JTP).EQ.0.AND.NFP(JTP,10).EQ.0) GO TO 1000
	IF(NFP(JTP,15).EQ.-1) THEN
		KF1=NFP(JTP,2)
		KF2=NFP(JTP,1)
		PQ21=PP(JTP,6)
		PQ22=PP(JTP,7)
		PQ11=PP(JTP,8)
		PQ12=PP(JTP,9)
		AM1=PP(JTP,15)
		AM2=PP(JTP,14)
	ELSE
		KF1=NFP(JTP,1)
		KF2=NFP(JTP,2)
		PQ21=PP(JTP,8)
		PQ22=PP(JTP,9)
		PQ11=PP(JTP,6)
		PQ12=PP(JTP,7)
		AM1=PP(JTP,14)
		AM2=PP(JTP,15)	
	ENDIF
C	********for NFP(JTP,15)=-1 NFP(JTP,1) IS IN -Z DIRECTION
	PB1=PQ11+PQ21
	PB2=PQ12+PQ22
	PB3=PP(JTP,3)
	PECM=PP(JTP,5)
	BTZ=PB3/PP(JTP,4)
	IF((ABS(PB1-PP(JTP,1)).GT.0.01.OR.
     &    ABS(PB2-PP(JTP,2)).GT.0.01).AND.IHPR2(10).NE.0)
     &	  WRITE(6,*) '  Pt of Q and QQ do not sum to the total'

	GO TO 300

200	IF(JTP.GT.IHNT2(3))  RETURN
	IF(NFT(JTP,5).NE.3.AND.NFT(JTP,3).NE.0
     &	   .AND.NTJ(JTP).EQ.0.AND.NFT(JTP,10).EQ.0) GO TO 1200
	IF(NFT(JTP,15).EQ.1) THEN
		KF1=NFT(JTP,1)
		KF2=NFT(JTP,2)
		PQ11=PT(JTP,6)
		PQ12=PT(JTP,7)
		PQ21=PT(JTP,8)
		PQ22=PT(JTP,9)
		AM1=PT(JTP,14)
		AM2=PT(JTP,15)
	ELSE
		KF1=NFT(JTP,2)
		KF2=NFT(JTP,1)
		PQ11=PT(JTP,8)
		PQ12=PT(JTP,9)
		PQ21=PT(JTP,6)
		PQ22=PT(JTP,7)
		AM1=PT(JTP,15)
		AM2=PT(JTP,14)
	ENDIF	
C	********for NFT(JTP,15)=1 NFT(JTP,1) IS IN +Z DIRECTION
	PB1=PQ11+PQ21
	PB2=PQ12+PQ22
	PB3=PT(JTP,3)
	PECM=PT(JTP,5)
	BTZ=PB3/PT(JTP,4)

	IF((ABS(PB1-PT(JTP,1)).GT.0.01.OR.
     &     ABS(PB2-PT(JTP,2)).GT.0.01).AND.IHPR2(10).NE.0)
     &     WRITE(6,*) '  Pt of Q and QQ do not sum to the total'

300	IF(PECM.LT.HIPR1(1)) THEN
	   IERROR=1
	   IF(IHPR2(10).EQ.0) RETURN
	   WRITE(6,*) ' ECM=',PECM,' energy of the string is too small'
	   RETURN
	ENDIF
	AMT=PECM**2+PB1**2+PB2**2
	AMT1=AM1**2+PQ11**2+PQ12**2
	AMT2=AM2**2+PQ21**2+PQ22**2
	PZCM=SQRT(ABS(AMT**2+AMT1**2+AMT2**2-2.0*AMT*AMT1
     &       -2.0*AMT*AMT2-2.0*AMT1*AMT2))/2.0/SQRT(AMT)
C		*******PZ of end-partons in c.m. frame of the string
	K(1,1)=2
	K(1,2)=KF1
	P(1,1)=PQ11
	P(1,2)=PQ12
	P(1,3)=PZCM
	P(1,4)=SQRT(AMT1+PZCM**2)
	P(1,5)=AM1
	K(2,1)=1
	K(2,2)=KF2
	P(2,1)=PQ21
	P(2,2)=PQ22
	P(2,3)=-PZCM
	P(2,4)=SQRT(AMT2+PZCM**2)
	P(2,5)=AM2

C BAC+++
C
C       Clear the starting point information in the Pythia arrays
C
	V(1,1)=0
	V(1,2)=0
	V(1,3)=0
	V(1,4)=0
	V(1,5)=0

	V(2,1)=0
	V(2,2)=0
	V(2,3)=0
	V(2,4)=0
	V(2,5)=0
C BAC---


	N=2
C*****
	CALL HIROBO(0.0,0.0,0.0,0.0,BTZ)
	JETOT=0
	IF((PQ21**2+PQ22**2).GT.(PQ11**2+PQ12**2)) THEN
		PMAX1=P(2,1)
		PMAX2=P(2,2)
		PMAX3=P(2,3)
	ELSE
		PMAX1=P(1,1)
		PMAX2=P(1,2)
		PMAX3=P(1,3)
	ENDIF
	IF(NTP.EQ.1) THEN
		PP(JTP,10)=PMAX1
		PP(JTP,11)=PMAX2
		PP(JTP,12)=PMAX3
	ELSE IF(NTP.EQ.2) THEN
		PT(JTP,10)=PMAX1
		PT(JTP,11)=PMAX2
		PT(JTP,12)=PMAX3
	ENDIF
C*******************attach produced jets to the leading partons****
	IF(NTP.EQ.1.AND.NPJ(JTP).NE.0) THEN
		JETOT=NPJ(JTP)
C		IF(NPJ(JTP).GE.2) CALL HIJSRT(JTP,1)
C			********sort jets in order of y
		IEX=0
		IF((ABS(KF1).GT.1000.AND.KF1.LT.0)
     &			.OR.(ABS(KF1).LT.1000.AND.KF1.GT.0)) IEX=1
		DO 520 I=N,2,-1
		DO 520 J=1,5
			II=NPJ(JTP)+I
			K(II,J)=K(I,J)
			P(II,J)=P(I,J)
			V(II,J)=V(I,J)
520		CONTINUE
		DO 540 I=1,NPJ(JTP)
			DO 542 J=1,5
				K(I+1,J)=0
				V(I+1,J)=0
542			CONTINUE				
			I0=I
			IF(IEX.EQ.1) I0=NPJ(JTP)-I+1
C				********reverse the order of jets
			KK1=KFPJ(JTP,I0)
			K(I+1,1)=2
			K(I+1,2)=KK1
			IF(KK1.NE.21 .AND. KK1.NE.0)  K(I+1,1)=
     &			  1+(ABS(KK1)+(2*IEX-1)*KK1)/2/ABS(KK1)
			P(I+1,1)=PJPX(JTP,I0)
			P(I+1,2)=PJPY(JTP,I0)
			P(I+1,3)=PJPZ(JTP,I0)
			P(I+1,4)=PJPE(JTP,I0)
			P(I+1,5)=PJPM(JTP,I0)
540		CONTINUE
		N=N+NPJ(JTP)
	ELSE IF(NTP.EQ.2.AND.NTJ(JTP).NE.0) THEN
		JETOT=NTJ(JTP)
c		IF(NTJ(JTP).GE.2)  CALL HIJSRT(JTP,2)
C			********sort jets in order of y
		IEX=1
		IF((ABS(KF2).GT.1000.AND.KF2.LT.0)
     &			.OR.(ABS(KF2).LT.1000.AND.KF2.GT.0)) IEX=0
		DO 560 I=N,2,-1
		DO 560 J=1,5
			II=NTJ(JTP)+I
			K(II,J)=K(I,J)
			P(II,J)=P(I,J)
			V(II,J)=V(I,J)
560		CONTINUE
		DO 580 I=1,NTJ(JTP)
			DO 582 J=1,5
				K(I+1,J)=0
				V(I+1,J)=0
582			CONTINUE				
			I0=I
			IF(IEX.EQ.1) I0=NTJ(JTP)-I+1
C				********reverse the order of jets
			KK1=KFTJ(JTP,I0)
			K(I+1,1)=2
			K(I+1,2)=KK1
			IF(KK1.NE.21 .AND. KK1.NE.0) K(I+1,1)=
     &		           1+(ABS(KK1)+(2*IEX-1)*KK1)/2/ABS(KK1)
			P(I+1,1)=PJTX(JTP,I0)
			P(I+1,2)=PJTY(JTP,I0)
			P(I+1,3)=PJTZ(JTP,I0)
			P(I+1,4)=PJTE(JTP,I0)
			P(I+1,5)=PJTM(JTP,I0)
580		CONTINUE
		N=N+NTJ(JTP)
	ENDIF
	IF(IHPR2(1).GT.0.AND.ATL_RAN(NSEED).LE.HIDAT(3)) THEN
	     HIDAT20=HIDAT(2)
	     HIPR150=HIPR1(5)
	     IF(IHPR2(8).EQ.0.AND.IHPR2(3).EQ.0.AND.IHPR2(9).EQ.0)
     &			HIDAT(2)=2.0
	     IF(HINT1(1).GE.1000.0.AND.JETOT.EQ.0)THEN
		HIDAT(2)=3.0
		HIPR1(5)=5.0
	     ENDIF
	     CALL ATTRAD(IERROR)
	     HIDAT(2)=HIDAT20
	     HIPR1(5)=HIPR150
	ELSE IF(JETOT.EQ.0.AND.IHPR2(1).GT.0.AND.
     &                       HINT1(1).GE.1000.0.AND.
     &		ATL_RAN(NSEED).LE.0.8) THEN
		HIDAT20=HIDAT(2)
		HIPR150=HIPR1(5)
		HIDAT(2)=3.0
		HIPR1(5)=5.0
	     IF(IHPR2(8).EQ.0.AND.IHPR2(3).EQ.0.AND.IHPR2(9).EQ.0)
     &			HIDAT(2)=2.0
		CALL ATTRAD(IERROR)
		HIDAT(2)=HIDAT20
		HIPR1(5)=HIPR150
	ENDIF
	IF(IERROR.NE.0) RETURN
C		******** conduct soft radiations
C****************************
C
C
C	CALL LULIST(1)
	CALL LUEXEC
	RETURN

1000	N=1
	K(1,1)=1
       	K(1,2)=NFP(JTP,3)
	DO 1100 JJ=1,5
       		P(1,JJ)=PP(JTP,JJ)
C BAC+++
C
C       Clear the starting point information in the Pythia arrays
C
       		V(1,JJ)=0
C BAC---
		
1100		CONTINUE
C			********proj remain as a nucleon or delta
	CALL LUEXEC
C	call lulist(1)
	RETURN
C
1200	N=1
	K(1,1)=1
	K(1,2)=NFT(JTP,3)
	DO 1300 JJ=1,5
		P(1,JJ)=PT(JTP,JJ)
C BAC+++
C
C       Clear the starting point information in the Pythia arrays
C
       		V(1,JJ)=0
C BAC---

1300	CONTINUE
C			********targ remain as a nucleon or delta
	CALL LUEXEC
C	call lulist(1)
	RETURN
	END

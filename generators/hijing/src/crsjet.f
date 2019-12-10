C
C
C	THIS PROGRAM IS TO CALCULATE THE JET CROSS SECTION
C	THE INTEGRATION IS DONE BY USING VEGAS
C
	SUBROUTINE CRSJET
	IMPLICIT REAL*8(A-H,O-Z)
	REAL HIPR1(100),HINT1(100)
        COMMON/HIPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
        SAVE  /HIPARNT/
	COMMON/NJET/N,IP_CRS
	SAVE  /NJET/
	COMMON/BVEG1/XL(10),XU(10),ACC,NDIM,NCALL,ITMX,NPRN
	SAVE  /BVEG1/
	COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI,NDO,IT
	SAVE  /BVEG2/
	COMMON/BVEG3/F,TI,TSI
	SAVE  /BVEG3/
	COMMON/SEEDVAX/NUM1
	SAVE  /SEEDVAX/
	EXTERNAL FJET,FJETRIG
C
c************************
c	NCALL give the number of inner-iteration, ITMX 
C       gives the limit of out-iteration. Nprn is an option
C       ( 1: print the integration process. 0: do not print)
C

C 	print *, 'ncall = ', ncall
C 	print *, 'itmx = ', itmx
C 	print *, 'nprn = ', nprn
C 	print *, 'acc = ', acc
C 	print *, 'xl(1), xu(1) = ', XL(1), XU(1)
C 	print *, 'xl(2), xu(2) = ', XL(2), XU(2)
C 	print *, 'xl(3), xu(3) = ', XL(3), XU(3)

C +++BAC
C
C  The following line inserted to improve the accuracy of the jet cross-section
C    integration for LHC purposes where integration region runs out to XT=1 where the
C    contribution is infinitesimal. Found to be necessary because large jet cross-section
C    errors were affecting total and inelastic cross-sections.
C
	ncall = 4000

C ---BAC

	NDIM=3
	IP_CRS=0
	CALL VEGAS(FJET,AVGI,SD,CHI2A)
	HINT1(14)=AVGI/2.5682
	IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
		IP_CRS=1
		CALL VEGAS(FJET,AVGI,SD,CHI2A)
		HINT1(15)=AVGI/2.5682
	ENDIF
	IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
		IP_CRS=2
		CALL VEGAS(FJET,AVGI,SD,CHI2A)
		HINT1(16)=AVGI/2.5682
	ENDIF
	IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
		IP_CRS=3
		CALL VEGAS(FJET,AVGI,SD,CHI2A)
		HINT1(17)=AVGI/2.5682
	ENDIF
C		********Total inclusive jet cross section(Pt>P0) 
C
	IF(IHPR2(3).NE.0) THEN
	   IP_CRS=0
	   CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	   HINT1(61)=AVGI/2.5682
	   IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
	      IP_CRS=1
	      CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	      HINT1(62)=AVGI/2.5682
	   ENDIF
	   IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
	      IP_CRS=2
	      CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	      HINT1(63)=AVGI/2.5682
	   ENDIF
	   IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
	      IP_CRS=3
	      CALL VEGAS(FJETRIG,AVGI,SD,CHI2A)
	      HINT1(64)=AVGI/2.5682
	   ENDIF
	ENDIF
C			********cross section of trigger jet
C
	RETURN
	END

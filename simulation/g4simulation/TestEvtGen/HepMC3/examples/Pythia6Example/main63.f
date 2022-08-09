C...A simple skeleton program, illustrating a typical Pythia run:
C...Z0 production at LEP 1.
C...Toy task: compare multiplicity distribution with matrix elements
C...and with parton showers (using same fragmentation parameters).
C...This code contains modifications for HepMC3 examples
C-----------------------------------------------------------------
      PROGRAM MAIN
C...Preamble: declarations.

C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
C...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...Parameters.
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
C...Generation and cross section statistics.
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
C...Random number generator information.
      COMMON/PYDATR/MRPY(6),RRPY(100)

C...HepMC3
      PARAMETER (NMXHEP=4000)
      COMMON /HEPEVT/  NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                 JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &                 VHEP(4,NMXHEP)
      INTEGER          NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION PHEP,VHEP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      include "Pythia6ToHepMC3.inc"
      INTEGER OUTID(2), HEPMC3STATUS

C-----------------------------------------------------------------

C...First section: initialization.

C...Main parameters of run: c.m. energy and number of events.
      ECM=91.2D0
      NEV=1000

C...Select gamma*/Z0 production process.
      MSEL=0
      MSUB(1)=1

C...Only allow Z0 decay to quarks (i.e. no leptonic final states).
      DO 100 IDC=MDCY(23,2),MDCY(23,2)+MDCY(23,3)-1
        IF(IABS(KFDP(IDC,1)).GE.6) MDME(IDC,1)=MIN(0,MDME(IDC,1))
  100 CONTINUE

C...Initialize.
      CALL PYINIT('CMS','e+','e-',ECM)

C...Check that Z0 decay channels set correctly.
C      CALL PYSTAT(2)

C...Book histograms.
      CALL PYBOOK(1,'charged multiplicity ME',100,-0.5D0,99.5D0)
      CALL PYBOOK(2,'charged multiplicity PS',100,-0.5D0,99.5D0)
C...Create output writers
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      OUTID(1)=hepmc3_new_writer(0,1,'ME.hepmc'//char(0))
      HEPMC3STATUS=hepmc3_new_weight(OUTID(1),'Default'//char(0))
      HEPMC3STATUS=hepmc3_new_weight(OUTID(1),'weme1'//char(0))
      HEPMC3STATUS=hepmc3_new_weight(OUTID(1),'weme2'//char(0))
      OUTID(2)=hepmc3_new_writer(0,1,'PS.hepmc'//char(0))
      HEPMC3STATUS=hepmc3_new_weight(OUTID(2),'Default'//char(0))
      HEPMC3STATUS=hepmc3_new_weight(OUTID(2),'weps1'//char(0))
      HEPMC3STATUS=hepmc3_new_weight(OUTID(2),'weps2'//char(0))
      NEVHEP=-123456
      HEPMC3STATUS=hepmc3_set_hepevt_address(NEVHEP)
C...Or one can set the pointer to some predefined block size      
C      HEPMC3STATUS=hepmc3_set_hepevt_address(NEVHEPL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------------------------------------------

C...Second section: event loop.

C...Outer loop over ME and PS options.
      DO 300 ICA=1,2
        IF(ICA.EQ.1) THEN
          MSTP(48)=1
          MSTJ(101)=2
        ELSE
          MSTP(48)=0
        ENDIF


C...Begin event loop.
        DO 200 IEV=1,NEV
          CALL PYEVNT

C...List first few events.
          IF(IEV.LE.2) CALL PYLIST(1)

C...Write output
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          CALL PYHEPC(1)
          NEVHEP=IEV
C...One can copy to some predefined block size
C          NEVHEPL=NEVHEP
C          NHEPL=NHEP
C           DO 500 J=1,NHEP
C          ISTHEPL(J)=ISTHEP(J)
C          IDHEPL(J)=IDHEP(J)
C          JMOHEPL(1,J)=JMOHEP(1,J)
C          JMOHEPL(2,J)=JMOHEP(2,J)
C          JDAHEPL(1,J)=JDAHEP(1,J)
C          JDAHEPL(2,J)=JDAHEP(2,J)
C          PHEPL(1,J)=PHEP(1,J)
C          PHEPL(2,J)=PHEP(2,J)
C          PHEPL(3,J)=PHEP(3,J)
C          PHEPL(4,J)=PHEP(4,J)
C          PHEPL(5,J)=PHEP(5,J)
C          VHEPL(1,J)=VHEP(1,J)
C          VHEPL(2,J)=VHEP(2,J)
C          VHEPL(3,J)=VHEP(3,J)
C          VHEPL(4,J)=VHEP(4,J)
C  500     CONTINUE
          HEPMC3STATUS=hepmc3_convert_event(OUTID(ICA))
C...Note: no explicit XS uncertainty
          HEPMC3STATUS=hepmc3_set_cross_section(OUTID(ICA),
     &    1.0E9*XSEC(0,3),
     &    1.0E9*XSEC(0,3)/sqrt(1.0*NGEN(0,3)),
     &    NGEN(0,3),0)
          HEPMC3STATUS=hepmc3_set_pdf_info(OUTID(ICA),
     &    MSTI(15),MSTI(16),PARI(33),PARI(34),PARI(23),
     &    MSTP(51),MSTP(52))
C...The values below are not always meaningful
          HEPMC3STATUS=hepmc3_set_attribute_int(OUTID(ICA),-1,
     &   'mpi'//char(0))
          HEPMC3STATUS=hepmc3_set_attribute_int(OUTID(ICA),MSUB(1),
     &   'signal_process_id'//char(0))
          HEPMC3STATUS=hepmc3_set_attribute_int(OUTID(ICA),MRPY(1),
     &   'random_states1'//char(0))
          HEPMC3STATUS=hepmc3_set_attribute_double(OUTID(ICA),-1.0D0,
     &   'alphaEM'//char(0))
          HEPMC3STATUS=hepmc3_set_attribute_double(OUTID(ICA),-1.0D0,
     &   'alphaQCD'//char(0))
          HEPMC3STATUS=hepmc3_set_attribute_double(OUTID(ICA),q2pdfeval,
     &   'event_scale'//char(0))
          HEPMC3STATUS=hepmc3_set_weight_by_index(OUTID(ICA),1.0D0,0)
          if (ICA.eq.1) then
          HEPMC3STATUS=hepmc3_set_weight_by_name(OUTID(1),
     &    1.1111D0,'weme1'//char(0))
          HEPMC3STATUS=hepmc3_set_weight_by_name(OUTID(1),
     &    1.2222D0,'weme2'//char(0))
          endif
          if (ICA.eq.2) then
          HEPMC3STATUS=hepmc3_set_weight_by_name(OUTID(2),
     &    1.5555D0,'weps1'//char(0))
          HEPMC3STATUS=hepmc3_set_weight_by_name(OUTID(2),
     &    1.6666D0,'weps2'//char(0))
          endif

C Note there should be PDF ids
          HEPMC3STATUS=hepmc3_write_event(OUTID(ICA))
          HEPMC3STATUS=hepmc3_clear_event(OUTID(ICA))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Extract and fill event properties.
          CALL PYEDIT(3)
          CALL PYFILL(ICA,DBLE(N),1D0)
C...End event loop.
  200   CONTINUE

C...End outer loop.
  300 CONTINUE
C...Delete output writers
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      HEPMC3STATUS=hepmc3_delete_writer(OUTID(1))
      HEPMC3STATUS=hepmc3_delete_writer(OUTID(2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.
      CALL PYSTAT(1)

C...Histograms.
      CALL PYHIST

      END

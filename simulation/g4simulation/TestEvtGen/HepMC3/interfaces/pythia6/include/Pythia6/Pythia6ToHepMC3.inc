CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Copyright (C) 2014-2020 The HepMC collaboration
C...These variables will be used for conversion. Block size is different
C...in Pythia6 and in HepMC3, so the most simple portable way is to have
C... a second block of same size as in HepMC3 and  copy the content of
C...block directly.
      INTEGER NMXHEPL
      PARAMETER (NMXHEPL=10000)
      COMMON /HEPEVTL/  NEVHEPL,NHEPL,ISTHEPL(NMXHEPL),IDHEPL(NMXHEPL),
     &           JMOHEPL(2,NMXHEPL),JDAHEPL(2,NMXHEPL),PHEPL(5,NMXHEPL),
     &                 VHEPL(4,NMXHEPL)
      INTEGER          NEVHEPL,NHEPL,ISTHEPL,IDHEPL,JMOHEPL,JDAHEPL
      DOUBLE PRECISION PHEPL,VHEPL


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Glue functions
      INTEGER hepmc3_new_writer,hepmc3_delete_writer
      INTEGER hepmc3_set_hepevt_address
      INTEGER hepmc3_convert_event,hepmc3_write_event,hepmc3_clear_event
      INTEGER hepmc3_set_attribute_int,hepmc3_set_attribute_double
      INTEGER hepmc3_set_pdf_info,hepmc3_set_cross_section
      INTEGER hepmc3_new_weight
      INTEGER hepmc3_set_weight_by_index
      INTEGER hepmc3_set_weight_by_name
      EXTERNAL hepmc3_new_writer,hepmc3_delete_writer
      EXTERNAL hepmc3_set_hepevt_address
      EXTERNAL hepmc3_convert_event,hepmc3_write_event
      EXTERNAL hepmc3_clear_event
      EXTERNAL hepmc3_set_attribute_int,hepmc3_set_attribute_double
      EXTERNAL hepmc3_set_pdf_info,hepmc3_set_cross_section
      EXTERNAL hepmc3_new_weight
      EXTERNAL hepmc3_set_weight_by_index
      EXTERNAL hepmc3_set_weight_by_name
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

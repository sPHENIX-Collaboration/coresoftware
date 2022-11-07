AUTOMAKE_OPTIONS = foreign 

lib_LTLIBRARIES = libpdbcalBase.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

libpdbcalBase_la_LDFLAGS = \
  `root-config --evelibs`

libpdbcalBase_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool \
  -lXMLIO

ROOT_DICTS = \
  PdbBankID_Dict.cc \
  PdbCalBank_Dict.cc \
  PdbCalChan_Dict.cc \
  PdbParameter_Dict.cc \
  PdbParameterError_Dict.cc \
  PdbParameterMap_Dict.cc \
  PdbParameterMapContainer_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  PdbBankID_Dict_rdict.pcm \
  PdbCalBank_Dict_rdict.pcm \
  PdbCalChan_Dict_rdict.pcm \
  PdbParameter_Dict_rdict.pcm \
  PdbParameterError_Dict_rdict.pcm \
  PdbParameterMap_Dict_rdict.pcm \
  PdbParameterMapContainer_Dict_rdict.pcm

#please add new classes sorted according to the roman alphabet
libpdbcalBase_la_SOURCES = \
  $(ROOT_DICTS) \
  PdbApplication.cc \
  PdbBankID.cc \
  PdbBankManager.cc \
  PdbCalBank.cc \
  PdbParameter.cc \
  PdbParameterError.cc \
  PdbParameterMap.cc \
  PdbParameterMapContainer.cc \
  PHGenericFactoryT.C \
  RunToTime.cc


#please add new includes sorted according to the roman alphabet
pkginclude_HEADERS = \
  Pdb.h \
  PdbApplication.h \
  PdbApplicationFactory.h \
  PdbBankID.h \
  PdbBankList.h \
  PdbBankListIterator.h \
  PdbBankManager.h \
  PdbBankManagerFactory.h \
  PdbCalBank.h \
  PdbCalBankIterator.h \
  PdbCalChan.h \
  PdbClassMap.h \
  PdbParameter.h \
  PdbParameterError.h \
  PdbParameterMap.h \
  PdbParameterMapContainer.h \
  PHGenericFactoryT.h \
  RunToTime.h

BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = \
  libpdbcalBase.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


clean-local:
	rm -f *_LinkDef.h *_Dict* *.pcm $(BUILT_SOURCES)

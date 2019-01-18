AUTOMAKE_OPTIONS = foreign 

lib_LTLIBRARIES = libpdbcalBase.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

libpdbcalBase_la_LDFLAGS = \
  `root-config --evelibs`

libpdbcalBase_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool \
  -lXMLIO

# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOT_DICTS = \
  PdbBankID_Dict.cc \
  PdbCalBank_Dict.cc \
  PdbCalChan_Dict.cc \
  PdbParameter_Dict.cc \
  PdbParameterError_Dict.cc \
  PdbParameterMap_Dict.cc \
  PdbParameterMapContainer_Dict.cc

# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# MAKEROOT6 is set in the configure.ac, true for root6
if MAKEROOT6
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    PdbBankID_Dict_rdict.pcm \
    PdbCalBank_Dict_rdict.pcm \
    PdbCalChan_Dict_rdict.pcm \
    PdbParameter_Dict_rdict.pcm \
    PdbParameterError_Dict_rdict.pcm \
    PdbParameterMap_Dict_rdict.pcm \
    PdbParameterMapContainer_Dict_rdict.pcm
else
  ROOT5_DICTS = \
    PdbApplication_Dict.cc \
    PdbApplicationFactory_Dict.cc \
    PdbBankList_Dict.cc \
    PdbBankListIterator_Dict.cc \
    PdbBankManager_Dict.cc \
    PdbBankManagerFactory_Dict.cc \
    PdbCalBankIterator_Dict.cc \
    RunToTime_Dict.cc
endif

#please add new classes sorted according to the roman alphabet
libpdbcalBase_la_SOURCES = \
  $(ROOT_DICTS) \
  $(ROOT5_DICTS) \
  PdbApplication.cc \
  PdbBankID.cc \
  PdbBankList.cc \
  PdbBankManager.cc \
  PdbBankManagerFactory.cc \
  PdbCalBank.cc \
  PdbCalBankIterator.cc \
  PdbCalChan.cc \
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
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;


clean-local:
	rm -f *_LinkDef.h *_Dict* *.pcm $(BUILT_SOURCES)

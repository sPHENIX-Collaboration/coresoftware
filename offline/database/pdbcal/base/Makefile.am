AUTOMAKE_OPTIONS = foreign 

lib_LTLIBRARIES = libpdbcalBase.la 

BUILT_SOURCES = \
  testexternals.C

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

ROOTIFY  = perl $(srcdir)/PdbRootify.pl

libpdbcalBase_la_LDFLAGS = \
  `root-config --evelibs`

libpdbcalBase_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool \
  -lXMLIO

#please add new classes sorted according to the roman alphabet
libpdbcalBase_la_SOURCES = \
  PdbApplication.cc \
  PdbApplication_dict.C \
  PdbApplicationFactory.cc \
  PdbApplicationFactory_dict.C \
  PdbBankID.cc \
  PdbBankID_dict.C \
  PdbBankListIterator_dict.C \
  PdbBankList.cc \
  PdbBankList_dict.C \
  PdbBankManager.cc \
  PdbBankManager_dict.C \
  PdbBankManagerFactory.cc \
  PdbBankManagerFactory_dict.C \
  PdbCalBank.cc \
  PdbCalBank_dict.C \
  PdbCalBankIterator.cc \
  PdbCalBankIterator_dict.C \
  PdbCalChan.cc \
  PdbCalChan_dict.C \
  PdbParameter.cc \
  PdbParameter_dict.C \
  PdbParameterError.cc \
  PdbParameterError_dict.C \
  PdbParameterMap.cc \
  PdbParameterMap_dict.C \
  PdbParameterMapContainer.cc \
  PdbParameterMapContainer_dict.C \
  PHGenericFactoryT.C \
  RunToTime.cc \
  RunToTime_dict.C


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

# this generates the standard LinkDef.h files
%_LinkDef.h : 
	 $(ROOTIFY) $@

# keep the LinkDef files in the build area for debugging
.PRECIOUS : %_LinkDef.h

# this generates the standard dictionary
%_dict.C : %.h %_LinkDef.h
	 rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#
# Some Root dependencies have to be put in 'by hand'.
# Do not remove or modify this comment. PdbNewClass depends on it.
#

PdbBankManagerFactory_dict.C : PHGenericFactoryT.h PdbBankManagerFactory.h PdbBankManagerFactory_LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

RunToTimeFactory_dict.C: PHGenericFactoryT.h RunToTimeFactory.h RunToTimeFactory_LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *_LinkDef.h *_dict* *.pcm

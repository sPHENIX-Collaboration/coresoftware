AUTOMAKE_OPTIONS = foreign 

lib_LTLIBRARIES = libpdbcalBase.la 

BUILT_SOURCES = \
  testexternals.C

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

ROOTIFY  = perl $(srcdir)/PdbRootify.pl

libpdbcalBase_la_LDFLAGS = \
  `root-config --evelibs`

libpdbcalBase_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool

#please add new classes sorted according to the roman alphabet
libpdbcalBase_la_SOURCES = \
  PdbApplication.cc \
  PdbApplication_dict.C \
  PdbApplicationFactory.cc \
  PdbApplicationFactory_dict.C \
  PdbBankID.cc \
  PdbBankID_dict.C \
  PdbBankListIterator.cc \
  PdbBankListIterator_dict.C \
  PdbBankList.cc \
  PdbBankList_dict.C \
  PdbBankManager.cc \
  PdbBankManager_dict.C \
  PdbBankManagerFactory.cc \
  PdbBankManagerFactory_dict.C \
  PdbCalBank.cc \
  PdbCalBank_dict.C \
  PdbCalBase.cc \
  PdbCalBase_dict.C \
  PdbCalChan.cc \
  PdbCalChan_dict.C \
  PdbCalHeader.cc \
  PdbCalHeader_dict.C \
  PdbCalParameters.cc \
  PdbCalParameters_dict.C \
  PdbCoordinate.cc \
  PdbCoordinate_dict.C \
  PdbDouble.cc \
  PdbDouble_dict.C \
  PdbFloatVector.cc \
  PdbFloatVector_dict.C \
  PdbIndex.cc \
  PdbIndex_dict.C \
  PdbInt.cc \
  PdbInt_dict.C \
  PdbMapIntVecFloat.cc \
  PdbMapIntVecFloat_dict.C \
  PdbMapIntInt.cc \
  PdbMapIntInt_dict.C \
  PdbParameter.cc \
  PdbParameter_dict.C \
  PdbParameterError.cc \
  PdbParameterError_dict.C \
  RunToTime.cc \
  RunToTime_dict.C \
  RunToTimeFactory.cc \
  RunToTimeFactory_dict.C


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
  PdbCalBase.h \
  PdbCalChan.h \
  PdbCalHeader.h \
  PdbCalParameters.h \
  PdbClassMap.h \
  PdbCoordinate.h \
  PdbDouble.h \
  PdbFloatVector.h \
  PdbIndex.h \
  PdbInt.h \
  PdbMapIntVecFloat.h \
  PdbMapIntInt.h \
  PdbParameter.h \
  PdbParameterError.h \
  PHGenericFactoryT.h \
  RunToTime.h \
  RunToTimeFactory.h

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
	 rootcint -f $@ -c $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

#
# Some Root dependencies have to be put in 'by hand'.
# Do not remove or modify this comment. PdbNewClass depends on it.
#

PdbBankManagerFactory_dict.C : PHGenericFactoryT.h PdbBankManagerFactory.h PdbBankManagerFactory_LinkDef.h
	rootcint -f $@ -c $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

RunToTimeFactory_dict.C: PHGenericFactoryT.h RunToTimeFactory.h RunToTimeFactory_LinkDef.h
	rootcint -f $@ -c $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *_LinkDef.h
	rm -f *_dict*

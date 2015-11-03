AUTOMAKE_OPTIONS = foreign 

lib_LTLIBRARIES = libPgCal.la libPgCalInstance.la

BUILT_SOURCES = \
  testexternals.C

BANKOMAT = perl $(srcdir)/PgBankomat.pl
ROOTIFY =  perl $(srcdir)/PgRootify.pl

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(OPT_PHENIX)/include \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_PHENIX)/lib \
  `root-config --libs`

libPgCal_la_LIBADD = \
  -lpdbcalBase \
  -lRDBCodbc \
  -lRDBC \
  -lodbc++

libPgCalInstance_la_SOURCES = PgPostInstantiator.cc

libPgCalInstance_la_LIBADD = \
  libPgCal.la

libPgCal_la_SOURCES = \
  PgPostADCChanBank.cc \
  PgPostADCChanBank_dict.C \
  PgPostApplication.cc  \
  PgPostApplication_dict.C \
  PgPostBankManager.cc \
  PgPostBankManager_dict.C \
  PgPostBankWrapper.cc \
  PgPostBankWrapper_dict.C \
  PgPostBankWrapperManager.cc \
  PgPostBankWrapperManager_dict.C \
  PgPostBankBackupManager.cc \
  PgPostBankBackupManager_dict.C \
  PgPostBankBackupLog.cc \
  PgPostBankBackupLog_dict.C \
  PgPostBankBackupStorage.cc \
  PgPostBankBackupStorage_dict.C \
  PgPostCalBank.cc \
  PgPostCalBank_dict.C \
  PgPostCalBankIterator.cc \
  PgPostCalBankIterator2.cc \
  PgPostCalParametersBank.cc \
  PgPostCalParametersBank_dict.C \
  PgPostCoordinateBank.cc \
  PgPostCoordinateBank_dict.C \
  PgPostDoubleBank.cc \
  PgPostDoubleBank_dict.C \
  PgPostFloatVectorBank.cc \
  PgPostFloatVectorBank_dict.C \
  PgPostHelper.cc \
  PgPostIndexBank_dict.C \
  PgPostIndexBank.cc \
  PgPostIntBank.cc \
  PgPostIntBank_dict.C \
  PgPostMapIntVecFloatBank.cc \
  PgPostMapIntVecFloatBank_dict.C \
  PgPostMapIntIntBank.cc \
  PgPostMapIntIntBank_dict.C \
  PgPostParameterBank.cc \
  PgPostParameterBank_dict.C \
  PgPostParameterErrorBank.cc \
  PgPostParameterErrorBank_dict.C \
  RunToTimePg.cc \
  RunToTimePg_dict.C

pkginclude_HEADERS = \
  PgPostApplication.hh \
  PgPostBankManager.hh \
  PgPostCalBank.hh \
  RunToTimePg.hh

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = \
  libPgCal.la \
  libPgCalInstance.la

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

# this generates the classes stored in our postgres DB
PgPost%Bank.hh PgPost%Bank.cc:
	$(BANKOMAT) $@

# this generates the standard dictionaries
%_dict.C : %.hh %_LinkDef.h
	 rootcint -f $@ -c $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -rf *_dict*

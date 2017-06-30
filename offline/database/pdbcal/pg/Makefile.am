AUTOMAKE_OPTIONS = foreign 

lib_LTLIBRARIES = \
  libPgCal.la \
  libPgCalInstance.la

BANKOMAT = perl $(srcdir)/PgBankomat.pl
ROOTIFY =  perl $(srcdir)/PgRootify.pl

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I$(OPT_SPHENIX)/include \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/lib \
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
  PgPostParameterBank.cc \
  PgPostParameterBank_dict.C \
  PgPostParameterErrorBank.cc \
  PgPostParameterErrorBank_dict.C \
  PgPostParameterMapBank.cc \
  PgPostParameterMapBank_dict.C \
  PgPostParameterMapContainerBank.cc \
  PgPostParameterMapContainerBank_dict.C \
  RunToTimePg.cc \
  RunToTimePg_dict.C

pkginclude_HEADERS = \
  PgPostApplication.h \
  PgPostBankManager.h \
  PgPostCalBank.h \
  RunToTimePg.h

BUILT_SOURCES = testexternals.C

noinst_PROGRAMS = \
  testexternals_PgCal \
  testexternals_PgCalInstance

testexternals_PgCal_SOURCES = testexternals.C
testexternals_PgCal_LDADD = libPgCal.la 

testexternals_PgCalInstance_SOURCES = testexternals.C
testexternals_PgCalInstance_LDADD = libPgCalInstance.la

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
PgPost%Bank.h PgPost%Bank.cc:
	$(BANKOMAT) $@

# this generates the standard dictionaries
%_dict.C : %.h %_LinkDef.h
	 rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -rf *_dict*
	rm -rf *_LinkDef.h
	rm -rf PgPost*Bank.*

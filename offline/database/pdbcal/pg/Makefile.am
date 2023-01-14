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
  -isystem`root-config --incdir`

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

# I/O dictionaries have to exist for root5 and root6. For ROOT6 we need
# pcm files in addition. If someone can figure out how to make a list
# so this list of dictionaries is transformed into a list of pcm files
# following a simple naming convention, please change this accordingly and
# let me know
ROOT_DICTS = \
  PgPostBankBackupStorage_Dict.cc \
  PgPostBankWrapper_Dict.cc \
  PgPostCalBank_Dict.cc \
  PgPostParameterBank_Dict.cc \
  PgPostParameterErrorBank_Dict.cc \
  PgPostParameterMapBank_Dict.cc \
  PgPostParameterMapContainerBank_Dict.cc
# for root6 we need pcm and dictionaries but only for
# i/o classes. For root5 we need only dictionaries but
# those for i/o and classes available on the cmd line
# this is a tweak to install files in $(libdir), automake refuses
# to install other files in libdir, this construct gets around this
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
  PgPostBankBackupStorage_Dict_rdict.pcm \
  PgPostBankWrapper_Dict_rdict.pcm \
  PgPostCalBank_Dict_rdict.pcm \
  PgPostParameterBank_Dict_rdict.pcm \
  PgPostParameterErrorBank_Dict_rdict.pcm \
  PgPostParameterMapBank_Dict_rdict.pcm \
  PgPostParameterMapContainerBank_Dict_rdict.pcm


libPgCal_la_SOURCES = \
  $(ROOT_DICTS) \
  PgPostApplication.cc  \
  PgPostBankManager.cc \
  PgPostBankWrapper.cc \
  PgPostBankWrapperManager.cc \
  PgPostBankBackupManager.cc \
  PgPostBankBackupLog.cc \
  PgPostBankBackupStorage.cc \
  PgPostCalBankIterator.cc \
  PgPostParameterBank.cc \
  PgPostParameterErrorBank.cc \
  PgPostParameterMapBank.cc \
  PgPostParameterMapContainerBank.cc \
  RunToTimePg.cc

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

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

# this generates the standard LinkDef.h files
%LinkDef.h : 
	$(ROOTIFY) $@

# keep the LinkDef files in the build area for debugging
.PRECIOUS : %LinkDef.h
# make install-data removes the generated PgPost%Bank.h files which then
# crashes the dictionary compilation. We could probably list all of them
# under noinclude_HEADERS but declaring them PRECIOUS is easier
.PRECIOUS : PgPost%Bank.h

# this generates the classes stored in our postgres DB
PgPost%Bank.h PgPost%Bank.cc:
	$(BANKOMAT) $@

clean-local:
	rm -rf *_dict.* *_Dict.* *LinkDef.h PgPost*Bank.* *.pcm $(BUILT_SOURCES)

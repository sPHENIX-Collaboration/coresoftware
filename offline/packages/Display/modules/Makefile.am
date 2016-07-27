
AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce

noinst_LTLIBRARIES = libpheve_modules.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib `root-config --evelibs`

AM_CXXFLAGS = `root-config --cflags`

INCLUDES = -I$(top_srcdir)/display -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I$(OPT_PHENIX)/include 

libpheve_modules_la_LIBADD = \
	-lboost_thread \
	-lphool

##############################################
# please add new classes in alphabetical order

noinst_HEADERS = \
	mPHEveModuleBase.h \
	mG4EveDisplay.h \
	mSvtxEveDisplay.h \
	mCaloEveDisplay.h \
	mJetEveDisplay.h


##############################################
# please add new classes in alphabetical order

libpheve_modules_la_SOURCES = \
	mG4EveDisplay.cxx \
	mSvtxEveDisplay.cxx \
	mCaloEveDisplay.cxx \
	mJetEveDisplay.cxx

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libpheve_modules.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

%Dict.C : %.h %LinkDef.h
	rootcint -f $@ -c -D__linux__ $(DEFAULT_INCLUDES) $(INCLUDES) $^

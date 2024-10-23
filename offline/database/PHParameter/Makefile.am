AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libphparameter_io.la \
  libphparameter.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

libphparameter_io_la_SOURCES = \
  PHParameters.cc \
  PHParametersContainer.cc

libphparameter_la_SOURCES = \
  PHParameterContainerInterface.cc \
  PHParameterInterface.cc \
  PHParameterUtils.cc

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

libphparameter_la_LIBADD = \
  libphparameter_io.la \
  -lboost_system \
  -lffamodules \
  -lpdbcalBase \
  -lphool

libphparameter_io_la_LIBADD = \
  -lpdbcalBase \
  -lphool


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  PHParameterContainerInterface.h \
  PHParameterInterface.h \
  PHParameters.h \
  PHParametersContainer.h \
  PHParameterUtils.h

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_phparameter \
  testexternals_phparameter_io

testexternals_phparameter_SOURCES = testexternals.cc
testexternals_phparameter_LDADD = libphparameter.la

testexternals_phparameter_io_SOURCES = testexternals.cc
testexternals_phparameter_io_LDADD = libphparameter_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)

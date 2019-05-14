AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(top_srcdir)/Seamstress \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include/eigen3

AM_CXXFLAGS = \
  -DNDEBUG

lib_LTLIBRARIES = \
  libFitNewton.la

# that is a bad kludge, but AC_INIT's package name is lower case
# which then gives the wrong install dir (helixhough instead of HelixHough)
# if anyone has a better idea, let me know
packageincludedir = $(includedir)/@PACKAGE_NAME@

packageinclude_HEADERS = \
  NewtonMinimizerGradHessian.h \
  FunctionGradHessian.h \
  GaussianGradHessian.h \
  GaussianIntegralGradHessian.h \
  ChiSquareGradHessian.h \
  GaussianRegGradHessian.h \
  SquareGradient.h

libFitNewton_la_SOURCES = \
  NewtonMinimizerGradHessian.cpp \
  GaussianGradHessian.cpp \
  GaussianIntegralGradHessian.cpp \
  ChiSquareGradHessian.cpp \
  GaussianRegGradHessian.cpp \
  SquareGradient.cpp


libFitNewton_la_LIBADD = \
  $(top_builddir)/Seamstress/libSeamstress.la

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals 

BUILT_SOURCES = testexternals.cc

testexternals_SOURCES = testexternals.cc
testexternals_LDADD = libFitNewton.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)

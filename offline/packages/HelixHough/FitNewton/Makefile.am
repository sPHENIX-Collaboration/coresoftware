AM_CXXFLAGS = -I$(OFFLINE_MAIN)/include/eigen3 -I$(top_srcdir)/Seamstress

lib_LTLIBRARIES = \
libFitNewton.la

pkginclude_HEADERS = \
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
$(Seamstress_LIBS) $(Eigen3_LIBS)

libFitNewton_la_LDFLAGS = -version-info 1:0:0

libFitNewton_la_CPPFLAGS = $(Seamstress_CFLAGS) $(Eigen3_CFLAGS) -DNDEBUG

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = fitnewton.pc

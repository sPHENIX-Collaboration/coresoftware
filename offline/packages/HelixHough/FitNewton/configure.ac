AC_INIT(FitNewton, 1.0.0, alan.dion@gmail.com)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(NewtonMinimizerGradHessian.cpp)

AM_INIT_AUTOMAKE

AC_PROG_CXX
AC_PROG_LIBTOOL

CXXFLAGS="-Werror -Wall ${CXXFLAGS}"

PKG_CHECK_MODULES([Seamstress], [seamstress])
PKG_CHECK_MODULES([Eigen3], [eigen3])

dnl test for root 6
if test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([Makefile fitnewton.pc])

AC_OUTPUT


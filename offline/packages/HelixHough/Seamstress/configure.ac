AC_INIT(Seamstress, 1.0, alan.dion@gmail.com)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(Seamstress.cpp)

AM_INIT_AUTOMAKE

AC_PROG_CXX
AC_PROG_LIBTOOL

CXXFLAGS="-Wall -Werror ${CXXFLAGS}"

dnl test for root 6
if test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([Makefile seamstress.pc])

AC_OUTPUT


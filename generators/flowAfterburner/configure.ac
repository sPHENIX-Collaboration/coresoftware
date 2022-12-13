AC_INIT(flowAfterburner,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror -pedantic -Wextra -Wno-deprecated-declarations"
fi

dnl test for root 6
if test `root-config --version | awk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader -Wno-inconsistent-missing-override "
AC_SUBST(CINTDEFS)
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_INIT(pdbcalbase,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

dnl order is important, if AC_PROG_LIBTOOL is called before AC_PROG_CXX(CC g++)
dnl it will choose the g++ compiler over CC on Solaris
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wno-deprecated-declarations"
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

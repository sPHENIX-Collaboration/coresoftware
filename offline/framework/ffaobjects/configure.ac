AC_INIT(ffaobjects,[3.0])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

LT_INIT([disable-static])

AC_PROG_CXX(CC g++)

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror -pedantic -Wextra -Wshadow"
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

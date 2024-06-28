AC_INIT(inttrawhitqa,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror -Wno-error=unused-variable -Wno-error=sign-compare -Wno-error=unused-but-set-variable"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

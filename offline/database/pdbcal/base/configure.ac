AC_CONFIG_SRCDIR([configure.ac])
AC_INIT(pdbcalbase,[1.00])

AM_INIT_AUTOMAKE

dnl order is important, if AC_PROG_LIBTOOL is called before AC_PROG_CXX(CC g++)
dnl it will choose the g++ compiler over CC on Solaris
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror"
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

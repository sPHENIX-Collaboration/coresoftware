AC_INIT(ffaobjects,[3.0])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

LT_INIT([disable-static])

AC_PROG_CXX(CC g++)

dnl test for root 6
if test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi
AM_CONDITIONAL([MAKEROOT6],[test `root-config --version | gawk '{print $1>=6.?"1":"0"}'` = 1])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

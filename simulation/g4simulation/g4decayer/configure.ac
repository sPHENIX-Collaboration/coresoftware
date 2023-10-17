AC_INIT(g4decayer, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
dnl leaving this here in case we want to play with different compiler 
dnl specific flags
dnl case $CXX in
dnl  clang++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl  ;;
dnl  *g++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl  ;;
dnl esac

if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"
fi
case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wno-undefined-var-template"
 ;;
esac

dnl   AM_CONDITIONAL(GCC_GE_48, test `g++ -dumpversion | awk '{print $1>=4.8?"1":"0"}'` = 1)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT



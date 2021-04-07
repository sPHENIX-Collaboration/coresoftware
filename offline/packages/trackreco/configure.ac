AC_INIT(trackreco,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl leaving this here in case we want to play with different compiler 
dnl specific flags
dnl case $CXX in
dnl  clang++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl  ;;
dnl  g++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl  ;;
dnl esac

if test $ac_cv_prog_gxx = yes; then
     CXXFLAGS="$CXXFLAGS -Wall"
fi

AM_CONDITIONAL([MAKE_ACTS],false)
dnl gcc 8.3 creates warning in boost header, needs -Wno-class-memaccess
dnl need to check for *g++ since $CXX contains full path to g++
case $CXX in
 clang++)
   CXXFLAGS="$CXXFLAGS -Werror -Wno-undefined-var-template -Wno-unused-private-field -Wno-range-loop-construct -Wno-unused-local-typedef"
   AM_CONDITIONAL([MAKE_ACTS],test `clang++ -dumpversion | awk '{print $1>=10.0?"1":"0"}'` = 1)
 ;;
 *analyzer)
   CXXFLAGS="$CXXFLAGS"
   AM_CONDITIONAL([MAKE_ACTS],test `clang++ -dumpversion | awk '{print $1>=10.0?"1":"0"}'` = 1)
 ;;
 *g++)
  if test `g++ -dumpversion | awk '{print $1>=8.3?"1":"0"}'` = 1; then
     CXXFLAGS="$CXXFLAGS -Werror -Wno-class-memaccess -Wno-unused-local-typedefs -Wno-sign-compare -Wno-switch -Wno-unused-function -Wno-unused-value"
  fi
dnl ACTS needs at least gcc 6
  AM_CONDITIONAL([MAKE_ACTS],test `g++ -dumpversion | awk '{print $1>=6.0?"1":"0"}'` = 1)
 ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

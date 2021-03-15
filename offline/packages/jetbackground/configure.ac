AC_INIT(jetbackground,[1.00])
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
  if test `g++ -dumpversion | gawk '{print $1>=8.0?"1":"0"}'` = 1; then
   CXXFLAGS="$CXXFLAGS -Wall -Wno-deprecated-declarations -Werror"
  else
   CXXFLAGS="$CXXFLAGS -Wall -Werror"
  fi
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT


AC_INIT(jetbase, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"
dnl leaving this here in case we want to play with different compiler 
dnl specific flags
dnl case $CXX in
dnl  *analyzer)
dnl    CXXFLAGS="$CXXFLAGS -Wno-deprecated-declarations"
dnl  ;;
dnl  clang++)
dnl  ;;
dnl  *g++)
dnl    CXXFLAGS="$CXXFLAGS -Wno-deprecated-declarations"
dnl  ;;
dnl esac


CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

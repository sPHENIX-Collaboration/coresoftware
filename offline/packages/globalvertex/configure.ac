AC_INIT(globalvertex, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror"
dnl leaving this here in case we want to play with different compiler 
dnl specific flags
dnl case $CXX in
dnl  clang++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
dnl  ;;
dnl  *g++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
dnl  ;;
dnl esac


CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

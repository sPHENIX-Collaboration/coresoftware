AC_INIT(g4main, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"

dnl  AM_CONDITIONAL(GCC_GE_48, test `g++ -dumpversion | awk '{print $1>=4.8?"1":"0"}'` = 1)

dnl case $CXX in
dnl   clang++)
dnl    CXXFLAGS="$CXXFLAGS -Wno-c11-extensions -Wno-deprecated-copy"
dnl  ;;
dnl esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

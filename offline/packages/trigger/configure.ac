AC_INIT(calotrigger,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall  -Werror -Wextra"
 ;;
 *g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra"
 ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT


AC_INIT(g4main, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl leaving this here in case we want to play with different compiler 
dnl specific flags
case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-return-stack-address"
 ;;
 g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
 ;;
esac


dnl  AM_CONDITIONAL(GCC_GE_48, test `g++ -dumpversion | awk '{print $1>=4.8?"1":"0"}'` = 1)


CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

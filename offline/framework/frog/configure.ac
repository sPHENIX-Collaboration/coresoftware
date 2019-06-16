AC_INIT(frog,[2.0])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)

LT_INIT([disable-static])

dnl   no point in suppressing warnings people should
dnl   at least see them, so here we go for g++: -Wall
dnl   treat warnings as errors: -Werror
case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall -Wno-overloaded-virtual -Wno-deprecated-declarations -Wno-undefined-bool-conversion -Werror"
 ;;
 g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
 ;;
esac

dnl test for root 6
if test `root-config --version | awk '{print $1>=6.?"1":"0"}'` = 1; then
CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)
fi
AM_CONDITIONAL([MAKEROOT6],[test `root-config --version | awk '{print $1>=6.?"1":"0"}'` = 1])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

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
dnl odbc++ uses auto_ptr
  CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-deprecated-declarations -Wextra"
 ;;
 *g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -pedantic -Wno-deprecated-declarations -Wextra"
 ;;
esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

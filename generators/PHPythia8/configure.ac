AC_INIT(PHPythia8,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should
dnl   at least see them, so here we go for g++: -Wall
dnl   make warnings fatal errors: -Werror
dnl leaving this here in case we want to play with different compiler 
dnl specific flags
 case $CXX in
  clang++)
   CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-c11-extensions -Wextra"
  ;;
  *g++)
   CXXFLAGS="$CXXFLAGS -Wall -Werror -pedantic -Wextra "
  ;;
 esac

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

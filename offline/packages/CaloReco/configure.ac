AC_INIT(caloreco,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AC_CHECK_HEADERS([phool/phool.h])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wextra -Wshadow"

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wno-c11-extensions"
 ;;
esac

AC_ARG_ENABLE(online,
        [  --enable-online	build using for online [default=no]],
        [case "${enableval}" in
                yes) online=true ;;
                no)  online=false ;;
                *) AC_MSG_ERROR(bad value ${enableval} for --enable-online) ;;
                esac],
        online=false)
AM_CONDITIONAL(USE_ONLINE, test "x$online" = xtrue)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

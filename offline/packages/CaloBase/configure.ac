AC_INIT(calobase,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wall  -Werror -Wextra -Wshadow"
 ;;
 *g++)
  CXXFLAGS="$CXXFLAGS -Wall -Werror -pedantic -Wextra -Wshadow"
 ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

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

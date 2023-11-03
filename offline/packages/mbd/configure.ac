AC_INIT(mbd,[0.20])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE
AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl leaving this here in case we want to play with different compiler 
dnl specific flags
dnl case $CXX in
dnl  clang++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl  ;;
dnl  *g++)
dnl   CXXFLAGS="$CXXFLAGS -Wall -Werror"
dnl  ;;
dnl esac

if test $ac_cv_prog_gxx = yes; then
     CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror"
fi

AC_ARG_ENABLE(online,
        [  --enable-online	build using for online [default=no]],
        [case "${enableval}" in
                yes) online=true ;;
                no)  online=false ;;
                *) AC_MSG_ERROR(bad value ${enableval} for --enable-online) ;;
                esac],
        online=false)
AM_CONDITIONAL(USE_ONLINE, test "x$online" = xtrue)

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

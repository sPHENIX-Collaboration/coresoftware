AC_INIT(intt, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror"
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader -Wno-unknown-warning-option"
AC_SUBST(CINTDEFS)

AC_ARG_ENABLE(online,
	[ --enable-online	build using for online [default=no]],
	[case "${enableval}" in
		yes) online=true ;;
		no)  online=false ;;
		*) AC_MSG_ERROR(bad value ${enableval} for --enable-online) ;;
		esac],
	online=false)
AM_CONDITIONAL(USE_ONLINE, test "x$online" = xtrue)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

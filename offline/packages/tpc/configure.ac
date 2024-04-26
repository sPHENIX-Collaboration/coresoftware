AC_INIT(tpc, [1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should 
dnl   at least see them, so here we go for g++: -Wall
if test $ac_cv_prog_gxx = yes; then
   CXXFLAGS="$CXXFLAGS -Wall -Wextra -Werror"
fi

dnl gcc 8.3 creates warning in boost header, needs -Wno-class-memaccess
dnl need to check for *g++ since $CXX contains full path to g++
case $CXX in
 clang++)
   CXXFLAGS="$CXXFLAGS -Werror -Wno-undefined-var-template -Wno-unused-private-field -Wno-range-loop-construct -Wno-unused-local-typedef -Wno-deprecated-copy -Wno-sign-compare -Wno-unused-parameter -Wno-unused-value -Wno-bitwise-instead\
-of-logical -Wno-unknown-warning-option -Wno-deprecated-declarations"
 ;;
 *analyzer)
   CXXFLAGS="$CXXFLAGS"
 ;;
 *g++)
  if test `g++ -dumpversion | awk -F\. '{print $1>=8?"1":"0"}'` = 1; then
     CXXFLAGS="$CXXFLAGS -Werror -Wno-class-memaccess -Wno-unused-local-typedefs -Wno-sign-compare -Wno-switch -Wno-unused-function -Wno-unused-value -Wno-maybe-uninitialized"
  fi
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

head	1.1;
access;
symbols;
locks; strict;
comment	@# @;


1.1
date	2018.07.03.19.38.36;	author vassalli;	state Exp;
branches;
next	;


desc
@@


1.1
log
@Initial revision
@
text
@AC_INIT(configure.in)

AM_INIT_AUTOMAKE(ClusterIso, 1.00)
AC_PROG_CXX(CC g++)
AC_ENABLE_STATIC(no)
LT_INIT

if test $ac_cv_prog_gxx = yes; then
  CXXFLAGS="$CXXFLAGS -Wall -Werror"
fi
CXXFLAGS="$CXXFLAGS"
ROOTLIBS=`root-config --libs`
AC_SUBST(ROOTLIBS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
@

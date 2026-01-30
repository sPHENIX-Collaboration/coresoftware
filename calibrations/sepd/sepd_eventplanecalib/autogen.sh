#!/bin/sh
srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

(cd "$srcdir" || exit 1; aclocal -I "${OFFLINE_MAIN}/share" &&
libtoolize --force && automake -a --add-missing && autoconf)

$srcdir/configure  "$@"

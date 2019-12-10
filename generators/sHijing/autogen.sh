#!/bin/sh
srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

(cd $srcdir; autoreconf --force --install)

$srcdir/configure  "$@"

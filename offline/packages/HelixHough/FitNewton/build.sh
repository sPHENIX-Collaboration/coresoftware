#!/bin/tcsh

echo "   autoreconf..."
autoreconf --force --install

echo "   configure..."
./configure --prefix=$HOUGH_INSTALL_PATH

echo "   make..."
make

echo "   make install..."
make install

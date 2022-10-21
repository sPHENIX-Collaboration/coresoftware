#!/bin/tcsh

g++ -m32 -O3 -pipe `pkg-config --cflags eigen2`  `pkg-config --cflags helix_hough`  `root-config --cflags` -c  test_with_vertex.cpp

g++ -m32 -o test_with_vertex test_with_vertex.o `pkg-config --libs helix_hough` `root-config --libs` /usr/lib/gcc/i386-redhat-linux6E/4.3.2/libstdc++.a


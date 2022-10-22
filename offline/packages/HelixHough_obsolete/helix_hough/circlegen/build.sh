#!/bin/tcsh

g++ -m32 -O3 -pipe  `root-config --cflags` -c circlegen.cpp

g++ -m32 -o circlegen circlegen.o  `root-config --libs` /usr/lib/gcc/i386-redhat-linux6E/4.3.2/libstdc++.a


all: genRho genE genD

genRho: genRho.cc QPileUp.cxx QPileUpToy.cxx
	`root-config --cxx` -o genRho `root-config --cflags --libs` $^


genE: genE.cc FieldMaps.cxx FieldMapsLaplace.cxx LaplaceSolution.cxx airy.o d1mach.o src.o
	`root-config --cxx` -o genE `root-config --cflags --libs` -lgfortran $^
#	`root-config --cxx` -o genE `root-config --cflags --libs` -lgfortran -L /usr/local/gfortran/lib/ $^


genD: genD.cc Langevin.cxx
	`root-config --cxx` -o genD `root-config --cflags --libs` $^

airy.o: airy.f d1mach.f src.f
	`root-config --f77` -c `root-config --cflags` $^

simRho: simRho.cc sHelix.cxx sChargeMap.cxx
	`root-config --cxx` -o simRho `root-config --cflags --libs` $^

clean:
	rm *.o genRho genE genD simRho

// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "ValidationControl.h"
#include <iostream>
#include <stdio.h>
int main(int /*argc*/, char** /*argv*/)
{
    FILE* Finput=fopen("testPythia3.input","w");
    fprintf(Finput,"\
#\n\
# Process: pp -> @ 14TeV\n\
#\n\
Beams:frameType = 3\n\
Beams:idA = 2212\n\
Beams:idB = 2212\n\
SoftQCD:singleDiffractive = on\n\
Beams:allowMomentumSpread = on\n\
\n");
    fclose(Finput);

    FILE* Fconfig=fopen("testPythia3.config","w");
    fprintf(Fconfig,"\
INPUT  pythia8 testPythia3.input\n\
TOOL   photos                   \n\
TOOL   output                   \n\
EVENTS 10\n\
\n");
    fclose(Fconfig);

    ValidationControl control;
    control.read_file("testPythia3.config");
    control.initialize();
    int counter=0;
    while( control.new_event() )
    {
        GenEvent HepMCEvt(Units::GEV,Units::MM);
        control.process(HepMCEvt);
        counter++;
    }
    control.finalize();
    return 1*(counter-10);
}

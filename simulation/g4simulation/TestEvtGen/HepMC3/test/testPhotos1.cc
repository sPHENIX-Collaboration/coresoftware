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
    FILE* Finput=fopen("testPhotos1.input","w");
    fprintf(Finput,"\
#\n\
# Process: ee -> Z -> tau+ tau- @ 92GeV\n\
#\n\
\n\
WeakSingleBoson:ffbar2gmZ = on\n\
\n\
Beams:idA =  11\n\
Beams:idB = -11\n\
Beams:eCM =  92.\n\
\n\
# Simplify event as much as possible\n\
HadronLevel:all = off\n\
HadronLevel:Hadronize = off\n\
SpaceShower:QEDshowerByL = off\n\
SpaceShower:QEDshowerByQ = off\n\
PartonLevel:ISR = off\n\
PartonLevel:FSR = off\n\
\n\
# Set Z properties\n\
23:onMode = off\n\
23:onIfAny = 15\n\
\n\
# Leave tau undecayed (tau decays are very slow in Pythia 8.185)\n\
15:mayDecay  = off\n\
-15:mayDecay = off\n\
\n");
    fclose(Finput);


    FILE* Fconfig=fopen("testPhotos1.config","w");
    fprintf(Fconfig,"\
INPUT  pythia8      testPhotos1.input\n\
TOOL   photos                   \n\
EVENTS 1000\n\
\n");
    fclose(Fconfig);

    ValidationControl control;
    control.read_file("testPhotos1.config");
    control.set_event_limit(1000);
    control.initialize();
    int counter=0;
    while( control.new_event() )
    {
        GenEvent HepMCEvt(Units::GEV,Units::MM);
        control.process(HepMCEvt);
        counter++;
    }
    control.finalize();
    return 1*(counter-1000);
}

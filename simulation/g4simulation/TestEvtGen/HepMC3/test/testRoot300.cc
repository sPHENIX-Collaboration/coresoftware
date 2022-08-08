// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderRoot.h"
using namespace HepMC3;
int main()
{
    ReaderRoot inputA("inputRoot300.root");
    if(inputA.failed()) return 1;
    int particles=0;
    while( !inputA.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        particles+=evt.particles().size();
        evt.clear();
    }
    inputA.close();
    return (particles!=1200);
}

// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3TestUtils.h"
#include <fstream>
#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::filebuf
using namespace HepMC3;
int main()
{

    std::filebuf isrA;
    isrA.open("inputI05.hepmc",std::ios::in );
    std::istream SisrA(&isrA);
    ReaderAsciiHepMC2 inputA(SisrA);
    if(inputA.failed()) return 1;
    std::filebuf osrA;
    osrA.open("frominputI05.hepmc",std::ios::out);
    std::ostream SosrA(&osrA);
    WriterAscii       outputA(SosrA);
    if(outputA.failed()) return 2;
    while( !inputA.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        outputA.write_event(evt);
        evt.clear();
    }
    inputA.close();
    outputA.close();

    std::filebuf isrB;
    isrB.open("frominputI05.hepmc",ios_base::in );
    std::istream SisrB(&isrB);
    ReaderAscii inputB(SisrB);
    if(inputB.failed()) return 3;
    std::filebuf osrB;
    osrB.open ("fromfrominputI05.hepmc",ios_base::out );
    std::ostream SosrB(&osrB);
    WriterAsciiHepMC2       outputB(SosrB);
    if(outputB.failed()) return 4;
    while( !inputB.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputB.read_event(evt);
        if( inputB.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        outputB.write_event(evt);
        evt.clear();
    }
    inputB.close();
    outputB.close();
    return COMPARE_ASCII_FILES("fromfrominputI05.hepmc","inputI05.hepmc");
}

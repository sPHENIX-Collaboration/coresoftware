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
using namespace HepMC3;
int main()
{
    auto tA = chrono::steady_clock::now();
    ReaderAsciiHepMC2 inputA("inputSkip1.hepmc");
    if(inputA.failed()) return 1;
    WriterAscii       outputA("frominputSkip1A.hepmc");
    if(outputA.failed()) return 2;

    int i=0;
    while( !inputA.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        if (i%10==0) outputA.write_event(evt);
        i++;
        evt.clear();
    }
    inputA.close();
    outputA.close();
    printf("Time taken A: %.2fms\n", chrono::duration <double, milli> (chrono::steady_clock::now()-tA).count());
    auto tB = chrono::steady_clock::now();

    ReaderAsciiHepMC2 inputB("inputSkip1.hepmc");
    if(inputB.failed()) return 1;

    WriterAscii       outputB("frominputSkip1B.hepmc");
    if(outputB.failed()) return 2;

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
        inputB.skip(9);
    }
    inputB.close();
    outputB.close();
    printf("Time taken B: %.2fms\n", chrono::duration <double, milli> (chrono::steady_clock::now()-tB).count());

    return COMPARE_ASCII_FILES("frominputSkip1A.hepmc","frominputSkip1B.hepmc");
}

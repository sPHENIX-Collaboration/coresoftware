// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMCCompatibility.h"
#include "HepMC3TestUtils.h"
int main()
{
    HepMC3::ReaderAsciiHepMC2 inputA("inputConvert1.hepmc");
    if(inputA.failed()) return 1;
    std::ofstream outputA("frominputConvert1.hepmc");
    outputA.precision(12);
    HepMC::write_HepMC_IO_block_begin(outputA);
    while( !inputA.failed() )
    {
        HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        HepMC::GenEvent* evt3=ConvertHepMCGenEvent_3to2(evt);
        if (!evt3) return 3;
        evt3->write(outputA);
        evt.clear();
        delete evt3;
    }
    inputA.close();
    HepMC::write_HepMC_IO_block_end( outputA );
    outputA.close();
    return COMPARE_ASCII_FILES("frominputConvert1.hepmc","inputConvert1.hepmc");
}

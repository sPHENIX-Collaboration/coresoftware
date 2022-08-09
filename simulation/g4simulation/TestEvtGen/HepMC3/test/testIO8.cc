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
    ReaderAsciiHepMC2 inputA("inputIO8.hepmc");
    if(inputA.failed()) return 1;
    WriterAscii       outputA("frominputIO8.hepmc");
    if(outputA.failed()) return 2;
    auto optionsA =  outputA.get_options();
    optionsA["float_printf_specifier"] = "g";
    outputA.set_options(optionsA);
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


    ReaderAscii inputB("frominputIO8.hepmc");
    if(inputB.failed()) return 3;
    WriterAsciiHepMC2       outputB("fromfrominputIO8.hepmc");
    if(outputB.failed()) return 4;
    auto optionsB =  outputB.get_options();
    optionsB["float_printf_specifier"] = "g";
    outputB.set_options(optionsB);
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
    return COMPARE_ASCII_FILES("fromfrominputIO8.hepmc","inputIO8.hepmc");
}

// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3TestUtils.h"
using namespace HepMC3;
int main()
{
    std::shared_ptr<Reader> inputA = deduce_reader("inputReaderFactory1.hepmc");
    if(inputA->failed()) return 1;
    WriterAscii       outputA("frominputReaderFactory1.hepmc");
    if(outputA.failed()) return 2;
    while( !inputA->failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA->read_event(evt);
        if( inputA->failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        outputA.write_event(evt);
        evt.clear();
    }
    inputA->close();
    outputA.close();


    std::shared_ptr<Reader> inputB = deduce_reader("frominputReaderFactory1.hepmc");
    if(inputB->failed()) return 3;
    WriterAsciiHepMC2       outputB("fromfrominputReaderFactory1.hepmc");
    if(outputB.failed()) return 4;
    while( !inputB->failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputB->read_event(evt);
        if( inputB->failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        outputB.write_event(evt);
        evt.clear();
    }
    inputB->close();
    outputB.close();
    return COMPARE_ASCII_FILES("fromfrominputReaderFactory1.hepmc","inputReaderFactory1.hepmc");
}

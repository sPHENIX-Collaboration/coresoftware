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
#include "HepMC3/ReaderMT.h"
#include <algorithm>
using namespace HepMC3;
int main()
{
    ReaderMT<ReaderAsciiHepMC2,3> inputA("inputIO10.hepmc");
    std::vector<GenEvent> inputA_events;
    if(inputA.failed()) return 1;
    while( !inputA.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        inputA_events.push_back(evt);
    }
    inputA.close();
    WriterAscii       outputA("frominputIO10.hepmc");
    if(outputA.failed()) return 2;
    std::sort(inputA_events.begin(), inputA_events.end(), [](const GenEvent& lhs, const GenEvent& rhs) {
        return lhs.event_number() < rhs.event_number();
    });
    for (auto& e: inputA_events) outputA.write_event(e);
    outputA.close();
    inputA_events.clear();

    ReaderMT<ReaderAscii,2> inputB("frominputIO10.hepmc");
    std::vector<GenEvent> inputB_events;
    if(inputB.failed()) return 3;
    while( !inputB.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputB.read_event(evt);
        if( inputB.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        inputB_events.push_back(evt);
    }
    inputB.close();
    WriterAsciiHepMC2       outputB("fromfrominputIO10.hepmc");
    if(outputB.failed()) return 4;
    std::sort(inputB_events.begin(), inputB_events.end(), [](const GenEvent& lhs, const GenEvent& rhs) {
        return lhs.event_number() < rhs.event_number();
    });
    for (auto& e: inputB_events) outputB.write_event(e);
    outputB.close();
    inputB_events.clear();

    return COMPARE_ASCII_FILES("fromfrominputIO10.hepmc","inputIO10.hepmc");
}

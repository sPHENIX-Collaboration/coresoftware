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
#include "HepMC3/ReaderHEPEVT.h"
#include "HepMC3/WriterHEPEVT.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3TestUtils.h"
using namespace HepMC3;
int main()
{
    std::shared_ptr<Reader> input = deduce_reader("inputReaderFactory2.hepmc");
    if(input->failed()) return 1;
    WriterAscii             outputA("frominputReaderFactory2.hepmc3");
    WriterAsciiHepMC2       outputB("frominputReaderFactory2.hepmc2");
    WriterHEPEVT            outputC("frominputReaderFactory2.hepevt");
    WriterRootTree          outputD("frominputReaderFactory2.root");
    if(outputA.failed()) return 2;
    if(outputB.failed()) return 3;
    if(outputC.failed()) return 4;
    if(outputD.failed()) return 5;
    while( !input->failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        input->read_event(evt);
        if( input->failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        outputA.write_event(evt);
        outputB.write_event(evt);
        outputC.write_event(evt);
        outputD.write_event(evt);
        evt.clear();
    }
    input->close();
    outputA.close();
    outputB.close();
    outputC.close();
    outputD.close();

    std::vector<std::shared_ptr<Reader> > inputv;
    inputv.push_back(deduce_reader("frominputReaderFactory2.hepmc3"));
    inputv.push_back(deduce_reader("frominputReaderFactory2.hepmc2"));
    inputv.push_back(deduce_reader("frominputReaderFactory2.hepevt"));
    inputv.push_back(deduce_reader("frominputReaderFactory2.root"));

    std::vector<WriterAsciiHepMC2*> outputv;

    outputv.push_back(new WriterAsciiHepMC2("AA.hepmc2"));
    outputv.push_back(new WriterAsciiHepMC2("BB.hepmc2"));
    outputv.push_back(new WriterAsciiHepMC2("CC.hepmc2"));
    outputv.push_back(new WriterAsciiHepMC2("DD.hepmc2"));

    for (size_t i=0; i<inputv.size(); i++)
        while( !inputv.at(i)->failed() )
        {
            GenEvent evt(Units::GEV,Units::MM);
            inputv.at(i)->read_event(evt);
            if( inputv.at(i)->failed() )  {
                printf("End of file reached. Exit.\n");
                break;
            }
            outputv.at(i)->write_event(evt);
            evt.clear();
        }
    for (size_t i=0; i<outputv.size(); i++) outputv.at(i)->close();

    return COMPARE_ASCII_FILES("AA.hepmc2","BB.hepmc2")+COMPARE_ASCII_FILES("BB.hepmc2","DD.hepmc2");
}

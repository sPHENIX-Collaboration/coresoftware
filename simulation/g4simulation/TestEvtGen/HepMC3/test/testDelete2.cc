// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/Print.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3TestUtils.h"
using namespace HepMC3;
int main()
{
    ReaderAsciiHepMC2 inputA("inputDelete2.hepmc");
    if(inputA.failed()) return 1;
    std::vector<std::shared_ptr<GenEvent> > evts;
    while( !inputA.failed() )
    {
        std::shared_ptr<GenEvent>  evt= std::make_shared<GenEvent>();
        inputA.read_event(*evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        evts.push_back(evt);
    }
    inputA.close();
    int i=0;
    int j=0;
    while(i==j)
    {
        i=rand()% evts.size();
        j=rand()% evts.size();
    }

    std::vector<GenVertexPtr> selectedVtx;
    for (std::vector<GenParticlePtr>::const_iterator p=evts.at(i)->particles().begin(); p!=evts.at(i)->particles().end(); ++p)
        for (std::vector<GenVertexPtr>::iterator v=evts.at(j)->vertices().begin(); v!=evts.at(j)->vertices().end(); ++v)
        {
            selectedVtx.push_back(*v);
        }

    WriterAscii       outputA("frominputDelete2.hepmc");
    if(outputA.failed()) return 2;
    for (size_t i=0; i<evts.size(); i++) outputA.write_event(*evts[i]);
    evts.clear();
    outputA.close();
    for (std::vector<GenVertexPtr>::const_iterator v=selectedVtx.begin(); v!=selectedVtx.end(); ++v)
    {
        Print::line(*v);
    }
    if (selectedVtx.front()->parent_event())
        Print::listing(*(selectedVtx.front()->parent_event()));
    ReaderAscii inputB("frominputDelete2.hepmc");
    if(inputB.failed()) return 3;
    WriterAsciiHepMC2       outputB("fromfrominputDelete2.hepmc");
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
    return COMPARE_ASCII_FILES("fromfrominputDelete2.hepmc","inputDelete2.hepmc");
}

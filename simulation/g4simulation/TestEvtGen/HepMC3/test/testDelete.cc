// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
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
    ReaderAsciiHepMC2 inputA("inputDelete.hepmc");
    if(inputA.failed()) return 1;
    std::vector<GenEvent> evts;
    while( !inputA.failed() )
    {
        GenEvent evt=GenEvent(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        evts.push_back(evt);
    }
    inputA.close();
//No alien particles should be detached from vertices or removed from events
    int i=0;
    int j=0;
    while(i==j)
    {
        i=rand()% evts.size();
        j=rand()% evts.size();
    }
    evts[i].remove_particles(evts[j].particles());

    for (GenParticlePtr p: evts.at(i).particles())
        evts[j].remove_particle(p);

    for (GenParticlePtr p: evts.at(i).particles()) {
        for (GenVertexPtr v: evts.at(j).vertices()) {
            (v)->remove_particle_in(p);
            (v)->remove_particle_out(p);
        }
    }

    WriterAscii       outputA("frominputDelete.hepmc");
    if(outputA.failed()) return 2;
    for (size_t i=0; i<evts.size(); i++) outputA.write_event(evts[i]);
    evts.clear();
    outputA.close();


    ReaderAscii inputB("frominputDelete.hepmc");
    if(inputB.failed()) return 3;
    WriterAsciiHepMC2       outputB("fromfrominputDelete.hepmc");
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
    return COMPARE_ASCII_FILES("fromfrominputDelete.hepmc","inputDelete.hepmc");
}

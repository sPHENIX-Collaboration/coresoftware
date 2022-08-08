// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/Print.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3TestUtils.h"
using namespace HepMC3;
int main()
{
    Setup::set_debug_level(60);
    ReaderAsciiHepMC2 inputA("inputSingleVertexHepMC2.hepmc");
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

    if (evts[0]->particles().size()==120&&evts[0]->vertices().size()==1) return EXIT_SUCCESS;
    return EXIT_FAILURE;
}

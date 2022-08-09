// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file rootIOTree_example_write.cc
 *  @brief Basic example of use of root I/O: writing events to file
 *
 *  @author Witold Pokorski/Andrii Verbytskyi
 *  @date   29/10/15
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/Print.h"

#include <iostream>

using namespace HepMC3;

/** Main */
int main(int argc, char **argv) {

    if( argc<3 ) {
        std::cout << "Usage: " << argv[0] << " <input_hepmc3_file> <output_root_file>" << std::endl;
        exit(-1);
    }

    ReaderAscii text_input (argv[1]);
    WriterRootTree  root_output(argv[2]);

    int events_parsed = 0;

    while( !text_input.failed() ) {

        GenEvent evt(Units::GEV,Units::MM);

        text_input.read_event(evt);

        if( text_input.failed() ) break;

        if( events_parsed == 0 ) {
            std::cout << "First event: " << std::endl;
            Print::listing(evt);
        }

        root_output.write_event(evt);
        ++events_parsed;

        if( events_parsed%1000 == 0 ) {
            std::cout << "Event: " << events_parsed << std::endl;
        }
    }

    text_input.close();
    root_output.close();

    std::cout << "Events parsed and written: " << events_parsed << std::endl;

    return 0;
}

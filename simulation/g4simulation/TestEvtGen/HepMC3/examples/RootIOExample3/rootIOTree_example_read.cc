// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file rootIOTree_example_read.cc
 *  @brief Basic example of use of root I/O with tree: reading events from file
 *
 *  @author Witold Pokorski/Andrii Verbytskyi
 *  @date   29/10/15
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/Print.h"

#include <iostream>

using namespace HepMC3;

/** Main */
int main(int argc, char **argv) {

    if( argc<3 ) {
        std::cout << "Usage: " << argv[0] << " <input_root_file> <output_hepmc3_file>" << std::endl;
        exit(-1);
    }

    ReaderRootTree  root_input (argv[1]);
    WriterAscii text_output(argv[2]);

    int events_parsed = 0;

    while( !root_input.failed() ) {

        GenEvent evt;

        root_input.read_event(evt);

        if( root_input.failed() ) break;

        if( events_parsed == 0 ) {
            std::cout << "First event: " << std::endl;
            Print::listing(evt);
        }

        text_output.write_event(evt);
        ++events_parsed;

        if( events_parsed%1000 == 0 ) {
            std::cout << "Event: " << events_parsed << std::endl;
        }
    }

    root_input.close();
    text_output.close();

    std::cout << "Events parsed and written: " << events_parsed << std::endl;

    return 0;
}

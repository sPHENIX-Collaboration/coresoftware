// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @example HepMC2_reader_example.cc
 *  @brief Example of use of HepMC2 adapter
 *
 *  Converts selected HepMC2 file to HepMC3 file
 *
 */
#include "HepMC3/GenEvent.h"

#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include <iostream>
#include <cstdlib> // atoi
using namespace HepMC3;

/** Main program */
int main(int argc, char **argv) {

    if( argc < 3 ) {
        std::cout<<"Usage: " << argv[0] <<" <input_hepmc2_file> <output_hepmc3_file> [<optional_events_limit>]" << std::endl;
        exit(-1);
    }

    // Open input and output files
    ReaderAsciiHepMC2 adapter(argv[1]);
    WriterAscii    output_file(argv[2]);
    int events_parsed = 0;
    int events_limit  = 0;

    if( argc >= 4 ) events_limit = atoi(argv[3]);

    while( !adapter.failed() ) {
        GenEvent evt(Units::GEV,Units::MM);

        // Read event from input file
        adapter.read_event(evt);

        // If reading failed - exit loop
        if( adapter.failed() ) break;

        // Save event to output file
        output_file.write_event(evt);

        if(events_parsed==0) {
            std::cout << " First event: " << std::endl;
            Print::listing(evt);
        }

        ++events_parsed;
        if( events_parsed%100 == 0 ) std::cout<<"Events parsed: "<<events_parsed<<std::endl;
        if( events_limit && events_parsed >= events_limit ) break;
    }

    adapter.close();
    output_file.close();

    return 0;
}

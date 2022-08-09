// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @example HepMC3_fileIO_example.cc
 *  @brief Test of file I/O
 *
 *  Parses HepMC3 file and saves it as a new HepMC3 file.
 *  The resulting file should be an exact copy of the input file
 *
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include <iostream>
using namespace HepMC3;

/** Main program */
int main(int argc, char **argv) {

    if( argc<3 ) {
        std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> <output_file>" << std::endl;
        exit(-1);
    }

    ReaderAscii input_file (argv[1]);
    WriterAscii output_file(argv[2]);

    int events_parsed = 0;

    while(!input_file.failed()) {
        GenEvent evt(Units::GEV,Units::MM);

        // Read event from input file
        input_file.read_event(evt);

        // If reading failed - exit loop
        if( input_file.failed() ) break;

        // Save event to output file
        output_file.write_event(evt);

        if(events_parsed==0) {
            std::cout << " First event: " << std::endl;
            Print::listing(evt);
            Print::content(evt);

            std::cout << " Testing attribute reading for the first event: " << std::endl;

            std::shared_ptr<GenCrossSection> cs = evt.attribute<GenCrossSection>("GenCrossSection");
            std::shared_ptr<GenHeavyIon>     hi = evt.attribute<GenHeavyIon>("GenHeavyIon");
            std::shared_ptr<GenPdfInfo>      pi = evt.attribute<GenPdfInfo>("GenPdfInfo");

            if(cs) {
                std::cout << " Has GenCrossSection:   ";
                Print::line(cs);
            }
            else std::cout << " No GenCrossSection " << std::endl;

            if(pi) {
                std::cout << " Has GenPdfInfo:        ";
                Print::line(pi);
            }
            else std::cout << " No GenPdfInfo " << std::endl;

            if(hi) {
                std::cout << " Has GenHeavyIon:       ";
                Print::line(hi);
            }
            else std::cout << " No GenHeavyIon " << std::endl;
        }

        ++events_parsed;
        if( events_parsed%100 == 0 ) std::cout<<"Events parsed: "<<events_parsed<<std::endl;
    }

    input_file.close();
    output_file.close();

    return 0;
}

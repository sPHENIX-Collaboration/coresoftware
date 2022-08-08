// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @example class_example_write.cc
 *  @brief Basic example of use of root I/O: writing events to file
 *
 *  @author Witold Pokorski
 *  @date   16/10/14
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"

#include "TFile.h"

#include <iostream>
#include <sstream>

using namespace HepMC3;

#include "MyClass.h"
#include "MyRunClass.h"

/** Main */
int main(int argc, char **argv) {

    if( argc<3 ) {
        std::cout << "Usage: " << argv[0] << " <input_hepmc3_file> <output_root_file>" <<std:: endl;
        exit(-1);
    }

    ReaderAscii text_input(argv[1]);

    TFile* fFile = new TFile(argv[2],"RECREATE");

    int events_parsed = 0;

    bool is_gen_run_info_written = false;

    while( !text_input.failed() ) {

        GenEvent evt(Units::GEV,Units::MM);

        text_input.read_event(evt);

        if( text_input.failed() ) break;

        if( events_parsed == 0 ) {
            std::cout << "First event: " << std::endl;
            Print::listing(evt);
        }

        if(!is_gen_run_info_written) {
            if(evt.run_info()) {
                GenRunInfo run_info(*evt.run_info());

                MyRunClass *my_run = new MyRunClass();

                my_run->SetRunInfo(&run_info);

                fFile->WriteObject(my_run,"MyRunClass");

                is_gen_run_info_written = true;
            }
        }

        MyClass* myclass = new MyClass();

        myclass->SetEvent(&evt);
        //
        std::ostringstream os;
        os << events_parsed;
        std::string stevt = "Event_" + os.str();
        const char* chevt = stevt.c_str();

        std::cout << "writing " << stevt << std::endl;

        fFile->WriteObject(myclass, chevt);

        ++events_parsed;

        if( events_parsed%100 == 0 ) {
            std::cout << "Event: " << events_parsed << std::endl;
        }
    }

    text_input.close();
    fFile->Close();
    std::cout << "Events parsed and written: " << events_parsed << std::endl;

    return 0;
}

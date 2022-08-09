// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @example class_example_read.cc
 *  @brief Basic example of use of root I/O: reading events from file
 *
 *  @author Witold Pokorski
 *  @date   16/10/14
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "MyClass.h"
#include "MyRunClass.h"

#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"

#include <iostream>

using namespace HepMC3;


/** Main */
int main(int argc, char **argv) {

    if( argc<3 ) {
        std::cout << "Usage: " << argv[0] << " <input_root_file> <output_hepmc3_file>" << std::endl;
        exit(-1);
    }


    TFile fo(argv[1]);
    WriterAscii text_output(argv[2]);

    MyClass* myevent;
    int events_parsed = 0;

    // Get GenRunInfo, if available
    MyRunClass *my_run = (MyRunClass*)fo.Get("MyRunClass");
    std::shared_ptr<GenRunInfo> run_info;

    if( my_run ) run_info.reset(my_run->GetRunInfo());


    fo.GetListOfKeys()->Print();

    TIter next(fo.GetListOfKeys());
    TKey *key;

    while ((key=(TKey*)next()))
    {
        const char *cl = key->GetClassName();

        if( strncmp(cl,"MyClass",7) != 0 ) continue;

        fo.GetObject(key->GetName(), myevent);

        std::cout << "Event: " << key->GetName() << std::endl;

        if( events_parsed == 0 ) {
            std::cout << "First event: " << std::endl;
            Print::listing(*(myevent->GetEvent()));
        }

        if( run_info ) {
            std::cout << "Setting run info" << std::endl;
            myevent->GetEvent()->set_run_info(run_info);
            run_info.reset();
        }

        text_output.write_event(*(myevent->GetEvent()));
        ++events_parsed;

        if( events_parsed%100 == 0 ) {
            std::cout << "Event: " << events_parsed << std::endl;
        }

        delete myevent->GetEvent();
        delete myevent;
    }

    text_output.close();

    std::cout << "Events parsed and written: " << events_parsed << std::endl;

    return 0;
}

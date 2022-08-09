// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3TestUtils.h"
using namespace HepMC3;
int main()
{
    ReaderAsciiHepMC2 inputA("inputIO7.hepmc");
    if(inputA.failed()) return 1;
    WriterAscii       outputA("frominputIO7.hepmc");
    if(outputA.failed()) return 2;
    size_t n = 0;
    while( !inputA.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        evt.set_run_info(std::make_shared<GenRunInfo>());
        std::vector<std::string> w_names;
        std::vector<double> w_values;
        for (size_t i=0; i<n+2; i++) {
            w_names.push_back(std::string("testname")+std::to_string(i));
            w_values.push_back(1.0+0.1*i);
        }
        evt.run_info()->add_attribute("testrunattribute",std::make_shared<IntAttribute>(10000+n));
        if (n%2==0)
            evt.run_info()->set_weight_names(w_names);//Try run info with names
        else
            evt.run_info()->set_weight_names(std::vector<std::string>());//Try run info w/o names
        evt.weights() = w_values;
        outputA.set_run_info(nullptr);//This instructs reader to take run info from the event
        outputA.write_event(evt);
        evt.clear();
        n++;
    }
    inputA.close();
    outputA.close();


    ReaderAscii inputB("frominputIO7.hepmc");
    if(inputB.failed()) return 3;
    WriterAscii outputB("fromfrominputIO7.hepmc");
    if(outputB.failed()) return 4;
    while( !inputB.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputB.read_event(evt);
        if( inputB.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        outputB.set_run_info(nullptr);//This instructs reader to take run info from the event
        outputB.write_event(evt);
        evt.clear();
    }
    inputB.close();
    outputB.close();
    return COMPARE_ASCII_FILES("fromfrominputIO7.hepmc","frominputIO7.hepmc");
}

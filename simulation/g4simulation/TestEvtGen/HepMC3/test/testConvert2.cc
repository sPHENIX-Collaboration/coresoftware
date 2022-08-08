// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMCCompatibility.h"
#include "HepMC3TestUtils.h"
int main()
{
    std::ifstream inputA( "inputConvert2.hepmc" );
    if( !inputA ) return 1;
    HepMC3::WriterAsciiHepMC2 outputA("frominputConvert2.hepmc");
    std::shared_ptr<HepMC3::GenRunInfo> run =std::make_shared<HepMC3::GenRunInfo>();
    while(inputA)
    {
        HepMC::GenEvent evt;
        evt.clear();
        evt.read( inputA );
        if( !evt.is_valid() )  break;
        HepMC3::GenEvent* evt3=ConvertHepMCGenEvent_2to3(evt,run);
        if (!evt3) return 4;
        outputA.write_event(*evt3);
        delete evt3;
    }
    inputA.close();
    outputA.close();
    return COMPARE_ASCII_FILES("frominputConvert2.hepmc","inputConvert2.hepmc");
}

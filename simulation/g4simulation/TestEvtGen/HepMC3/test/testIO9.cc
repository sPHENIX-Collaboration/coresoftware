// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderGZ.h"
#include "HepMC3/WriterGZ.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3TestUtils.h"
using namespace HepMC3;
int main()
{
    ReaderAsciiHepMC2 inputA("inputIO9.hepmc");
    if(inputA.failed()) return 1;
    const std::vector<Compression> supported = HepMC3::supported_compression_types;
    std::vector<std::shared_ptr<Writer> > writersGZ;
    for (auto w: supported) {
        switch (w) {
        case Compression::z: {
            writersGZ.push_back(std::shared_ptr<Writer>(new WriterGZ<WriterAsciiHepMC2,Compression::z>("frominputIO9.hepmc."+std::to_string(w))));
            break;
        }
        case Compression::lzma: {
            writersGZ.push_back(std::shared_ptr<Writer>(new WriterGZ<WriterAsciiHepMC2,Compression::lzma>("frominputIO9.hepmc."+std::to_string(w))));
            break;
        }
        case Compression::bz2: {
            writersGZ.push_back(std::shared_ptr<Writer>(new WriterGZ<WriterAsciiHepMC2,Compression::bz2>("frominputIO9.hepmc."+std::to_string(w))));
            break;
        }
        default: {
            return 9;
        }
        }
        if(writersGZ.back()->failed()) return 10+writersGZ.size();

    }
    while( !inputA.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        for (auto w: writersGZ) w->write_event(evt);
        evt.clear();
    }
    inputA.close();
    for (auto w: writersGZ) w->close();
    writersGZ.clear();

    int result = 0;
    for (auto w: supported) {
        ReaderGZ<ReaderAsciiHepMC2> inputB("frominputIO9.hepmc."+std::to_string(w));
        if(inputB.failed()) return 20;

        WriterAsciiHepMC2 outputB("fromfrominputIO9" + std::to_string(w) + ".hepmc");
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
        result += COMPARE_ASCII_FILES("fromfrominputIO9" + std::to_string(w) + ".hepmc","inputIO9.hepmc");
    }
    return result;
}

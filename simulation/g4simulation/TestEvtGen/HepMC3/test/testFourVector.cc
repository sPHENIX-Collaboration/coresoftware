// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#include <iostream>
#include <fstream>
#include <vector>
#include "HepMC3/FourVector.h"
#include "HepMC3/PrintStreams.h"
#undef NDEBUG
#include <assert.h>
using namespace HepMC3;
int main()
{
    std::vector<FourVector>  vectors_to_test{
        FourVector(0.0, 0.0, 0.0, 0.0),
        FourVector(1.0, 2.0, 0.0, 0.0),
        FourVector(0.0, 0.0, 0.0, 1.0),
        FourVector(0.0, 0.0, 0.0,-1.0),
        FourVector(0.0, 0.0, 1.0, 0.0),
        FourVector(0.0, 0.0,-1.0, 0.0),
        FourVector(1.0, 2.0, 3.0, 4.0),
        FourVector(1.0, 2.0, 3.0, -4.0),
        FourVector(1.0, 2.0, -3.0, 4.0),
        FourVector(1.0, 2.0, -3.0, -4.0),
        FourVector(1.0, 2.0, -3.0, 40.0),
        FourVector(1.0, 2.0, -3.0, -40.0)
    };
    std::vector<double>   correct_eta{
        0.0,
        0.0,
        0.0,
        0.0,
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        /*
                std::numeric_limits<double>::max(),
                -std::numeric_limits<double>::max(),
        */
        std::log((std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) + 3.0) / (std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) - 3.0))*0.5,
        std::log((std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) + 3.0) / (std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) - 3.0))*0.5,
        std::log((std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) - 3.0) / (std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) + 3.0))*0.5,
        std::log((std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) - 3.0) / (std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) + 3.0))*0.5,
        std::log((std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) - 3.0) / (std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) + 3.0))*0.5,
        std::log((std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) - 3.0) / (std::sqrt(1.0*1.0 + 2.0*2.0 + 3.0*3.0) + 3.0))*0.5
    };

    std::vector<double>   correct_rap{
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        std::log((4.0 + 3.0) / (4.0 - 3.0))*0.5,
        std::log((-4.0 + 3.0) / (-4.0 - 3.0))*0.5,
        std::log((4.0 - 3.0) / (4.0 + 3.0))*0.5,
        std::log((-4.0 - 3.0) / (-4.0 + 3.0))*0.5,
        std::log((40.0 - 3.0) / (40.0 + 3.0))*0.5,
        std::log((-40.0 - 3.0) / (-40.0 + 3.0))*0.5

    };
    cout.setf(ios_base::scientific);
    cout.precision(10);
    cout.width(15);
    cout.setf(std::ios_base::showpos);
    for (size_t i = 0; i < vectors_to_test.size(); i++) {
        std::cout << "         eta() = " << vectors_to_test.at(i).eta() << "         rap()=" << vectors_to_test.at(i).rap() << " for " << vectors_to_test.at(i) << std::endl;
        std::cout << " Correct eta() = " << correct_eta.at(i)         << " Correct rap()=" << correct_rap.at(i) << std::endl << std::endl;
    }
    for (size_t i=0; i<vectors_to_test.size(); i++)
    {
        std::cout << "Testing " << vectors_to_test.at(i) << std::endl;
        assert(vectors_to_test.at(i).eta() == correct_eta.at(i) );
        assert(vectors_to_test.at(i).rap() == correct_rap.at(i) );

    }
    return 0;
}

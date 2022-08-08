// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include <iostream>
#include <fstream>
#include <vector>

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/Print.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#include "HepMC3TestUtils.h"
using namespace HepMC3;
int main()
{
    //
    // In this example we will place the following event into HepMC "by hand"
    //
    //     name status pdg_id  parent Px       Py    Pz       Energy      Mass
    //  1  !p+!    3   2212    0,0    0.000    0.000 7000.000 7000.000    0.938
    //  2  !p+!    3   2212    0,0    0.000    0.000-7000.000 7000.000    0.938
    //=========================================================================
    //  3  !d!     3      1    1,1    0.750   -1.569   32.191   32.238    0.000
    //  4  !u~!    3     -2    2,2   -3.047  -19.000  -54.629   57.920    0.000
    //  5  !W-!    3    -24    1,2    1.517   -20.68  -20.605   85.925   80.799
    //  6  !gamma! 1     22    1,2   -3.813    0.113   -1.833    4.233    0.000
    //  7  !d!     1      1    5,5   -2.445   28.816    6.082   29.552    0.010
    //  8  !u~!    1     -2    5,5    3.962  -49.498  -26.687   56.373    0.006
    //  9  !gamma! 3     22    3,4    0.000    0.000    0.000    0.000    0.000


    // declare several WriterAscii instances for comparison
    WriterAscii xout1("testLoops1.out");
    // output in old format
    WriterAsciiHepMC2 xout2( "testLoops2.out" );

    // build the graph, which will look like
    //                       p7                   #
    // p1                   /                     #
    //   \v1__p3      p5---v4                     #
    //         \_v3_/       \                     #
    //         /   |\        p8                   #
    //    v2__p4   | \                            #
    //   /  \     /  p6                           #
    // p2    \p9_/                                #
    //
    // define a flow pattern as  p1 -> p3 -> p6
    //                       and p2 -> p4 -> p5
    //

    // First create the event container, with Signal Process 20, event number 1
    //
    GenEvent evt(Units::GEV,Units::MM);
    evt.set_event_number(1);
    evt.add_attribute("signal_process_id", std::make_shared<IntAttribute>(20));
    // create vertex 1
    GenVertexPtr v1 = std::make_shared<GenVertex>();
    evt.add_vertex( v1 );
    GenParticlePtr p1 = std::make_shared<GenParticle>( FourVector(0,0,7000,7000),2212, 3 );
    v1->add_particle_in( p1 );
    p1->add_attribute("flow1", std::make_shared<IntAttribute>(231));
    p1->add_attribute("flow1", std::make_shared<IntAttribute>(231));
    p1->add_attribute("theta", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI));
    p1->add_attribute("phi", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI*2));

    GenVertexPtr v2 = std::make_shared<GenVertex>();
    evt.add_vertex( v2 );
    GenParticlePtr p2 = std::make_shared<GenParticle>(  FourVector(0,0,-7000,7000),2212, 3 );
    v2->add_particle_in( p2 );
    p2->add_attribute("flow1", std::make_shared<IntAttribute>(243));
    p2->add_attribute("theta", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI));
    p2->add_attribute("phi", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI*2));
    //
    // create the outgoing particles of v1 and v2
    GenParticlePtr p3 = std::make_shared<GenParticle>( FourVector(.750,-1.569,32.191,32.238),1, 3 );
    v1->add_particle_out( p3 );
    p3->add_attribute("flow1", std::make_shared<IntAttribute>(231));
    p3->add_attribute("theta", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI));
    p3->add_attribute("phi", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI*2));

    GenParticlePtr p4 = std::make_shared<GenParticle>( FourVector(-3.047,-19.,-54.629,57.920),-2, 3 );
    v2->add_particle_out( p4 );
    p4->add_attribute("flow1", std::make_shared<IntAttribute>(243));
    p4->add_attribute("theta", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI));
    p4->add_attribute("phi", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI*2));
    //
    // create v3
    GenVertexPtr v3 = std::make_shared<GenVertex>();
    evt.add_vertex( v3 );
    v3->add_particle_in( p3 );
    v3->add_particle_in( p4 );
    GenParticlePtr p6 = std::make_shared<GenParticle>(  FourVector(-3.813,0.113,-1.833,4.233 ),22, 1 );
    evt.add_particle( p6 );
    p6->add_attribute("flow1", std::make_shared<IntAttribute>(231));
    p6->add_attribute("theta", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI));
    p6->add_attribute("phi", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI*2));
    v3->add_particle_out( p6 );
    GenParticlePtr p5 = std::make_shared<GenParticle>( FourVector(1.517,-20.68,-20.605,85.925),-24, 3 );
    v3->add_particle_out( p5 );
    p5->add_attribute("flow1", std::make_shared<IntAttribute>(243));
    p5->add_attribute("theta", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI));
    p5->add_attribute("phi", std::make_shared<DoubleAttribute>(std::rand()/double(RAND_MAX)*M_PI*2));

    //
    // create v4
    GenVertexPtr v4 = std::make_shared<GenVertex>(FourVector(0.12,-0.3,0.05,0.004));
    evt.add_vertex( v4 );
    v4->add_particle_in( p5 );
    GenParticlePtr p7(new GenParticle( FourVector(-2.445,28.816,6.082,29.552), 1,1 ));
    v4->add_particle_out( p7 );
    GenParticlePtr p8(new GenParticle( FourVector(3.962,-49.498,-26.687,56.373), -2,1 ));
    v4->add_particle_out( p8 );


    GenParticlePtr ploop = std::make_shared<GenParticle>(  FourVector(0.0,0.0,0.0,0.0 ),21, 3 );
    v3->add_particle_out( ploop );
    v2->add_particle_in( ploop );

    //
    // tell the event which vertex is the signal process vertex
    //evt.set_signal_process_vertex( v3 );
    evt.add_attribute("signal_process_vertex", std::make_shared<IntAttribute>(v3->id()));
    // the event is complete, we now print it out
    Print::content(evt);
    //we now print it out in old format
    Print::listing(evt,8);
    // print each particle so we can see the polarization
    for ( ConstGenParticlePtr ip: evt.particles()) {
        Print::line(ip,true);
    }

    // write event
    xout1.write_event(evt);
    // write event in old format
    xout2.write_event(evt);

    // now clean-up by deleteing all objects from memory
    //
    // deleting the event deletes all contained vertices, and all particles
    // contained in those vertices
    evt.clear();
    xout1.close();
    xout2.close();

    int Nxin1=0;
    ReaderAscii xin1("testLoops1.out");
    if(xin1.failed()) {
        xin1.close();
        return 102;
    }
    while( !xin1.failed() )
    {
        xin1.read_event(evt);
        if( xin1.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        evt.clear();
        Nxin1++;
    }
    xin1.close();

    int Nxin2=0;
    ReaderAsciiHepMC2 xin2("testLoops2.out");
    if(xin2.failed()) {
        xin2.close();
        return 103;
    }
    while( !xin2.failed() )
    {
        xin2.read_event(evt);
        if( xin2.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        evt.clear();
        Nxin2++;
    }
    xin2.close();
    int ret = 10*std::abs(Nxin1-1)+std::abs(Nxin2-1);
    return ret;
}

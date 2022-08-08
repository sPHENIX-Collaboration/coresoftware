// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
///We set some non-default value
#define HEPMC3_HEPEVT_NMXHEP  4000
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/HEPEVT_Wrapper_Runtime.h"
#include "HepMC3/HEPEVT_Wrapper_Runtime_Static.h"
#include "HepMC3/HEPEVT_Wrapper_Template.h"
#include "HepMC3/HEPEVT_Wrapper.h"
#include "HepMC3TestUtils.h"
using namespace HepMC3;

GenEvent generate1() {
    GenEvent evt;
    std::shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>();
    evt.set_run_info(run);
    GenParticlePtr b2 = std::make_shared<GenParticle>( FourVector( 0.0,    0.0,   7000.0,  7000.0  ),2212,  3 );
    GenParticlePtr b1 = std::make_shared<GenParticle>( FourVector( 0.750, -1.569,   32.191,  32.238),   1,  3 );
    GenParticlePtr b3 = std::make_shared<GenParticle>( FourVector( 0.750, -1.569,   32.191,  -32.238),   1,  3 );
    GenVertexPtr v1 = std::make_shared<GenVertex>();
    v1->add_particle_in (b1);
    v1->add_particle_in(b2);
    v1->add_particle_out(b3);
    evt.add_vertex(v1);
    for (size_t z= 0; z < 5; z++) {
        std::vector<GenParticlePtr> particles = evt.particles();
        for (auto p: particles) {
            if (p->end_vertex()) continue;
            GenParticlePtr p2 = std::make_shared<GenParticle>( FourVector( 0.0,    0.0,   7000.0+0.01*evt.particles().size(),  7000.0  ),2212,  3 );
            GenParticlePtr p1 = std::make_shared<GenParticle>( FourVector( 0.750, -1.569,   32.191+0.01*evt.particles().size(),  32.238),   1,  3 );
            GenVertexPtr v = std::make_shared<GenVertex>();
            v->add_particle_in (p);
            v->add_particle_out(p1);
            v->add_particle_out(p2);
            evt.add_vertex(v);
        }
    }
    return evt;
}



int main()
{
    struct HEPEVT_Templated<HEPMC3_HEPEVT_NMXHEP,double>  X;
    GenEvent evt1 = generate1();
    HEPEVT_Wrapper_Runtime  test1;
    test1.set_max_number_entries(HEPMC3_HEPEVT_NMXHEP);
    test1.set_hepevt_address((char*)&X);
    test1.GenEvent_to_HEPEVT(&evt1);

    HEPEVT_Wrapper_Template<HEPMC3_HEPEVT_NMXHEP,double>  test2;
    GenEvent evt2;
    test2.set_hepevt_address((char*)&X);
    test2.HEPEVT_to_GenEvent(&evt2);

    HEPEVT_Wrapper_Template<20000,double>  test3;
    GenEvent evt3;
    test3.allocate_internal_storage();
    test3.copy_to_internal_storage((char*)&X, HEPMC3_HEPEVT_NMXHEP);
    test3.HEPEVT_to_GenEvent(&evt3);

    GenEvent evt4;
    HEPEVT_Wrapper_Runtime_Static::set_max_number_entries(HEPMC3_HEPEVT_NMXHEP);
    HEPEVT_Wrapper_Runtime_Static::set_hepevt_address((char*)&X);
    HEPEVT_Wrapper_Runtime_Static::print_hepevt();
    HEPEVT_Wrapper_Runtime_Static::HEPEVT_to_GenEvent(&evt4);


    GenEvent evt5;
    HEPEVT_Wrapper::set_hepevt_address((char*)&X);
    HEPEVT_Wrapper::HEPEVT_to_GenEvent(&evt5);


    GenEvent evt6;
    HEPEVT_Wrapper_Runtime  test6;
    test6.set_max_number_entries(20000);
    test6.allocate_internal_storage();
    test6.copy_to_internal_storage((char*)&X,HEPMC3_HEPEVT_NMXHEP);
    test6.HEPEVT_to_GenEvent(&evt6);


    std::shared_ptr<WriterAscii> w1 = std::make_shared<WriterAscii>("testHEPEVTWrapper1output1.txt");
    w1->write_event(evt1);
    w1->close();
    std::shared_ptr<WriterAscii> w2= std::make_shared<WriterAscii>("testHEPEVTWrapper1output2.txt");
    w2->write_event(evt2);
    w2->close();
    std::shared_ptr<WriterAscii> w3= std::make_shared<WriterAscii>("testHEPEVTWrapper1output3.txt");
    w3->write_event(evt3);
    w3->close();
    std::shared_ptr<WriterAscii> w4= std::make_shared<WriterAscii>("testHEPEVTWrapper1output4.txt");
    w4->write_event(evt4);
    w4->close();
    std::shared_ptr<WriterAscii> w5= std::make_shared<WriterAscii>("testHEPEVTWrapper1output5.txt");
    w5->write_event(evt5);
    w5->close();
    std::shared_ptr<WriterAscii> w6= std::make_shared<WriterAscii>("testHEPEVTWrapper1output6.txt");
    w6->write_event(evt6);
    w6->close();
    return COMPARE_ASCII_FILES("testHEPEVTWrapper1output1.txt","testHEPEVTWrapper1output2.txt") +
           COMPARE_ASCII_FILES("testHEPEVTWrapper1output2.txt","testHEPEVTWrapper1output3.txt") +
           COMPARE_ASCII_FILES("testHEPEVTWrapper1output3.txt","testHEPEVTWrapper1output4.txt") +
           COMPARE_ASCII_FILES("testHEPEVTWrapper1output4.txt","testHEPEVTWrapper1output5.txt") +
           COMPARE_ASCII_FILES("testHEPEVTWrapper1output5.txt","testHEPEVTWrapper1output6.txt")
           ;
}

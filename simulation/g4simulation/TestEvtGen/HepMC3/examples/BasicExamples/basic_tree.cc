// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/// @example basic_tree.cc
/// @brief Basic example of building HepMC3 tree by hand
///
///  Based on HepMC2/examples/example_BuildEventFromScratch.cc

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"

using namespace HepMC3;


/** Main program */
int main() {
    //
    // In this example we will place the following event into HepMC "by hand"
    //
    //     name status pdg_id  parent Px       Py    Pz       Energy      Mass
    //  1  !p+!    3   2212    0,0    0.000    0.000 7000.000 7000.000    0.938
    //  3  !p+!    3   2212    0,0    0.000    0.000-7000.000 7000.000    0.938
    //=========================================================================
    //  2  !d!     3      1    1,1    0.750   -1.569   32.191   32.238    0.000
    //  4  !u~!    3     -2    2,2   -3.047  -19.000  -54.629   57.920    0.000
    //  5  !W-!    3    -24    1,2    1.517   -20.68  -20.605   85.925   80.799
    //  6  !gamma! 1     22    1,2   -3.813    0.113   -1.833    4.233    0.000
    //  7  !d!     1      1    5,5   -2.445   28.816    6.082   29.552    0.010
    //  8  !u~!    1     -2    5,5    3.962  -49.498  -26.687   56.373    0.006

    // now we build the graph, which will looks like
    //                       p7                         #
    // p1                   /                           #
    //   \v1__p2      p5---v4                           #
    //         \_v3_/       \                           #
    //         /    \        p8                         #
    //    v2__p4     \                                  #
    //   /            p6                                #
    // p3                                               #
    //                                                  #
    GenEvent evt(Units::GEV,Units::MM);

    //                                                               px      py        pz       e     pdgid status
    GenParticlePtr p1 = std::make_shared<GenParticle>( FourVector( 0.0,    0.0,   7000.0,  7000.0  ),2212,  3 );
    GenParticlePtr p2 = std::make_shared<GenParticle>( FourVector( 0.750, -1.569,   32.191,  32.238),   1,  3 );
    GenParticlePtr p3 = std::make_shared<GenParticle>( FourVector( 0.0,    0.0,  -7000.0,  7000.0  ),2212,  3 );
    GenParticlePtr p4 = std::make_shared<GenParticle>( FourVector(-3.047,-19.0,    -54.629,  57.920),  -2,  3 );

    GenVertexPtr v1 = std::make_shared<GenVertex>();
    v1->add_particle_in (p1);
    v1->add_particle_out(p2);
    evt.add_vertex(v1);

    // Set vertex status if needed
    v1->set_status(4);

    GenVertexPtr v2 = std::make_shared<GenVertex>();
    v2->add_particle_in (p3);
    v2->add_particle_out(p4);
    evt.add_vertex(v2);

    GenVertexPtr v3 = std::make_shared<GenVertex>();
    v3->add_particle_in(p2);
    v3->add_particle_in(p4);
    evt.add_vertex(v3);

    GenParticlePtr p5 = std::make_shared<GenParticle>( FourVector(-3.813,  0.113, -1.833, 4.233),  22, 1 );
    GenParticlePtr p6 = std::make_shared<GenParticle>( FourVector( 1.517,-20.68, -20.605,85.925), -24, 3 );

    v3->add_particle_out(p5);
    v3->add_particle_out(p6);

    GenVertexPtr v4 =std:: make_shared<GenVertex>();
    v4->add_particle_in (p6);
    evt.add_vertex(v4);

    GenParticlePtr p7 = std::make_shared<GenParticle>( FourVector(-2.445, 28.816,  6.082,29.552),  1, 1 );
    GenParticlePtr p8 = std::make_shared<GenParticle>( FourVector( 3.962,-49.498,-26.687,56.373), -2, 1 );

    v4->add_particle_out(p7);
    v4->add_particle_out(p8);

    //
    // Example of adding event attributes
    //
    std::shared_ptr<GenPdfInfo> pdf_info = std::make_shared<GenPdfInfo>();
    evt.add_attribute("GenPdfInfo",pdf_info);

    pdf_info->set(1,2,3.4,5.6,7.8,9.0,1.2,3,4);

    std::shared_ptr<GenHeavyIon> heavy_ion = std::make_shared<GenHeavyIon>();
    evt.add_attribute("GenHeavyIon",heavy_ion);

    heavy_ion->set( 1,2,3,4,5,6,7,8,9,0.1,2.3,4.5,6.7);

    std::shared_ptr<GenCrossSection> cross_section = std::make_shared<GenCrossSection>();
    evt.add_attribute("GenCrossSection",cross_section);

    cross_section->set_cross_section(1.2,3.4);

    //
    // Example of manipulating the attributes
    //

    std::cout << std::endl << " Manipulating attributes:" << std::endl;

    // get attribute
    std::shared_ptr<GenCrossSection> cs = evt.attribute<GenCrossSection>("GenCrossSection");

    // if attribute exists - do something with it
    if(cs) {
        cs->set_cross_section(-1.0,0.0);
        Print::line(cs);
    }
    else std::cout << "Problem accessing attribute!" <<std::endl;

    // remove attribute
    evt.remove_attribute("GenCrossSection");
    evt.remove_attribute("GenCrossSection"); // This call will do nothing

    // now this should be null
    cs = evt.attribute<GenCrossSection>("GenCrossSection");

    if(!cs)std::cout << "Successfully removed attribute" <<std::endl;
    else   std::cout << "Problem removing attribute!" <<std::endl;

    //
    // Example of adding attributes and finding particles with attributes
    //

    std::shared_ptr<Attribute> tool1           = std::make_shared<IntAttribute>(1);
    std::shared_ptr<Attribute> tool999         = std::make_shared<IntAttribute>(999);
    std::shared_ptr<Attribute> test_attribute  = std::make_shared<StringAttribute>("test attribute");
    std::shared_ptr<Attribute> test_attribute2 = std::make_shared<StringAttribute>("test attribute2");

    p2->add_attribute( "tool" ,  tool1           );
    p2->add_attribute( "other" , test_attribute  );

    p4->add_attribute( "tool" ,  tool1           );

    p6->add_attribute( "tool" ,  tool999         );
    p6->add_attribute( "other" , test_attribute2 );

    v3->add_attribute( "vtx_att" , test_attribute );
    v4->add_attribute( "vtx_att" , test_attribute2 );
/* TODO: Make this code portable

    std::cout << std::endl << "Find all particles with attribute 'tool' "<< std::endl;
    std::cout << "(should return particles 2,4,6):" << std::endl;

    /// @todo can we add some utility funcs to simplify creation of Features from Attributes and check they exist.
    /// Features and Attributes are quite similar concepts anyway, can they be unified (but Features can also be
    ///  non-attribute-like e.g. pT, rapidity or any quantity it is possible to obtain from a particle)

    for(ConstGenParticlePtr p: applyFilter(Selector::ATTRIBUTE("tool"), evt.particles())){
      Print::line(p);
    }

    std::cout <<std::endl << "Find all particles with attribute 'tool' equal 1 "<< std::endl;
    std::cout << "(should return particles 2,4):" <<std::endl;

    for(ConstGenParticlePtr p: applyFilter(Selector::ATTRIBUTE("tool") && Selector::ATTRIBUTE("tool") == tool1, evt.particles())){
      Print::line(p);
    }

    std::cout << std::endl << "Find all particles with a string attribute 'other' equal 'test attribute' "<< std::endl;
    std::cout << "(should return particle 2):" << std::endl;


    for(ConstGenParticlePtr p: applyFilter(Selector::ATTRIBUTE("other") && Selector::ATTRIBUTE("other") == "test_attribute", evt.particles())){
      Print::line(p);
    }
*/

    std::cout << std::endl << "Offsetting event position by 5,5,5,5" << std::endl;

    evt.shift_position_by( FourVector(5,5,5,5) );

    Print::listing(evt);

    std::cout << std::endl << "Printing full content of the GenEvent object " << std::endl
                 << "(including particles and vertices in one-line format):" << std::endl << std::endl;

    Print::content(evt);

   std::cout <<std::endl << "Now: removing particle with id 6 and printing again:" <<std::endl <<std::endl;
    evt.remove_particle(p6);

    Print::listing(evt);
    Print::content(evt);

   std::cout <<std::endl << "Now: removing beam particles, leaving an empty event" <<std::endl <<std::endl;
    evt.remove_particles( evt.beams() );

    Print::listing(evt);
    Print::content(evt);
    return 0;
}

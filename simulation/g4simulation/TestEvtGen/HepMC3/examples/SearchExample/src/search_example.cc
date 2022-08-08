// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Relatives.h"

HepMC3::GenEvent* generate_event(const size_t  n,const  size_t  iterations)
{
    HepMC3::GenEvent* e=new HepMC3::GenEvent();
    HepMC3::GenVertexPtr v0=std::make_shared<HepMC3::GenVertex>();
    e->add_vertex(v0);
    size_t it=0;
    for (;;)
    {
        if (it>iterations)
        {
            for (auto v: e->vertices())
            {
                if (v->particles_out().size()!=0) continue;
                for (size_t i=0; i<n; i++) v->add_particle_out(std::make_shared<HepMC3::GenParticle>());
            }
            break;
        }
        auto vertices=e->vertices();
        for (auto v: vertices)
        {
            if (v->particles_out().size()!=0) continue;
            for (size_t i=0; i<n; i++) v->add_particle_out(std::make_shared<HepMC3::GenParticle>());
            for (auto p: v->particles_out())
            {
                HepMC3::GenVertexPtr vx=std::make_shared<HepMC3::GenVertex>();
                vx->add_particle_in(p);
                e->add_vertex(vx);
            }
        }
        it++;
    }
    return e;
}

int main()
{
    std::cout<<"search_example: start"<<std::endl;
    auto start0 = std::chrono::system_clock::now();
    for (int i=0; i<10000; i++)
    {
        HepMC3::GenEvent* evt=generate_event(3,4);
        delete evt;
    }
    auto end0 = std::chrono::system_clock::now();
    std::cout<<"search_example: generation of events                                  "<<std::chrono::duration_cast<std::chrono::milliseconds>(end0 - start0).count()<<" ms"<<std::endl;

    auto start1 = std::chrono::system_clock::now();
    size_t  np1=0;
    for (int i=0; i<10000; i++)
    {
        HepMC3::GenEvent* evt=generate_event(3,4);
        for (auto p: evt->particles()) np1+=HepMC3::descendant_particles(p).size();
        delete evt;
    }
    auto end1 = std::chrono::system_clock::now();
    std::cout<<"search_example: generation of events and descendants_of_same_type()  "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count()<<" ms"
             << " total number of decandants: "<<np1<<std::endl;

    auto start2 = std::chrono::system_clock::now();
    size_t  np2=0;
    for (int i=0; i<10000; i++)
    {
        HepMC3::GenEvent* evt=generate_event(3,4);
        for (auto p: evt->particles()) np2+=(HepMC3::Relatives::DESCENDANTS(p)).size();
        delete evt;
    }
    auto end2 = std::chrono::system_clock::now();
    std::cout<<"search_example: generation of events and Relatives::DESCENDANTS()    "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count()<<" ms"
             << " total number of decandants: "<<np2<<std::endl;

    auto start3 = std::chrono::system_clock::now();
    size_t  np3=0;
    for (int i=0; i<10000; i++)
    {
        HepMC3::GenEvent* evt=generate_event(3,4);
        for (auto p: evt->particles()) np3+=HepMC3::ancestor_particles(p).size();
        delete evt;
    }
    auto end3 = std::chrono::system_clock::now();
    std::cout<<"search_example: generation of events and ancestors_of_same_type()    "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count()<<" ms"
             << " total number of ancestors: "<<np3<<std::endl;


    auto start4 = std::chrono::system_clock::now();
    size_t  np4=0;
    for (int i=0; i<10000; i++)
    {
        HepMC3::GenEvent* evt=generate_event(3,4);
        for (auto p: evt->particles()) np4+=(HepMC3::Relatives::ANCESTORS(p)).size();
        delete evt;
    }
    auto end4 = std::chrono::system_clock::now();
    std::cout<<"search_example: generation of events and Relatives::ANCESTORS()      "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(end4 - start4).count()<<" ms"
             << " total number of ancestors: "<<np4<<std::endl;

    auto start1o = std::chrono::system_clock::now();
    size_t  np1o=0;
    for (int i=0; i<10000; i++)
    {
        HepMC3::GenEvent* evt=generate_event(3,4);
        for (auto p: evt->particles()) np1o+=HepMC3::descendant_vertices(p).size();
        delete evt;
    }
    auto end1o = std::chrono::system_clock::now();
    std::cout<<"search_example: generation of events and descendants_of_other_type() "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(end1o - start1o).count()<<" ms"
             << " total number of decandants: "<<np1o<<std::endl;


    auto start3o = std::chrono::system_clock::now();
    size_t  np3o=0;
    for (int i=0; i<10000; i++)
    {
        HepMC3::GenEvent* evt=generate_event(3,4);
        for (auto p: evt->particles()) np3o+=HepMC3::ancestor_vertices(p).size();
        delete evt;
    }
    auto end3o = std::chrono::system_clock::now();
    std::cout<<"search_example: generation of events and ancestors_of_other_type()   "
             <<std::chrono::duration_cast<std::chrono::milliseconds>(end3o - start3o).count()<<" ms"
             << " total number of decandants: "<<np3o<<std::endl;

    std::cout<<"search_example: end"<<std::endl;
    return EXIT_SUCCESS;
}

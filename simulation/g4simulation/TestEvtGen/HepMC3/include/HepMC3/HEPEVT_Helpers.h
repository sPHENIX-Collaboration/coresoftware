// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPEVT_HELPERS_H
#define HEPEVT_HELPERS_H
/**
 *  @file HEPEVT_Helpers.h
 *  @brief Helper functions used to manipulate with HEPEVT block
 *
 */
#if defined(WIN32)&&(!defined(HEPMC3_NO_EXPORTS))
#if defined(HepMC3_EXPORTS)
#define HEPMC3_EXPORT_API __declspec(dllexport)
#else
#define HEPMC3_EXPORT_API __declspec(dllimport)
#endif
#else
#define HEPMC3_EXPORT_API
#endif
#include <algorithm>
#include <map>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
namespace HepMC3
{

/** @struct HEPEVT_Templated
 *  @brief  C structure representing Fortran common block HEPEVT
 * T. Sjöstrand et al., "A proposed standard event record",
 *  in `Z physics at LEP 1', eds. G. Altarelli, R. Kleiss and C. Verzegnassi,
 * Geneva, Switzerland, September 4-5, 1989, CERN 89-08 (Geneva, 1989), Vol. 3, p. 327
 * Disk representation is given by Fortran WRITE/READ format.
 */
template <int max_particles, typename momentum_type = double>
struct HEPEVT_Templated
{
    int        nevhep;             //!< Event number
    int        nhep;               //!< Number of entries in the event
    int        isthep[max_particles];     //!< Status code
    int        idhep [max_particles];     //!< PDG ID
    int        jmohep[max_particles][2];  //!< Position of 1st and 2nd (or last!) mother
    int        jdahep[max_particles][2];  //!< Position of 1nd and 2nd (or last!) daughter
    momentum_type phep  [max_particles][5];  //!< Momentum: px, py, pz, e, m
    momentum_type vhep  [max_particles][4];  //!< Time-space position: x, y, z, t
};

/** @struct HEPEVT_Pointers
 *  @brief  C structure representing Fortran common block HEPEVT
 * T. Sjöstrand et al., "A proposed standard event record",
 *  in `Z physics at LEP 1', eds. G. Altarelli, R. Kleiss and C. Verzegnassi,
 * Geneva, Switzerland, September 4-5, 1989, CERN 89-08 (Geneva, 1989), Vol. 3, p. 327
 * Disk representation is given by Fortran WRITE/READ format.
 */
template<typename momentum_type = double>
struct HEPEVT_Pointers
{
    int*        nevhep;             //!< Pointer to Event number
    int*        nhep;               //!< Pointer to Number of entries in the event
    int*        isthep;     //!< Pointer to Status code
    int*        idhep;      //!< Pointer to PDG ID
    int*        jmohep;     //!< Pointer to position of 1st and 2nd (or last!) mother
    int*        jdahep;     //!< Pointer to position of 1nd and 2nd (or last!) daughter
    momentum_type* phep;    //!< Pointer to momentum: px, py, pz, e, m
    momentum_type* vhep;    //!< Pointer to time-space position: x, y, z, t
};

/** @brief comparison of two particles */
struct GenParticlePtr_greater
{
    /** @brief comparison of two particles */
    bool operator()(ConstGenParticlePtr lx, ConstGenParticlePtr rx) const;
};
/** @brief  Order vertices with equal paths. */
struct pair_GenVertexPtr_int_greater
{
    /** @brief  Order vertices with equal paths. If the paths are equal, order in other quantities.
     * We cannot use id, as it can be assigned in different way*/
    bool operator()(const std::pair<ConstGenVertexPtr, int>& lx, const std::pair<ConstGenVertexPtr, int>& rx) const;
};
/** @brief Calculates the path to the top (beam) particles */
void calculate_longest_path_to_top(ConstGenVertexPtr v, std::map<ConstGenVertexPtr, int>& pathl);


/** @brief  Converts HEPEVT into GenEvent. */
template <class T>
bool HEPEVT_to_GenEvent_nonstatic(GenEvent* evt, T* A)
{
    if ( !evt ) { std::cerr << "HEPEVT_to_GenEvent_nonstatic  - passed null event." << std::endl; return false;}
    evt->set_event_number(A->event_number());
    std::map<GenParticlePtr, int > hepevt_particles;
    std::map<int, GenParticlePtr > particles_index;
    std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > > hepevt_vertices;
    std::map<int, GenVertexPtr > vertex_index;
    int ne=0;
    ne=A->number_entries();
    for ( int i = 1; i <= ne; i++ )
    {
        GenParticlePtr p = std::make_shared<GenParticle>();
        p->set_momentum(FourVector(A->px(i), A->py(i), A->pz(i), A->e(i)));
        p->set_status(A->status(i));
        p->set_pid(A->id(i)); //Confusing!
        p->set_generated_mass(A->m(i));
        hepevt_particles[p] = i;
        particles_index[i] = p;
        GenVertexPtr v = std::make_shared<GenVertex>();
        v->set_position(FourVector(A->x(i), A->y(i), A->z(i), A->t(i)));
        v->add_particle_out(p);
        std::set<int> in;
        std::set<int> out;
        out.insert(i);
        vertex_index[i] = v;
        hepevt_vertices[v] = std::pair<std::set<int>, std::set<int> >(in, out);
    }
    /* The part above is always correct as it is a raw information without any topology.*/

    /* In this way we trust mother information. The "Trust daughters" is not implemented.*/
    for (std::map<GenParticlePtr, int >::iterator it1 = hepevt_particles.begin(); it1 != hepevt_particles.end(); ++it1)
        for (std::map<GenParticlePtr, int >::iterator it2 = hepevt_particles.begin(); it2 != hepevt_particles.end(); ++it2) {
            if   (A->first_parent(it2->second) <= it1->second && it1->second <= A->last_parent(it2->second)) hepevt_vertices[it2->first->production_vertex()].first.insert(it1->second);
        }
    /* Now all incoming sets are correct for all vertices. But we have duplicates.*/

    /* Disconnect all particles from the vertices*/
    for ( int i = 1; i <= A->number_entries(); i++ ) vertex_index[i]->remove_particle_out(particles_index[i]);

    /*Fill container with vertices with unique sets of incoming particles. Merge the outcoming particle sets.*/
    std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > > final_vertices_map;
    for (std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > >::iterator vs = hepevt_vertices.begin(); vs != hepevt_vertices.end(); ++vs)
    {
        if ((final_vertices_map.size() == 0) || (vs->second.first.size() == 0 && vs->second.second.size() != 0)) { final_vertices_map.insert(*vs);  continue; } /*Always insert particles out of nowhere*/
        std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > >::iterator  v2;
        for (v2 = final_vertices_map.begin(); v2 != final_vertices_map.end(); ++v2) if (vs->second.first == v2->second.first) {v2->second.second.insert(vs->second.second.begin(), vs->second.second.end()); break;}
        if (v2 == final_vertices_map.end()) final_vertices_map.insert(*vs);
    }

    std::vector<GenParticlePtr> final_particles;
    std::set<int> used;
    for (std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > >:: iterator it = final_vertices_map.begin(); it != final_vertices_map.end(); ++it)
    {
        GenVertexPtr v = it->first;
        std::set<int> in = it->second.first;
        std::set<int> out = it->second.second;
        used.insert(in.begin(), in.end());
        used.insert(out.begin(), out.end());
        for (std::set<int>::iterator el = in.begin(); el != in.end(); ++el) v->add_particle_in(particles_index[*el]);
        if (in.size() !=0 ) for (std::set<int>::iterator el = out.begin(); el != out.end(); ++el) v->add_particle_out(particles_index[*el]);
    }
    for (std::set<int>::iterator el = used.begin(); el != used.end(); ++el) final_particles.push_back(particles_index[*el]);
    /* One can put here a check on the number of particles/vertices*/

    evt->add_tree(final_particles);

    return true;
}
/** @brief  Converts GenEvent into HEPEVT. */
template <class T>
bool GenEvent_to_HEPEVT_nonstatic(const GenEvent* evt, T* A)
{
    /// This writes an event out to the HEPEVT common block. The daughters
    /// field is NOT filled, because it is possible to contruct graphs
    /// for which the mothers and daughters cannot both be make sequential.
    /// This is consistent with how pythia fills HEPEVT (daughters are not
    /// necessarily filled properly) and how IO_HEPEVT reads HEPEVT.
    if ( !evt ) return false;

    /*AV Sorting the vertices by the lengths of their longest incoming paths assures the mothers will not go before the daughters*/
    /* Calculate all paths*/
    std::map<ConstGenVertexPtr, int> longest_paths;
    for (ConstGenVertexPtr v: evt->vertices()) calculate_longest_path_to_top(v, longest_paths);
    /* Sort paths*/
    std::vector<std::pair<ConstGenVertexPtr, int> > sorted_paths;
    std::copy(longest_paths.begin(), longest_paths.end(), std::back_inserter(sorted_paths));
    std::sort(sorted_paths.begin(), sorted_paths.end(), pair_GenVertexPtr_int_greater());

    std::vector<ConstGenParticlePtr> sorted_particles;
    std::vector<ConstGenParticlePtr> stable_particles;
    /*For a valid "Trust mothers" HEPEVT record we must  keep mothers together*/
    for (std::pair<ConstGenVertexPtr, int> it: sorted_paths)
    {
        std::vector<ConstGenParticlePtr> Q = it.first->particles_in();
        std::sort(Q.begin(), Q.end(), GenParticlePtr_greater());
        std::copy(Q.begin(), Q.end(), std::back_inserter(sorted_particles));
        /*For each vertex put all outgoing particles w/o end vertex. Ordering of particles to produces reproduceable record*/
        for (ConstGenParticlePtr pp: it.first->particles_out())
            if (!(pp->end_vertex())) stable_particles.push_back(pp);
    }
    std::sort(stable_particles.begin(), stable_particles.end(), GenParticlePtr_greater());
    std::copy(stable_particles.begin(), stable_particles.end(), std::back_inserter(sorted_particles));

    int particle_counter;
    particle_counter = std::min(int(sorted_particles.size()), A->max_number_entries());
    /* fill the HEPEVT event record (MD code)*/
    A->set_event_number(evt->event_number());
    A->set_number_entries(particle_counter);
    for ( int i = 1; i <= particle_counter; ++i )
    {
        A->set_status(i, sorted_particles[i-1]->status());
        A->set_id(i, sorted_particles[i-1]->pid());
        FourVector m = sorted_particles[i-1]->momentum();
        A->set_momentum(i, m.px(), m.py(), m.pz(), m.e());
        A->set_mass(i, sorted_particles[i-1]->generated_mass());
        if ( sorted_particles[i-1]->production_vertex()  &&
                sorted_particles[i-1]->production_vertex()->particles_in().size())
        {
            FourVector p = sorted_particles[i-1]->production_vertex()->position();
            A->set_position(i, p.x(), p.y(), p.z(), p.t() );
            std::vector<int> mothers;
            mothers.clear();

            for (ConstGenParticlePtr it: sorted_particles[i-1]->production_vertex()->particles_in())
                for ( int j = 1; j <= particle_counter; ++j )
                    if (sorted_particles[j-1] == (it))
                        mothers.push_back(j);
            std::sort(mothers.begin(), mothers.end());
            if (mothers.size() == 0)
                mothers.push_back(0);
            if (mothers.size() == 1) mothers.push_back(mothers[0]);

            A->set_parents(i, mothers.front(), mothers.back());
        }
        else
        {
            A->set_position(i, 0, 0, 0, 0);
            A->set_parents(i, 0, 0);
        }
        A->set_children(i, 0, 0);
    }
    return true;
}


/** @brief  Converts  HEPEVT into GenEvent. */
template <class T>
bool HEPEVT_to_GenEvent_static(GenEvent* evt)
{
    if ( !evt ) { std::cerr << "HEPEVT_to_GenEvent_static  - passed null event." << std::endl; return false;}
    evt->set_event_number(T::event_number());
    std::map<GenParticlePtr, int > hepevt_particles;
    std::map<int, GenParticlePtr > particles_index;
    std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > > hepevt_vertices;
    std::map<int, GenVertexPtr > vertex_index;
    int ne=0;
    ne=T::number_entries();
    for ( int i = 1; i <= ne; i++ )
    {
        GenParticlePtr p = std::make_shared<GenParticle>();
        p->set_momentum(FourVector(T::px(i), T::py(i), T::pz(i), T::e(i)));
        p->set_status(T::status(i));
        p->set_pid(T::id(i)); //Confusing!
        p->set_generated_mass(T::m(i));
        hepevt_particles[p] = i;
        particles_index[i] = p;
        GenVertexPtr v = std::make_shared<GenVertex>();
        v->set_position(FourVector(T::x(i), T::y(i), T::z(i), T::t(i)));
        v->add_particle_out(p);
        std::set<int> in;
        std::set<int> out;
        out.insert(i);
        vertex_index[i] = v;
        hepevt_vertices[v] = std::pair<std::set<int>, std::set<int> >(in, out);
    }
    /* The part above is always correct as it is a raw information without any topology.*/

    /* In this way we trust mother information. The "Trust daughters" is not implemented.*/
    for (std::map<GenParticlePtr, int >::iterator it1 = hepevt_particles.begin(); it1 != hepevt_particles.end(); ++it1)
        for (std::map<GenParticlePtr, int >::iterator it2 = hepevt_particles.begin(); it2 != hepevt_particles.end(); ++it2) {
            if   (T::first_parent(it2->second) <= it1->second && it1->second <= T::last_parent(it2->second)) hepevt_vertices[it2->first->production_vertex()].first.insert(it1->second);
        }
    /* Now all incoming sets are correct for all vertices. But we have duplicates.*/

    /* Disconnect all particles from the vertices*/
    for ( int i = 1; i <= T::number_entries(); i++ ) vertex_index[i]->remove_particle_out(particles_index[i]);

    /*Fill container with vertices with unique sets of incoming particles. Merge the outcoming particle sets.*/
    std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > > final_vertices_map;
    for (std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > >::iterator vs = hepevt_vertices.begin(); vs != hepevt_vertices.end(); ++vs)
    {
        if ((final_vertices_map.size() == 0) || (vs->second.first.size() == 0 && vs->second.second.size() != 0)) { final_vertices_map.insert(*vs);  continue; } /*Always insert particles out of nowhere*/
        std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > >::iterator  v2;
        for (v2 = final_vertices_map.begin(); v2 != final_vertices_map.end(); ++v2) if (vs->second.first == v2->second.first) {v2->second.second.insert(vs->second.second.begin(), vs->second.second.end()); break;}
        if (v2 == final_vertices_map.end()) final_vertices_map.insert(*vs);
    }

    std::vector<GenParticlePtr> final_particles;
    std::set<int> used;
    for (std::map<GenVertexPtr, std::pair<std::set<int>, std::set<int> > >:: iterator it = final_vertices_map.begin(); it != final_vertices_map.end(); ++it)
    {
        GenVertexPtr v = it->first;
        std::set<int> in = it->second.first;
        std::set<int> out = it->second.second;
        used.insert(in.begin(), in.end());
        used.insert(out.begin(), out.end());
        for (std::set<int>::iterator el = in.begin(); el != in.end(); ++el) v->add_particle_in(particles_index[*el]);
        if (in.size() !=0 ) for (std::set<int>::iterator el = out.begin(); el != out.end(); ++el) v->add_particle_out(particles_index[*el]);
    }
    for (std::set<int>::iterator el = used.begin(); el != used.end(); ++el) final_particles.push_back(particles_index[*el]);
    /* One can put here a check on the number of particles/vertices*/

    evt->add_tree(final_particles);

    return true;
}

/** @brief  Converts GenEvent into HEPEVT. */
template <class T>
bool GenEvent_to_HEPEVT_static(const GenEvent* evt)
{
    /// This writes an event out to the HEPEVT common block. The daughters
    /// field is NOT filled, because it is possible to contruct graphs
    /// for which the mothers and daughters cannot both be make sequential.
    /// This is consistent with how pythia fills HEPEVT (daughters are not
    /// necessarily filled properly) and how IO_HEPEVT reads HEPEVT.
    if ( !evt ) return false;

    /*AV Sorting the vertices by the lengths of their longest incoming paths assures the mothers will not go before the daughters*/
    /* Calculate all paths*/
    std::map<ConstGenVertexPtr, int> longest_paths;
    for (ConstGenVertexPtr v: evt->vertices()) calculate_longest_path_to_top(v, longest_paths);
    /* Sort paths*/
    std::vector<std::pair<ConstGenVertexPtr, int> > sorted_paths;
    std::copy(longest_paths.begin(), longest_paths.end(), std::back_inserter(sorted_paths));
    std::sort(sorted_paths.begin(), sorted_paths.end(), pair_GenVertexPtr_int_greater());

    std::vector<ConstGenParticlePtr> sorted_particles;
    std::vector<ConstGenParticlePtr> stable_particles;
    /*For a valid "Trust mothers" HEPEVT record we must  keep mothers together*/
    for (std::pair<ConstGenVertexPtr, int> it: sorted_paths)
    {
        std::vector<ConstGenParticlePtr> Q = it.first->particles_in();
        std::sort(Q.begin(), Q.end(), GenParticlePtr_greater());
        std::copy(Q.begin(), Q.end(), std::back_inserter(sorted_particles));
        /*For each vertex put all outgoing particles w/o end vertex. Ordering of particles to produces reproduceable record*/
        for (ConstGenParticlePtr pp: it.first->particles_out())
            if (!(pp->end_vertex())) stable_particles.push_back(pp);
    }
    std::sort(stable_particles.begin(), stable_particles.end(), GenParticlePtr_greater());
    std::copy(stable_particles.begin(), stable_particles.end(), std::back_inserter(sorted_particles));

    int particle_counter;
    particle_counter = std::min(int(sorted_particles.size()), T::max_number_entries());
    /* fill the HEPEVT event record (MD code)*/
    T::set_event_number(evt->event_number());
    T::set_number_entries(particle_counter);
    for ( int i = 1; i <= particle_counter; ++i )
    {
        T::set_status(i, sorted_particles[i-1]->status());
        T::set_id(i, sorted_particles[i-1]->pid());
        FourVector m = sorted_particles[i-1]->momentum();
        T::set_momentum(i, m.px(), m.py(), m.pz(), m.e());
        T::set_mass(i, sorted_particles[i-1]->generated_mass());
        if ( sorted_particles[i-1]->production_vertex()  &&
                sorted_particles[i-1]->production_vertex()->particles_in().size())
        {
            FourVector p = sorted_particles[i-1]->production_vertex()->position();
            T::set_position(i, p.x(), p.y(), p.z(), p.t() );
            std::vector<int> mothers;
            mothers.clear();

            for (ConstGenParticlePtr it: sorted_particles[i-1]->production_vertex()->particles_in())
                for ( int j = 1; j <= particle_counter; ++j )
                    if (sorted_particles[j-1] == (it))
                        mothers.push_back(j);
            std::sort(mothers.begin(), mothers.end());
            if (mothers.size() == 0)
                mothers.push_back(0);
            if (mothers.size() == 1) mothers.push_back(mothers[0]);

            T::set_parents(i, mothers.front(), mothers.back());
        }
        else
        {
            T::set_position(i, 0, 0, 0, 0);
            T::set_parents(i, 0, 0);
        }
        T::set_children(i, 0, 0);
    }
    return true;
}

}
#endif

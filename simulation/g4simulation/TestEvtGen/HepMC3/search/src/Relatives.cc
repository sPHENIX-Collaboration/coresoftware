// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file Relatives.cc
/// @brief Implementation of \b Relatives class
///
#include "HepMC3/Relatives.h"

namespace HepMC3 {
const Parents     Relatives::PARENTS;
const Children    Relatives::CHILDREN;
#ifdef _MSC_VER
const Ancestors   Relatives::ANCESTORS;
const Descendants Relatives::DESCENDANTS;
#else
thread_local const Ancestors   Relatives::ANCESTORS;
thread_local const Descendants Relatives::DESCENDANTS;
#endif
}

namespace HepMC3 {
/// @brief Returns children of vertex, i.e. outgoing particles.
std::vector<HepMC3::GenParticlePtr>      children(HepMC3::GenVertexPtr O) {
    if (O) return O->particles_out();
    return  std::vector<HepMC3::GenParticlePtr>();
}
/// @brief Returns children of const vertex, i.e. outgoing particles.
std::vector<HepMC3::ConstGenParticlePtr> children(HepMC3::ConstGenVertexPtr O) {
    if (O) return O->particles_out();
    return  std::vector<HepMC3::ConstGenParticlePtr>();
}
/// @brief Returns children of particle, i.e. the end vertex.
std::vector<HepMC3::GenVertexPtr>        children(HepMC3::GenParticlePtr O) {
    std::vector<HepMC3::GenVertexPtr> result;
    if (O->end_vertex()) result.push_back(O->end_vertex());
    return result;
}
/// @brief Returns children of const particle, i.e. the end vertex.
std::vector<HepMC3::ConstGenVertexPtr>   children(HepMC3::ConstGenParticlePtr O) {
    std::vector<HepMC3::ConstGenVertexPtr> result;
    if (O->end_vertex()) result.push_back(O->end_vertex());
    return result;
}
/// @brief Returns grandchildren of particle, i.e. the outgoing particles of the end vertex.
std::vector<HepMC3::GenParticlePtr>      grandchildren(HepMC3::GenParticlePtr O) {
    if (O) if (O->end_vertex()) return O->end_vertex()->particles_out();
    return std::vector<HepMC3::GenParticlePtr> ();
}
/// @brief Returns grandchildren of const particle, i.e. the outgoing particles of the end vertex.
std::vector<HepMC3::ConstGenParticlePtr> grandchildren(HepMC3::ConstGenParticlePtr O) {
    if (O) if (O->end_vertex()) return O->end_vertex()->particles_out();
    return std::vector<HepMC3::ConstGenParticlePtr> ();
}
/// @brief Returns grandchildren of vertex, i.e. the end vertices of the outgoing particles.
std::vector<HepMC3::GenVertexPtr>        grandchildren(HepMC3::GenVertexPtr O) {
    std::vector<HepMC3::GenVertexPtr> result;
    if (O) for (auto o: O->particles_out()) if (o->end_vertex()) result.push_back(o->end_vertex());
    return result;
}
/// @brief Returns grandchildren of const vertex, i.e. the end vertices of the outgoing particles.
std::vector<HepMC3::ConstGenVertexPtr>   grandchildren(HepMC3::ConstGenVertexPtr O) {
    std::vector<HepMC3::ConstGenVertexPtr> result;
    if (O)  for (auto o:O->particles_out()) if (o->end_vertex()) result.push_back(o->end_vertex());
    return result;
}
/// @brief Returns parents of vertex, i.e. incoming particles.
std::vector<HepMC3::GenParticlePtr>      parents(HepMC3::GenVertexPtr O) {
    if (O) return O->particles_in();
    return  std::vector<GenParticlePtr>();
}
/// @brief Returns parents of const vertex, i.e. incoming particles.
std::vector<HepMC3::ConstGenParticlePtr> parents(HepMC3::ConstGenVertexPtr O) {
    if (O) return O->particles_in();
    return  std::vector<HepMC3::ConstGenParticlePtr>();
}
/// @brief Returns parents of particle, i.e. production vertex.
std::vector<HepMC3::GenVertexPtr>        parents(HepMC3::GenParticlePtr O) {
    std::vector<HepMC3::GenVertexPtr> result;
    if (O->production_vertex()) result.push_back(O->production_vertex());
    return result;
}
/// @brief Returns parents of const particle, i.e. production vertex.
std::vector<HepMC3::ConstGenVertexPtr>   parents(HepMC3::ConstGenParticlePtr O) {
    std::vector<HepMC3::ConstGenVertexPtr> result;
    if (O->production_vertex()) result.push_back(O->production_vertex());
    return result;
}
/// @brief Returns grandparents of particle, i.e. incoming particles of production vertex.
std::vector<HepMC3::GenParticlePtr>      grandparents(HepMC3::GenParticlePtr O) {
    if (O) if (O->production_vertex()) return O->production_vertex()->particles_in();
    return std::vector<HepMC3::GenParticlePtr> ();
}
/// @brief Returns grandparents of const particle, i.e. incoming particles of production vertex.
std::vector<HepMC3::ConstGenParticlePtr> grandparents(HepMC3::ConstGenParticlePtr O) {
    if (O) if (O->production_vertex()) return O->production_vertex()->particles_in();
    return std::vector<HepMC3::ConstGenParticlePtr> ();
}
/// @brief Returns grandparents of vertex, i.e. production vertices of incoming particles.
std::vector<HepMC3::GenVertexPtr>        grandparents(HepMC3::GenVertexPtr O) {
    std::vector<HepMC3::GenVertexPtr> result;
    if (O) for (auto o: O->particles_in()) if (o->production_vertex()) result.push_back(o->production_vertex());
    return result;
}
/// @brief Returns grandparents of const vertex, i.e. production vertices of incoming particles.
std::vector<HepMC3::ConstGenVertexPtr>   grandparents(HepMC3::ConstGenVertexPtr O) {
    std::vector<HepMC3::ConstGenVertexPtr> result;
    if (O)  for (auto o:O->particles_in()) if (o->end_vertex()) result.push_back(o->production_vertex());
    return result;
}
/// @brief Returns descendands of the same type, i.e. vertices for vertex and particles for particle
template <class O>  std::vector<O> descendants_of_same_type(O obj)
{
    std::vector<O>  result = grandchildren(obj);
    size_t gc = 0;
    for (;;)
    {
        std::vector<O>  temp;
        for (; gc < result.size(); gc++)
        {
            auto  temp0 = grandchildren(result[gc]);
            temp.insert(temp.end(), temp0.begin(), temp0.end());
        }
        for (auto p2: temp) if (std::find(result.begin(), result.end(), p2) == result.end()) result.push_back(p2);
        if (gc >= result.size()) break;
    }
    return result;
}
/// @brief Returns descendands of the other type, i.e. vertices for  particle and particles for vertex
template <class O, class R>  std::vector<R> descendants_of_other_type(O obj)
{
    std::vector<R> localchildren = children(obj);
    std::vector<R>  result = localchildren;
    for (auto c: localchildren)
    {
        std::vector<R> desc = descendants_of_same_type(c);
        for (auto d: desc) if (std::find(result.begin(), result.end(), d) == result.end()) result.push_back(d);
    }
    return result;
}
/// @brief Returns ancestors of the same type, i.e. vertices for vertex and particles for particle
template <class O>  std::vector<O> ancestors_of_same_type(O obj)
{
    std::vector<O>  result = grandparents(obj);
    size_t gc = 0;
    for (;;)
    {
        std::vector<O>  temp;
        for (; gc < result.size(); gc++)
        {
            auto  temp0 = grandparents(result[gc]);
            temp.insert(temp.end(), temp0.begin(), temp0.end());
        }
        for (auto p2: temp) if (std::find(result.begin(), result.end(), p2) == result.end()) result.push_back(p2);
        if (gc >= result.size()) break;
    }
    return result;
}
/// @brief Returns ancestors of the other type, i.e. vertices for  particle and particles for vertex
template <class O, class R>  std::vector<R> ancestors_of_other_type(O obj)
{
    std::vector<R> localparents = parents(obj);
    std::vector<R>  result = localparents;
    for (auto c: localparents)
    {
        std::vector<R> desc = ancestors_of_same_type(c);
        for (auto d: desc) if (std::find(result.begin(), result.end(), d) == result.end()) result.push_back(d);
    }
    return result;
}

std::vector<HepMC3::ConstGenParticlePtr> descendant_particles(HepMC3::ConstGenVertexPtr obj) {
    return  descendants_of_other_type<HepMC3::ConstGenVertexPtr, HepMC3::ConstGenParticlePtr>(obj);
}
std::vector<HepMC3::GenParticlePtr> descendant_particles(HepMC3::GenVertexPtr obj) {
    return descendants_of_other_type<HepMC3::GenVertexPtr, HepMC3::GenParticlePtr>(obj);
}

std::vector<ConstGenVertexPtr> descendant_vertices(HepMC3::ConstGenParticlePtr obj) {
    return descendants_of_other_type<HepMC3::ConstGenParticlePtr, HepMC3::ConstGenVertexPtr>(obj);
}
std::vector<HepMC3::GenVertexPtr> descendant_vertices(HepMC3::GenParticlePtr obj) {
    return descendants_of_other_type<HepMC3::GenParticlePtr, HepMC3::GenVertexPtr>(obj);
}

std::vector<HepMC3::ConstGenParticlePtr> ancestor_particles(HepMC3::ConstGenVertexPtr obj) {
    return  ancestors_of_other_type<HepMC3::ConstGenVertexPtr, HepMC3::ConstGenParticlePtr>(obj);
}
std::vector<HepMC3::GenParticlePtr> ancestor_particles(HepMC3::GenVertexPtr obj) {
    return ancestors_of_other_type<HepMC3::GenVertexPtr, HepMC3::GenParticlePtr>(obj);
}

std::vector<HepMC3::ConstGenVertexPtr> ancestor_vertices(HepMC3::ConstGenParticlePtr obj) {
    return ancestors_of_other_type<HepMC3::ConstGenParticlePtr, HepMC3::ConstGenVertexPtr>(obj);
}
std::vector<HepMC3::GenVertexPtr> ancestor_vertices(HepMC3::GenParticlePtr obj) {
    return ancestors_of_other_type<HepMC3::GenParticlePtr, HepMC3::GenVertexPtr>(obj);
}


std::vector<HepMC3::ConstGenParticlePtr> descendant_particles(HepMC3::ConstGenParticlePtr obj)  { return descendants_of_same_type(obj); }
std::vector<HepMC3::GenParticlePtr>      descendant_particles(HepMC3::GenParticlePtr obj)       { return descendants_of_same_type(obj); }
std::vector<HepMC3::ConstGenVertexPtr>   descendant_vertices(HepMC3::ConstGenVertexPtr obj)     { return descendants_of_same_type(obj); }
std::vector<HepMC3::GenVertexPtr>        descendant_vertices(HepMC3::GenVertexPtr obj)          { return descendants_of_same_type(obj); }
std::vector<HepMC3::ConstGenParticlePtr> ancestor_particles(HepMC3::ConstGenParticlePtr obj)    { return ancestors_of_same_type(obj); }
std::vector<HepMC3::GenParticlePtr>      ancestor_particles(HepMC3::GenParticlePtr obj)         { return ancestors_of_same_type(obj); }
std::vector<HepMC3::ConstGenVertexPtr>   ancestor_vertices(HepMC3::ConstGenVertexPtr obj)       { return ancestors_of_same_type(obj); }
std::vector<HepMC3::GenVertexPtr>        ancestor_vertices(HepMC3::GenVertexPtr obj)            { return ancestors_of_same_type(obj); }
std::vector<HepMC3::GenParticlePtr>      children_particles(HepMC3::GenVertexPtr O)             { return children(O); }
std::vector<HepMC3::ConstGenParticlePtr> children_particles(HepMC3::ConstGenVertexPtr O)        { return children(O); }
std::vector<HepMC3::GenVertexPtr>        children_vertices(HepMC3::GenParticlePtr O)            { return children(O); }
std::vector<HepMC3::ConstGenVertexPtr>   children_vertices(HepMC3::ConstGenParticlePtr O)       { return children(O); }
std::vector<HepMC3::GenParticlePtr>      grandchildren_particles(HepMC3::GenParticlePtr O)      { return grandchildren(O); }
std::vector<HepMC3::ConstGenParticlePtr> grandchildren_particles(HepMC3::ConstGenParticlePtr O) { return grandchildren(O); }
std::vector<HepMC3::GenVertexPtr>        grandchildren_vertices(HepMC3::GenVertexPtr O)         { return grandchildren(O); }
std::vector<HepMC3::ConstGenVertexPtr>   grandchildren_vertices(HepMC3::ConstGenVertexPtr O)    { return grandchildren(O); }
std::vector<HepMC3::GenParticlePtr>      parent_particles(HepMC3::GenVertexPtr O)               { return parents(O); }
std::vector<HepMC3::ConstGenParticlePtr> parent_particles(HepMC3::ConstGenVertexPtr O)          { return parents(O); }
std::vector<HepMC3::GenVertexPtr>        parent_vertices(HepMC3::GenParticlePtr O)              { return parents(O); }
std::vector<HepMC3::ConstGenVertexPtr>   parent_vertices(HepMC3::ConstGenParticlePtr O)         { return parents(O); }
std::vector<HepMC3::GenParticlePtr>      grandparent_particles(HepMC3::GenParticlePtr O)        { return grandparents(O); }
std::vector<HepMC3::ConstGenParticlePtr> grandparent_particles(HepMC3::ConstGenParticlePtr O)   { return grandparents(O); }
std::vector<HepMC3::GenVertexPtr>        grandparent_vertices(HepMC3::GenVertexPtr O)           { return grandparents(O); }
std::vector<HepMC3::ConstGenVertexPtr>   grandparent_vertices(HepMC3::ConstGenVertexPtr O)      { return grandparents(O); }


} // namespace HepMC3

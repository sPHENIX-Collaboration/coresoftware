// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file Filter.h
/// @brief Defines Filter operations for combingin Filters
///
#ifndef HEPMC3_FILTER_H
#define HEPMC3_FILTER_H

#include <vector>
#include <functional>
#include "HepMC3/GenParticle.h"

namespace HepMC3 {
/// @brief type of Filter
using Filter = std::function<bool(ConstGenParticlePtr)>;

/// @brief Apply a Filter to a list of GenParticles
/// Returns a vector of GenParticles that satisfy the Filter
inline std::vector<GenParticlePtr> applyFilter(const Filter &filter, const std::vector<GenParticlePtr> &particles) {
    std::vector<GenParticlePtr> result;
    for (GenParticlePtr p: particles) {
        if (filter(p)) result.push_back(p);
    }
    return result;
}

/// @brief Apply a Filter to a list of ConstGenParticles
/// Returns a vector of ConstGenParticles that satisfy the Filter
inline std::vector<ConstGenParticlePtr> applyFilter(const Filter &filter, const std::vector<ConstGenParticlePtr> &particles) {
    std::vector<ConstGenParticlePtr> result;
    for (ConstGenParticlePtr p: particles) {
        if (filter(p)) result.push_back(p);
    }
    return result;
}

/// @brief A Filter that will accept all particles
/// This might be needed if a signature requires a default Filter
inline bool ACCEPT_ALL(ConstGenParticlePtr /*dummy*/) {
    return true;
}

/// @brief The logical AND of two Filters is itself a Filter
inline Filter operator && (const Filter & lhs, const Filter &rhs) {
    return [lhs, rhs](ConstGenParticlePtr p)->bool{return lhs(p) && rhs(p); };
}

/// @brief The logical OR of two Filters is itself a Filter
inline Filter operator || (const Filter & lhs, const Filter &rhs) {
    return [lhs, rhs](ConstGenParticlePtr p)->bool{return lhs(p) || rhs(p); };
}

/// @brief The negation of a Filter is itself a Filter
inline Filter operator !(const Filter &rhs) {
    return [rhs](ConstGenParticlePtr p)->bool{return !(rhs(p));};
}

}
#endif

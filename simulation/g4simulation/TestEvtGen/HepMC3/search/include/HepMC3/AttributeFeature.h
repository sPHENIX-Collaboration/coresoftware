// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file AttributeFeature.h
/// @brief Defines AttributeFeature for obtaining Filters to search by Attribute.
/// @class HepMC3::AttributeFeature
/// @brief AttributeFeature

#ifndef HEPMC3_ATTRIBUTE_FEATURE_H
#define HEPMC3_ATTRIBUTE_FEATURE_H

#include <memory>
#include <string>
#include "HepMC3/Attribute.h"
#include "HepMC3/Filter.h"

namespace HepMC3 {

class AttributeFeature {
public:
    /// @brief constructor
    AttributeFeature(const std::string &name): m_name(name) {}

    /// @brief existence
    Filter exists() const {
        std::string name = m_name;
        return [name](ConstGenParticlePtr p)->bool{return p->attribute_as_string(name).length() != 0;};
    }

    /// @brief evaluate
    bool operator()(ConstGenParticlePtr p) const {
        return p->attribute_as_string(m_name).length() != 0;
    }

    /// @brief equality operator
    Filter operator == (const Attribute &rhs) const {
        std::string name = m_name;
        std::string other;
        rhs.to_string(other);
        return [other, name](ConstGenParticlePtr p)->bool{return p->attribute_as_string(name).compare(other) == 0;};
    }

    /// @brief equality operator
    Filter operator == (std::shared_ptr<const Attribute> rhs) const {
        std::string name = m_name;
        std::string other;
        rhs->to_string(other);
        return [other, name](ConstGenParticlePtr p)->bool{return p->attribute_as_string(name).compare(other) == 0;};
    }

    /// @brief equality operator
    Filter operator == (std::string rhs) const {
        const std::string &name = m_name;
        return [name, rhs](ConstGenParticlePtr p)->bool{return p->attribute_as_string(name).compare(rhs) == 0;};
    }

private:
    std::string m_name;  ///< holds name
};
}
#endif

// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_SEARCH_FILTEATTRIBUTE_H
#define HEPMC3_SEARCH_FILTEATTRIBUTE_H
///
/// @file FilterAttribute.h
/// @brief Definition of \b class ATTRIBUTE
///
/// @class HepMC3::ATTRIBUTE
/// @brief Filter for the attributes
///
/// Used to construct filters that can check if an attribute exists
/// or to compare against other attribute.
///
/// @ingroup searchengine
#include <string>
#include <memory>
#include "HepMC3/Filter.h"
#include "HepMC3/Attribute.h"

namespace HepMC3 {
/** Deprecated */
using std::string;

class ATTRIBUTE : public Filter {
//
// Constructors
//
public:
    /// @brief Default constructor
    ///
    /// Provides the name of the attribute used in by the filter
    ATTRIBUTE(const std::string &name):Filter(ATTRIBUTE_EXISTS, name) {}

//
// Operators
//
public:
    /// @brief Compare if this attribute is equal to other attribute
    Filter& operator==(std::shared_ptr<Attribute> &at) {
        m_attribute = ATTRIBUTE_IS_EQUAL;
        at->to_string(m_attribute_str);
        return *this;
    }

    /// @brief Compare if this attribute is not equal to other attribute
    Filter& operator!=(std::shared_ptr<Attribute> &at) {
        m_bool_value = !m_bool_value;
        m_attribute  = ATTRIBUTE_IS_EQUAL;
        at->to_string(m_attribute_str);
        return *this;
    }

    /// @brief Compare if string version of this attribute is equal value
    Filter& operator==(const std::string &value) {
        m_attribute     = ATTRIBUTE_IS_EQUAL;
        m_attribute_str = value;
        return *this;
    }

    /// @brief Compare if string version of this attribute is not equal value
    Filter& operator!=(const std::string &value) {
        m_bool_value    = !m_bool_value;
        m_attribute     = ATTRIBUTE_IS_EQUAL;
        m_attribute_str = value;
        return *this;
    }

    /// @brief Negate logic of the result (eg. check if attribute does not exist)
    Filter& operator!() {
        m_bool_value = !m_bool_value;
        return *this;
    }
};

} // namespace HepMC3

#endif

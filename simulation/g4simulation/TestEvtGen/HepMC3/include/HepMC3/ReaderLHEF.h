// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERLHEF_H
#define HEPMC3_READERLHEF_H
/**
 *  @file  ReaderLHEF.h
 *  @brief Definition of \b class ReaderLHEF
 *
 *  @class HepMC3::ReaderLHEF
 *  @brief GenEvent I/O parsing and serialization for LHEF files
 *
 *
 *  @ingroup IO
 *
 */
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"
#include <deque>
#include <string>
#include <fstream>
#include <istream>
#include "HepMC3/LHEFAttributes.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include <iomanip>


namespace HepMC3
{
class ReaderLHEF : public Reader
{
public:
    /// The ctor to read from stream
    ReaderLHEF(std::istream &);
    /** @brief Constructor */
    ReaderLHEF(const std::string& filename);
    /// The ctor to read from temp stream
    ReaderLHEF(std::shared_ptr<std::istream> s_stream);
    /// @brief skip events
    bool skip(const int)  override;
    /** @brief Reading event */
    bool read_event(GenEvent& ev)  override;
    /** @brief Close */
    void close()  override;
    /** @brief State */
    bool failed()  override;
    /** @brief Destructor */
    ~ReaderLHEF() ;
private:
    void init();                       ///< Init helper
    std::shared_ptr<std::istream> m_shared_stream; ///< Holds temporary stream
    std::shared_ptr<LHEF::Reader> m_reader;            ///< The actual reader
    std::shared_ptr<HEPRUPAttribute> m_hepr; ///< Holder of attributes
    int m_neve;                         ///< Event counter
    bool m_failed;                      ///< State of reader
    std::deque<GenEvent> m_storage; ///<storage used for subevents.
};
}
#endif

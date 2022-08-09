// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERHEPEVT_H
#define HEPMC3_READERHEPEVT_H
/**
 *  @file  ReaderHEPEVT.h
 *  @brief Definition of \b class ReaderHEPEVT
 *
 *  @class HepMC3::ReaderHEPEVT
 *  @brief GenEvent I/O parsing and serialization for HEPEVT files
 *
 *
 *  @ingroup IO
 *
 */
#include <set>
#include <string>
#include <fstream>
#include <istream>
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/HEPEVT_Wrapper_Template.h"

namespace HepMC3
{

class ReaderHEPEVT : public Reader
{
//
// Constructors
//
public:
    /** @brief Default constructor */
    ReaderHEPEVT(const std::string &filename);
    /// The ctor to read from stream
    ReaderHEPEVT(std::istream &);
    /// The ctor to read from temp stream
    ReaderHEPEVT(std::shared_ptr<std::istream> s_stream);
//
// Functions
//
public:
    /** @brief Find and read event header line  from file
    *
    */
    virtual bool read_hepevt_event_header();
    /** @brief read particle from file
    *
    * @param[in] i Particle id
    */
    virtual bool read_hepevt_particle(int i);

    /// @brief skip events
    bool skip(const int)  override;


    /** @brief Read event from file*/
    bool read_event(GenEvent &evt)  override;


    /** @brief Close file stream */
    void close()  override;

    /** @brief Get stream error state */
    bool failed()  override;

    char* hepevtbuffer; //!< Pointer to HEPEVT Fortran common block/C struct
private:
    std::ifstream m_file; //!< Input file
    std::shared_ptr<std::istream> m_shared_stream; //!< For ctor when reading from temp stream
    std::istream* m_stream; //!< For ctor when reading from stream
    bool m_isstream; //!< toggles usage of m_file or m_stream
    HEPEVT_Wrapper_Template<100000> m_hepevt_interface; //!< Templated HEPEVT interface
};

} // namespace HepMC3

#endif

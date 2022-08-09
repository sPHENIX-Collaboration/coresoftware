// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_WRITERHEPEVT_H
#define HEPMC3_WRITERHEPEVT_H
/**
 *  @file  WriterHEPEVT.h
 *  @brief Definition of \b class WriterHEPEVT
 *
 *  @class HepMC3::WriterHEPEVT
 *  @brief GenEvent I/O serialization for HEPEVT files
 *
 *
 *  @ingroup IO
 *
 */
#include <fstream>
#include "HepMC3/Writer.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/HEPEVT_Wrapper_Template.h"
namespace HepMC3
{

class WriterHEPEVT : public Writer
{
//
// Constructors
//
public:
    /** @brief Default constructor
     *  @warning If file exists, it will be overwritten
     */
    WriterHEPEVT(const std::string &filename,
                 std::shared_ptr<GenRunInfo> run = nullptr);

    /// @brief Constructor from ostream
    WriterHEPEVT(std::ostream& stream,
                 std::shared_ptr<GenRunInfo> run = nullptr);
    /// @brief Constructor from temp ostream
    WriterHEPEVT(std::shared_ptr<std::ostream> s_stream,
                 std::shared_ptr<GenRunInfo> run = nullptr);
//
// Functions
//
public:

    /** @brief Write particle to file
     *
     *  @param[in] index Particle to be serialized
     *  @param[in] iflong Format of record
     */

    virtual void write_hepevt_particle( int index, bool iflong = true );
    /** @brief Write event header to file
     *
     */
    virtual void write_hepevt_event_header();

    /** @brief Write event to file
     *
     *  @param[in] evt Event to be serialized
     */
    void write_event(const GenEvent &evt)  override;

    /** @brief Close file stream */
    void close()  override;

    /** @brief Get stream error state flag */
    bool failed()  override;
    /** @brief  set flag if vertex positions are available.
     *  Effectively this adds or removes key "vertices_positions_are_absent"
     *  to/from the m_options.*/
    void set_vertices_positions_present(bool iflong);

    /** @brief  get flag if vertex positions are available.
     * The flag is deduced from m_options. If the m_options have the key
     * "vertices_positions_are_absent" the result if false. True otherwise. */
    bool get_vertices_positions_present() const;

protected:
    std::ofstream m_file; //!< Output file
    std::shared_ptr<std::ostream> m_shared_stream;///< Output temp. stream
    std::ostream* m_stream; //!< Output stream
    char* hepevtbuffer;   //!< Pointer to HEPEVT Fortran common block/C struct
    int   m_events_count; //!< Events count. Needed to generate unique object name
    HEPEVT_Wrapper_Template<100000> m_hepevt_interface; //!< Templated HEPEVT interface
};

} // namespace HepMC3
#endif

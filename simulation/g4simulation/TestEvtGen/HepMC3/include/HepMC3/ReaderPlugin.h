// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERPLUGIN_H
#define HEPMC3_READERPLUGIN_H
/**
 *  @file  ReaderPlugin.h
 *  @brief Definition of \b class ReaderPlugin
 *
 *  @class HepMC3::ReaderPlugin
 *  @brief GenEvent I/O parsing and serialization using external plugin
 *
 *
 *  @ingroup IO
 *
 */
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"
namespace HepMC3
{
class ReaderPlugin : public Reader
{
public:
    /** @brief Constructor  to read from stream*/
    ReaderPlugin(std::istream & stream,const std::string &libname, const std::string &newreader);
    /** @brief Constructor to read from file*/
    ReaderPlugin(const std::string& filename,const std::string &libname, const std::string &newreader);
    /** @brief Reading event */
    bool read_event(GenEvent& ev)  override {if (!m_reader) return false; return m_reader->read_event(ev);};
    /** @brief Close */
    void close() override { if (!m_reader) return; m_reader->close(); };
    /** @brief State */
    bool failed() override {if (!m_reader) return true; return m_reader->failed();};
    /** @brief Get the global GenRunInfo object. */
    std::shared_ptr<GenRunInfo> run_info() const  override  { return m_reader?m_reader->run_info():nullptr; }
    /** @brief  Set options */
    void set_options(const std::map<std::string, std::string>& options)  override { if (!m_reader) return; else m_reader->set_options(options); }
    /** @brief  Get options  */
    std::map<std::string, std::string> get_options()  const  override { return m_reader?m_reader->get_options(): std::map<std::string, std::string>();  }
    /** @brief Destructor */
    ~ReaderPlugin()  override;
//protected:
    /// Set the global GenRunInfo object.
    void set_run_info(std::shared_ptr<GenRunInfo> run) override { if (!m_reader) return; else m_reader->set_run_info(run); }
private:
    Reader* m_reader; ///< The actual reader
    void*  dll_handle; ///< library handler
};
}
#endif

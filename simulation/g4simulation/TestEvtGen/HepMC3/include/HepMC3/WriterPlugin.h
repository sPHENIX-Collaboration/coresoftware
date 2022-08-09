// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_WRITERPLUGIN_H
#define HEPMC3_WRITERPLUGIN_H
/**
 *  @file  WriterPlugin.h
 *  @brief Definition of \b class WriterPlugin
 *
 *  @class HepMC3::WriterPlugin
 *  @brief GenEvent I/O parsing and serialization using external plugin
 *
 *
 *  @ingroup IO
 *
 */
#include "HepMC3/Writer.h"
#include "HepMC3/GenEvent.h"
namespace HepMC3
{
class WriterPlugin : public Writer
{
public:

    /** @brief Constructor  to read from stream */
    WriterPlugin(std::ostream & stream,const std::string &libname, const std::string &newwriter, std::shared_ptr<HepMC3::GenRunInfo> run = std::shared_ptr<GenRunInfo>());

    /** @brief Constructor to read from file */
    WriterPlugin(const std::string& filename,const std::string &libname, const std::string &newwriter, std::shared_ptr<HepMC3::GenRunInfo> run = std::shared_ptr<GenRunInfo>());

    /** @brief Reading event */
    void write_event(const GenEvent& ev)  override {if (!m_writer) return; return m_writer->write_event(ev);};
    /** @brief Close */
    void close() override { if (!m_writer) return; m_writer->close();};
    /** @brief State */
    bool failed() override {if (!m_writer) return true; return m_writer->failed();};
    /** @brief Get the global GenRunInfo object. */
    std::shared_ptr<GenRunInfo> run_info()  const override { return m_writer?m_writer->run_info():nullptr; }
    /** @brief  Set options */
    void set_options(const std::map<std::string, std::string>& options)  override { if (!m_writer) return; else m_writer->set_options(options); }
    /** @brief  Get options  */
    std::map<std::string, std::string> get_options()  const  override { return m_writer?m_writer->get_options(): std::map<std::string, std::string>();  }
    /// Set the global GenRunInfo object.
    void set_run_info(std::shared_ptr<GenRunInfo> run) override { if (!m_writer) return; else m_writer->set_run_info(run); }
    /** @brief Destructor */
    ~WriterPlugin()  override;
private:
    Writer* m_writer; ///< The actual writer
    void* dll_handle; ///< library handler
};
}
#endif

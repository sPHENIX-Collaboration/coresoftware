// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "OutputValidationTool.h"

OutputValidationTool::OutputValidationTool( const std::string &filename ):m_filename(filename),m_timer("IO time") {
}

void OutputValidationTool::initialize() {
    HEPMC2CODE( m_file = new IO_GenEvent((std::string("outputHepMC2")+m_filename).c_str(), std::ios::out); )
    HEPMC3CODE( std::shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>(); m_file = new WriterAsciiHepMC2(std::string("outputHepMC3")+m_filename,run);)
}

int OutputValidationTool::process(GenEvent &hepmc) {
    m_timer.start();
    HEPMC2CODE(HepMC::GenEvent* hepmcevt=&hepmc; (*m_file)<<hepmcevt;)
    HEPMC3CODE(m_file->write_event(hepmc);)
    return 0;
}

void OutputValidationTool::finalize() {
    delete m_file;
}

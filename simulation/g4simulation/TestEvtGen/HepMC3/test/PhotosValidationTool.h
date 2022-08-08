// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef PHOTOS_VALIDATION_TOOL_H
#define PHOTOS_VALIDATION_TOOL_H

#ifdef HEPMC2
#include "Photos/PhotosHepMCEvent.h"
#include "HepMC/GenEvent.h"
#else
#include "Photos/PhotosHepMC3Event.h"
#include "HepMC3/GenEvent.h"
#endif // ifdef HEPMC2

#include "ValidationTool.h"
#include "Timer.h"

#include "Photos/Photos.h"
#include "Photos/Log.h"
/// @class PhotosValidationTool
/// @brief Interface for validatio to Photos
class PhotosValidationTool : public ValidationTool {
public:
    PhotosValidationTool();

public:
    const  std::string name()    { return "Photos++"; }
    bool   tool_modifies_event() { return true;       }
    Timer* timer()               { return &m_timer;   }

    void initialize();
    int  process(GenEvent &hepmc);
    void finalize();

private:
    static const int MAX_PHOTONS_TO_KEEP_TRACK_OF = 4;  ///< Number of tracked photons
    int    m_photons_added[MAX_PHOTONS_TO_KEEP_TRACK_OF]; ///< Added photons
    int    m_more_photons_added;                          ///< More added photons
    Timer  m_timer; ///< Timer
};

#endif

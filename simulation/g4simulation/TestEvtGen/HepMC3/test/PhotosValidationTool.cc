// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "PhotosValidationTool.h"

#include <cstring> // memset
#include <cstdio> // printf

PhotosValidationTool::PhotosValidationTool():m_more_photons_added(0),m_timer("Photos++ processing time") {
    memset(m_photons_added,0,sizeof(int)*MAX_PHOTONS_TO_KEEP_TRACK_OF);
}

void PhotosValidationTool::initialize() {
    Photospp::Photos::initialize();
    Photospp::Photos::setInfraredCutOff(0.001/200);

    HEPMC2CODE( Photospp::Photos::createHistoryEntries(false,3); )
    HEPMC3CODE( Photospp::Photos::createHistoryEntries(false,3);  )
}

int PhotosValidationTool::process(GenEvent &hepmc) {

    HEPMC2CODE( int buf = -hepmc.particles_size(); )
    HEPMC3CODE(
        int buf=-hepmc.particles().size();
    )
    // Time only processing
    m_timer.start();

    // Process by Photos++
    HEPMC2CODE( Photospp::PhotosHepMCEvent  *p_event = new Photospp::PhotosHepMCEvent (&hepmc); )
    HEPMC3CODE( Photospp::PhotosHepMC3Event *p_event = new Photospp::PhotosHepMC3Event(&hepmc); )

    p_event->process();
    delete p_event;

    m_timer.stop();

    // Check number of photons created
    HEPMC2CODE( buf += hepmc.particles_size(); )

    HEPMC3CODE(
        buf +=hepmc.particles().size();
    )

    if(buf<MAX_PHOTONS_TO_KEEP_TRACK_OF) ++m_photons_added[buf];
    else                                 ++m_more_photons_added;

    return 0;
}

void PhotosValidationTool::finalize() {
    Photospp::Log::Summary();

    int sum = m_more_photons_added;
    for(int i=0; i<MAX_PHOTONS_TO_KEEP_TRACK_OF; ++i) sum += m_photons_added[i];

    if( sum == 0 ) sum = 1;

    printf("---------------------------------------------------\n");
    printf(" Number of photons added by Photos++ (per event):\n");
    printf("---------------------------------------------------\n");
    for(int i=0; i<MAX_PHOTONS_TO_KEEP_TRACK_OF; ++i) {
        printf("%5i: %7i events (%6.2f%%)\n",i,m_photons_added[i],m_photons_added[i]*100./sum );
    }
    printf(" more: %7i events (%6.2f%%)\n",m_more_photons_added,m_more_photons_added*100./sum );
    printf("total: %7i events\n",sum );
    printf("---------------------------------------------------\n");
}

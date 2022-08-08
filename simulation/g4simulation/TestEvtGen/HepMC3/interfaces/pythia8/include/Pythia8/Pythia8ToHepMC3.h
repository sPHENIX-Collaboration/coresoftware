// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
// Part of code was adopted from Pythia8-HepMC interface by Mikhail Kirsanov.
#ifndef Pythia8_Pythia8ToHepMC3_H
#define Pythia8_Pythia8ToHepMC3_H
#ifdef  _MSC_VER
#pragma message("HepMC3 interface is available in the latest versions of Pythia8, see https://pythia.org. This interface will be removed in the future HepMC3 versions.")
#else
#warning "HepMC3 interface is available in the latest versions of Pythia8, see https://pythia.org. This interface will be removed in the future HepMC3 versions."
#endif
#include <vector>
#include "Pythia8/Pythia.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenEvent.h"

namespace HepMC3 {


class Pythia8ToHepMC3 {

public:

    // Constructor and destructor
    Pythia8ToHepMC3(): m_internal_event_number(0),
        m_print_inconsistency(true),
        m_free_parton_warnings(true),
        m_crash_on_problem(false),
        m_convert_gluon_to_0(false),
        m_store_pdf(true),
        m_store_proc(true),
        m_store_xsec(true),
        m_store_weights(true) {}

    virtual ~Pythia8ToHepMC3() {}

    // The recommended method to convert Pythia events into HepMC ones
    bool fill_next_event( Pythia8::Pythia& pythia, GenEvent* evt, int ievnum = -1 )
    {
        return fill_next_event( pythia.event, evt, ievnum, &pythia.info, &pythia.settings);
    }

    // Alternative method to convert Pythia events into HepMC ones
#if defined(PYTHIA_VERSION_INTEGER) && (PYTHIA_VERSION_INTEGER>8299)
    bool fill_next_event( Pythia8::Event& pyev, GenEvent* evt,
                          int ievnum = -1, const Pythia8::Info* pyinfo = 0,
                          Pythia8::Settings* pyset = 0)
#else
    bool fill_next_event( Pythia8::Event& pyev, GenEvent* evt,
                          int ievnum = -1, Pythia8::Info* pyinfo = 0,
                          Pythia8::Settings* pyset = 0)
#endif
    {

        // 1. Error if no event passed.
        if (!evt) {
            std::cerr << "Pythia8ToHepMC3::fill_next_event error - passed null event."
                      << std::endl;
            return false;
        }

        // Event number counter.
        if ( ievnum >= 0 ) {
            evt->set_event_number(ievnum);
            m_internal_event_number = ievnum;
        }
        else {
            evt->set_event_number(m_internal_event_number);
            ++m_internal_event_number;
        }

        evt->set_units(Units::GEV,Units::MM);

        // 2. Fill particle information
        std::vector<GenParticlePtr> hepevt_particles;
        hepevt_particles.reserve( pyev.size() );

        for(int i=0; i<pyev.size(); ++i) {
            hepevt_particles.push_back( std::make_shared<GenParticle>( FourVector( pyev[i].px(), pyev[i].py(),
                                        pyev[i].pz(), pyev[i].e() ),
                                        pyev[i].id(), pyev[i].statusHepMC() )
                                      );
            hepevt_particles[i]->set_generated_mass( pyev[i].m() );
        }

        // 3. Fill vertex information and find beam particles.
        std::vector<GenVertexPtr> vertex_cache;
        std::vector<GenParticlePtr> beam_particles;
        for(int  i=0; i<pyev.size(); ++i) {

            std::vector<int> mothers = pyev[i].motherList();

            if(mothers.size()) {
                GenVertexPtr prod_vtx = hepevt_particles[mothers[0]]->end_vertex();

                if(!prod_vtx) {
                    prod_vtx = std::make_shared<GenVertex>();
                    vertex_cache.push_back(prod_vtx);

                    for(unsigned int j=0; j<mothers.size(); ++j) {
                        prod_vtx->add_particle_in( hepevt_particles[mothers[j]] );
                    }
                }
                FourVector prod_pos( pyev[i].xProd(), pyev[i].yProd(),pyev[i].zProd(), pyev[i].tProd() );

                // Update vertex position if necessary
                if(!prod_pos.is_zero() && prod_vtx->position().is_zero()) prod_vtx->set_position( prod_pos );

                prod_vtx->add_particle_out( hepevt_particles[i] );
            }
            else beam_particles.push_back(hepevt_particles[i]);
        }

        // Reserve memory for the event
        evt->reserve( hepevt_particles.size(), vertex_cache.size() );

        // Add particles and vertices in topological order
        if (beam_particles.size()!=2) {
            std::cerr << "There are  " << beam_particles.size() <<"!=2 particles without mothers"<< std::endl;
            if ( m_crash_on_problem ) exit(1);
        }
        evt->add_tree( beam_particles );
        //Attributes should be set after adding the particles to event
        for(int  i=0; i<pyev.size(); ++i) {
            /* TODO: Set polarization */
            // Colour flow uses index 1 and 2.
            int colType = pyev[i].colType();
            if (colType ==  -1 ||colType ==  1 || colType == 2)
            {
                int flow1=0, flow2=0;
                if (colType ==  1 || colType == 2) flow1=pyev[i].col();
                if (colType == -1 || colType == 2) flow2=pyev[i].acol();
                hepevt_particles[i]->add_attribute("flow1",std::make_shared<IntAttribute>(flow1));
                hepevt_particles[i]->add_attribute("flow2",std::make_shared<IntAttribute>(flow2));
            }
        }

        // If hadronization switched on then no final coloured particles.
        bool doHadr = (pyset == 0) ? m_free_parton_warnings : pyset->flag("HadronLevel:all") && pyset->flag("HadronLevel:Hadronize");

        // 4. Check for particles which come from nowhere, i.e. are without
        // mothers or daughters. These need to be attached to a vertex, or else
        // they will never become part of the event.
        for (int  i = 1; i < pyev.size(); ++i) {

            // Check for particles not added to the event
            // NOTE: We have to check if this step makes any sense in HepMC event standard
            if ( !hepevt_particles[i]->in_event() ) {
                std::cerr << "hanging particle " << i << std::endl;
                GenVertexPtr prod_vtx = std::make_shared<GenVertex>();
                prod_vtx->add_particle_out( hepevt_particles[i] );
                evt->add_vertex(prod_vtx);
            }

            // Also check for free partons (= gluons and quarks; not diquarks?).
            if ( doHadr && m_free_parton_warnings ) {
                if ( hepevt_particles[i]->pid() == 21 && !hepevt_particles[i]->end_vertex() ) {
                    std::cerr << "gluon without end vertex " << i << std::endl;
                    if ( m_crash_on_problem ) exit(1);
                }
                if ( std::abs(hepevt_particles[i]->pid()) <= 6 && !hepevt_particles[i]->end_vertex() ) {
                    std::cerr << "quark without end vertex " << i << std::endl;
                    if ( m_crash_on_problem ) exit(1);
                }
            }
        }


        // 5. Store PDF, weight, cross section and other event information.
        // Flavours of incoming partons.
        if (m_store_pdf && pyinfo != 0) {
            int id1pdf = pyinfo->id1pdf();
            int id2pdf = pyinfo->id2pdf();
            if ( m_convert_gluon_to_0 ) {
                if (id1pdf == 21) id1pdf = 0;
                if (id2pdf == 21) id2pdf = 0;
            }

            GenPdfInfoPtr pdfinfo = std::make_shared<GenPdfInfo>();
            pdfinfo->set(id1pdf, id2pdf, pyinfo->x1pdf(), pyinfo->x2pdf(), pyinfo->QFac(), pyinfo->pdf1(), pyinfo->pdf2() );
            // Store PDF information.
            evt->set_pdf_info( pdfinfo );
        }

        // Store process code, scale, alpha_em, alpha_s.
        if (m_store_proc && pyinfo != 0) {
            evt->add_attribute("mpi",std::make_shared<IntAttribute>(pyinfo->nMPI()));
            evt->add_attribute("signal_process_id",std::make_shared<IntAttribute>( pyinfo->code()));
            evt->add_attribute("event_scale",std::make_shared<DoubleAttribute>(pyinfo->QRen()));
            evt->add_attribute("alphaQCD",std::make_shared<DoubleAttribute>(pyinfo->alphaS()));
            evt->add_attribute("alphaQED",std::make_shared<DoubleAttribute>(pyinfo->alphaEM()));
        }

        // Store cross-section information in pb.
        if (m_store_xsec && pyinfo != 0) {
            GenCrossSectionPtr xsec = std::make_shared<GenCrossSection>();
            xsec->set_cross_section( pyinfo->sigmaGen() * 1e9, pyinfo->sigmaErr() * 1e9);
            evt->set_cross_section(xsec);
        }

        // Store event weights.
        if (m_store_weights && pyinfo != 0) {
            evt->weights().clear();
            for (int iweight=0; iweight < pyinfo->nWeights(); ++iweight) {
                evt->weights().push_back(pyinfo->weight(iweight));
            }
        }

        // Done.
        return true;
    }

    // Read out values for some switches
    bool print_inconsistency()  const {
        return m_print_inconsistency;
    }
    bool free_parton_warnings() const {
        return m_free_parton_warnings;
    }
    bool crash_on_problem()     const {
        return m_crash_on_problem;
    }
    bool convert_gluon_to_0()   const {
        return m_convert_gluon_to_0;
    }
    bool store_pdf()            const {
        return m_store_pdf;
    }
    bool store_proc()           const {
        return m_store_proc;
    }
    bool store_xsec()           const {
        return m_store_xsec;
    }
    bool store_weights()        const {
        return m_store_weights;
    }

    // Set values for some switches
    void set_print_inconsistency(bool b = true)  {
        m_print_inconsistency  = b;
    }
    void set_free_parton_warnings(bool b = true) {
        m_free_parton_warnings = b;
    }
    void set_crash_on_problem(bool b = false)    {
        m_crash_on_problem     = b;
    }
    void set_convert_gluon_to_0(bool b = false)  {
        m_convert_gluon_to_0   = b;
    }
    void set_store_pdf(bool b = true)            {
        m_store_pdf            = b;
    }
    void set_store_proc(bool b = true)           {
        m_store_proc           = b;
    }
    void set_store_xsec(bool b = true)           {
        m_store_xsec           = b;
    }
    void set_store_weights(bool b = true)        {
        m_store_weights        = b;
    }

private:

    // Following methods are not implemented for this class
    virtual bool fill_next_event( GenEvent*  )  {
        return 0;
    }
    virtual void write_event( const GenEvent* ) {}

    // Use of copy constructor is not allowed
    Pythia8ToHepMC3( const Pythia8ToHepMC3& ) {}

    // Data members
    int  m_internal_event_number;
    bool m_print_inconsistency;
    bool m_free_parton_warnings;
    bool m_crash_on_problem;
    bool m_convert_gluon_to_0;
    bool m_store_pdf;
    bool m_store_proc;
    bool m_store_xsec;
    bool m_store_weights;
};

} // namespace HepMC3
#endif

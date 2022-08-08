// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file ReaderLHEF.cc
 *  @brief Implementation of \b class ReaderLHEF
 *
 */
#include "HepMC3/ReaderLHEF.h"
namespace HepMC3
{
ReaderLHEF::ReaderLHEF(const std::string& filename)
{
    m_reader = std::make_shared<LHEF::Reader>(filename);
    init();
}
ReaderLHEF::ReaderLHEF(std::istream & stream)
{
    m_reader = std::make_shared<LHEF::Reader>(stream);
    init();
}

ReaderLHEF::ReaderLHEF(std::shared_ptr<std::istream> s_stream)
    : m_shared_stream(s_stream)
{
    m_reader = std::make_shared<LHEF::Reader>(*(m_shared_stream.get()));
    init();
}

bool ReaderLHEF::skip(const int n)
{
    GenEvent evt;
    for (int nn = n; nn > 0; --nn)
    {
        if (!read_event(evt)) return false;
        evt.clear();
    }
    return !failed();
}


void ReaderLHEF::init()
{
    m_neve = 0;
    m_failed = false;
    // Create a HEPRUP attribute and initialize it from the reader.
    m_hepr = std::make_shared<HEPRUPAttribute>();
    m_hepr->heprup = m_reader->heprup;
    // There may be some XML tags in the LHE file which are
    // non-standard, but we can save them as well.
    m_hepr->tags = LHEF::XMLTag::findXMLTags(m_reader->headerBlock + m_reader->initComments);
    // This code is ugly and should be replaced.
    size_t nweights = 0;
    for (auto t1: m_hepr->tags) {
        if (t1->name != "header") continue;
        for (auto t2: t1->tags) {
            if (t2->name != "initrwgt") continue;
            for (auto t3: t2->tags) {
                if (t3->name != "weightgroup") continue;
                for (auto t4: t3->tags) if (t4->name == "weight") nweights++;
                break;
            }
            break;
        }
        break;
    }
    //
    // Now we want to create a GenRunInfo object for the HepMC file, and
    // we add the LHEF attribute to that.
    set_run_info(std::make_shared<GenRunInfo>());
    run_info()->add_attribute("HEPRUP", m_hepr);

    // This is just a test to make sure we can add other attributes as
    // well.
    run_info()->add_attribute("NPRUP",
                              std::make_shared<FloatAttribute>(m_hepr->heprup.NPRUP));

    // We want to be able to convey the different event weights to
    // HepMC. In particular we need to add the names of the weights to
    // the GenRunInfo object.

    std::vector<std::string> weightnames;
    for ( int i = 0, N = m_hepr->heprup.weightinfo.size(); i < N; ++i ) weightnames.push_back(m_hepr->heprup.weightNameHepMC(i));
    if (nweights == 0) nweights=1;
    for ( size_t i = weightnames.size(); i < nweights; ++i ) weightnames.push_back(std::to_string(i));
    run_info()->set_weight_names(weightnames);

    // We also want to convey the information about which generators was
    // used.
    for ( int i = 0, N = m_hepr->heprup.generators.size(); i < N; ++i )
    {
        GenRunInfo::ToolInfo tool;
        tool.name =  m_hepr->heprup.generators[i].name;
        tool.version =  m_hepr->heprup.generators[i].version;
        tool.description =  m_hepr->heprup.generators[i].contents;
        run_info()->tools().push_back(tool);
    }
}
/// @brief Destructor
ReaderLHEF::~ReaderLHEF() {close();}

bool ReaderLHEF::read_event(GenEvent& ev)
{
    if (m_storage.size() > 0)
    {
        ev = m_storage.front();
        m_storage.pop_front();
        return m_failed;
    }
    m_failed = !(m_reader->readEvent());
    if (m_failed) return m_failed;
    // To each GenEvent we want to add an attribute corresponding to
    // the HEPEUP. Also here there may be additional non-standard
    // information outside the LHEF <event> tags, which we may want to
    // add.
    std::shared_ptr<HEPEUPAttribute> hepe = std::make_shared<HEPEUPAttribute>();
    if ( m_reader->outsideBlock.length() )
        hepe->tags =  LHEF::XMLTag::findXMLTags(m_reader->outsideBlock);

    hepe->hepeup = m_reader->hepeup;
    std::vector<LHEF::HEPEUP*> input;
    if (m_reader->hepeup.subevents.size() > 0) input.insert(input.end(), hepe->hepeup.subevents.begin(), hepe->hepeup.subevents.end());
    else { input.push_back(&m_reader->hepeup);}
    int first_group_event = m_neve;
    m_neve++;
    for (auto ahepeup: input)
    {
        GenEvent evt;
        evt.set_event_number(first_group_event);
        evt.add_attribute("AlphaQCD", std::make_shared<DoubleAttribute>(ahepeup->AQCDUP));
        evt.add_attribute("AlphaEM", std::make_shared<DoubleAttribute>(ahepeup->AQEDUP));
        evt.add_attribute("NUP",  std::make_shared<IntAttribute>(ahepeup->NUP));
        evt.add_attribute("IDPRUP", std::make_shared<LongAttribute>(ahepeup->IDPRUP));
        // Now add the Particles from the LHE event to HepMC
        std::vector<GenParticlePtr> particles;
        std::map< std::pair<int, int>, GenVertexPtr> vertices;
        for ( int i = 0; i < ahepeup->NUP; ++i )
        {
            FourVector mom((ahepeup->PUP)[i][0], (ahepeup->PUP)[i][1], (ahepeup->PUP)[i][2], (ahepeup->PUP)[i][3]);
            particles.push_back(std::make_shared<GenParticle>(mom, ahepeup->IDUP[i], ahepeup->ISTUP[i]));
            if ( i < 2 ) continue;
            std::pair<int, int> vertex_index(ahepeup->MOTHUP[i].first, ahepeup->MOTHUP[i].second);
            if (vertices.find(vertex_index) == vertices.end()) vertices[vertex_index] = std::make_shared<GenVertex>();
            vertices[vertex_index]->add_particle_out(particles.back());
        }
        for ( auto v: vertices )
        {
            std::pair<int, int> vertex_index = v.first;
            GenVertexPtr          vertex = v.second;
            for (int i = vertex_index.first-1; i < vertex_index.second; ++i)
                if ( i >= 0 && i < (int)particles.size())
                    vertex->add_particle_in(particles[i]);
        }
        std::pair<int, int> vertex_index(0, 0);
        if (vertices.find(vertex_index) == vertices.end()) vertices[vertex_index] = std::make_shared<GenVertex>();
        for (size_t i = 0; i < particles.size(); ++i)
            if (!particles[i]->end_vertex() && !particles[i]->production_vertex())
            {
                if ( i < 2 ) vertices[vertex_index]->add_particle_in(particles[i]);
                else vertices[vertex_index]->add_particle_out(particles[i]);
            }
        for ( auto v: vertices ) evt.add_vertex(v.second);
        if (particles.size() > 1)
        {
            particles[0]->set_status(4);
            particles[1]->set_status(4);
            evt.set_beam_particles(particles[0], particles[1]);
        }
        // And we also want to add the weights.


        std::vector<double> wts;
        for ( int i = 0, N = ahepeup->weights.size(); i < N; ++i )
        {
            wts.push_back(ahepeup->weights[i].first);
        }
        evt.weights() = wts;
        m_storage.push_back(evt);
    }
    ev = m_storage.front();
    m_storage.pop_front();
    return m_failed;
}
/// @brief Return status of the stream
bool ReaderLHEF::failed() { return m_failed;}

/// @brief Close file stream
void ReaderLHEF::close() { }
} // namespace HepMC3

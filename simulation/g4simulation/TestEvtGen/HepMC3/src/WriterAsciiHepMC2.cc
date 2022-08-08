// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file WriterAsciiHepMC2.cc
/// @brief Implementation of \b class WriterAsciiHepMC2
///
#include <cstring>

#include "HepMC3/WriterAsciiHepMC2.h"

#include "HepMC3/Version.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"

namespace HepMC3
{


WriterAsciiHepMC2::WriterAsciiHepMC2(const std::string &filename, std::shared_ptr<GenRunInfo> run)
    : m_file(filename),
      m_stream(&m_file),
      m_precision(16),
      m_buffer(nullptr),
      m_cursor(nullptr),
      m_buffer_size(256*1024),
      m_particle_counter(0)
{
    HEPMC3_WARNING("WriterAsciiHepMC2::WriterAsciiHepMC2: HepMC2 IO_GenEvent format is outdated. Please use HepMC3 Asciiv3 format instead.")
    set_run_info(run);
    if ( !run_info() ) set_run_info(std::make_shared<GenRunInfo>());
    if ( !m_file.is_open() )
    {
        HEPMC3_ERROR("WriterAsciiHepMC2: could not open output file: " << filename )
    }
    else
    {
        const std::string header = "HepMC::Version " + version() + "\nHepMC::IO_GenEvent-START_EVENT_LISTING\n";
        m_file.write(header.data(), header.length());
    }
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
}

WriterAsciiHepMC2::WriterAsciiHepMC2(std::ostream &stream, std::shared_ptr<GenRunInfo> run)
    : m_file(),
      m_stream(&stream),
      m_precision(16),
      m_buffer(nullptr),
      m_cursor(nullptr),
      m_buffer_size(256*1024),
      m_particle_counter(0)
{
    HEPMC3_WARNING("WriterAsciiHepMC2::WriterAsciiHepMC2: HepMC2 IO_GenEvent format is outdated. Please use HepMC3 Asciiv3 format instead.")
    set_run_info(run);
    if ( !run_info() ) set_run_info(std::make_shared<GenRunInfo>());
    const std::string header = "HepMC::Version " + version() + "\nHepMC::IO_GenEvent-START_EVENT_LISTING\n";
    m_stream->write(header.data(), header.length());
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
}

WriterAsciiHepMC2::WriterAsciiHepMC2(std::shared_ptr<std::ostream> s_stream, std::shared_ptr<GenRunInfo> run)
    : m_file(),
      m_shared_stream(s_stream),
      m_stream(s_stream.get()),
      m_precision(16),
      m_buffer(nullptr),
      m_cursor(nullptr),
      m_buffer_size(256*1024),
      m_particle_counter(0)
{
    HEPMC3_WARNING("WriterAsciiHepMC2::WriterAsciiHepMC2: HepMC2 IO_GenEvent format is outdated. Please use HepMC3 Asciiv3 format instead.")
    set_run_info(run);
    if ( !run_info() ) set_run_info(std::make_shared<GenRunInfo>());
    const std::string header = "HepMC::Version " + version() + "\nHepMC::IO_GenEvent-START_EVENT_LISTING\n";
    m_stream->write(header.data(), header.length());
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
}


WriterAsciiHepMC2::~WriterAsciiHepMC2()
{
    close();
    if ( m_buffer ) delete[] m_buffer;
}


void WriterAsciiHepMC2::write_event(const GenEvent &evt)
{
    allocate_buffer();
    if ( !m_buffer ) return;
    auto float_printf_specifier_option = m_options.find("float_printf_specifier");
    std::string  letter=(float_printf_specifier_option != m_options.end())?float_printf_specifier_option->second.substr(0,2):"e";
    if (letter != "e" && letter != "E" && letter != "G" && letter != "g" && letter != "f" && letter != "F" ) letter = "e";
    m_float_printf_specifier = " %." + std::to_string(m_precision) + letter;
    // Make sure nothing was left from previous event
    flush();

    if ( !run_info() ) set_run_info(evt.run_info());
    if ( evt.run_info() && run_info() != evt.run_info() ) set_run_info(evt.run_info());


    std::shared_ptr<DoubleAttribute> A_event_scale = evt.attribute<DoubleAttribute>("event_scale");
    std::shared_ptr<DoubleAttribute> A_alphaQED = evt.attribute<DoubleAttribute>("alphaQED");
    std::shared_ptr<DoubleAttribute> A_alphaQCD = evt.attribute<DoubleAttribute>("alphaQCD");
    std::shared_ptr<IntAttribute> A_signal_process_id = evt.attribute<IntAttribute>("signal_process_id");
    std::shared_ptr<IntAttribute> A_mpi = evt.attribute<IntAttribute>("mpi");
    std::shared_ptr<IntAttribute> A_signal_process_vertex = evt.attribute<IntAttribute>("signal_process_vertex");

    double event_scale = A_event_scale?(A_event_scale->value()):0.0;
    double alphaQED = A_alphaQED?(A_alphaQED->value()):0.0;
    double alphaQCD = A_alphaQCD?(A_alphaQCD->value()):0.0;
    int signal_process_id = A_signal_process_id?(A_signal_process_id->value()):0;
    int mpi = A_mpi?(A_mpi->value()):0;
    int signal_process_vertex = A_signal_process_vertex?(A_signal_process_vertex->value()):0;

    std::vector<long> m_random_states;
    std::shared_ptr<VectorLongIntAttribute> random_states_a = evt.attribute<VectorLongIntAttribute>("random_states");
    if (random_states_a) {
        m_random_states = random_states_a->value();
    } else {
        m_random_states.reserve(100);
        for (int i = 0; i < 100; i++)
        {
            std::shared_ptr<LongAttribute> rs = evt.attribute<LongAttribute>("random_states"+std::to_string((long long unsigned int)i));
            if (!rs) break;
            m_random_states.push_back(rs->value());
        }
    }
    // Write event info
    //Find beam particles
    std::vector<int> beams;
    beams.reserve(2);
    int idbeam = 0;
    for (ConstGenVertexPtr v: evt.vertices())
    {
        for (ConstGenParticlePtr p: v->particles_in())
        {
            if (!p->production_vertex())                { if (p->status() == 4) beams.push_back(idbeam); idbeam++; }
            else if (p->production_vertex()->id() == 0) { if (p->status() == 4) beams.push_back(idbeam); idbeam++; }
        }
        for (ConstGenParticlePtr p: v->particles_out()) { if (p->status() == 4) beams.push_back(idbeam); idbeam++; }
    }
    //
    int idbeam1 = 10000;
    int idbeam2 = 10000;
    if (beams.size() > 0) idbeam1 += beams[0] + 1;
    if (beams.size() > 1) idbeam2 += beams[1] + 1;
    m_cursor += sprintf(m_cursor, "E %d %d %e %e %e %d %d %zu %i %i",
                        evt.event_number(),
                        mpi,
                        event_scale,
                        alphaQCD,
                        alphaQED,
                        signal_process_id,
                        signal_process_vertex,
                        evt.vertices().size(),
                        idbeam1, idbeam2);

    // This should be the largest single add to the buffer. Its size 11+4*11+3*22+2*11+10=153
    flush();
    m_cursor += sprintf(m_cursor, " %zu", m_random_states.size());
    for (size_t q = 0; q < m_random_states.size(); q++)
    {
        m_cursor += sprintf(m_cursor, " %i", (int)q);
        flush();
    }
    flush();
    m_cursor += sprintf(m_cursor, " %zu", evt.weights().size());
    if ( evt.weights().size() )
    {
        for (double w: evt.weights()) {
            m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), w);
            flush();
        }
        m_cursor += sprintf(m_cursor, "\n");
        flush();
        m_cursor += sprintf(m_cursor, "N %zu", evt.weights().size());
        const std::vector<std::string> names = run_info()->weight_names();
        for (size_t q = 0; q < evt.weights().size(); q++)
        {
            if (q < names.size())
                write_string(" \""+names[q]+"\"");
            else
                write_string(" \""+std::to_string(q)+"\"");
            flush();
        }
    }
    m_cursor += sprintf(m_cursor, "\n");
    flush();
    // Write units
    m_cursor += sprintf(m_cursor, "U %s %s\n", Units::name(evt.momentum_unit()).c_str(), Units::name(evt.length_unit()).c_str());
    flush();
    std::shared_ptr<GenCrossSection> cs = evt.attribute<GenCrossSection>("GenCrossSection");
    if (cs) {
        m_cursor += sprintf(m_cursor, ("C" + m_float_printf_specifier + m_float_printf_specifier + "\n").c_str(), cs->xsec(), cs->xsec_err() );
        flush();
    }

    std::shared_ptr<GenHeavyIon> hi = evt.attribute<GenHeavyIon>("GenHeavyIon");
    if (hi) {
        m_cursor += sprintf(m_cursor, "H %i %i %i %i %i %i %i %i %i %e %e %e %e\n",
                            hi->Ncoll_hard,
                            hi->Npart_proj,
                            hi->Npart_targ,
                            hi->Ncoll,
                            hi->spectator_neutrons,
                            hi->spectator_protons,
                            hi->N_Nwounded_collisions,
                            hi->Nwounded_N_collisions,
                            hi->Nwounded_Nwounded_collisions,
                            hi->impact_parameter,
                            hi->event_plane_angle,
                            hi->eccentricity,
                            hi->sigma_inel_NN);
        flush();
    }

    std::shared_ptr<GenPdfInfo> pi = evt.attribute<GenPdfInfo>("GenPdfInfo");
    if (pi) {
        std::string st;
        // We use it here because the HepMC3 GenPdfInfo has the same format as in HepMC2 IO_GenEvent and get error handeling for free.
        bool status = pi->to_string(st);
        if ( !status )
        {
            HEPMC3_WARNING("WriterAsciiHepMC2::write_event: problem serializing GenPdfInfo attribute")
        } else {
            m_cursor += sprintf(m_cursor, "F ");
            flush();
            write_string(escape(st));
            m_cursor += sprintf(m_cursor, "\n");
            flush();
        }
    }


    m_particle_counter = 0;
    for (ConstGenVertexPtr v: evt.vertices() )
    {
        int production_vertex = 0;
        production_vertex = v->id();
        write_vertex(v);
        for (ConstGenParticlePtr p: v->particles_in())
        {
            if (!p->production_vertex()) write_particle( p, production_vertex );
            else
            {
                if (p->production_vertex()->id() == 0) write_particle( p, production_vertex );
            }
        }
        for (ConstGenParticlePtr p: v->particles_out())
            write_particle(p, production_vertex);
    }

    // Flush rest of the buffer to file
    forced_flush();
}


void WriterAsciiHepMC2::allocate_buffer()
{
    if ( m_buffer ) return;
    while ( m_buffer == nullptr && m_buffer_size >= 512 ) {
        try {
            m_buffer = new char[ m_buffer_size ]();
        } catch (const std::bad_alloc& e) {
            delete[] m_buffer;
            m_buffer_size /= 2;
            HEPMC3_WARNING("WriterAsciiHepMC2::allocate_buffer:" << e.what() << " buffer size too large. Dividing by 2. New size: " << m_buffer_size)
        }
    }

    if ( !m_buffer )
    {
        HEPMC3_ERROR("WriterAsciiHepMC2::allocate_buffer: could not allocate buffer!")
        return;
    }

    m_cursor = m_buffer;
}


std::string WriterAsciiHepMC2::escape(const std::string& s) const
{
    std::string ret;
    ret.reserve(s.length()*2);
    for ( std::string::const_iterator it = s.begin(); it != s.end(); ++it )
    {
        switch ( *it )
        {
        case '\\':
            ret += "\\\\";
            break;
        case '\n':
            ret += "\\|";
            break;
        default:
            ret += *it;
        }
    }
    return ret;
}

void WriterAsciiHepMC2::write_vertex(ConstGenVertexPtr v)
{
    std::vector<double> weights;
    std::shared_ptr<VectorDoubleAttribute> weights_a = v->attribute<VectorDoubleAttribute>("weights");
    if (weights_a) {
        weights = weights_a->value();
    } else {
        weights.reserve(100);
        for (int i = 0; i < 100; i++)
        {
            std::shared_ptr<DoubleAttribute> rs = v->attribute<DoubleAttribute>("weight"+std::to_string((long long unsigned int)i));
            if (!rs) break;
            weights.push_back(rs->value());
        }
    }
    flush();
    m_cursor += sprintf(m_cursor, "V %i %i", v->id(), v->status());
    int orph = 0;
    for (ConstGenParticlePtr p: v->particles_in())
    {
        if (!p->production_vertex()) orph++;
        else
        {
            if (p->production_vertex()->id() == 0) orph++;
        }
    }
    const FourVector &pos = v->position();
    if (pos.is_zero())
    {
        m_cursor += sprintf(m_cursor, " 0 0 0 0");
    }
    else
    {
        m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), pos.x());
        m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), pos.y());
        m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), pos.z());
        m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), pos.t());
    }
    m_cursor += sprintf(m_cursor, " %i %zu %zu", orph, v->particles_out().size(), weights.size());
    flush();
    for (size_t i = 0; i < weights.size(); i++) { m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), weights[i]); flush(); }
    m_cursor += sprintf(m_cursor, "\n");
    flush();
}


inline void WriterAsciiHepMC2::flush()
{
    // The maximum size of single add to the buffer should not be larger than 256. This is a safe value as
    // we will not allow precision larger than 24 anyway
    if ( m_buffer + m_buffer_size < m_cursor + 512 )
    {
        std::ptrdiff_t length = m_cursor - m_buffer;
        m_stream->write(m_buffer, length);
        m_cursor = m_buffer;
    }
}


inline void WriterAsciiHepMC2::forced_flush()
{
    std::ptrdiff_t length = m_cursor - m_buffer;
    m_stream->write(m_buffer, length);
    m_cursor = m_buffer;
}


void WriterAsciiHepMC2::write_run_info() {}

void WriterAsciiHepMC2::write_particle(ConstGenParticlePtr p, int /*second_field*/)
{
    flush();
    m_cursor += sprintf(m_cursor, "P %i", int(10001+m_particle_counter));
    m_particle_counter++;
    m_cursor += sprintf(m_cursor, " %i", p->pid() );
    m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), p->momentum().px() );
    m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), p->momentum().py());
    m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), p->momentum().pz() );
    m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), p->momentum().e() );
    m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), p->generated_mass() );
    m_cursor += sprintf(m_cursor, " %i", p->status() );
    flush();
    int ev = 0;
    if (p->end_vertex())
        if (p->end_vertex()->id() != 0)
            ev = p->end_vertex()->id();

    std::shared_ptr<DoubleAttribute> A_theta = p->attribute<DoubleAttribute>("theta");
    std:: shared_ptr<DoubleAttribute> A_phi = p->attribute<DoubleAttribute>("phi");
    if (A_theta) m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), A_theta->value());
    else m_cursor += sprintf(m_cursor, " 0");
    if (A_phi) m_cursor += sprintf(m_cursor, m_float_printf_specifier.c_str(), A_phi->value());
    else m_cursor += sprintf(m_cursor, " 0");
    m_cursor += sprintf(m_cursor, " %i", ev);
    flush();
    std::shared_ptr<VectorIntAttribute> A_flows = p->attribute<VectorIntAttribute>("flows");
    if (A_flows)
    {
        std::vector<int> flowsv = A_flows->value();
        std::string flowss = " " + std::to_string(flowsv.size());
        for (size_t k = 0; k < flowsv.size(); k++) { flowss += ( " " + std::to_string(k+1) + " " + std::to_string(flowsv.at(k))); }
        flowss += "\n";
        write_string(flowss);
    } else {
        std::shared_ptr<IntAttribute> A_flow1 = p->attribute<IntAttribute>("flow1");
        std::shared_ptr<IntAttribute> A_flow2 = p->attribute<IntAttribute>("flow2");
        std::shared_ptr<IntAttribute> A_flow3 = p->attribute<IntAttribute>("flow3");
        int flowsize = 0;
        if (A_flow1) flowsize++;
        if (A_flow2) flowsize++;
        if (A_flow3) flowsize++;
        std::string flowss = " " + std::to_string(flowsize);
        if (A_flow1) flowss += ( " 1 " + std::to_string(A_flow1->value()));
        if (A_flow2) flowss += ( " 2 " + std::to_string(A_flow2->value()));
        if (A_flow3) flowss += ( " 3 " + std::to_string(A_flow3->value()));
        flowss += "\n";
        write_string(flowss);
    }
    flush();
}


inline void WriterAsciiHepMC2::write_string(const std::string &str)
{
    // First let's check if string will fit into the buffer
    if ( m_buffer + m_buffer_size > m_cursor + str.length() )
    {
        strncpy(m_cursor, str.data(), str.length());
        m_cursor += str.length();
        flush();
    }
    // If not, flush the buffer and write the string directly
    else
    {
        forced_flush();
        m_stream->write(str.data(), str.length());
    }
}


void WriterAsciiHepMC2::close()
{
    std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
    if (ofs && !ofs->is_open()) return;
    forced_flush();
    const std::string footer("HepMC::IO_GenEvent-END_EVENT_LISTING\n\n");
    if (m_stream) m_stream->write(footer.data(),footer.length());
    if (ofs) ofs->close();
}
bool WriterAsciiHepMC2::failed() { return (bool)m_file.rdstate(); }

void WriterAsciiHepMC2::set_precision(const int& prec ) {
    if (prec < 2 || prec > 24) return;
    m_precision = prec;
}

int WriterAsciiHepMC2::precision() const {
    return m_precision;
}

void WriterAsciiHepMC2::set_buffer_size(const size_t& size ) {
    if (m_buffer) return;
    if (size < 1024) return;
    m_buffer_size = size;
}
} // namespace HepMC3

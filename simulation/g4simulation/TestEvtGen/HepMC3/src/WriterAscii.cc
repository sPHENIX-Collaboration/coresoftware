// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file WriterAscii.cc
/// @brief Implementation of \b class WriterAscii
///
#include <cstring>
#include <algorithm>//min max for VS2017

#include "HepMC3/WriterAscii.h"

#include "HepMC3/Version.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"

namespace HepMC3 {


WriterAscii::WriterAscii(const std::string &filename, std::shared_ptr<GenRunInfo> run)
    : m_file(filename),
      m_stream(&m_file),
      m_precision(16),
      m_buffer(nullptr),
      m_cursor(nullptr),
      m_buffer_size(256*1024)
{
    set_run_info(run);
    if ( !m_file.is_open() ) {
        HEPMC3_ERROR("WriterAscii: could not open output file: " << filename)
    } else {
        const std::string header = "HepMC::Version " + version() + "\nHepMC::Asciiv3-START_EVENT_LISTING\n";
        m_file.write(header.data(), header.length());
        if ( run_info() ) write_run_info();
    }
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";
}


WriterAscii::WriterAscii(std::ostream &stream, std::shared_ptr<GenRunInfo> run)
    : m_file(),
      m_stream(&stream),
      m_precision(16),
      m_buffer(nullptr),
      m_cursor(nullptr),
      m_buffer_size(256*1024)
{
    set_run_info(run);
    const std::string header = "HepMC::Version " + version() + "\nHepMC::Asciiv3-START_EVENT_LISTING\n";
    m_stream->write(header.data(), header.length());
    if ( run_info() ) write_run_info();
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";
}

WriterAscii::WriterAscii(std::shared_ptr<std::ostream> s_stream, std::shared_ptr<GenRunInfo> run)
    : m_file(),
      m_shared_stream(s_stream),
      m_stream(s_stream.get()),
      m_precision(16),
      m_buffer(nullptr),
      m_cursor(nullptr),
      m_buffer_size(256*1024)
{
    set_run_info(run);
    const std::string header = "HepMC::Version " + version() + "\nHepMC::Asciiv3-START_EVENT_LISTING\n";
    m_stream->write(header.data(), header.length());
    if ( run_info() ) write_run_info();
    m_float_printf_specifier = " %." + std::to_string(m_precision) + "e";
    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";
}

WriterAscii::~WriterAscii() {
    close();
    if ( m_buffer ) delete[] m_buffer;
}


void WriterAscii::write_event(const GenEvent &evt) {
    allocate_buffer();
    if ( !m_buffer ) return;
    auto float_printf_specifier_option = m_options.find("float_printf_specifier");
    std::string  letter=(float_printf_specifier_option != m_options.end())?float_printf_specifier_option->second.substr(0,2):"e";
    if (letter != "e" && letter != "E" && letter != "G" && letter != "g" && letter != "f" && letter != "F" ) letter = "e";
    m_float_printf_specifier = " %." + std::to_string(m_precision) + letter;


    m_particle_printf_specifier = "P %i %i %i"
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier
                                  + m_float_printf_specifier + " %i\n";
    m_vertex_short_printf_specifier = "V %i %i [%s]\n";
    m_vertex_long_printf_specifier = "V %i %i [%s] @"+ m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n";

    // Make sure nothing was left from previous event
    flush();

    if ( !run_info() ) {
        set_run_info(evt.run_info());
        write_run_info();
    } else {
        if ( evt.run_info() && (run_info() != evt.run_info()) ) {
            HEPMC3_WARNING("WriterAscii::write_event: GenEvents contain "
                           "different GenRunInfo objects from - only the "
                           "first such object will be serialized.")
        }
    }

    // Write event info
    flush();
    std::string especifier =  "E " + std::to_string(evt.event_number()) + " "
                              + std::to_string(evt.vertices().size()) + " "
                              + std::to_string(evt.particles().size());
    // Write event position if not zero
    const FourVector &pos = evt.event_pos();
    if ( !pos.is_zero() ) {
        especifier += ( " @"  + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + m_float_printf_specifier + "\n" );
        m_cursor += sprintf(m_cursor, especifier.c_str(), pos.x(), pos.y(), pos.z(), pos.t());
    } else {
        m_cursor += sprintf(m_cursor, "%s\n", especifier.c_str());
    }
    flush();

    // Write units
    m_cursor += sprintf(m_cursor, "U %s %s\n", Units::name(evt.momentum_unit()).c_str(), Units::name(evt.length_unit()).c_str());
    flush();

    // Write weight values if present
    if ( evt.weights().size() ) {
        m_cursor += sprintf(m_cursor, "W");
        for (auto w: evt.weights())
        {
            m_cursor += sprintf(m_cursor, " %.*e", std::min(3*m_precision, 22), w);
            flush();
        }
        m_cursor += sprintf(m_cursor, "\n");
        flush();
    }

    // Write attributes
    for ( auto vt1: evt.attributes() ) {
        for ( auto vt2: vt1.second ) {
            std::string st;
            bool status = vt2.second->to_string(st);

            if ( !status ) {
                HEPMC3_WARNING("WriterAscii::write_event: problem serializing attribute: " << vt1.first)
            }
            else {
                m_cursor += sprintf(m_cursor, "A %i ", vt2.first);
                write_string(escape(vt1.first));
                flush();
                m_cursor += sprintf(m_cursor, " ");
                write_string(escape(st));
                m_cursor += sprintf(m_cursor, "\n");
                flush();
            }
        }
    }


    // Print particles
    std::map<int, bool> alreadywritten;
    for (ConstGenParticlePtr p: evt.particles()) {
        // Check to see if we need to write a vertex first
        ConstGenVertexPtr v = p->production_vertex();
        int parent_object = 0;

        if (v) {
            // Check if we need this vertex at all
            // Yes, use vertex as parent object
            if ( v->particles_in().size() > 1 || !v->data().is_zero() ) parent_object = v->id();
            // No, use particle as parent object
            // Add check for attributes of this vertex
            else if ( v->particles_in().size() == 1 )                   parent_object = v->particles_in().front()->id();
            else if ( v->particles_in().size() == 0 ) HEPMC3_DEBUG(30, "WriterAscii::write_event - found a vertex without incoming particles: " << v->id());
            // Usage of map instead of simple counter helps to deal with events with random ids of vertices.
            if (alreadywritten.count(v->id()) == 0 && parent_object < 0)
            { write_vertex(v); alreadywritten[v->id()] = true; }
        }

        write_particle(p, parent_object);
    }
    alreadywritten.clear();

    // Flush rest of the buffer to file
    forced_flush();
}


void WriterAscii::allocate_buffer() {
    if ( m_buffer ) return;
    while ( m_buffer == nullptr && m_buffer_size >= 512 ) {
        try {
            m_buffer = new char[ m_buffer_size ]();
        } catch (const std::bad_alloc& e) {
            delete[] m_buffer;
            m_buffer_size /= 2;
            HEPMC3_WARNING("WriterAscii::allocate_buffer:" << e.what() << " buffer size too large. Dividing by 2. New size: " << m_buffer_size)
        }
    }

    if ( !m_buffer ) {
        HEPMC3_ERROR("WriterAscii::allocate_buffer: could not allocate buffer!")
        return;
    }
    m_cursor = m_buffer;
}


std::string WriterAscii::escape(const std::string& s) const {
    std::string ret;
    ret.reserve(s.length()*2);
    for ( std::string::const_iterator it = s.begin(); it != s.end(); ++it ) {
        switch ( *it ) {
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

void WriterAscii::write_vertex(ConstGenVertexPtr v) {
    flush();
    std::string vlist;
    std::vector<int> pids;
    pids.reserve(v->particles_in().size());
    for (ConstGenParticlePtr p: v->particles_in()) pids.push_back(p->id());
    //We order pids to be able to compare ascii files
    std::sort(pids.begin(), pids.end());
    for (auto p: pids) vlist.append( std::to_string(p).append(",") );
    if ( pids.size() ) vlist.pop_back();
    const FourVector &pos = v->position();
    if ( !pos.is_zero() ) {
        m_cursor += sprintf(m_cursor, m_vertex_long_printf_specifier.c_str(),  v->id(), v->status(), vlist.c_str(), pos.x(), pos.y(), pos.z(), pos.t() );
    } else {
        m_cursor += sprintf(m_cursor, m_vertex_short_printf_specifier.c_str(), v->id(), v->status(), vlist.c_str());
    }
    flush();
}


inline void WriterAscii::flush() {
    // The maximum size of single add to the buffer (other than by
    // using WriterAscii::write_string) should not be larger than 256. This is a safe value as
    // we will not allow precision larger than 24 anyway
    if ( m_buffer + m_buffer_size < m_cursor + 512 ) {
        std::ptrdiff_t length = m_cursor - m_buffer;
        m_stream->write(m_buffer, length);
        m_cursor = m_buffer;
    }
}


inline void WriterAscii::forced_flush() {
    std::ptrdiff_t length = m_cursor - m_buffer;
    m_stream->write(m_buffer, length);
    m_cursor = m_buffer;
}


void WriterAscii::write_run_info() {
    allocate_buffer();

    // If no run info object set, create a dummy one.
    if ( !run_info() ) set_run_info(std::make_shared<GenRunInfo>());

    const std::vector<std::string> names = run_info()->weight_names();

    if ( !names.empty() ) {
        std::string out = names[0];
        for ( int i = 1, N = names.size(); i < N; ++i )
            out += "\n" + names[i];
        m_cursor += sprintf(m_cursor, "W ");
        flush();
        write_string(escape(out));
        m_cursor += sprintf(m_cursor, "\n");
    }

    for (int i = 0, N = run_info()->tools().size(); i < N; ++i) {
        std::string out = "T " + run_info()->tools()[i].name + "\n"
                          + run_info()->tools()[i].version + "\n"
                          + run_info()->tools()[i].description;
        write_string(escape(out));
        m_cursor += sprintf(m_cursor, "\n");
    }


    for ( auto att: run_info()->attributes() ) {
        std::string st;
        if ( !att.second->to_string(st) ) {
            HEPMC3_WARNING("WriterAscii::write_run_info: problem serializing attribute: " << att.first)
        }
        else {
            m_cursor += sprintf(m_cursor, "A ");
            write_string(att.first);
            flush();
            m_cursor += sprintf(m_cursor, " ");
            write_string(escape(st));
            m_cursor += sprintf(m_cursor, "\n");
            flush();
        }
    }
}

void WriterAscii::write_particle(ConstGenParticlePtr p, int second_field) {
    flush();
    m_cursor += sprintf(m_cursor, m_particle_printf_specifier.c_str(), p->id(), second_field, p->pid(), p->momentum().px(), p->momentum().py(), p->momentum().pz(), p->momentum().e(), p->generated_mass(), p->status());
    flush();
}


inline void WriterAscii::write_string(const std::string &str) {
    // First let's check if string will fit into the buffer
    if ( m_buffer + m_buffer_size > m_cursor + str.length() ) {
        strncpy(m_cursor, str.data(), str.length());
        m_cursor += str.length();
        flush();
    }
    // If not, flush the buffer and write the string directly
    else {
        forced_flush();
        m_stream->write(str.data(), str.length());
    }
}


void WriterAscii::close() {
    std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
    if (ofs && !ofs->is_open()) return;
    forced_flush();
    const std::string footer("HepMC::Asciiv3-END_EVENT_LISTING\n\n");
    if (m_stream) m_stream->write(footer.data(),footer.length());
    if (ofs) ofs->close();
}
bool WriterAscii::failed() { return (bool)m_file.rdstate(); }

void WriterAscii::set_precision(const int& prec ) {
    if (prec < 2 || prec > 24) return;
    m_precision = prec;
}

int WriterAscii::precision() const {
    return m_precision;
}

void WriterAscii::set_buffer_size(const size_t& size ) {
    if (m_buffer) return;
    if (size < 1024) return;
    m_buffer_size = size;
}


} // namespace HepMC3

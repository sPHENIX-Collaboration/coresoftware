#include "WriterDOT.h"
namespace HepMC3
{
WriterDOT::WriterDOT(const std::string &filename,std::shared_ptr<GenRunInfo> /*run*/): m_file(filename),
    m_stream(&m_file),
    m_style(0),
    m_buffer(nullptr),
    m_cursor(nullptr),
    m_buffer_size( 256*1024 )
{
    if ( !m_file.is_open() ) {
        HEPMC3_ERROR( "WriterDOT: could not open output file: "<<filename )
    }
}

WriterDOT::WriterDOT(std::ostream &stream, std::shared_ptr<GenRunInfo> /*run*/)
    : m_file(),
      m_stream(&stream),
      m_style(0),
      m_buffer(nullptr),
      m_cursor(nullptr),
      m_buffer_size( 256*1024 )
{}


void WriterDOT::close() {
    std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
    if (ofs && !ofs->is_open()) return;
    forced_flush();
    if (ofs) ofs->close();
}
/// @brief Detects if particle is parton. Might be used to draw partons different from hadrons
bool is_parton(const int& pd )
{
    bool parton=false;

    if (pd==81||pd==82||pd<25) parton=true;
    if (
        (pd/1000==1||pd/1000==2||pd/1000==3||pd/1000==4||pd/1000==5)
        &&(pd%1000/100==1||pd%1000/100==2||pd%1000/100==3||pd%1000/100==4)
        &&(pd%100==1||pd%100==3)
    )
        parton = true;
    return parton;
}
void WriterDOT::write_event(const GenEvent &evt)
{
    allocate_buffer();
    if ( !m_buffer ) return;
    flush();
    m_cursor += sprintf(m_cursor, "digraph graphname%d {\n",evt.event_number());
    m_cursor += sprintf(m_cursor, "v0[label=\"Machine\"];\n");
    for(auto v: evt.vertices() ) {
        if (m_style != 0)
        {
            if (m_style == 1) //paint decay and fragmentation vertices in green
            {
                if (v->status() == 2) m_cursor += sprintf(m_cursor, "node [color=\"green\"];\n");
                else  m_cursor += sprintf(m_cursor, "node [color=\"black\"];\n");
            }
        }
        m_cursor += sprintf(m_cursor, "node [shape=ellipse];\n");
        m_cursor += sprintf(m_cursor, "v%d[label=\"%d\"];\n", -v->id(),v->id());
        flush();
    }
    for(auto p: evt.beams() ) {
        if (!p->end_vertex()) continue;
        m_cursor += sprintf(m_cursor, "node [shape=point];\n");
        m_cursor += sprintf(m_cursor, "v0 -> v%d [label=\"%d(%d)\"];\n", -p->end_vertex()->id(),p->id(),p->pid());
    }

    for(auto v: evt.vertices() ) {
        for(auto p: v->particles_out() ) {
            {
                if (m_style != 0)
                {
                    if (m_style == 1) //paint suspected partons and 81/82 in red
                    {
                        if (is_parton(std::abs(p->pid()))&&p->status()!=1) m_cursor += sprintf(m_cursor, "edge [color=\"red\"];\n");
                        else        m_cursor +=sprintf(m_cursor, "edge [color=\"black\"];\n");
                    }
                }
                if (!p->end_vertex())
                {
                    m_cursor += sprintf(m_cursor, "node [shape=point];\n");
                    m_cursor += sprintf(m_cursor, "v%d -> o%d [label=\"%d(%d)\"];\n", -v->id(), p->id(), p->id(), p->pid());
                    flush();
                    continue;
                }
                m_cursor += sprintf(m_cursor, "node [shape=ellipse];\n");
                m_cursor += sprintf(m_cursor, "v%d -> v%d [label=\"%d(%d)\"];\n", -v->id(), -p->end_vertex()->id(), p->id(), p->pid());
                flush();
            }
        }
    }
    m_cursor += sprintf(m_cursor, "labelloc=\"t\";\nlabel=\"Event %d; Vertices %lu; Particles %lu;\";\n", evt.event_number(), evt.vertices().size(), evt.particles().size());
    m_cursor += sprintf(m_cursor,"}\n\n");
    forced_flush();
}
void WriterDOT::allocate_buffer() {
    if ( m_buffer ) return;
    while( m_buffer == nullptr && m_buffer_size >= 256 ) {
        try {
            m_buffer = new char[ m_buffer_size ]();
        }     catch (const std::bad_alloc& e) {
            delete[] m_buffer;
            m_buffer_size /= 2;
            HEPMC3_WARNING( "WriterDOT::allocate_buffer: buffer size too large. Dividing by 2. New size: " << m_buffer_size << e.what())
        }
    }

    if ( !m_buffer ) {
        HEPMC3_ERROR( "WriterDOT::allocate_buffer: could not allocate buffer!" )
        return;
    }
    m_cursor = m_buffer;
}
inline void WriterDOT::flush() {
    // The maximum size of single add to the buffer (other than by
    // using WriterDOT::write) is 32 bytes. This is a safe value as
    // we will not allow precision larger than 24 anyway
    if ( m_buffer + m_buffer_size < m_cursor + 256 ) {
        m_stream->write( m_buffer, m_cursor - m_buffer );
        m_cursor = m_buffer;
    }
}
inline void WriterDOT::forced_flush() {
    m_stream->write( m_buffer, m_cursor - m_buffer );
    m_cursor = m_buffer;
}

} // namespace HepMC3

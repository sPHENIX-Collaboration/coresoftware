#include "WriterHEPEVTZEUS.h"
#include "HepMC3/HEPEVT_Wrapper.h"
namespace HepMC3
{
WriterHEPEVTZEUS::WriterHEPEVTZEUS(const std::string &filename):WriterHEPEVT(filename) {}
void WriterHEPEVTZEUS::write_hepevt_event_header()
{
    char buf[512];//Note: the format is fixed, so no reason for complicatied tratment
    char* cursor = &(buf[0]);
    cursor += sprintf(cursor, " E % 12i% 12i% 12i\n", m_hepevt_interface.event_number(), 0, m_hepevt_interface.number_entries());
    unsigned long length = cursor - &(buf[0]);
    m_stream->write( buf, length );
}
void WriterHEPEVTZEUS::write_hepevt_particle( int index, bool iflong)
{
    if (!iflong) printf("INFO: the parameter is ignored as HEPEVTZEUS always uses long format\n");
    char buf[512];//Note: the format is fixed, so no reason for complicatied tratment
    char* cursor = &(buf[0]);
    cursor += sprintf(cursor,"% 12i% 8i", m_hepevt_interface.status(index), m_hepevt_interface.id(index));
    cursor += sprintf(cursor,"% 8i% 8i", m_hepevt_interface.first_parent(index), m_hepevt_interface.last_parent(index));
    cursor += sprintf(cursor,"% 8i% 8i", m_hepevt_interface.first_child(index), m_hepevt_interface.last_child(index));
    cursor += sprintf(cursor,      "% 19.11E% 19.11E% 19.11E% 19.11E% 19.11E\n", m_hepevt_interface.px(index), m_hepevt_interface.py(index), m_hepevt_interface.pz(index), m_hepevt_interface.e(index), m_hepevt_interface.m(index));
    cursor += sprintf(cursor, "%-52s% 19.11E% 19.11E% 19.11E% 19.11E% 19.11E\n", " ", m_hepevt_interface.x(index), m_hepevt_interface.y(index), m_hepevt_interface.z(index), m_hepevt_interface.t(index), 0.0);
    std::ptrdiff_t length = cursor - &(buf[0]);
    m_stream->write( buf, length );
}
}// namespace HepMC3

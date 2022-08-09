#ifndef HEPMC3_WRITERHEPEVTZEUS_H
#define HEPMC3_WRITERHEPEVTZEUS_H
///
/// @file  WriterHEPEVTZEUS.h
/// @brief Definition of class \b WriterHEPEVTZEUS
///
/// @class HepMC3::WriterHEPEVTZEUS
/// @brief GenEvent I/O output to files readable by ZEUS software
///
/// @ingroup Examples
///
#include "HepMC3/WriterHEPEVT.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/Data/GenEventData.h"
namespace HepMC3
{
class WriterHEPEVTZEUS : public  WriterHEPEVT
{
public:
    /** @brief Constructor */
    WriterHEPEVTZEUS(const std::string &filename);
    /** @brief Write the header */
    void write_hepevt_event_header() override;
    /** @brief Write particles */
    void write_hepevt_particle( int index, bool iflong = true )  override;
};
}
#endif

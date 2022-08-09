#ifndef HEPMC3_WRITERDOT_H
#define HEPMC3_WRITERDOT_H
///
/// @file  WriterDOT.h
/// @brief Definition of class \b WriterDOT
///
/// @class HepMC3::WriterDOT
/// @brief GenEvent I/O output to dot files that should be processed by graphviz or other software
///
/// @ingroup Examples
///
#include <string>
#include <fstream>
#include "HepMC3/Writer.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Data/GenEventData.h"
namespace HepMC3
{
class WriterDOT : public Writer
{
public:
    /// @brief Constructor
    /// @warning If file already exists, it will be cleared before writing
    WriterDOT(const std::string &filename,std::shared_ptr<GenRunInfo> run = std::shared_ptr<GenRunInfo>());
    /// @brief Constructor from ostream
    WriterDOT(std::ostream& stream,std::shared_ptr<GenRunInfo> run =std:: shared_ptr<GenRunInfo>());
    /// @brief Write event to file
    ///
    /// @param[in] evt Event to be serialized
    void write_event(const GenEvent &evt);
    /// @brief Return status of the stream
    bool failed() override {
        return (bool)m_file.rdstate();
    }
    /// @brief Close file stream
    void close() override;
    /// @brief Close file stream
    void set_style(const int& istyle) {
        m_style=istyle;
    };

private:
    void allocate_buffer(); //!< allocates buffer for output
    void flush(); //!< flushes output buffer
    void forced_flush(); //!< flushes output buffer
    std::ofstream m_file; //!< Output file
    std::ostream* m_stream; //!< Output stream
    int m_style; //!< style of dot file
    char* m_buffer;  //!< Stream buffer
    char* m_cursor;  //!< Cursor inside stream buffer
    unsigned long m_buffer_size; //!< Buffer size
};
}
#endif

#ifndef HEPMC3_ANALYSISEXAMPLE_H
#define HEPMC3_ANALYSISEXAMPLE_H
///
/// @file  AnalysisExample.h
/// @brief Definition of class \b AnalysisExample
///
/// @class HepMC3::AnalysisExample
/// @brief Example analysis. Produces a rapidity distribution of final state particles.
///
/// @ingroup Examples
///
#include <string>
#include <fstream>
#include "HepMC3/Writer.h"
#include "HepMC3/Version.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
namespace HepMC3
{
class AnalysisExample : public Writer
{
public:
    /// @brief Constructor
    /// @warning If file already exists, it will be cleared before writing
    AnalysisExample(const std::string &filename,std::shared_ptr<GenRunInfo> run);
    /// @brief Constructor from ostream
    AnalysisExample(std::ostream& stream,std::shared_ptr<GenRunInfo> run);
    /// @brief Write event to file
    ///
    /// @param[in] evt Event to be serialized
    void write_event(const GenEvent &evt)  override;
    /// @brief Return status of the stream
    bool failed() override {
        return (bool)m_file.rdstate();
    }
    /// @brief Close file stream
    void close() override;
    /// @brief destructor
    ~AnalysisExample() { close(); }

    double m_sum_of_weights=0;  //!< Sum of event weights
    double m_sum_of_weights2=0; //!< Sum of event weights**2
    std::map<std::string, std::vector<double> > m_bins;  //!< Binings
    std::map<std::string, std::vector<double> > m_vals;  //!< Values
    std::map<std::string, std::vector<double> > m_errs;  //!< Uncertainties
private:
    std::ofstream m_file; //!< Output file
    std::ostream* m_stream; //!< Output stream
};
}
#endif

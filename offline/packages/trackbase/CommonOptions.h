#ifndef __COMMONOPTIONS_H__
#define __COMMONOPTIONS_H__

#include <Acts/Utilities/EnumBitwiseOperators.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <ActsExamples/Framework/RandomNumbers.hpp>
#include <ActsExamples/Framework/Sequencer.hpp>

#include <string>

#include <boost/program_options.hpp>

namespace ActsExamples {

enum class OutputFormat : uint8_t {
  DirectoryOnly = 0,
  Root = 1,
  Csv = 2,
  Obj = 4,
  Json = 8,
  Cbor = 16,
  Txt = 32,
  All = std::numeric_limits<uint8_t>::max()
};

namespace Options {

/// Construct the options description with minimal default options.
///
/// @param caption Optional help text caption
boost::program_options::options_description makeDefaultOptions(
    const std::string& caption = std::string());

/// Add common geometry-related options.
void addGeometryOptions(boost::program_options::options_description& opt);

/// Add common material-related options.
void addMaterialOptions(boost::program_options::options_description& opt);

/// Parse options and return the resulting variables map.
///
/// Automatically prints the help text if requested.
///
/// @returns Empty variables map if help text was shown.
boost::program_options::variables_map parse(
    const boost::program_options::options_description& opt, int argc,
    char* argv[]);



}  // namespace Options
}  // namespace ActsExamples

#endif  // __COMMONOPTIONS_H__
#ifndef _MAGNETICFIELDOPTIONS_H
#define _MAGNETICFIELDOPTIONS_H

#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <ActsExamples/MagneticField/MagneticField.hpp>

namespace boost::program_options
{
  class options_description;
  class variables_map;
}  // namespace boost::program_options

namespace ActsExamples::Options
{
  using Description = ::boost::program_options::options_description;
  using Variables = ::boost::program_options::variables_map;
}  // namespace ActsExamples::Options

namespace ActsExamples {

namespace Options {

/// Add magnetic field options with a `bf-` prefix.
void addMagneticFieldOptions(Description& desc);

/// Read and create the magnetic field from the given user variables.
std::shared_ptr<Acts::MagneticFieldProvider> readMagneticField(
    const Variables& vars);

}  // namespace Options
}  // namespace ActsExamples

#endif  // _MAGNETICFIELDOPTIONS_H

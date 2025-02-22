#ifndef _MAGNETICFIELDOPTIONS_H
#define _MAGNETICFIELDOPTIONS_H

#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <ActsExamples/MagneticField/MagneticField.hpp>
#include <ActsExamples/Utilities/OptionsFwd.hpp>

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
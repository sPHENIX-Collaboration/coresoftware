#ifndef PHGARFIELD_ROOTOUTPUT_H
#define PHGARFIELD_ROOTOUTPUT_H
#include <TFile.h>
#include <ostream>

// Final file Write may return zero after all detached objects were explicitly
// written. The caller must independently reopen and validate the required data.
inline bool PHGarfieldFlushAndCloseRootOutput(TFile& file, std::ostream& error)
{
  file.Write();
  file.Flush();
  file.Close();
  if (file.TestBit(TFile::kWriteError))
  {
    error << "Forward-transport ROOT write/flush/close failed\n";
    return false;
  }
  return true;
}
#endif

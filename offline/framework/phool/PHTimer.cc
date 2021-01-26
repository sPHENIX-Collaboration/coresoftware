// $Id: PHTimer.C,v 1.2 2010/12/05 02:30:28 bbannier Exp $

/*!
\file    PHTimer.cxx
\brief   high precision timer
\author  Sean Kelly, Hugo Pereira
\version $Revision: 1.2 $
\date    $Date: 2010/12/05 02:30:28 $
*/

#include "PHTimer.h"

#include <climits>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>

//______________________________________________________________________
// static members
PHTimer::Frequency PHTimer::_frequency = PHTimer::Frequency();
const double PHTimer::_twopower32 = pow(2, 32);

//______________________________________________________________________
void PHTimer::Frequency::set_cpu_freq(const std::string &path)
{
  // Set the default to 2 GHz
  _frequency = 2e9;

  // Open the cpuinfo file
  std::ifstream cpuProcFile(path);
  if (!cpuProcFile.is_open())
    throw std::runtime_error(std::string("cpu info. unavailable"));
  else
  {
    // Now parse it looking for the string "cpu MHz"
    char readLine[1024];
    std::string searchString("cpu MHz");
    while ((cpuProcFile.rdstate() & std::ios::failbit) == 0)
    {
      // Read into the raw char array and then construct a std::string
      // (std::string) to do the searching
      cpuProcFile.getline(readLine, 1024);
      std::string readLineString(readLine);
      if (readLineString.find(searchString) != std::string::npos)
      {
        // Now look for the :, the clock frequency will follow it
        size_t semicolonPosition = readLineString.find(':', 0);
        if (semicolonPosition == std::string::npos) throw std::runtime_error(std::string("wrong format for cpu info file"));
        std::string frequencyString(readLineString.substr(semicolonPosition + 1));

        // Make a string stream for the conversion to floating number
        double freqMHz = 0;
        std::istringstream frequencySstream(frequencyString);

        frequencySstream >> freqMHz;
        _frequency = freqMHz * 1e6;
      }
    }
  }
}

//______________________________________________________________________
double PHTimer::get_difference(const PHTimer::time_struct& t0, const PHTimer::time_struct& t1)
{
  unsigned long diff_high = t0._high - t1._high;
  unsigned long diff_low;
  if (t0._low < t1._low)
  {
    --diff_high;
    diff_low = (UINT_MAX - t1._low) + t0._low + 1;
  }
  else
    diff_low = t0._low - t1._low;

  return (_twopower32 * diff_high + diff_low) * _frequency.period();
}

//_______________________________________________________
void PHTimer::PRINT(std::ostream& os, const std::string& message)
{
  const int max_col = 80;
  if (!message.size())
  {
    os << std::string(max_col, '-') << std::endl;
    return;
  }
  int fill = max_col - message.size() - 2;
  int pre = static_cast<int>(std::floor(fill / 2.0));
  int post = fill - pre;
  os << std::string(pre, '-') << " ";
  os << message << " ";
  os << std::string(post, '-') << std::endl;
}

#include "mvtx_utils.h"




float mvtx_utils::getStrobeLength(const int& runNumber)
{
  float strobeWidth = std::numeric_limits<float>::quiet_NaN();

  std::string executable_command = "psql -h sphnxdaqdbreplica daq --csv -c \"SELECT strobe FROM mvtx_strobe WHERE hostname = \'mvtx0\' AND runnumber = ";
  executable_command += std::to_string(runNumber);
  executable_command += ";\" | tail -n 1";

  std::array<char, 128> buffer = {};
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(executable_command.c_str(), "r"), pclose);
  if (!pipe)
  {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
  {
    result += buffer.data();
  }
  try
  {
    strobeWidth = stof(result);
  }
  catch (std::invalid_argument const& ex)
  {
    std::cout << "mvtx_utils::getStrobeLength() Run number " << runNumber << " has no strobe length in the DAQ database, returning NAN" << std::endl;
  }
  
  return strobeWidth;
}
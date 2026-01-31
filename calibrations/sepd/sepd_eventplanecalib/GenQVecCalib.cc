#include "QVecCalib.h"

#include <iostream>

int main(int argc, const char* const argv[])
{
  const std::vector<std::string> args(argv, argv + argc);

  if (args.size() < 4 || args.size() > 7)
  {
    std::cout << "Usage: " << args[0] << " <input_file> <input_hist_file> <input_Q_calib> [pass] [events] [output_directory]" << std::endl;
    return 1;  // Indicate error
  }

  const std::string &input_file = args[1];
  const std::string &input_hist = args[2];
  const std::string &input_Q_calib = args[3];
  const std::string &pass_str = (args.size() >= 5) ? args[4] : "ComputeRecentering"; // Default to the first pass
  std::string output_dir = (args.size() >= 7) ? args[6] : ".";

  const std::map<std::string, QVecCalib::Pass> pass_map = {
      {"ComputeRecentering", QVecCalib::Pass::ComputeRecentering},
      {"ApplyRecentering", QVecCalib::Pass::ApplyRecentering},
      {"ApplyFlattening", QVecCalib::Pass::ApplyFlattening}
  };

  QVecCalib::Pass pass = QVecCalib::Pass::ComputeRecentering;
  if (pass_map.contains(pass_str))
  {
    pass = pass_map.at(pass_str);
  }
  else
  {
    std::cout << "Error: Invalid pass specified: " << pass_str << std::endl;
    std::cout << "Available passes are: ComputeRecentering, ApplyRecentering, ApplyFlattening" << std::endl;
    return 1;
  }

  try
  {
    long long events = (args.size() >= 6) ? std::stoll(args[5]) : 0;
    QVecCalib analysis(input_file, input_hist, input_Q_calib, static_cast<int>(pass), events, output_dir);
    analysis.run();
  }
  catch (const std::invalid_argument& e)
  {
    std::cout << "Error: events must be an integer" << std::endl;
    return 1;
  }
  catch (const std::exception& e)
  {
    std::cout << "An exception occurred: " << e.what() << std::endl;
    return 1;
  }

  std::cout << "======================================" << std::endl;
  std::cout << "done" << std::endl;
  return 0;
}

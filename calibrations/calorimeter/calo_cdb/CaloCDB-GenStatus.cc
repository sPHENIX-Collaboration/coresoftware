#include "GenStatus.h"

#include <iostream>
#include <string>
#include <vector>

int main(int argc, const char* const argv[])
{
  const std::vector<std::string> args(argv, argv + argc);

  if (args.size() < 2 || args.size() > 4)
  {
    std::cerr << "usage: " << args[0] << " <input_list> [output_directory] [calofitting_prefix]" << std::endl;
    std::cerr << "  input_list: path to the input list file" << std::endl;
    std::cerr << "  output_directory: (optional) path to the output directory (default: 'output')" << std::endl;
    std::cerr << "  calofitting_prefix: (optional) dsttype prefix of the calofitting QA files (default: 'HIST_DST_CALOFITTING')" << std::endl;
    return 1;  // Indicate error
  }

  const std::string& input_list_path = args[1];
  std::string output_dir_path = "output";
  std::string calofitting_prefix = "HIST_DST_CALOFITTING";

  if (args.size() >= 3)
  {
    output_dir_path = args[2];
  }
  if (args.size() >= 4)
  {
    calofitting_prefix = args[3];
  }

  GenStatus genStat(calofitting_prefix);
  genStat.process(input_list_path, output_dir_path);

  std::cout << "======================================" << std::endl;
  std::cout << "done" << std::endl;
  return 0;
}

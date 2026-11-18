#include "FilterDatasets.h"

#include <phool/recoConsts.h>

#include <iostream>
#include <string>
#include <vector>

int main(int argc, const char* const argv[])
{
  const std::vector<std::string> args(argv, argv + argc);

  if (args.size() < 2 || args.size() > 5)
  {
    std::cout << "usage: " << args[0] << " <input_csv> [output_directory] [debug] [calofitting_prefix]" << std::endl;
    std::cout << "input: input csv" << std::endl;
    std::cout << "output: Output directory. Default: ." << std::endl;
    std::cout << "debug: debug mode. Default false" << std::endl;
    std::cout << "calofitting_prefix: dsttype prefix of the calofitting QA files. Default: HIST_DST_CALOFITTING" << std::endl;
    return 1;  // Indicate error
  }

  const std::string& input_csv = args[1];
  std::string output_dir_path = ".";
  bool debug = false;
  std::string calofitting_prefix = "HIST_DST_CALOFITTING";

  if (args.size() >= 3)
  {
    output_dir_path = args[2];
  }
  if (args.size() >= 4)
  {
    debug = std::stoi(args[3]);
  }
  if (args.size() >= 5)
  {
    calofitting_prefix = args[4];
  }

  recoConsts* rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  FilterDatasets filter(debug, calofitting_prefix);
  filter.process(input_csv, output_dir_path);

  std::cout << "======================================" << std::endl;
  std::cout << "done" << std::endl;
  return 0;
}

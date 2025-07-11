#include "filter-datasets.h"

#include <iostream>

int main(int argc, const char* const argv[])
{
  const std::vector<std::string> args(argv, argv + argc);

  if (args.size() < 2 || args.size() > 4)
  {
    std::cout << "usage: " << args[0] << " <input_csv> [output_directory] [debug]" << std::endl;
    std::cout << "input: input csv" << std::endl;
    std::cout << "output: Output directory. Default: ." << std::endl;
    std::cout << "debug: debug mode. Default false" << std::endl;
    return 1;  // Indicate error
  }

  const std::string& input_csv = args[1];
  std::string output_dir_path = ".";
  Bool_t debug = false;

  if (args.size() >= 3)
  {
    output_dir_path = args[2];
  }
  if (args.size() >= 4)
  {
    debug = std::stoi(args[3]);
  }

  FilterDatasets filter(debug);
  filter.process(input_csv, output_dir_path);

  std::cout << "======================================" << std::endl;
  std::cout << "done" << std::endl;
  return 0;
}

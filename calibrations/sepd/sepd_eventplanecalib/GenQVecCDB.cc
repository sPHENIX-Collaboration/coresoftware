#include "QVecCDB.h"

#include <iostream>
#include <format>

int main(int argc, const char* const argv[])
{
  const std::vector<std::string> args(argv, argv + argc);

  if (args.size() < 3 || args.size() > 5)
  {
    std::cout << "Usage: " << args[0] << " input_file runnumber [output_dir] [cdb_tag]" << std::endl;
    return 1;
  }

  const std::string &input_file = args[1];
  const std::string output_dir = (args.size() >= 4) ? args[3] : ".";
  const std::string cdb_tag = (args.size() >= 5) ? args[4] : "new_newcdbtag_v008";

  try
  {
    int runnumber = std::stoi(args[2]);

    std::cout << std::format("{:#<20}\n", "");
    std::cout << std::format("Analysis Params\n");
    std::cout << std::format("Input File: {}\n", input_file);
    std::cout << std::format("Run: {}\n", runnumber);
    std::cout << std::format("Output Dir: {}\n", output_dir);
    std::cout << std::format("CDB Tag: {}\n", cdb_tag);
    std::cout << std::format("{:#<20}\n", "");

    QVecCDB analysis(input_file, runnumber, output_dir, cdb_tag);
    analysis.run();
  }
  catch (const std::invalid_argument& e)
  {
    std::cout << "Error: runnumber must be an integer: " << args[2] << std::endl;
    return 1;
  }
  catch (const std::out_of_range& e)
  {
    std::cout << "Error: runnumber is out of range for an integer." << std::endl;
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

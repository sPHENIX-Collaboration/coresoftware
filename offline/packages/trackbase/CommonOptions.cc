
#include "CommonOptions.h"

#include <Acts/Utilities/Helpers.hpp>
#include <ActsExamples/Utilities/Options.hpp>

#include <exception>
#include <fstream>
#include <regex>
#include <system_error>

using namespace boost::program_options;

boost::program_options::options_description
ActsExamples::Options::makeDefaultOptions(const std::string& caption) {

  options_description opt(caption);

  opt.add_options()("help,h", "Produce help message");
  opt.add_options()(
      "loglevel,l", value<size_t>()->default_value(2),
      "The output log level. Please set the wished number (0 = VERBOSE, 1 = "
      "DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).");
  opt.add_options()(
      "response-file", value<std::string>()->default_value(""),
      "Configuration file (response file) replacing command line options.");

  return opt;
}


void ActsExamples::Options::addGeometryOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("geo-surface-loglevel", value<size_t>()->default_value(3),
                    "The outoput log level for the surface building.")(
      "geo-layer-loglevel", value<size_t>()->default_value(3),
      "The output log level for the layer building.")(
      "geo-volume-loglevel", value<size_t>()->default_value(3),
      "The output log level "
      "for the volume "
      "building.");
}

void ActsExamples::Options::addMaterialOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()(
      "mat-input-type", value<std::string>()->default_value("build"),
      "The way material is loaded: 'none', 'build', 'proto', 'file'.")(
      "mat-input-file", value<std::string>()->default_value(""),
      "Name of the material map input file, supported: '.json', '.cbor' or "
      "'.root'.")("mat-output-file", value<std::string>()->default_value(""),
                  "Name of the material map output file (without extension).")(
      "mat-output-sensitives", value<bool>()->default_value(true),
      "Write material information of sensitive surfaces.")(
      "mat-output-approaches", value<bool>()->default_value(true),
      "Write material information of approach surfaces.")(
      "mat-output-representing", value<bool>()->default_value(true),
      "Write material information of representing surfaces.")(
      "mat-output-boundaries", value<bool>()->default_value(true),
      "Write material information of boundary surfaces.")(
      "mat-output-volumes", value<bool>()->default_value(true),
      "Write material information of volumes.")(
      "mat-output-dense-volumes", value<bool>()->default_value(false),
      "Write material information of dense volumes.")(
      "mat-output-allmaterial", value<bool>()->default_value(false),
      "Add protoMaterial to all surfaces and volume for the mapping.");
}


boost::program_options::variables_map ActsExamples::Options::parse(
    const boost::program_options::options_description& opt, int argc,
    char* argv[]) noexcept(false) {
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).run(), vm);
  notify(vm);

  if (vm.count("response-file") != 0u and
      not vm["response-file"].template as<std::string>().empty()) {
    // Load the file and tokenize it
    std::ifstream ifs(vm["response-file"].as<std::string>().c_str());
    if (!ifs) {
      throw(std::system_error(std::error_code(),
                              "Could not open response file."));
    }
    // Read the whole file into a string
    std::stringstream ss;
    ss << ifs.rdbuf();
    std::string rString = ss.str();
    std::vector<std::string> args;
    const std::regex rgx("[ \t\r\n\f]");
    std::sregex_token_iterator iter(rString.begin(), rString.end(), rgx, -1);
    std::sregex_token_iterator end;
    for (; iter != end; ++iter) {
      if (std::string(*iter).empty()) {
        continue;
      }
      args.push_back(*iter);
    }
    // Parse the file and store the options
    store(command_line_parser(args).options(opt).run(), vm);
  }

  // Automatically handle help
  if (vm.count("help") != 0u) {
    std::cout << opt << std::endl;
    vm.clear();
  }
  return vm;
}

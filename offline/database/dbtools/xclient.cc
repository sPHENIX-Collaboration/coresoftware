#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <xpload/xpload.h>

using namespace std;

class ArgParser
{
 public:
  ArgParser(int &argc, char **argv);
  string get_value(const string &option) const;
 private:
   vector<string> args;
};


int main(int argc, char **argv)
{
  ArgParser arg_parser(argc, argv);

  string cfg = arg_parser.get_value("-c");
  string tag = arg_parser.get_value("-t");
  string domain = arg_parser.get_value("-d");
  uint64_t timestamp = std::stoul(arg_parser.get_value("-s"));

  xpload::Configurator config(cfg);

  xpload::Result result = xpload::fetch(tag, domain, timestamp, config);

  if (result.paths.empty())
  {
    cout << "No paths found\n";
  }
  else
  {
    cout << "Found paths:\n";

    for (const string path : result.paths)
      cout << path << '\n';
  }

  return EXIT_SUCCESS;
} 


ArgParser::ArgParser(int &argc, char **argv)
{
  for (int i=1; i < argc; ++i)
    args.push_back( string(argv[i]) );
}

string ArgParser::get_value(const string &option) const
{
  auto itr = find(args.begin(), args.end(), option);

  string value{""};

  if (itr != args.end() && ++itr != args.end()) {
    value = *itr;
  }

  // Default values
  if (value.empty() && option == "-c") return "test";
  if (value.empty() && option == "-t") return "example_tag_1";
  if (value.empty() && option == "-d") return "CEMC";
  if (value.empty() && option == "-s") return std::to_string(UINT64_MAX);

  return value;
}

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <xpload/xpload.h>


int main()
{
  using namespace std;

  string tag = "example_tag_1";
  uint64_t timestamp = 9999999999;

  vector<string> paths = xpload::fetch(tag, timestamp);

  if (paths.empty())
  {
    cout << "No paths found\n";
  }
  else
  {
    cout << "Found paths:\n";

    for (const string& path : paths) 
      cout << path << '\n';
  }

  return EXIT_SUCCESS;
} 

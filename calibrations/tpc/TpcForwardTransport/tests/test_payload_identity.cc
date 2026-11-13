#include "../PHGarfieldForwardMap.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
  void check(bool valid, const std::string& message)
  {
    if (!valid)
    {
      throw std::runtime_error(message);
    }
  }

}  // namespace

int main()
{
  char pattern[] = "payload-identity-XXXXXX";
  const char* directory = mkdtemp(pattern);
  if (!directory)
  {
    return 1;
  }
  const std::filesystem::path root(directory);
  int status = 0;
  try
  {
    const auto digest = [](const std::filesystem::path& path)
    { return PHGarfieldForwardMap::PayloadMd5(path.string()); };
    check(digest(root).empty(), "Empty gas directory accepted");
    std::ofstream(root / "MERGED_E200.gas") << "high\n";
    std::ofstream(root / "MERGED_E010.gas") << "low\n";
    const auto original = digest(root);
    check(original == "2e34e02643f3b293096c14d106848c3f", "Gas inventory ordering or checksum format changed");
    std::ofstream(root / "ignored.gas") << "not an active gas payload";
    check(digest(root) == original, "Inactive gas files changed identity");
    std::ofstream(root / "MERGED_E010.gas") << "changed\n";
    check(digest(root) != original, "Gas content change was not detected");
    std::ofstream(root / "ordinary") << "abc";
    check(digest(root / "ordinary") == "900150983cd24fb0d6963f7d28e17f72", "Regular-file MD5 changed");
    check(digest(root / "absent").empty(), "Missing payload accepted");
    check(PHGarfieldForwardMap::PayloadMd5("").empty(), "Empty payload accepted");
    check(PHGarfieldForwardMap::Fingerprint("") == "cbf29ce484222325" &&
              PHGarfieldForwardMap::Fingerprint("hello") == "a430d84680aabd0b",
          "FNV-1a format changed");
  }
  catch (const std::exception& exception)
  {
    std::cerr << exception.what() << '\n';
    status = 1;
  }
  std::filesystem::remove_all(root);
  return status;
}

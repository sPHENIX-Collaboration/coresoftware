#include "../PHGarfieldAdaptiveMap.h"
#include "../PHGarfieldCompactMap.h"
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

  void check_map_fingerprints()
  {
    // Fixed primitive identities captured before sharing the hash code.
    // These raw records exercise signed indices and signed zero, not a fit.
    PHGarfieldAdaptiveMap::Data data;
    data.geometry = {72, 32, 28, 21.6, 76.4, .28, 102.605, .00755};
    data.policy = {10., 12, 127, 2048, 5000000, 2000000, 2ULL * 1024 * 1024 * 1024};
    data.first_root = 1;
    data.end_root = 2;
    data.maximum_drift_time_ns = 123.5;
    data.points.push_back({7, {.25, 44., -5.}, 1, 1, {-0., .003, 12.}});
    PHGarfieldAdaptiveMap::Cell cell;
    cell.id = 3;
    cell.root_id = 1;
    cell.lo = {.1, 42., -8.};
    cell.hi = {.5, 47., -2.};
    cell.corner_id = {-1, 7, 8, 9, 10, 11, 12, 13};
    cell.child_id = {-1, 4};
    cell.depth = 2;
    cell.directional_error_um = {.2, .3, .4};
    cell.directional_reach_changes = {1, 2, 3};
    data.cells.push_back(cell);
    data.roots.push_back({1, 3});
    check(PHGarfieldAdaptiveMap::ContentFingerprint(data) == "16c2b1f7bb9efb53", "Adaptive primitive identity changed");

    PHGarfieldCompactMap::Data compact;
    compact.geometry = data.geometry;
    compact.source = {"source", "config", "content", "policy", "converter"};
    compact.nodes = {{3, 0, 7}, {1, 2, 0}};
    compact.corners.push_back({0, 1, 2, 3, 4, 5, 6, 7});
    compact.values.push_back({-0., .003, 12.});
    compact.root_lookup = {3, 4};
    check(PHGarfieldCompactMap::ContentFingerprint(compact) == "8b53c8b73c2f6a08", "Compact primitive identity changed");
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
    check_map_fingerprints();
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

#include "../PHGarfieldFieldMapValidation.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  using Domain = PHGarfieldFieldMapValidation::CartesianDomain;
  using Entry = std::array<float, 6>;
  using Axes = std::array<std::vector<float>, 3>;

  const std::array<const char*, 6> branch_names = {{"x", "y", "z", "bx", "by", "bz"}};

  void check(const bool condition, const std::string& message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  class TestDirectory
  {
   public:
    TestDirectory()
    {
      char pattern[] = "field-map-validation-XXXXXX";
      const char* directory = mkdtemp(pattern);
      check(directory != nullptr, "Could not create temporary field-map test directory");
      m_path = directory;
    }

    ~TestDirectory()
    {
      std::error_code error;
      std::filesystem::remove_all(m_path, error);
    }

    std::filesystem::path file(const std::string& name) const
    {
      return m_path / name;
    }

   private:
    std::filesystem::path m_path;
  };

  std::vector<Entry> grid(const Axes& axes)
  {
    std::vector<Entry> entries;
    for (const float x : axes[0])
    {
      for (const float y : axes[1])
      {
        for (const float z : axes[2])
        {
          entries.push_back({{x, y, z, 0.1F * x, -0.2F * y, 1.4F + 0.01F * z}});
        }
      }
    }
    return entries;
  }

  void write_map(const std::filesystem::path& filename,
                 const std::vector<Entry>& entries,
                 const int missing_branch = -1,
                 const int double_branch = -1,
                 const Long64_t declared_entries = -1)
  {
    TFile output(filename.c_str(), "RECREATE");
    check(!output.IsZombie(), "Could not write test map " + filename.string());
    TTree tree("fieldmap", "Self-contained Cartesian magnetic-field test grid");
    Entry value{};
    double double_value = 0.0;
    for (std::size_t branch = 0; branch < branch_names.size(); ++branch)
    {
      if (static_cast<int>(branch) == missing_branch)
      {
        continue;
      }
      const std::string leaf = std::string(branch_names[branch]) +
                               (static_cast<int>(branch) == double_branch ? "/D" : "/F");
      if (static_cast<int>(branch) == double_branch)
      {
        tree.Branch(branch_names[branch], &double_value, leaf.c_str());
      }
      else
      {
        tree.Branch(branch_names[branch], &value[branch], leaf.c_str());
      }
    }
    for (const auto& entry : entries)
    {
      value = entry;
      if (double_branch >= 0)
      {
        double_value = entry[double_branch];
      }
      check(tree.Fill() > 0, "Could not fill test map " + filename.string());
    }
    if (declared_entries >= 0)
    {
      tree.SetEntries(declared_entries);
    }
    check(tree.Write() > 0, "Could not save test map " + filename.string());
  }

  void expect_invalid(const std::string& filename, const std::string& diagnostic)
  {
    // A failed inspection must clear a previously usable domain.
    Domain domain;
    domain.minimum_cm = {{-100.0, -100.0, -100.0}};
    domain.maximum_cm = {{100.0, 100.0, 100.0}};
    std::ostringstream errors;
    check(!PHGarfieldFieldMapValidation::InspectCartesianMap(filename, domain, errors),
          "Accepted invalid field map " + filename);
    check(errors.str().find(diagnostic) != std::string::npos,
          "Unexpected diagnostic for " + filename + ": " + errors.str());
    for (std::size_t axis = 0; axis < domain.minimum_cm.size(); ++axis)
    {
      check(std::isnan(domain.minimum_cm[axis]) && std::isnan(domain.maximum_cm[axis]),
            "Failed inspection retained usable bounds for " + filename);
    }
  }

  void expect_valid(const std::filesystem::path& filename, const Axes& axes)
  {
    Domain domain;
    std::ostringstream errors;
    check(PHGarfieldFieldMapValidation::InspectCartesianMap(filename.string(), domain, errors),
          "Rejected complete Cartesian field map " + filename.string() + ": " + errors.str());
    check(errors.str().empty(), "Successful inspection produced an error diagnostic");
    for (std::size_t axis = 0; axis < axes.size(); ++axis)
    {
      const auto extrema = std::minmax_element(axes[axis].begin(), axes[axis].end());
      check(domain.minimum_cm[axis] == static_cast<double>(*extrema.first) &&
                domain.maximum_cm[axis] == static_cast<double>(*extrema.second),
            "Inspection changed the map-frame coordinate extrema");
    }
  }

  void run_tests(const TestDirectory& directory)
  {
    const Axes axes = {{{-4.0F, -1.0F, 2.0F}, {7.0F, 9.0F}, {-11.0F, -6.0F}}};
    const auto entries = grid(axes);
    auto unordered = entries;
    std::reverse(unordered.begin(), unordered.end());
    std::rotate(unordered.begin(), unordered.begin() + 5, unordered.end());
    const auto valid_map = directory.file("shifted-unordered.root");
    write_map(valid_map, unordered);
    expect_valid(valid_map, axes);

    expect_invalid("", "filename is empty");
    expect_invalid(directory.file("absent.root").string(), "Could not inspect");
    const auto not_root = directory.file("not-root.root");
    {
      std::ofstream output(not_root);
      output << "This file is not a ROOT field map.\n";
    }
    expect_invalid(not_root.string(), "Could not inspect");

    const auto no_tree = directory.file("no-tree.root");
    {
      TFile output(no_tree.c_str(), "RECREATE");
      TTree tree("other", "Not the fieldmap tree");
      check(tree.Write() > 0, "Could not save absent-tree fixture");
    }
    expect_invalid(no_tree.string(), "no fieldmap tree");

    for (std::size_t branch = 0; branch < branch_names.size(); ++branch)
    {
      const auto missing = directory.file("missing-" + std::string(branch_names[branch]) + ".root");
      write_map(missing, entries, static_cast<int>(branch));
      expect_invalid(missing.string(), "missing branch " + std::string(branch_names[branch]));

      const auto wrong_type = directory.file("double-" + std::string(branch_names[branch]) + ".root");
      write_map(wrong_type, entries, -1, static_cast<int>(branch));
      expect_invalid(wrong_type.string(), "Failed to bind");

      for (const float invalid : {std::numeric_limits<float>::quiet_NaN(),
                                  std::numeric_limits<float>::infinity()})
      {
        auto nonfinite_entries = entries;
        nonfinite_entries[3][branch] = invalid;
        const auto nonfinite = directory.file("nonfinite.root");
        write_map(nonfinite, nonfinite_entries);
        expect_invalid(nonfinite.string(), "Invalid magnetic-field grid entry 3");
      }
    }

    const auto empty = directory.file("empty.root");
    write_map(empty, {});
    expect_invalid(empty.string(), "Unsupported magnetic-field grid size 0");

    const auto oversized = directory.file("oversized.root");
    write_map(oversized, {}, -1, -1, 100000001);
    expect_invalid(oversized.string(), "Unsupported magnetic-field grid size 100000001");

    const auto unreadable_entry = directory.file("unreadable-entry.root");
    write_map(unreadable_entry, {}, -1, -1, 8);
    expect_invalid(unreadable_entry.string(), "Invalid magnetic-field grid entry 0");

    auto missing_entries = entries;
    missing_entries.pop_back();
    const auto missing_point = directory.file("missing-point.root");
    write_map(missing_point, missing_entries);
    expect_invalid(missing_point.string(), "not a complete Cartesian lattice");

    auto duplicate_entries = entries;
    duplicate_entries.back() = duplicate_entries.front();
    const auto duplicate_point = directory.file("duplicate-point.root");
    write_map(duplicate_point, duplicate_entries);
    expect_invalid(duplicate_point.string(), "duplicate grid point");

    for (std::size_t axis = 0; axis < axes.size(); ++axis)
    {
      auto nonuniform_axes = axes;
      nonuniform_axes[axis] = {-1.0F, 0.0F, 2.0F};
      const auto nonuniform = directory.file("nonuniform.root");
      write_map(nonuniform, grid(nonuniform_axes));
      expect_invalid(nonuniform.string(), "do not form uniform axes");

      auto singleton_axes = axes;
      singleton_axes[axis] = {2.0F};
      const auto singleton = directory.file("singleton.root");
      write_map(singleton, grid(singleton_axes));
      expect_invalid(singleton.string(), "do not form uniform axes");
    }

    // Preserve the existing 1e-6 * max(1, |spacing|) uniformity tolerance.
    for (const float spacing : {0.25F, 2.0F})
    {
      const float scale = std::max(1.0F, spacing);
      auto tolerant_axes = axes;
      tolerant_axes[0] = {0.0F, spacing + 0.5e-6F * scale, 2.0F * spacing};
      const auto tolerated = directory.file("tolerated-spacing.root");
      write_map(tolerated, grid(tolerant_axes));
      expect_valid(tolerated, tolerant_axes);

      tolerant_axes[0][1] = spacing + 2.0e-6F * scale;
      const auto outside_tolerance = directory.file("outside-tolerance.root");
      write_map(outside_tolerance, grid(tolerant_axes));
      expect_invalid(outside_tolerance.string(), "do not form uniform axes");
    }

    // Failed probes do not affect subsequent successful inspections.
    expect_valid(valid_map, axes);
  }
}  // namespace

int main()
{
  try
  {
    const TestDirectory directory;
    run_tests(directory);
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}

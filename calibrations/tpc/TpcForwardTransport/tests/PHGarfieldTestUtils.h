#ifndef PHGARFIELDTESTUTILS_H
#define PHGARFIELDTESTUTILS_H

#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>

#include <stdexcept>
#include <string>

namespace PHGarfieldTest
{
  inline void check(bool condition, const std::string& message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  inline void check_persisted_multibasket_tree(TFile& file, TDirectory* directory, const char* name)
  {
    auto* tree = dynamic_cast<TTree*>(directory->Get(name));
    check(tree && tree->GetEntries() > 0, "missing nonempty streamed primitive tree");
    auto* branches = tree->GetListOfBranches();
    for (int branch_index = 0; branch_index < branches->GetEntriesFast(); ++branch_index)
    {
      auto* branch = dynamic_cast<TBranch*>(branches->At(branch_index));
      check(branch && branch->GetEntries() == tree->GetEntries(), "streamed primitive branch count mismatch");
      const auto* starts = branch->GetBasketEntry();
      const auto* sizes = branch->GetBasketBytes();
      const int last = branch->GetWriteBasket();
      check(starts && sizes && last >= 0 && last < branch->GetMaxBaskets(), "missing streamed basket metadata");
      Long64_t covered = 0;
      unsigned int persisted = 0;
      // Check persisted coverage before GetEntry can materialize any basket.
      // Inline-only baskets would let the old detached writer falsely pass a
      // value-only roundtrip at small size but cannot scale to a full map.
      for (int basket_index = 0; basket_index <= last; ++basket_index)
      {
        const auto begin = starts[basket_index];
        const auto end = basket_index < last ? starts[basket_index + 1] : tree->GetEntries();
        if (end == begin)
        {
          continue;
        }
        const auto seek = branch->GetBasketSeek(basket_index);
        check(begin == covered && end > begin && end <= tree->GetEntries(), "persisted baskets do not tile every primitive row");
        check(seek > 0 && sizes[basket_index] > 0 && seek <= file.GetSize() && sizes[basket_index] <= file.GetSize() - seek,
              "primitive payload was retained inline instead of streamed into the declared output");
        covered = end;
        ++persisted;
      }
      check(covered == tree->GetEntries() && persisted >= 2, "fixture did not exercise multiple persisted baskets per branch");
      for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
      {
        check(branch->GetEntry(entry) > 0, "streamed primitive branch row is unreadable");
      }
    }
  }

}  // namespace PHGarfieldTest

#endif

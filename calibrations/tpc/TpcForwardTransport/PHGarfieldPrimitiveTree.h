#ifndef PHGARFIELDPRIMITIVETREE_H
#define PHGARFIELDPRIMITIVETREE_H

// Internal primitive ROOT storage shared by adaptive and compact maps.
#include <TBasket.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <ostream>
#include <string>
#include <vector>

namespace PHGarfieldPrimitiveTree
{
  inline bool require(bool valid, std::ostream& error, const char* source, const char* message)
  {
    if (!valid)
    {
      error << source << " LUT: " << message << '\n';
    }
    return valid;
  }
  struct Field
  {
    const char* name;
    void* address;
    const char* type;
    int count;
    const char* leaf;
  };
  template <class T, class Describe>
  bool Write(TDirectory* directory, const char* name, const std::vector<T>& records,
             std::ostream& error, Describe describe, const char* source)
  {
    auto* file = directory ? directory->GetFile() : nullptr;
    if (!require(file && file->IsOpen() && file->IsWritable() && !file->TestBit(TFile::kWriteError),
                 error, source, "primitive tree output is not a writable ROOT file"))
    {
      return false;
    }
    TDirectory::TContext context(directory);
    TTree tree(name, name);
    tree.SetDirectory(directory);
    T value{};
    // A detached tree retains every basket inline. Full maps can exceed ROOT's
    // per-object byte-count limit even when their total allocation budget fits.
    // Stream baskets into the declared destination and keep cluster/header work
    // bounded independently of the file's total number of primitive rows.
    constexpr int BasketBytes = 32 * 1024;
    tree.SetAutoFlush(16384);
    tree.SetAutoSave(0);
    tree.SetMaxVirtualSize(8 * 1024 * 1024);
    for (const auto& field : describe(value))
    {
      if (!require(tree.Branch(field.name, field.address, field.leaf, BasketBytes) != nullptr,
                   error, source, "cannot create primitive branch"))
      {
        return false;
      }
    }
    for (const auto& record : records)
    {
      value = record;
      if (tree.Fill() < 0 || file->TestBit(TFile::kWriteError))
      {
        error << source << " LUT: primitive tree Fill failed for " << name << '\n';
        return false;
      }
    }
    if (!require(tree.FlushBaskets() >= 0 && !file->TestBit(TFile::kWriteError) && tree.Write() > 0 &&
                     !file->TestBit(TFile::kWriteError),
                 error, source, "primitive tree flush/write failed"))
    {
      return false;
    }
    // The scoped tree unregisters itself before its output file can close.
    return true;
  }

  template <class T, class Describe>
  bool Read(TTree* tree, std::vector<T>& records, std::ostream& error, Describe describe, const char* source)
  {
    T value{};
    const auto definitions = describe(value);
    if (!require(tree && tree->GetEntries() == static_cast<Long64_t>(records.size()) &&
                     tree->GetListOfBranches()->GetEntriesFast() == static_cast<int>(definitions.size()),
                 error, source, "wrong tree rows/branches"))
    {
      return false;
    }
    struct Reset
    {
      TTree* tree;
      ~Reset() { tree->ResetBranchAddresses(); }
    } reset{tree};
    std::vector<TBranch*> branches;
    for (const auto& field : definitions)
    {
      auto* branch = tree->GetBranch(field.name);
      auto* leaf = branch ? branch->GetLeaf(field.name) : nullptr;
      if (!require(branch && leaf && !leaf->GetLeafCount() && leaf->GetLenStatic() == field.count &&
                       field.type == std::string(leaf->GetTypeName()) && branch->GetEntries() == tree->GetEntries(),
                   error, source, "wrong primitive branch type/shape/count"))
      {
        return false;
      }
      // Inspect metadata before GetEntry/GetBasket can materialize or replace a
      // basket. A positive aggregate read is not proof of persisted coverage.
      if (!require(std::string(branch->GetFileName()).empty(), error, source, "external branch storage is not declared"))
      {
        return false;
      }
      const auto entries = tree->GetEntries();
      if (entries)
      {
        const auto maximum = branch->GetMaxBaskets();
        const auto last = branch->GetWriteBasket();
        const auto* starts = branch->GetBasketEntry();
        const auto* sizes = branch->GetBasketBytes();
        auto* resident = branch->GetListOfBaskets();
        auto* file = tree->GetCurrentFile();
        if (!require(last >= 0 && last < maximum && starts && sizes && file, error, source, "invalid basket index metadata"))
        {
          return false;
        }
        Long64_t covered = 0;
        for (int index = 0; index <= last; ++index)
        {
          const auto begin = starts[index];
          const auto end = index < last ? starts[index + 1] : entries;
          if (end == begin)
          {
            continue;
          }
          if (!require(begin == covered && begin >= 0 && end > begin && end <= entries, error, source, "basket ranges do not tile declared entries"))
          {
            return false;
          }
          const auto seek = branch->GetBasketSeek(index);
          const bool persisted = seek > 0 && sizes[index] > 0 && seek <= file->GetSize() && sizes[index] <= file->GetSize() - seek;
          const auto* basket = resident ? dynamic_cast<TBasket*>(resident->At(index)) : nullptr;
          const bool inline_complete = basket && basket->GetNevBuf() == end - begin;
          if (!require(persisted || inline_complete, error, source, "declared branch rows have no persisted/inline basket"))
          {
            return false;
          }
          covered = end;
        }
        if (!require(covered == entries, error, source, "basket coverage is incomplete"))
        {
          return false;
        }
      }
      branch->SetAddress(field.address);
      branches.push_back(branch);
    }
    for (std::size_t i = 0; i < records.size(); ++i)
    {
      value = T{};
      for (auto* branch : branches)
      {
        if (!require(branch->GetEntry(i) > 0, error, source, "unreadable primitive branch entry"))
        {
          return false;
        }
      }
      records[i] = value;
    }
    return true;
  }
}  // namespace PHGarfieldPrimitiveTree

#endif

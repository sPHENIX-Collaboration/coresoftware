#ifndef INTT_ARBORIST_H
#define INTT_ARBORIST_H

#include "InttMap.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

class InttArborist : public SubsysReco
{
 public:
  InttArborist();
  ~InttArborist();

  typedef Short_t small_t;   // for branches that store many values per event
  typedef uint64_t large_t;  // for branches that store one value per event

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int EndRun(int const) override;

  int CreateOutputFile(std::string const&);
  int WriteOutputFile();  // called by EndRun

 private:
  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;

  const std::string m_intt_raw_node_name = "INTTRAWHIT";
  const std::string m_gl1_raw_node_name = "GL1RAWHIT";
  const std::string m_tree_name = "intt_tree";

  typedef std::map<InttMap::RawData_s, unsigned int, InttMap::RawDataComparator> clone_map_t;
  clone_map_t m_clone_map;

  typedef std::map<std::string, std::vector<small_t>*> small_map_t;
  small_map_t m_small_branches;

  typedef std::map<std::string, large_t> large_map_t;
  large_map_t m_large_branches;
  // see the definition (.cc) file for branch names
};

#endif  // INTT_ARBORIST_H

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef METADATA_METADATATREEWRITER_H
#define METADATA_METADATATREEWRITER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

// Fun4All module that puts the provenance tree (MetaDataTree::Make) of the job into an output
// file shared through PHTFileServer, so it ends up next to the histograms/trees of the analysis
// modules writing to the same file. In the macro:
//
//   R__LOAD_LIBRARY(libmetadata.so)
//   #include <metadata/MetaDataTreeWriter.h>
//   ...
//   se->registerSubsystem(new MetaDataTreeWriter("MetaDataTreeWriter", outfile));
//
// To share a file with an analysis module, give the same name the module passes to
// PHTFileServer::open(). Do not use the name of a file a module opens itself with new TFile: two
// TFile objects writing one file corrupt it. Any other name gets a file of its own, holding only
// the provenance tree. The tree is made in End(), when the
// CDB payloads and input files of the whole job are known, and written to the file by
// PHTFileServer when the last module using the file calls PHTFileServer::write().
class MetaDataTreeWriter : public SubsysReco
{
 public:
  explicit MetaDataTreeWriter(const std::string &name = "MetaDataTreeWriter", const std::string &filename = "");
  ~MetaDataTreeWriter() override = default;

  // opens (or joins) the output file with PHTFileServer
  int Init(PHCompositeNode *topNode) override;

  // makes the provenance tree in the output file and hands the file back to PHTFileServer
  int End(PHCompositeNode *topNode) override;

  void set_filename(const std::string &filename) { m_FileName = filename; }
  const std::string &get_filename() const { return m_FileName; }

  // name of the tree (default "provenance", which check_provenance etc. expect)
  void set_treename(const std::string &treename) { m_TreeName = treename; }
  const std::string &get_treename() const { return m_TreeName; }

 private:
  std::string m_FileName;
  std::string m_TreeName{"provenance"};
  bool m_Opened{false};
};

#endif  // METADATA_METADATATREEWRITER_H

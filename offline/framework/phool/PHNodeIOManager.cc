//  Implementation of class PHNodeIOManager
//  Author: Matthias Messer

#include "PHNodeIOManager.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHNodeIterator.h"
#include "phooldefs.h"

#include <TBranch.h>  // for TBranch
#include <TBranchElement.h>
#include <TBranchObject.h>
#include <TClass.h>
#include <TDirectory.h>  // for TDirectory
#include <TFile.h>
#include <TLeafObject.h>
#include <TObjArray.h>  // for TObjArray
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/algorithm/string.hpp>
#pragma GCC diagnostic pop

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

PHNodeIOManager::PHNodeIOManager(const std::string& f,
                                 const PHAccessType a)
{
  isFunctionalFlag = setFile(f, "titled by PHOOL", a) ? 1 : 0;
}

PHNodeIOManager::PHNodeIOManager(const std::string& f, const std::string& title,
                                 const PHAccessType a)
{
  isFunctionalFlag = setFile(f, title, a) ? 1 : 0;
}

PHNodeIOManager::PHNodeIOManager(const std::string& f, const PHAccessType a,
                                 const PHTreeType treeindex)
{
  if (treeindex != PHEventTree)
  {
    std::ostringstream temp;
    temp << TreeName << treeindex;  // create e.g. T1
    TreeName = temp.str();
  }
  isFunctionalFlag = setFile(f, "titled by PHOOL", a) ? 1 : 0;
}

PHNodeIOManager::~PHNodeIOManager()
{
  closeFile();
  delete file;
}

void PHNodeIOManager::closeFile()
{
  if (file)
  {
    if (accessMode == PHWrite || accessMode == PHUpdate)
    {
      file->Write();
    }
    file->Close();
  }
}

bool PHNodeIOManager::setFile(const std::string& f, const std::string& title,
                              const PHAccessType a)
{
  filename = f;
  accessMode = a;
  if (file)
  {
    if (file->IsOpen())
    {
      closeFile();
    }
    delete file;
    file = nullptr;
  }
  std::string currdir = gDirectory->GetPath();
  gROOT->cd();
  switch (accessMode)
  {
  case PHWrite:
    file = TFile::Open(filename.c_str(), "RECREATE", title.c_str());
    if (!file)
    {
      return false;
    }
    file->SetCompressionLevel(CompressionLevel);
    tree = new TTree(TreeName.c_str(), title.c_str());
    tree->SetMaxTreeSize(900000000000LL);  // set max size to ~900 GB
    gROOT->cd(currdir.c_str());
    return true;
    break;
  case PHReadOnly:
    file = TFile::Open(filename.c_str());
    tree = nullptr;
    if (!file)
    {
      return false;
    }
    selectObjectToRead("*", true);
    gROOT->cd(currdir.c_str());
    return true;
    break;
  case PHUpdate:
    file = TFile::Open(filename.c_str(), "UPDATE", title.c_str());
    if (!file)
    {
      return false;
    }
    file->SetCompressionLevel(CompressionLevel);
    tree = new TTree(TreeName.c_str(), title.c_str());
    gROOT->cd(currdir.c_str());
    return true;
    break;
  }

  return false;
}

bool PHNodeIOManager::write(PHCompositeNode* topNode)
{
  // The write function of the PHCompositeNode topNode will
  // recursively call the write functions of its subnodes, thus
  // constructing the path-string which is then stored as name of the
  // Root-branch corresponding to the data of each PHRootIODataNode.
  topNode->write(this);

  // Now all PHRootIODataNodes should have called the write function
  // of this I/O-manager and thus created their branch. The tree can
  // be filled.
  if (file && tree)
  {
    tree->Fill();
    eventNumber++;
    return true;
  }

  return false;
}

bool PHNodeIOManager::write(TObject** data, const std::string& path, int buffersize, int splitlevel)
{
  if (file && tree)
  {
    TBranch* thisBranch = tree->GetBranch(path.c_str());
    if (!thisBranch)
    {
      // the buffersize and splitlevel are set on the first call
      // when the branch is created, the values come from the caller
      // which is the node which writes itself
      tree->Branch(path.c_str(), (*data)->ClassName(),
                   data, buffersize, splitlevel);
    }
    else
    {
      thisBranch->SetAddress(data);
    }
    return true;
  }

  return false;
}

bool PHNodeIOManager::read(size_t requestedEvent)
{
  if (readEventFromFile(requestedEvent))
  {
    return true;
  }
  else
  {
    return false;
  }
}

PHCompositeNode*
PHNodeIOManager::read(PHCompositeNode* topNode, size_t requestedEvent)
{
  // No tree means we have not yet looked at the file,
  // so we'll reconstruct the node tree now.
  if (!tree)
  {
    topNode = reconstructNodeTree(topNode);
  }

  // If everything worked, there should be a tree now.
  if (tree && readEventFromFile(requestedEvent))
  {
    return topNode;
  }
  else
  {
    return nullptr;
  }
}

void PHNodeIOManager::print() const
{
  if (file)
  {
    if (accessMode == PHReadOnly)
    {
      std::cout << "PHNodeIOManager reading  " << filename << std::endl;
    }
    else
    {
      std::cout << "PHNodeIOManager writing  " << filename << std::endl;
    }
  }
  if (file && tree)
  {
    tree->Print();
  }
  std::cout << "\n\nList of selected objects to read:" << std::endl;
  std::map<std::string, bool>::const_iterator classiter;
  for (classiter = objectToRead.begin(); classiter != objectToRead.end(); ++classiter)
  {
    std::cout << classiter->first << " is set to " << classiter->second << std::endl;
  }
}

std::string
PHNodeIOManager::getBranchClassName(TBranch* branch)
{
  // OK. Here all the game is to find out the name of the type
  // contained in this branch.  In ROOT pre-3.01/05 versions, all
  // branches we used were of the same type = TBranchObject, so that
  // was easy.  Since version 3.01/05 ROOT introduced new branch style
  // with some TBranchElement objects. So far so good.  The problem is
  // that I did not find a common way to grab the typename of the
  // object contained in those branches, so I hereby use some durty if
  // { } else if { } ...

#if ROOT_VERSION_CODE >= ROOT_VERSION(3, 01, 5)
  TBranchElement* be = dynamic_cast<TBranchElement*>(branch);

  if (be)
  {
    // TBranchElement has a nice GetClassName() method for us :
    return be->GetClassName();
  }
#endif

  TBranchObject* bo = dynamic_cast<TBranchObject*>(branch);
  if (bo)
  {
    // For this one we need to go down a little before getting the
    // name...
    TLeafObject* leaf = static_cast<TLeafObject*>(branch->GetLeaf(branch->GetName()));
    assert(leaf != nullptr);
    return leaf->GetTypeName();
  }
  std::cout << PHWHERE << "Fatal error, dynamic cast of TBranchObject failed" << std::endl;
  gSystem->Exit(1);
  exit(1);  // the compiler does not know gSystem->Exit() quits, needs exit to avoid warning
}

bool PHNodeIOManager::readEventFromFile(size_t requestedEvent)
{
  // Se non c'e niente, non possiamo fare niente.  Logisch, n'est ce
  // pas?
  if (!tree)
  {
    PHMessage("PHNodeIOManager::readEventFromFile", PHError,
              "Tree not initialized.");
    return false;
  }

  int bytesRead;

  // Due to the current implementation of TBuffer>>(Long_t) we need
  // to cd() in the current file before trying to fetch any event,
  // otherwise mixing of reading 2.25/03 DST with writing some
  // 3.01/05 trees will fail.
  std::string currdir = gDirectory->GetPath();
  TFile* file_ptr = gFile;  // save current gFile
  file->cd();

  if (requestedEvent)
  {
    if ((bytesRead = tree->GetEvent(requestedEvent)))
    {
      eventNumber = requestedEvent + 1;
    }
  }
  else
  {
    bytesRead = tree->GetEvent(eventNumber++);
  }

  gFile = file_ptr;  // recover gFile
  gROOT->cd(currdir.c_str());

  if (!bytesRead)
  {
    return false;
  }
  if (bytesRead == -1)
  {
    std::cout << PHWHERE << "Error: Input TTree corrupt, exiting now" << std::endl;
    exit(1);
  }
  return true;
}

int PHNodeIOManager::readSpecific(size_t requestedEvent, const std::string& objectName)
{
  // objectName should be one of the valid branch name of the "T" TTree, and
  // should be one of the branches selected by selectObjectToRead() method.
  // No wildcard allowed for the moment.
  std::map<std::string, TBranch*>::const_iterator p = fBranches.find(objectName);

  if (p != fBranches.end())
  {
    TBranch* branch = p->second;
    if (branch)
    {
      return branch->GetEvent(requestedEvent);
    }
  }
  else
  {
    PHMessage("PHNodeIOManager::readSpecific", PHError,
              "Unknown object name");
  }
  return 0;
}

PHCompositeNode*
PHNodeIOManager::reconstructNodeTree(PHCompositeNode* topNode)
{
  if (!file)
  {
    if (filename.empty())
    {
      std::cout << PHWHERE << "filename was never set" << std::endl;
    }
    else
    {
      std::cout << PHWHERE << "TFile " << filename << " NULL pointer" << std::endl;
    }
    return nullptr;
  }

  tree = static_cast<TTree*>(file->Get(TreeName.c_str()));

  if (!tree)
  {
    std::cout << PHWHERE << "PHNodeIOManager::reconstructNodeTree : Root Tree "
         << TreeName << " not found in file " << file->GetName() << std::endl;
    return nullptr;
  }

  // ROOT sucks, we need a unique name for the tree so we can open multiple
  // files. So we take the memory location of the file pointer which
  // should be unique within this process to create it
  std::ostringstream nname;
  nname << TreeName << file;

  tree->SetName(nname.str().c_str());

  // Select the branches according to objectToRead
  std::map<std::string, bool>::const_iterator it;

  if (tree->GetNbranches() > 0)
  {
    for (it = objectToRead.begin(); it != objectToRead.end(); ++it)
    {
      tree->SetBranchStatus((it->first).c_str(),
                            static_cast<bool>(it->second));
    }
  }
  // The file contains a TTree with a list of the TBranchObjects
  // attached to it.
  TObjArray* branchArray = tree->GetListOfBranches();

  // We need these in the loops down below...
  size_t i, j;

  // If a topNode was provided, we can feed the iterator with it.
  if (!topNode)
  {
    topNode = new PHCompositeNode("TOP");  // create topNode if we got a null pointer
  }
  PHNodeIterator nodeIter(topNode);

  // Loop over all branches in the tree. Each branch-name contains the
  // full 'path' of composite-nodes in the original node tree. We
  // split the name and reconstruct the tree.
  std::string delimeters = phooldefs::branchpathdelim + phooldefs::legacypathdelims;  // add old backslash for backward compat
  for (i = 0; i < (size_t)(branchArray->GetEntriesFast()); i++)
  {
    std::string branchname = (*branchArray)[i]->GetName();
    std::vector<std::string> splitvec;
    boost::split(splitvec, branchname, boost::is_any_of(delimeters));
    for (size_t ia = 1; ia < splitvec.size() - 1; ia++)  // -1 so we skip the node name
    {
      if (!nodeIter.cd(splitvec[ia]))
      {
        nodeIter.addNode(new PHCompositeNode(splitvec[ia]));
        nodeIter.cd(splitvec[ia]);
      }
    }
    TBranch* thisBranch = (TBranch*) ((*branchArray)[i]);

    // Skip non-selected branches
    if (thisBranch->TestBit(kDoNotProcess))
    {
      continue;
    }

    std::string branchClassName = getBranchClassName(thisBranch);
    std::string branchName = thisBranch->GetName();
    fBranches[branchName] = thisBranch;

    assert(gROOT != nullptr);
    TClass* thisClass = gROOT->GetClass(branchClassName.c_str());

    if (!thisClass)
    {
      std::cout << PHWHERE << std::endl;
      std::cout << "Missing Class: " << branchClassName << std::endl;
      std::cout << "Did you forget to load the shared library which contains "
           << branchClassName << "?" << std::endl;
    }
    // it does not make sense to continue - the code coredumps
    // later if a class is not loaded
    assert(thisClass != nullptr);

    PHIODataNode<TObject>* newIODataNode =
        static_cast<PHIODataNode<TObject>*>(nodeIter.findFirst("PHIODataNode", (*splitvec.rbegin()).c_str()));
    if (!newIODataNode)
    {
      TObject* newTObject = static_cast<TObject*>(thisClass->New());
      newIODataNode = new PHIODataNode<TObject>(newTObject, (*splitvec.rbegin()).c_str());
      nodeIter.addNode(newIODataNode);
    }
    else
    {
      TObject* oldobject = newIODataNode->getData();
      std::string oldclass = oldobject->ClassName();
      if (oldclass != branchClassName)
      {
        std::cout << "You only have to worry if you get this message when reading parallel files"
             << std::endl
             << "if you get this when opening the 2nd, 3rd,... file" << std::endl
             << "It looks like your objects are not of the same version in these files" << std::endl;
        std::cout << PHWHERE << "Found object " << oldobject->ClassName()
             << " in node tree but the  file "
             << filename << " contains a " << branchClassName
             << " object. The object will be replaced without harming you" << std::endl;
        std::cout << "CAVEAT: If you use local copies of pointers to data nodes" << std::endl
             << "instead of searching the node tree you are in trouble now" << std::endl;
        delete newIODataNode;
        TObject* newTObject = static_cast<TObject*>(thisClass->New());
        newIODataNode = new PHIODataNode<TObject>(newTObject, (*splitvec.rbegin()).c_str());
        nodeIter.addNode(newIODataNode);
      }
    }

    if (thisClass->InheritsFrom("PHObject"))
    {
      newIODataNode->setObjectType("PHObject");
    }
    else
    {
      std::cout << PHWHERE << branchClassName.c_str()
           << " inherits neither from PHTable nor from PHObject"
           << " setting type to PHObject" << std::endl;
      newIODataNode->setObjectType("PHObject");
    }
    thisBranch->SetAddress(&(newIODataNode->data));
    for (j = 1; j < splitvec.size() - 1; j++)
    {
      nodeIter.cd("..");
    }
  }
  return topNode;
}

void PHNodeIOManager::selectObjectToRead(const std::string& objectName, bool readit)
{
  objectToRead[objectName] = readit;

  // If tree is already open, loop over map and set branch status
  if (tree)
  {
    std::map<std::string, bool>::const_iterator it;

    for (it = objectToRead.begin(); it != objectToRead.end(); ++it)
    {
      tree->SetBranchStatus((it->first).c_str(),
                            static_cast<bool>(it->second));
    }
  }
  return;
}

bool PHNodeIOManager::isSelected(const std::string& objectName)
{
  std::map<std::string, TBranch*>::const_iterator p = fBranches.find(objectName);

  if (p != fBranches.end())
  {
    return true;
  }

  return false;
}

bool PHNodeIOManager::SetCompressionLevel(const int level)
{
  if (level < 0)
  {
    return false;
  }
  CompressionLevel = level;
  if (file)
  {
    file->SetCompressionLevel(CompressionLevel);
  }

  return true;
}

double
PHNodeIOManager::GetBytesWritten()
{
  if (file) return file->GetBytesWritten();
  return 0.;
}

std::map<std::string, TBranch*>*
PHNodeIOManager::GetBranchMap()
{
  FillBranchMap();
  return &fBranches;
}

int PHNodeIOManager::FillBranchMap()
{
  if (fBranches.empty())
  {
    TTree* treetmp = static_cast<TTree*>(file->Get(TreeName.c_str()));
    if (treetmp)
    {
      TObjArray* branchArray = treetmp->GetListOfBranches();
      for (size_t i = 0; i < (size_t)(branchArray->GetEntriesFast()); i++)
      {
        TBranch* thisBranch = (TBranch*) ((*branchArray)[i]);
        std::string branchName = (*branchArray)[i]->GetName();
        fBranches[branchName] = thisBranch;
      }
    }
    else
    {
      std::cout << PHWHERE << " No Root Tree " << TreeName
           << " on file " << filename << std::endl;
      return -1;
    }
  }
  return 0;
}

bool PHNodeIOManager::NodeExist(const std::string& nodename)
{
  if (fBranches.empty())
  {
    FillBranchMap();
  }
  std::string delimeters = phooldefs::branchpathdelim + phooldefs::legacypathdelims;  // add old backslash for backward compat
  for (auto & fBranche : fBranches)
  {
    std::vector<std::string> splitvec;
    boost::split(splitvec, fBranche.first, boost::is_any_of(delimeters));
    if (splitvec.back() == nodename)
    {
      return true;
    }
  }
  return false;
}

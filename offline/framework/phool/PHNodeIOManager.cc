//  Implementation of class PHNodeIOManager
//  Author: Matthias Messer

#include "PHNodeIOManager.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "PHNodeIterator.h"
#include "phooldefs.h"

#include <RVersion.h>
#include <TBranch.h>                                         // for TBranch
#include <TBranchObject.h>
#include <TClass.h>
#include <TDirectory.h>                                      // for TDirectory
#include <TFile.h>
#include <TLeafObject.h>
#include <TObject.h>
#include <TObjArray.h>                                       // for TObjArray
#include <TROOT.h>
#include <TTree.h>

// ROOT version taken from RVersion.h
#if ROOT_VERSION_CODE >= ROOT_VERSION(3, 01, 5)
#include <TBranchElement.h>
#endif

#include <boost/algorithm/string.hpp>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

PHNodeIOManager::PHNodeIOManager()
  : file(nullptr)
  , tree(nullptr)
  , TreeName("T")
  , accessMode(PHReadOnly)
  , CompressionLevel(3)
  , isFunctionalFlag(0)
{
}

PHNodeIOManager::PHNodeIOManager(const string& f,
                                 const PHAccessType a)
  : file(nullptr)
  , tree(nullptr)
  , TreeName("T")
  , CompressionLevel(3)
{
  isFunctionalFlag = setFile(f, "titled by PHOOL", a) ? 1 : 0;
}

PHNodeIOManager::PHNodeIOManager(const string& f, const string& title,
                                 const PHAccessType a)
  : file(nullptr)
  , tree(nullptr)
  , TreeName("T")
  , CompressionLevel(3)
{
  isFunctionalFlag = setFile(f, title, a) ? 1 : 0;
}

PHNodeIOManager::PHNodeIOManager(const string& f, const PHAccessType a,
                                 const PHTreeType treeindex)
  : file(nullptr)
  , tree(nullptr)
  , TreeName("T")
  , CompressionLevel(3)
{
  if (treeindex != PHEventTree)
  {
    ostringstream temp;
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

bool PHNodeIOManager::setFile(const string& f, const string& title,
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
    file = 0;
  }
  string currdir = gDirectory->GetPath();
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
    tree = 0;
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

bool PHNodeIOManager::write(TObject** data, const string& path, int buffersize, int splitlevel)
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
    return 0;
  }
}

void PHNodeIOManager::print() const
{
  if (file)
  {
    if (accessMode == PHReadOnly)
    {
      cout << "PHNodeIOManager reading  " << filename << endl;
    }
    else
    {
      cout << "PHNodeIOManager writing  " << filename << endl;
    }
  }
  if (file && tree)
  {
    tree->Print();
  }
  cout << "\n\nList of selected objects to read:" << endl;
  map<string, bool>::const_iterator classiter;
  for (classiter = objectToRead.begin(); classiter != objectToRead.end(); ++classiter)
  {
    cout << classiter->first << " is set to " << classiter->second << endl;
  }
}

string
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
    assert(leaf != 0);
    return leaf->GetTypeName();
  }
  cout << PHWHERE << "Fatal error, dynamic cast of TBranchObject failed" << endl;
  exit(1);
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
  string currdir = gDirectory->GetPath();
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
    cout << PHWHERE << "Error: Input TTree corrupt, exiting now" << endl;
    exit(1);
  }
  return true;
}

int PHNodeIOManager::readSpecific(size_t requestedEvent, const char* objectName)
{
  // objectName should be one of the valid branch name of the "T" TTree, and
  // should be one of the branches selected by selectObjectToRead() method.
  // No wildcard allowed for the moment.
  string name = objectName;
  map<string, TBranch*>::const_iterator p = fBranches.find(name);

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
      cout << PHWHERE << "filename was never set" << endl;
    }
    else
    {
      cout << PHWHERE << "TFile " << filename << " NULL pointer" << endl;
    }
    return nullptr;
  }

  tree = static_cast<TTree*>(file->Get(TreeName.c_str()));

  if (!tree)
  {
    cout << PHWHERE << "PHNodeIOManager::reconstructNodeTree : Root Tree "
         << TreeName << " not found in file " << file->GetName() << endl;
    return nullptr;
  }

  // ROOT sucks, we need a unique name for the tree so we can open multiple
  // files. So we take the memory location of the file pointer which
  // should be unique within this process to create it
  ostringstream nname;
  nname << TreeName << file;

  tree->SetName(nname.str().c_str());

  // Select the branches according to objectToRead
  map<string, bool>::const_iterator it;

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
  string delimeters = phooldefs::branchpathdelim + phooldefs::legacypathdelims;  // add old backslash for backward compat
  for (i = 0; i < (size_t)(branchArray->GetEntriesFast()); i++)
  {
    string branchname = (*branchArray)[i]->GetName();
    vector<string> splitvec;
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

    string branchClassName = getBranchClassName(thisBranch);
    string branchName = thisBranch->GetName();
    fBranches[branchName] = thisBranch;

    assert(gROOT != 0);
    TClass* thisClass = gROOT->GetClass(branchClassName.c_str());

    if (!thisClass)
    {
      cout << PHWHERE << endl;
      cout << "Missing Class: " << branchClassName.c_str() << endl;
      cout << "Did you forget to load the shared library which contains "
           << branchClassName.c_str() << "?" << endl;
    }
    // it does not make sense to continue - the code coredumps
    // later if a class is not loaded
    assert(thisClass != 0);

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
      string oldclass = oldobject->ClassName();
      if (oldclass != branchClassName)
      {
        cout << "You only have to worry if you get this message when reading parallel files"
             << endl
             << "if you get this when opening the 2nd, 3rd,... file" << endl
             << "It looks like your objects are not of the same version in these files" << endl;
        cout << PHWHERE << "Found object " << oldobject->ClassName()
             << " in node tree but the  file "
             << filename << " contains a " << branchClassName
             << " object. The object will be replaced without harming you" << endl;
        cout << "CAVEAT: If you use local copies of pointers to data nodes" << endl
             << "instead of searching the node tree you are in trouble now" << endl;
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
      cout << PHWHERE << branchClassName.c_str()
           << " inherits neither from PHTable nor from PHObject"
           << " setting type to PHObject" << endl;
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

void PHNodeIOManager::selectObjectToRead(const char* objectName, bool readit)
{
  objectToRead[objectName] = readit;

  // If tree is already open, loop over map and set branch status
  if (tree)
  {
    map<string, bool>::const_iterator it;

    for (it = objectToRead.begin(); it != objectToRead.end(); ++it)
    {
      tree->SetBranchStatus((it->first).c_str(),
                            static_cast<bool>(it->second));
    }
  }
  return;
}

bool PHNodeIOManager::isSelected(const char* objectName)
{
  string name = objectName;
  map<string, TBranch*>::const_iterator p = fBranches.find(name);

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

map<string, TBranch*>*
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
        string branchName = (*branchArray)[i]->GetName();
        fBranches[branchName] = thisBranch;
      }
    }
    else
    {
      cout << PHWHERE << " No Root Tree " << TreeName
           << " on file " << filename << endl;
      return -1;
    }
  }
  return 0;
}

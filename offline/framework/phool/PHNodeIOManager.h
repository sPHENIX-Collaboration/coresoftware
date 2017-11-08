#ifndef PHNODEIOMANAGER_H__
#define PHNODEIOMANAGER_H__

//  Declaration of class PHNodeIOManager
//  Purpose: manages file IO for PHIODataNodes
//  Author: Matthias Messer

#include "PHIOManager.h"
#include <string>
#include <vector>
#include <map>


class TObject;
class TFile;
class TTree;
class TBranch;
class PHObject;

class PHNodeIOManager : public PHIOManager { 
public: 
   PHNodeIOManager();
   PHNodeIOManager(const std::string&, const PHAccessType = PHReadOnly);
   PHNodeIOManager(const std::string&, const std::string&, const PHAccessType = PHReadOnly);
   PHNodeIOManager(const std::string&, const PHAccessType , const PHTreeType);
   virtual ~PHNodeIOManager(); 

public:
   virtual void closeFile() ;
   virtual PHBoolean write(PHCompositeNode *) ;   
   virtual void print() const ;

   PHBoolean setFile(const std::string&, const std::string&, const PHAccessType = PHReadOnly) ;
   PHCompositeNode * read(PHCompositeNode * = 0, size_t = 0); 
   PHBoolean read(size_t requestedEvent) ;
   int readSpecific(size_t requestedEvent, const char* objectName) ;
   void selectObjectToRead(const char* objectName, PHBoolean readit) ;
   PHBoolean isSelected(const char* objectName) ;
   int isFunctional() const {return isFunctionalFlag;}
   PHBoolean SetCompressionLevel(const int level);
   double GetBytesWritten();
   std::map<std::string,TBranch*> *GetBranchMap();

public:
   PHBoolean write(TObject**, const std::string&);
private:
   int FillBranchMap();
   PHCompositeNode * reconstructNodeTree(PHCompositeNode *);

   PHBoolean readEventFromFile(size_t requestedEvent);

   //! node = node + readObjectBuffer
   void integrateNodes(PHCompositeNode * node);

   std::string getBranchClassName(TBranch*) ;

  TFile *file;
  TTree *tree;
  std::string TreeName;
  int   bufSize;
  int   split;
  int   accessMode;
  int   CompressionLevel;
  std::map<std::string,TBranch*> fBranches ;
  std::map<std::string,PHBoolean> objectToRead ;

  int isFunctionalFlag;  // flag to tell if that object initialized properly

  //! buffer input objects prior to integrate to the current nodes. This is only used in accessMode = PHReadAndIntegrate
  //! it link TObject * object_on_nodes_to_be_updated -> TObject * object_read_from_DST_file
  std::map<PHObject * , PHObject *> readObjectBuffer;

}; 

#endif /* __PHNODEIOMANAGER_H__ */

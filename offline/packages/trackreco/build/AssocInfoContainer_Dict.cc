// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME AssocInfoContainer_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "../AssocInfoContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_AssocInfoContainer(void *p = 0);
   static void *newArray_AssocInfoContainer(Long_t size, void *p);
   static void delete_AssocInfoContainer(void *p);
   static void deleteArray_AssocInfoContainer(void *p);
   static void destruct_AssocInfoContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AssocInfoContainer*)
   {
      ::AssocInfoContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AssocInfoContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("AssocInfoContainer", ::AssocInfoContainer::Class_Version(), "", 18,
                  typeid(::AssocInfoContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::AssocInfoContainer::Dictionary, isa_proxy, 4,
                  sizeof(::AssocInfoContainer) );
      instance.SetNew(&new_AssocInfoContainer);
      instance.SetNewArray(&newArray_AssocInfoContainer);
      instance.SetDelete(&delete_AssocInfoContainer);
      instance.SetDeleteArray(&deleteArray_AssocInfoContainer);
      instance.SetDestructor(&destruct_AssocInfoContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AssocInfoContainer*)
   {
      return GenerateInitInstanceLocal((::AssocInfoContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AssocInfoContainer*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr AssocInfoContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *AssocInfoContainer::Class_Name()
{
   return "AssocInfoContainer";
}

//______________________________________________________________________________
const char *AssocInfoContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AssocInfoContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int AssocInfoContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AssocInfoContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AssocInfoContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AssocInfoContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AssocInfoContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AssocInfoContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void AssocInfoContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class AssocInfoContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AssocInfoContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(AssocInfoContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AssocInfoContainer(void *p) {
      return  p ? new(p) ::AssocInfoContainer : new ::AssocInfoContainer;
   }
   static void *newArray_AssocInfoContainer(Long_t nElements, void *p) {
      return p ? new(p) ::AssocInfoContainer[nElements] : new ::AssocInfoContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_AssocInfoContainer(void *p) {
      delete ((::AssocInfoContainer*)p);
   }
   static void deleteArray_AssocInfoContainer(void *p) {
      delete [] ((::AssocInfoContainer*)p);
   }
   static void destruct_AssocInfoContainer(void *p) {
      typedef ::AssocInfoContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AssocInfoContainer

namespace ROOT {
   static TClass *multimaplEunsignedsPlongcOunsignedsPintgR_Dictionary();
   static void multimaplEunsignedsPlongcOunsignedsPintgR_TClassManip(TClass*);
   static void *new_multimaplEunsignedsPlongcOunsignedsPintgR(void *p = 0);
   static void *newArray_multimaplEunsignedsPlongcOunsignedsPintgR(Long_t size, void *p);
   static void delete_multimaplEunsignedsPlongcOunsignedsPintgR(void *p);
   static void deleteArray_multimaplEunsignedsPlongcOunsignedsPintgR(void *p);
   static void destruct_multimaplEunsignedsPlongcOunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const multimap<unsigned long,unsigned int>*)
   {
      multimap<unsigned long,unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(multimap<unsigned long,unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("multimap<unsigned long,unsigned int>", -2, "map", 95,
                  typeid(multimap<unsigned long,unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &multimaplEunsignedsPlongcOunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(multimap<unsigned long,unsigned int>) );
      instance.SetNew(&new_multimaplEunsignedsPlongcOunsignedsPintgR);
      instance.SetNewArray(&newArray_multimaplEunsignedsPlongcOunsignedsPintgR);
      instance.SetDelete(&delete_multimaplEunsignedsPlongcOunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_multimaplEunsignedsPlongcOunsignedsPintgR);
      instance.SetDestructor(&destruct_multimaplEunsignedsPlongcOunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< multimap<unsigned long,unsigned int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const multimap<unsigned long,unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *multimaplEunsignedsPlongcOunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const multimap<unsigned long,unsigned int>*)0x0)->GetClass();
      multimaplEunsignedsPlongcOunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void multimaplEunsignedsPlongcOunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_multimaplEunsignedsPlongcOunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) multimap<unsigned long,unsigned int> : new multimap<unsigned long,unsigned int>;
   }
   static void *newArray_multimaplEunsignedsPlongcOunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) multimap<unsigned long,unsigned int>[nElements] : new multimap<unsigned long,unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_multimaplEunsignedsPlongcOunsignedsPintgR(void *p) {
      delete ((multimap<unsigned long,unsigned int>*)p);
   }
   static void deleteArray_multimaplEunsignedsPlongcOunsignedsPintgR(void *p) {
      delete [] ((multimap<unsigned long,unsigned int>*)p);
   }
   static void destruct_multimaplEunsignedsPlongcOunsignedsPintgR(void *p) {
      typedef multimap<unsigned long,unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class multimap<unsigned long,unsigned int>

namespace {
  void TriggerDictionaryInitialization_AssocInfoContainer_Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "AssocInfoContainer_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class AssocInfoContainer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "AssocInfoContainer_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef TRACKRECO_ASSOCINFOCONTAINER_H
#define TRACKRECO_ASSOCINFOCONTAINER_H

#include <phool/PHObject.h>

#include <trackbase/TrkrDefs.h>

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair
#include <vector>                // for vector

class AssocInfoContainer : public PHObject
{
 public:
  typedef std::multimap<TrkrDefs::cluskey, unsigned int> ClusterTrackMap;

  AssocInfoContainer();
  virtual ~AssocInfoContainer();

  void Reset();
  void identify(std::ostream& os = std::cout) const;

  void SetClusterTrackAssoc(const TrkrDefs::cluskey& cluster_id, const unsigned int& track_id)
  {
    _map_cluster_id_track_id.insert(ClusterTrackMap::value_type(cluster_id, track_id));
  }

  std::vector<unsigned int> GetTracksFromCluster(const TrkrDefs::cluskey& cluster_id) const
  {
    std::vector<unsigned int> ret;
    for (auto iter = _map_cluster_id_track_id.lower_bound(cluster_id);
         iter != _map_cluster_id_track_id.upper_bound(cluster_id); ++iter)
    {
      ret.push_back(iter->second);
    }
    return ret;
  }

 private:
  ClusterTrackMap _map_cluster_id_track_id;

  ClassDef(AssocInfoContainer, 1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"AssocInfoContainer", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("AssocInfoContainer_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_AssocInfoContainer_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_AssocInfoContainer_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_AssocInfoContainer_Dict() {
  TriggerDictionaryInitialization_AssocInfoContainer_Dict_Impl();
}

// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HelixHoughBin_Dict

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
#include "../HelixHoughBin.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_HelixHoughBin(void *p);
   static void deleteArray_HelixHoughBin(void *p);
   static void destruct_HelixHoughBin(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HelixHoughBin*)
   {
      ::HelixHoughBin *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HelixHoughBin >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HelixHoughBin", ::HelixHoughBin::Class_Version(), "", 25,
                  typeid(::HelixHoughBin), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HelixHoughBin::Dictionary, isa_proxy, 4,
                  sizeof(::HelixHoughBin) );
      instance.SetDelete(&delete_HelixHoughBin);
      instance.SetDeleteArray(&deleteArray_HelixHoughBin);
      instance.SetDestructor(&destruct_HelixHoughBin);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HelixHoughBin*)
   {
      return GenerateInitInstanceLocal((::HelixHoughBin*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HelixHoughBin*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HelixHoughBin::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HelixHoughBin::Class_Name()
{
   return "HelixHoughBin";
}

//______________________________________________________________________________
const char *HelixHoughBin::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HelixHoughBin::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HelixHoughBin::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HelixHoughBin::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HelixHoughBin::Streamer(TBuffer &R__b)
{
   // Stream an object of class HelixHoughBin.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HelixHoughBin::Class(),this);
   } else {
      R__b.WriteClassBuffer(HelixHoughBin::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HelixHoughBin(void *p) {
      delete ((::HelixHoughBin*)p);
   }
   static void deleteArray_HelixHoughBin(void *p) {
      delete [] ((::HelixHoughBin*)p);
   }
   static void destruct_HelixHoughBin(void *p) {
      typedef ::HelixHoughBin current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HelixHoughBin

namespace {
  void TriggerDictionaryInitialization_HelixHoughBin_Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "HelixHoughBin_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class HelixHoughBin;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "HelixHoughBin_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef G4HOUGH_HELIXHOUGHBIN_H
#define G4HOUGH_HELIXHOUGHBIN_H


#include "HelixHoughSpace.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <set>

class HelixHoughSpace;

// bin has hierarchy, 5 bins for helix parameters

class HelixHoughBin : public PHObject {

public:

  typedef std::set<unsigned int> ClusterSet;
  typedef std::set<unsigned int>::const_iterator ConstClusterIter;
  typedef std::set<unsigned int>::iterator       ClusterIter;


  virtual ~HelixHoughBin() {}

  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "HelixHoughBin base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const {return 0;}
  virtual PHObject* CloneMe() const {return nullptr;}

  virtual void init() {}

//     get_cluster_IDs() {}  
  virtual void add_cluster_ID(unsigned int cluster_ID) {}
  virtual unsigned int get_count() const {return UINT_MAX;}
  virtual void clear_clusters() {}
  virtual bool empty_clusters() {return true;}
  virtual size_t           erase_cluster(unsigned int cluster_id)       {return 0;}
  virtual ConstClusterIter begin_clusters() const                       {return ClusterSet().end();}
  virtual ConstClusterIter find_cluster(unsigned int cluster_id) const  {return ClusterSet().end();}
  virtual ConstClusterIter end_clusters() const                         {return ClusterSet().end();}
  virtual ClusterIter      begin_clusters()                             {return ClusterSet().end();}
  virtual ClusterIter      find_cluster(unsigned int cluster_id)        {return ClusterSet().end();}
  virtual ClusterIter      end_clusters()                               {return ClusterSet().end();}

  virtual unsigned int get_bin(unsigned int zoomlevel) const    {return UINT_MAX;}
  virtual void set_bin(unsigned int zoomlevel, unsigned int bin){}

  virtual unsigned int get_zoomlevel() const {return UINT_MAX;}
  virtual void set_zoomlevel(unsigned int zoomlevel) {}

  virtual unsigned int get_kappa_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_kappa_bin(unsigned int zoomlevel, unsigned int kappa_bin) {}
  virtual unsigned int get_phi_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_phi_bin(unsigned int zoomlevel, unsigned int phi_bin) {}
  virtual unsigned int get_phi_high_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_phi_high_bin(unsigned int zoomlevel) {}
  virtual void set_phi_high_bin(unsigned int zoomlevel, unsigned int phi_high_bin) {}
  virtual unsigned int get_phi_low_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_phi_low_bin(unsigned int zoomlevel) {}
  virtual void set_phi_low_bin(unsigned int zoomlevel, unsigned int phi_low_bin) {}
  virtual unsigned int get_d_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_d_bin(unsigned int zoomlevel, unsigned int d_bin) {}
  virtual unsigned int get_dzdl_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_dzdl_bin(unsigned int zoomlevel, unsigned int dzdl_bin) {}
  virtual unsigned int get_dzdl_high_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_dzdl_high_bin(unsigned int zoomlevel) {}
  virtual void set_dzdl_high_bin(unsigned int zoomlevel, unsigned int dzdl_high_bin) {}
  virtual unsigned int get_dzdl_low_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_dzdl_low_bin(unsigned int zoomlevel) {}
  virtual void set_dzdl_low_bin(unsigned int zoomlevel, unsigned int dzdl_low_bin) {}
  virtual unsigned int  get_z0_bin(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual void set_z0_bin(unsigned int zoomlevel, unsigned int z0_bin) {}

  virtual  void set_hough_space(HelixHoughSpace* hough_space) {};
  virtual  void set_bins(unsigned int zoomlevel, unsigned int bin) {};

  virtual unsigned int get_global_bin(unsigned int zoomlevel) {return UINT_MAX;}
  virtual void set_global_bin(unsigned int zoomlevel) {}

  virtual unsigned int get_neighbors_global_bin(unsigned int zoomlevel, unsigned int var, unsigned int bit_sign) {return UINT_MAX;}

  virtual float get_kappa_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_phi_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_d_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_dzdl_center(unsigned int zoomlevel) {return 999.;}
  virtual float get_z0_center(unsigned int zoomlevel) {return 999.;}

protected:
  HelixHoughBin() {}

  ClassDef(HelixHoughBin,1);
};

#endif


#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HelixHoughBin", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("HelixHoughBin_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_HelixHoughBin_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_HelixHoughBin_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_HelixHoughBin_Dict() {
  TriggerDictionaryInitialization_HelixHoughBin_Dict_Impl();
}

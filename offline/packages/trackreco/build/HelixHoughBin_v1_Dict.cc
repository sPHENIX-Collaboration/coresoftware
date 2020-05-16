// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HelixHoughBin_v1_Dict

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
#include "../HelixHoughBin_v1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_HelixHoughBin_v1(void *p);
   static void deleteArray_HelixHoughBin_v1(void *p);
   static void destruct_HelixHoughBin_v1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HelixHoughBin_v1*)
   {
      ::HelixHoughBin_v1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HelixHoughBin_v1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HelixHoughBin_v1", ::HelixHoughBin_v1::Class_Version(), "", 17,
                  typeid(::HelixHoughBin_v1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HelixHoughBin_v1::Dictionary, isa_proxy, 4,
                  sizeof(::HelixHoughBin_v1) );
      instance.SetDelete(&delete_HelixHoughBin_v1);
      instance.SetDeleteArray(&deleteArray_HelixHoughBin_v1);
      instance.SetDestructor(&destruct_HelixHoughBin_v1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HelixHoughBin_v1*)
   {
      return GenerateInitInstanceLocal((::HelixHoughBin_v1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HelixHoughBin_v1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HelixHoughBin_v1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HelixHoughBin_v1::Class_Name()
{
   return "HelixHoughBin_v1";
}

//______________________________________________________________________________
const char *HelixHoughBin_v1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin_v1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HelixHoughBin_v1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin_v1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HelixHoughBin_v1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin_v1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HelixHoughBin_v1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughBin_v1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HelixHoughBin_v1::Streamer(TBuffer &R__b)
{
   // Stream an object of class HelixHoughBin_v1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HelixHoughBin_v1::Class(),this);
   } else {
      R__b.WriteClassBuffer(HelixHoughBin_v1::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HelixHoughBin_v1(void *p) {
      delete ((::HelixHoughBin_v1*)p);
   }
   static void deleteArray_HelixHoughBin_v1(void *p) {
      delete [] ((::HelixHoughBin_v1*)p);
   }
   static void destruct_HelixHoughBin_v1(void *p) {
      typedef ::HelixHoughBin_v1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HelixHoughBin_v1

namespace ROOT {
   static TClass *setlEunsignedsPintgR_Dictionary();
   static void setlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_setlEunsignedsPintgR(void *p = 0);
   static void *newArray_setlEunsignedsPintgR(Long_t size, void *p);
   static void delete_setlEunsignedsPintgR(void *p);
   static void deleteArray_setlEunsignedsPintgR(void *p);
   static void destruct_setlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<unsigned int>*)
   {
      set<unsigned int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("set<unsigned int>", -2, "set", 90,
                  typeid(set<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(set<unsigned int>) );
      instance.SetNew(&new_setlEunsignedsPintgR);
      instance.SetNewArray(&newArray_setlEunsignedsPintgR);
      instance.SetDelete(&delete_setlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_setlEunsignedsPintgR);
      instance.SetDestructor(&destruct_setlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<unsigned int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const set<unsigned int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const set<unsigned int>*)0x0)->GetClass();
      setlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<unsigned int> : new set<unsigned int>;
   }
   static void *newArray_setlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) set<unsigned int>[nElements] : new set<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEunsignedsPintgR(void *p) {
      delete ((set<unsigned int>*)p);
   }
   static void deleteArray_setlEunsignedsPintgR(void *p) {
      delete [] ((set<unsigned int>*)p);
   }
   static void destruct_setlEunsignedsPintgR(void *p) {
      typedef set<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class set<unsigned int>

namespace {
  void TriggerDictionaryInitialization_HelixHoughBin_v1_Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "HelixHoughBin_v1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class HelixHoughBin_v1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "HelixHoughBin_v1_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef G4HOUGH_HELIXHOUGHBINV1_H
#define G4HOUGH_HELIXHOUGHBINV1_H

#include "HelixHoughBin.h"
#include "HelixHoughSpace.h"

#include <stddef.h>           // for size_t
#include <iostream>           // for cout, ostream

class PHObject;

class HelixHoughBin_v1 : public HelixHoughBin {

public:

  HelixHoughBin_v1(unsigned int bin);
  virtual ~HelixHoughBin_v1() {}

  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const;
  void Reset() {};
  int  isValid() const {return 1;}
  PHObject* CloneMe() const {return new HelixHoughBin_v1(*this);}

  void init();

//     get_cluster_IDs() {}  
  void 		   add_cluster_ID(unsigned int cluster_ID) 	{_cluster_IDs.insert(cluster_ID);} 
  unsigned int 	   get_count() const				{return _cluster_IDs.size();}
  void		   clear_clusters()				{_cluster_IDs.clear();}
  bool 		   empty_clusters()				{return _cluster_IDs.empty();}
  size_t       	   erase_cluster(unsigned int cluster_id)       {return _cluster_IDs.erase(cluster_id);}

  ConstClusterIter begin_clusters() const                    	{return _cluster_IDs.begin();}
  ConstClusterIter find_cluster(unsigned int cluster_id) const  {return _cluster_IDs.find(cluster_id);}
  ConstClusterIter end_clusters() const                      	{return _cluster_IDs.end();}
  ClusterIter      begin_clusters()                          	{return _cluster_IDs.begin();}
  ClusterIter      find_cluster(unsigned int cluster_id)        {return _cluster_IDs.find(cluster_id);}
  ClusterIter      end_clusters()                            	{return _cluster_IDs.end();}

  unsigned int get_bin(unsigned int zoomlevel) const 	{return _bin[zoomlevel];}
  void set_bin(unsigned int zoomlevel, unsigned int bin){_bin[zoomlevel]=bin;}

  unsigned int get_zoomlevel() const 	{return _zoomlevel;}
  void set_zoomlevel(unsigned int zoomlevel) {_zoomlevel = zoomlevel;}

  unsigned int get_kappa_bin(unsigned int zoomlevel) const {return _kappa_bins[zoomlevel];}
  void set_kappa_bin(unsigned int zoomlevel, unsigned int kappa_bin) {_kappa_bins[zoomlevel] = kappa_bin;}
  unsigned int get_phi_bin(unsigned int zoomlevel) const {return _phi_bins[zoomlevel];}
  void set_phi_bin(unsigned int zoomlevel, unsigned int phi_bin) {_phi_bins[zoomlevel] = phi_bin;}
  unsigned int get_phi_high_bin(unsigned int zoomlevel) const {return _phi_high_bins[zoomlevel];}
  void set_phi_high_bin(unsigned int zoomlevel, unsigned int phi_high_bin) {_phi_high_bins[zoomlevel] = phi_high_bin;}
  unsigned int get_phi_low_bin(unsigned int zoomlevel) const {return _phi_low_bins[zoomlevel];}
  void set_phi_low_bin(unsigned int zoomlevel, unsigned int phi_low_bin) {_phi_low_bins[zoomlevel] = phi_low_bin;}
  unsigned int get_d_bin(unsigned int zoomlevel) const {return _d_bins[zoomlevel];}
  void set_d_bin(unsigned int zoomlevel, unsigned int d_bin) {_d_bins[zoomlevel] = d_bin;}
  unsigned int get_dzdl_bin(unsigned int zoomlevel) const {return _dzdl_bins[zoomlevel];}
  void set_dzdl_bin(unsigned int zoomlevel, unsigned int dzdl_bin) {_dzdl_bins[zoomlevel] = dzdl_bin;}
  unsigned int get_dzdl_high_bin(unsigned int zoomlevel) const {return _dzdl_high_bins[zoomlevel];}
  void set_dzdl_high_bin(unsigned int zoomlevel, unsigned int dzdl_high_bin) {_dzdl_high_bins[zoomlevel] = dzdl_high_bin;}
  unsigned int get_dzdl_low_bin(unsigned int zoomlevel) const {return _dzdl_low_bins[zoomlevel];}
  void set_dzdl_low_bin(unsigned int zoomlevel, unsigned int dzdl_low_bin) {_dzdl_low_bins[zoomlevel] = dzdl_low_bin;}
  unsigned int get_z0_bin(unsigned int zoomlevel) const {return _z0_bins[zoomlevel];}
  void set_z0_bin(unsigned int zoomlevel, unsigned int z0_bin) {_z0_bins[zoomlevel] = z0_bin;}


  void set_hough_space(HelixHoughSpace* hough_space);   
  void set_bins(unsigned int zoomlevel, unsigned int bin);

  unsigned int get_global_bin(unsigned int zoomlevel);
  void set_global_bin(unsigned int zoomlevel);
  unsigned int get_neighbors_global_bin(unsigned int zoomlevel, unsigned int var, unsigned int bit_sign);

  float get_kappa_center(unsigned int zoomlevel);
  float get_phi_center(unsigned int zoomlevel);
  float get_d_center(unsigned int zoomlevel);
  float get_dzdl_center(unsigned int zoomlevel);
  float get_z0_center(unsigned int zoomlevel);

private:

  ClusterSet _cluster_IDs;// hits who voted for this bin

  unsigned int _global_bin;
  unsigned int _bin[ZOOMLEVEL_MAX];
  unsigned int _kappa_bins[ZOOMLEVEL_MAX];//={3,2,4} 3rd in most coarse bins, 4th in the narrowest bins 
  unsigned int _phi_bins[ZOOMLEVEL_MAX];
  unsigned int _phi_high_bins[ZOOMLEVEL_MAX];
  unsigned int _phi_low_bins[ZOOMLEVEL_MAX];
  unsigned int _d_bins[ZOOMLEVEL_MAX];
  unsigned int _dzdl_bins[ZOOMLEVEL_MAX];
  unsigned int _dzdl_high_bins[ZOOMLEVEL_MAX];
  unsigned int _dzdl_low_bins[ZOOMLEVEL_MAX];
  unsigned int _z0_bins[ZOOMLEVEL_MAX];
  unsigned int _zoomlevel;


  HelixHoughSpace* _hough_space;

  ClassDef(HelixHoughBin_v1,1);
};

#endif


#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HelixHoughBin_v1", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("HelixHoughBin_v1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_HelixHoughBin_v1_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_HelixHoughBin_v1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_HelixHoughBin_v1_Dict() {
  TriggerDictionaryInitialization_HelixHoughBin_v1_Dict_Impl();
}

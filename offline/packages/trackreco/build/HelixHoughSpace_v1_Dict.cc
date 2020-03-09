// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HelixHoughSpace_v1_Dict

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
#include "../HelixHoughSpace_v1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_HelixHoughSpace_v1(void *p = 0);
   static void *newArray_HelixHoughSpace_v1(Long_t size, void *p);
   static void delete_HelixHoughSpace_v1(void *p);
   static void deleteArray_HelixHoughSpace_v1(void *p);
   static void destruct_HelixHoughSpace_v1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HelixHoughSpace_v1*)
   {
      ::HelixHoughSpace_v1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HelixHoughSpace_v1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HelixHoughSpace_v1", ::HelixHoughSpace_v1::Class_Version(), "", 16,
                  typeid(::HelixHoughSpace_v1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HelixHoughSpace_v1::Dictionary, isa_proxy, 4,
                  sizeof(::HelixHoughSpace_v1) );
      instance.SetNew(&new_HelixHoughSpace_v1);
      instance.SetNewArray(&newArray_HelixHoughSpace_v1);
      instance.SetDelete(&delete_HelixHoughSpace_v1);
      instance.SetDeleteArray(&deleteArray_HelixHoughSpace_v1);
      instance.SetDestructor(&destruct_HelixHoughSpace_v1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HelixHoughSpace_v1*)
   {
      return GenerateInitInstanceLocal((::HelixHoughSpace_v1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HelixHoughSpace_v1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HelixHoughSpace_v1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HelixHoughSpace_v1::Class_Name()
{
   return "HelixHoughSpace_v1";
}

//______________________________________________________________________________
const char *HelixHoughSpace_v1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace_v1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HelixHoughSpace_v1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace_v1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HelixHoughSpace_v1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace_v1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HelixHoughSpace_v1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace_v1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HelixHoughSpace_v1::Streamer(TBuffer &R__b)
{
   // Stream an object of class HelixHoughSpace_v1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HelixHoughSpace_v1::Class(),this);
   } else {
      R__b.WriteClassBuffer(HelixHoughSpace_v1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HelixHoughSpace_v1(void *p) {
      return  p ? new(p) ::HelixHoughSpace_v1 : new ::HelixHoughSpace_v1;
   }
   static void *newArray_HelixHoughSpace_v1(Long_t nElements, void *p) {
      return p ? new(p) ::HelixHoughSpace_v1[nElements] : new ::HelixHoughSpace_v1[nElements];
   }
   // Wrapper around operator delete
   static void delete_HelixHoughSpace_v1(void *p) {
      delete ((::HelixHoughSpace_v1*)p);
   }
   static void deleteArray_HelixHoughSpace_v1(void *p) {
      delete [] ((::HelixHoughSpace_v1*)p);
   }
   static void destruct_HelixHoughSpace_v1(void *p) {
      typedef ::HelixHoughSpace_v1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HelixHoughSpace_v1

namespace {
  void TriggerDictionaryInitialization_HelixHoughSpace_v1_Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "HelixHoughSpace_v1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class HelixHoughSpace_v1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "HelixHoughSpace_v1_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef G4HOUGH_HELIXHOUGHSPACEV1_H
#define G4HOUGH_HELIXHOUGHSPACEV1_H

#include "HelixHoughSpace.h"

#include <iostream>           // for cout, ostream
#include <vector>             // for vector

class PHObject;

class HelixHoughSpace_v1 : public HelixHoughSpace {

public:

  HelixHoughSpace_v1();
  HelixHoughSpace_v1(const HelixHoughSpace_v1& hough_space);
  virtual ~HelixHoughSpace_v1() {};


  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const {};
  void Reset() {}
  int  isValid() const {return 1;}
  PHObject* CloneMe() const {return new HelixHoughSpace_v1(*this);}

  void add_one_zoom(std::vector<unsigned int>& one_zoom);
  unsigned int get_max_zoom();
  void print_zoom_profile();
  void print_para_range();

  void set_kappa_min(float kappa_min) 	{_para_min[0] = kappa_min;}
  float get_kappa_min() const		{return _para_min[0];}
  void set_kappa_max(float kappa_max) 	{_para_max[0] = kappa_max;}
  float get_kappa_max() const 		{return _para_max[0];}
  void set_phi_min(float phi_min) 	{_para_min[1] = phi_min;}
  float get_phi_min() const		{return _para_min[1];}
  void set_phi_max(float phi_max) 	{_para_max[1] = phi_max;}
  float get_phi_max() const 		{return _para_max[1];}
  void set_d_min(float d_min) 		{_para_min[2] = d_min;}
  float get_d_min() const 		{return _para_min[2];}
  void set_d_max(float d_max)		{_para_max[2] = d_max;}
  float get_d_max() const 		{return _para_max[2];}
  void set_dzdl_min(float dzdl_min)	{_para_min[3] = dzdl_min;}
  float get_dzdl_min() const 		{return _para_min[3];}
  void set_dzdl_max(float dzdl_max) 	{_para_max[3] = dzdl_max;}
  float get_dzdl_max() const 		{return _para_max[3];}
  void set_z0_min(float z0_min)         {_para_min[4] = z0_min;}
  float get_z0_min() const 		{return _para_min[4];}
  void set_z0_max(float z0_max)         {_para_max[4] = z0_max;}
  float get_z0_max() const 		{return _para_max[4];}

  unsigned int get_n_kappa_bins(unsigned int zoomlevel) const {return _zoom_profile[zoomlevel][0];}
  unsigned int get_n_phi_bins(unsigned int zoomlevel) const {return _zoom_profile[zoomlevel][1];}
  unsigned int get_n_d_bins(unsigned int zoomlevel) const {return _zoom_profile[zoomlevel][2];}
  unsigned int get_n_dzdl_bins(unsigned int zoomlevel) const {return _zoom_profile[zoomlevel][3];}
  unsigned int get_n_z0_bins(unsigned int zoomlevel) const {return _zoom_profile[zoomlevel][4];}

  float get_kappa_bin_size(unsigned int zoomlevel) const;
  float get_phi_bin_size(unsigned int zoomlevel) const;
  float get_d_bin_size(unsigned int zoomlevel) const;
  float get_dzdl_bin_size(unsigned int zoomlevel) const;
  float get_z0_bin_size(unsigned int zoomlevel) const;
/*
  float get_kappa_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ik) const;
  float get_phi_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ip) const;
  float get_d_center(unsigned int zoomlevel, std::vector<unsigned int>& v_id) const;
  float get_dzdl_center(unsigned int zoomlevel, std::vector<unsigned int>& v_il) const;
  float get_z0_center(unsigned int zoomlevel, std::vector<unsigned int>& v_iz) const;
*/
  unsigned int get_kappa_bin(unsigned int zoomlevel, float kappa) const;
  unsigned int get_phi_bin(unsigned int zoomlevel, float phi) const;
  unsigned int get_d_bin(unsigned int zoomlevel, float d) const;
  unsigned int get_dzdl_bin(unsigned int zoomlevel, float dzdl) const;
  unsigned int get_z0_bin(unsigned int zoomlevel, float z0) const;

  unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins) const;

private:

  float _para_min[5];
  float _para_max[5];

  unsigned int _zoom_profile[ZOOMLEVEL_MAX][5];
  unsigned int _max_zoom;

  ClassDef(HelixHoughSpace_v1,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HelixHoughSpace_v1", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("HelixHoughSpace_v1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_HelixHoughSpace_v1_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_HelixHoughSpace_v1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_HelixHoughSpace_v1_Dict() {
  TriggerDictionaryInitialization_HelixHoughSpace_v1_Dict_Impl();
}

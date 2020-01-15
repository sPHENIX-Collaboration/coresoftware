// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HelixHoughFuncs_v1_Dict

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
#include "../HelixHoughFuncs_v1.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_HelixHoughFuncs_v1(void *p = 0);
   static void *newArray_HelixHoughFuncs_v1(Long_t size, void *p);
   static void delete_HelixHoughFuncs_v1(void *p);
   static void deleteArray_HelixHoughFuncs_v1(void *p);
   static void destruct_HelixHoughFuncs_v1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HelixHoughFuncs_v1*)
   {
      ::HelixHoughFuncs_v1 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HelixHoughFuncs_v1 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HelixHoughFuncs_v1", ::HelixHoughFuncs_v1::Class_Version(), "", 17,
                  typeid(::HelixHoughFuncs_v1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HelixHoughFuncs_v1::Dictionary, isa_proxy, 4,
                  sizeof(::HelixHoughFuncs_v1) );
      instance.SetNew(&new_HelixHoughFuncs_v1);
      instance.SetNewArray(&newArray_HelixHoughFuncs_v1);
      instance.SetDelete(&delete_HelixHoughFuncs_v1);
      instance.SetDeleteArray(&deleteArray_HelixHoughFuncs_v1);
      instance.SetDestructor(&destruct_HelixHoughFuncs_v1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HelixHoughFuncs_v1*)
   {
      return GenerateInitInstanceLocal((::HelixHoughFuncs_v1*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HelixHoughFuncs_v1*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HelixHoughFuncs_v1::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HelixHoughFuncs_v1::Class_Name()
{
   return "HelixHoughFuncs_v1";
}

//______________________________________________________________________________
const char *HelixHoughFuncs_v1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs_v1*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HelixHoughFuncs_v1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs_v1*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HelixHoughFuncs_v1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs_v1*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HelixHoughFuncs_v1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs_v1*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HelixHoughFuncs_v1::Streamer(TBuffer &R__b)
{
   // Stream an object of class HelixHoughFuncs_v1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HelixHoughFuncs_v1::Class(),this);
   } else {
      R__b.WriteClassBuffer(HelixHoughFuncs_v1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HelixHoughFuncs_v1(void *p) {
      return  p ? new(p) ::HelixHoughFuncs_v1 : new ::HelixHoughFuncs_v1;
   }
   static void *newArray_HelixHoughFuncs_v1(Long_t nElements, void *p) {
      return p ? new(p) ::HelixHoughFuncs_v1[nElements] : new ::HelixHoughFuncs_v1[nElements];
   }
   // Wrapper around operator delete
   static void delete_HelixHoughFuncs_v1(void *p) {
      delete ((::HelixHoughFuncs_v1*)p);
   }
   static void deleteArray_HelixHoughFuncs_v1(void *p) {
      delete [] ((::HelixHoughFuncs_v1*)p);
   }
   static void destruct_HelixHoughFuncs_v1(void *p) {
      typedef ::HelixHoughFuncs_v1 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HelixHoughFuncs_v1

namespace {
  void TriggerDictionaryInitialization_HelixHoughFuncs_v1_Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "HelixHoughFuncs_v1_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class HelixHoughFuncs_v1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "HelixHoughFuncs_v1_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef G4HOUGH_HELIXHOUGHFUNCSV1_H
#define G4HOUGH_HELIXHOUGHFUNCSV1_H

#include "HelixHoughFuncs.h"

#include <iostream>           // for cout, ostream
#include <vector>             // for vector

class HelixHoughSpace;
class PHObject;

class HelixHoughFuncs_v1 : public HelixHoughFuncs {

public:

  HelixHoughFuncs_v1();
  HelixHoughFuncs_v1(const HelixHoughFuncs_v1& hough_funcs);
  virtual ~HelixHoughFuncs_v1() {};


  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const {};
  void Reset() {}
  int  isValid() const {return 1;}
  PHObject* CloneMe() const {return new HelixHoughFuncs_v1(*this);}

  void set_current_zoom(unsigned int cur_zoom) { _cur_zoom = cur_zoom;}
  void set_hough_space(HelixHoughSpace* hough_space); 


  void calculate_dzdl_range(float* hitpos3d, std::vector<float>& z0_range, std::vector<float>& kappa_phi_d_ranges, float* dzdl_range);
  void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, float* phi_r_range, float* phi_l_range);
  void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_next_range);
  void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_prev_range, float* phi_next_range);

/*
  unsigned int get_kappa_min() const	{return _para_min[0];}
  unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins);
*/
private:
   
/*
  float _para_min[5];
  float _para_max[5];

  unsigned int _zoom_profile[ZOOMLEVEL_MAX][5];
  unsigned int _max_zoom;
*/
  HelixHoughSpace* _hough_space;
  unsigned int _cur_zoom;

  ClassDef(HelixHoughFuncs_v1,1)
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HelixHoughFuncs_v1", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("HelixHoughFuncs_v1_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_HelixHoughFuncs_v1_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_HelixHoughFuncs_v1_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_HelixHoughFuncs_v1_Dict() {
  TriggerDictionaryInitialization_HelixHoughFuncs_v1_Dict_Impl();
}

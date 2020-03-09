// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HelixHoughFuncs_Dict

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
#include "../HelixHoughFuncs.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_HelixHoughFuncs(void *p);
   static void deleteArray_HelixHoughFuncs(void *p);
   static void destruct_HelixHoughFuncs(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HelixHoughFuncs*)
   {
      ::HelixHoughFuncs *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HelixHoughFuncs >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HelixHoughFuncs", ::HelixHoughFuncs::Class_Version(), "", 18,
                  typeid(::HelixHoughFuncs), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HelixHoughFuncs::Dictionary, isa_proxy, 4,
                  sizeof(::HelixHoughFuncs) );
      instance.SetDelete(&delete_HelixHoughFuncs);
      instance.SetDeleteArray(&deleteArray_HelixHoughFuncs);
      instance.SetDestructor(&destruct_HelixHoughFuncs);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HelixHoughFuncs*)
   {
      return GenerateInitInstanceLocal((::HelixHoughFuncs*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HelixHoughFuncs*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HelixHoughFuncs::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HelixHoughFuncs::Class_Name()
{
   return "HelixHoughFuncs";
}

//______________________________________________________________________________
const char *HelixHoughFuncs::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HelixHoughFuncs::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HelixHoughFuncs::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HelixHoughFuncs::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughFuncs*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HelixHoughFuncs::Streamer(TBuffer &R__b)
{
   // Stream an object of class HelixHoughFuncs.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HelixHoughFuncs::Class(),this);
   } else {
      R__b.WriteClassBuffer(HelixHoughFuncs::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HelixHoughFuncs(void *p) {
      delete ((::HelixHoughFuncs*)p);
   }
   static void deleteArray_HelixHoughFuncs(void *p) {
      delete [] ((::HelixHoughFuncs*)p);
   }
   static void destruct_HelixHoughFuncs(void *p) {
      typedef ::HelixHoughFuncs current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HelixHoughFuncs

namespace {
  void TriggerDictionaryInitialization_HelixHoughFuncs_Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "HelixHoughFuncs_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class HelixHoughFuncs;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "HelixHoughFuncs_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef G4HOUGH_HELIXHOUGHFUNCS_H
#define G4HOUGH_HELIXHOUGHFUNCS_H

#include "HelixHoughSpace_v1.h"

#include <phool/PHObject.h>

#include <climits>
#include <cmath>

class HelixHoughSpace;

class HelixHoughFuncs : public PHObject {

public :
  virtual ~HelixHoughFuncs() {}


  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "HelixHoughFuncs base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int  isValid() const 			{return 0;}
  virtual PHObject* CloneMe() const 	{return nullptr;}

  // Define Hough space for helical tracks 
  virtual void set_current_zoom(unsigned int cur_zoom)		{}
  virtual void set_hough_space(HelixHoughSpace* hough_space)    {}
  virtual void calculate_dzdl_range(float* hitpos3d, std::vector<float>& z0_range, std::vector<float>& kappa_phi_d_ranges, float* dzdl_range) {};
  virtual void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, float* phi_r_range, float* phi_l_range) {};
  virtual void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_next_range){};
  virtual void calculate_phi_range(float* hitpos2d, std::vector<float>& kappa_d_ranges, int helicity, float* phi_range, float* phi_prev_range, float* phi_next_range){};

/*
  virtual void set_z0_max(float z0_max)         {}
  virtual unsigned int get_z0_max() const       {return UINT_MAX;}

  virtual unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins){return -9999;}
*/
protected:
  HelixHoughFuncs(){};
  ClassDef(HelixHoughFuncs,1);
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HelixHoughFuncs", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("HelixHoughFuncs_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_HelixHoughFuncs_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_HelixHoughFuncs_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_HelixHoughFuncs_Dict() {
  TriggerDictionaryInitialization_HelixHoughFuncs_Dict_Impl();
}

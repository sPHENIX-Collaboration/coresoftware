// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HelixHoughSpace_Dict

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
#include "../HelixHoughSpace.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_HelixHoughSpace(void *p);
   static void deleteArray_HelixHoughSpace(void *p);
   static void destruct_HelixHoughSpace(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HelixHoughSpace*)
   {
      ::HelixHoughSpace *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HelixHoughSpace >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HelixHoughSpace", ::HelixHoughSpace::Class_Version(), "", 18,
                  typeid(::HelixHoughSpace), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HelixHoughSpace::Dictionary, isa_proxy, 4,
                  sizeof(::HelixHoughSpace) );
      instance.SetDelete(&delete_HelixHoughSpace);
      instance.SetDeleteArray(&deleteArray_HelixHoughSpace);
      instance.SetDestructor(&destruct_HelixHoughSpace);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HelixHoughSpace*)
   {
      return GenerateInitInstanceLocal((::HelixHoughSpace*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HelixHoughSpace*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HelixHoughSpace::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HelixHoughSpace::Class_Name()
{
   return "HelixHoughSpace";
}

//______________________________________________________________________________
const char *HelixHoughSpace::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HelixHoughSpace::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HelixHoughSpace::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HelixHoughSpace::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HelixHoughSpace*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HelixHoughSpace::Streamer(TBuffer &R__b)
{
   // Stream an object of class HelixHoughSpace.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HelixHoughSpace::Class(),this);
   } else {
      R__b.WriteClassBuffer(HelixHoughSpace::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HelixHoughSpace(void *p) {
      delete ((::HelixHoughSpace*)p);
   }
   static void deleteArray_HelixHoughSpace(void *p) {
      delete [] ((::HelixHoughSpace*)p);
   }
   static void destruct_HelixHoughSpace(void *p) {
      typedef ::HelixHoughSpace current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HelixHoughSpace

namespace {
  void TriggerDictionaryInitialization_HelixHoughSpace_Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "HelixHoughSpace_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class HelixHoughSpace;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "HelixHoughSpace_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#ifndef G4HOUGH_HELIXHOUGHSPACE_H
#define G4HOUGH_HELIXHOUGHSPACE_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>

#ifndef ZOOMLEVEL_MAX
#define ZOOMLEVEL_MAX           5
#endif

class HelixHoughSpace : public PHObject {

public :
  virtual ~HelixHoughSpace() {}


  // The "standard PHObject response" functions...
  virtual void identify(std::ostream &os=std::cout) const {
    os << "HelixHough base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int  isValid() const 			{return 0;}
  virtual PHObject* CloneMe() const 	{return nullptr;}

  // Define Hough space for helical tracks 
  virtual void add_one_zoom(std::vector<unsigned int>& one_zoom) {};
  virtual unsigned int get_max_zoom() 		{return UINT_MAX;}
  virtual void print_zoom_profile()		{}
  virtual void print_para_range()		{}

  virtual void set_kappa_min(float kappa_min)   {}
  virtual float get_kappa_min() const    	{return NAN;}
  virtual void set_kappa_max(float kappa_max)   {}
  virtual float get_kappa_max() const    	{return NAN;}
  virtual void set_phi_min(float phi_min)       {}
  virtual float get_phi_min() const      	{return NAN;} 
  virtual void set_phi_max(float phi_max)       {}
  virtual float get_phi_max() const      	{return NAN;}
  virtual void set_d_min(float d_min)           {}
  virtual float get_d_min() const        	{return NAN;}
  virtual void set_d_max(float d_max)           {}
  virtual float get_d_max() const        	{return NAN;}
  virtual void set_dzdl_min(float dzdl_min)     {}
  virtual float get_dzdl_min() const     	{return NAN;}
  virtual void set_dzdl_max(float dzdl_max)     {}
  virtual float get_dzdl_max() const     	{return NAN;}
  virtual void set_z0_min(float z0_min)         {}
  virtual float get_z0_min() const       	{return NAN;}
  virtual void set_z0_max(float z0_max)         {}
  virtual float get_z0_max() const       	{return NAN;}


  virtual unsigned int get_n_kappa_bins(unsigned int zoomlevel) const {return UINT_MAX;}
  virtual unsigned int get_n_phi_bins(unsigned int zoomlevel) const   {return UINT_MAX;}
  virtual unsigned int get_n_d_bins(unsigned int zoomlevel) const     {return UINT_MAX;}
  virtual unsigned int get_n_dzdl_bins(unsigned int zoomlevel) const  {return UINT_MAX;}
  virtual unsigned int get_n_z0_bins(unsigned int zoomlevel) const    {return UINT_MAX;}

  virtual float get_kappa_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_phi_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_d_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_dzdl_bin_size(unsigned int zoomlevel) const      {return NAN;}
  virtual float get_z0_bin_size(unsigned int zoomlevel) const      {return NAN;}

/*
  virtual float get_kappa_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ik) const {return NAN;}
  virtual float get_phi_center(unsigned int zoomlevel, std::vector<unsigned int>& v_ip) const   {return NAN;}
  virtual float get_d_center(unsigned int zoomlevel, std::vector<unsigned int>& v_id) const     {return NAN;}
  virtual float get_dzdl_center(unsigned int zoomlevel, std::vector<unsigned int>& v_il) const  {return NAN;}
  virtual float get_z0_center(unsigned int zoomlevel, std::vector<unsigned int>& v_iz) const    {return NAN;}
*/
  virtual unsigned int get_kappa_bin(unsigned int zoomlevel, float kappa) const {return UINT_MAX;} 
  virtual unsigned int get_phi_bin(unsigned int zoomlevel, float phi) const   {return UINT_MAX;}
  virtual unsigned int get_d_bin(unsigned int zoomlevel, float d) const     {return UINT_MAX;} 
  virtual unsigned int get_dzdl_bin(unsigned int zoomlevel, float dzdl) const  {return UINT_MAX;}
  virtual unsigned int get_z0_bin(unsigned int zoomlevel, float z0) const    {return UINT_MAX;}


  virtual unsigned int get_bin(unsigned int zoomlevel, unsigned int* bins) const {return UINT_MAX;}

protected:
  HelixHoughSpace(){};
  ClassDef(HelixHoughSpace,1);
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HelixHoughSpace", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("HelixHoughSpace_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_HelixHoughSpace_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_HelixHoughSpace_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_HelixHoughSpace_Dict() {
  TriggerDictionaryInitialization_HelixHoughSpace_Dict_Impl();
}

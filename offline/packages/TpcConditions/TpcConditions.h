#ifndef TPCCONDITIONS_H
#define TPCCONDITIONS_H

#include <phool/PHObject.h>

#include <iostream>

//================================
//  A simple container used by TpcConditionsReco
//  to transmit the conditions for the current event
//  to the consumers of that information.
//
//  This does not require schema evolution since
//  it is not intended for output records.
//
//                 TKH 9/4/2026
//================================

class TpcConditions : public PHObject
{
 public:
  TpcConditions();
  ~TpcConditions() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "TpcConditions object" << std::endl;
  }

  void Reset() override {}
  int isValid() const override { return 0; }

  // float is used throught since telemetric data is of limited precision.
  void set_LoadCurrent(float Value) { m_LoadCurrent = Value; }
  void set_LoadNorth  (float Value) { m_LoadNorth   = Value; }
  void set_LoadSouth  (float Value) { m_LoadSouth   = Value; }
  void set_LoadSR1    (float Value) { m_LoadSR1     = Value; }
  void set_LoadSR2    (float Value) { m_LoadSR2     = Value; }
  void set_LoadSR3    (float Value) { m_LoadSR3     = Value; }
  void set_LoadNR1    (float Value) { m_LoadNR1     = Value; }
  void set_LoadNR2    (float Value) { m_LoadNR2     = Value; }
  void set_LoadNR3    (float Value) { m_LoadNR3     = Value; }
  void set_Temperature(float Value) { m_Temperature = Value; }
  void set_Pressure(float Value) { m_Pressure = Value; }
  void set_FieldOK(bool Value) { m_FieldOK = Value; }
  void set_GainOK(bool Value) { m_GainOK = Value; }

  float get_LoadCurrent() const { return m_LoadCurrent; }
  float get_LoadNorth  () const { return m_LoadNorth  ; }
  float get_LoadSouth  () const { return m_LoadSouth  ; }
  float get_LoadSR1    () const { return m_LoadSR1    ; }
  float get_LoadSR2    () const { return m_LoadSR2    ; }
  float get_LoadSR3    () const { return m_LoadSR3    ; }
  float get_LoadNR1    () const { return m_LoadNR1    ; }
  float get_LoadNR2    () const { return m_LoadNR2    ; }
  float get_LoadNR3    () const { return m_LoadNR3    ; }
  float get_Temperature() const { return m_Temperature; }
  float get_Pressure() const { return m_Pressure; }
  bool get_FieldOK() const { return m_FieldOK; }
  bool get_GainOK() const { return m_GainOK; }

 protected:
  float m_LoadCurrent;
  float m_LoadNorth  ;
  float m_LoadSouth  ;
  float m_LoadSR1    ;
  float m_LoadSR2    ;
  float m_LoadSR3    ;
  float m_LoadNR1    ;
  float m_LoadNR2    ;
  float m_LoadNR3    ;
  float m_Temperature;
  float m_Pressure;
  bool m_FieldOK;
  bool m_GainOK;
};

#endif

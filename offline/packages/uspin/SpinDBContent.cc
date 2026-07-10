#include "SpinDBContent.h"

#include <iostream>

void SpinDBContent::identify(std::ostream &os) const
{
  os << "virtual SpinDBContent object" << std::endl;
}

int SpinDBContent::GetPolarizationBlue(int /*bunch*/, float &value, float &error) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  return -1;
}

int SpinDBContent::GetPolarizationBlue(int /*bunch*/, float &value, float &error, float &syserr) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  syserr = GetErrorValue();
  return -1;
}

int SpinDBContent::GetPolarizationBlue(int /*bunch*/, double &value, double &error) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  return -1;
}

int SpinDBContent::GetPolarizationBlue(int /*bunch*/, double &value, double &error, double &syserr) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  syserr = GetErrorValue();
  return -1;
}

int SpinDBContent::GetPolarizationYellow(int /*bunch*/, float &value, float &error) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  return -1;
}
int SpinDBContent::GetPolarizationYellow(int /*bunch*/, float &value, float &error, float &syserr) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  syserr = GetErrorValue();
  return -1;
}

int SpinDBContent::GetPolarizationYellow(int /*bunch*/, double &value, double &error) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  return -1;
}

int SpinDBContent::GetPolarizationYellow(int /*bunch*/, double &value, double &error, double &syserr) const
{
  value = GetErrorValue();
  error = GetErrorValue();
  syserr = GetErrorValue();
  return -1;
}

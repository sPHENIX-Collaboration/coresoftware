#include "phool.h"

#include <iostream>
#include <string>

void PHMessage(const std::string& functionName, int messageType, const std::string& message)
{
  switch (messageType)
  {
  case (PHError):
    std::cerr << functionName << std::endl;
    std::cerr << "\tERROR" << std::endl;
    std::cerr << "\t" << message << std::endl;
    break;
  case (PHWarning):
    std::cout << functionName << std::endl;
    std::cout << "\tWARNING" << std::endl;
    std::cout << "\t" << message << std::endl;
    break;
  case (PHHullo):
    std::cout << functionName << std::endl;
    std::cout << "\tHULLO HULLO HULLO" << std::endl;
    std::cout << "\t" << message << std::endl;
    break;
  }
}

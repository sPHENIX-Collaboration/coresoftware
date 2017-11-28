#include "phool.h"

#include <iostream>
#include <string>

using namespace std;

void PHMessage(const std::string& functionName, int messageType, const std::string& message)
{
  switch (messageType)
  {
  case (PHError):
    cerr << functionName << endl;
    cerr << "\tERROR" << endl;
    cerr << "\t" << message << endl;
    break;
  case (PHWarning):
    cout << functionName << endl;
    cout << "\tWARNING" << endl;
    cout << "\t" << message << endl;
    break;
  case (PHHullo):
    cout << functionName << endl;
    cout << "\tHULLO HULLO HULLO" << endl;
    cout << "\t" << message << endl;
    break;
  }
}

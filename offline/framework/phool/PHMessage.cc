#include "PHString.h"
#include "phool.h"

#include <iostream>

using namespace std;

void PHMessage(const PHString& functionName, int messageType, const PHString& message)
{
   switch (messageType) {
   case (PHError):
      cerr << functionName << endl;
      cerr << "\tERROR" << endl;
      cerr << "\t" << message << endl;
      break;
   case(PHWarning):
      cout << functionName << endl;
      cout << "\tWARNING" << endl;
      cout << "\t" << message << endl;
      break;
   case(PHHullo):
      cout << functionName << endl;
      cout << "\tHULLO HULLO HULLO" << endl;
      cout << "\t" << message << endl;
      break;
   }
}


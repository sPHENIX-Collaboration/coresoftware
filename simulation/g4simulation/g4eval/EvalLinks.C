
#include "EvalLinks.h"

ClassImp(EvalLinks)

using namespace std;

EvalLinks::EvalLinks(const std::string left_name, std::string right_name) {
}

void EvalLinks::identify(std::ostream& os) const {
  os << "EvalLinks" << endl;
  return;
}

// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HIJINGFLIPAFTERBURNER_H
#define HIJINGFLIPAFTERBURNER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
namespace HepMC
{
  class GenEvent;
}

class HIJINGFlipAfterburner : public SubsysReco
{
public:
    HIJINGFlipAfterburner(const std::string &name = "HIJINGFlipAfterburner");

    ~HIJINGFlipAfterburner() override;

    int process_event(PHCompositeNode *topNode) override;


private:
    void flipZDirection(HepMC::GenEvent *event);
    bool doFlip = true;
};

#endif // HIJINGFLIPAFTERBURNER_H

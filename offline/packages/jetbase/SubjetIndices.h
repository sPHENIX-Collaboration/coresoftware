#ifndef JETBASE_SUBJETINDICES
#define JETBASE_SUBJETINDICES

#include <phool/PHObject.h>


// ---------------------------------------------------------------------------------------
// SubjetIndices class -- used to reference the selected constituents of various jets
// ---------------------------------------------------------------------------------------
class SubjetIndices : public PHObject
{
public:
    SubjetIndices() = default;
    ~SubjetIndices() override = default;

    virtual unsigned int nsubjets() const { return 0; };
    virtual unsigned int index_begin(unsigned int /*which_jet=0*/) const;
    virtual unsigned int index_end  (unsigned int /*which_jet=0*/) const;

    virtual void add_index_pair (unsigned int /**/, unsigned int /**/) {};

private:
  ClassDefOverride(SubjetIndices, 1);
};


#endif // JETBASE_SUBJETINDICES_H

#ifndef G4TPC_DISTORTEDTRACKCONTAINER_H
#define G4TPC_DISTORTEDTRACKCONTAINER_H

#include "DistortionStruct.h"
#include <phool/PHObject.h>                             // for PHObject

 /// track container
  class DistortedTrackContainer: public PHObject
  {

    public:

    /// constructor
    explicit DistortedTrackContainer() = default;

    // destructor
    ~DistortedTrackContainer() {}

    /// copy constructor
    explicit DistortedTrackContainer(const DistortedTrackContainer &) = delete;

    /// assignment operator
    DistortedTrackContainer& operator = ( const DistortedTrackContainer& ) = delete;

    /// reset
    //virtual void Reset();
    void Reset() {}

    /// distortions
    const DistortionStruct::List& distortions() const
    { return _distortions; }
    
    /// add distortion
    void addDistortion( const DistortionStruct& distortion )
    { _distortions.push_back( distortion ); }

    /// reset distortions list
    void clearDistortions()
    { _distortions.clear(); }

    private:

    /// event struct
    DistortionStruct::List _distortions;

    ClassDef(DistortedTrackContainer,1)

  };

#endif // G4TPC_DISTORTEDTRACKCONTAINER_H

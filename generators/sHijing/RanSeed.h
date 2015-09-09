
#ifndef RanSeed_h
#define RanSeed_h

extern "C" { void* ranseed_address_(void); }

/**
@class RanSeed
@brief  Class definition for RanSeed, which is used
      to modify the Hijing RANSEED common.
*/
class RanSeed
{
public:
    RanSeed();
    ~RanSeed();
    
    int&	nseed	(void);
    
private: 

    struct RANSEED;
    friend struct RANSEED;

    struct RANSEED
    {
	int	nseed;
    };

    static RANSEED* _ranseed;
};

// set pointer to zero at start
RanSeed::RANSEED* RanSeed::_ranseed =0;

inline
RanSeed::RanSeed()
{}

inline
RanSeed::~RanSeed()
{}

inline int&
RanSeed::nseed	(void)
{
    if (!_ranseed) _ranseed = static_cast<RANSEED*>(ranseed_address_());
    return _ranseed->nseed;
}

#endif


#ifndef HijCrdn_h
#define HijCrdn_h

extern "C" { void* hijcrdn_address_(void); }
/**
@brief Class definition for HijCrdn, which is used
      to modify the Hijing HIJCRDN common.

*/
class HijCrdn {
public:
    HijCrdn();
    ~HijCrdn();
    
    float&  	yp	(int i, int j);
    float&  	yt	(int i, int j);
    
    void	init	(void);

    // return common array lengths
    inline int	lenI() const {return _lenI;}
    inline int	lenJ() const {return _lenJ;}

private: 

    // Lengths of array in HiMain2 common
    static const int _lenI	= 3;
    static const int _lenJ	= 300;

    struct HIJCRDN;
    friend struct HIJCRDN;

    struct HIJCRDN
    {
	float  	yp	[_lenJ][_lenI];
	float  	yt	[_lenJ][_lenI];
    };

//    int  _dummy;
    float  _realdummy;

    static HIJCRDN* _hijcrdn;
};

// set pointer to zero at start
HijCrdn::HIJCRDN* HijCrdn::_hijcrdn =0;

inline void
HijCrdn::init(void)
{ if (!_hijcrdn) _hijcrdn = static_cast<HIJCRDN*>(hijcrdn_address_()); }

// Constructor
inline
HijCrdn::HijCrdn()
	:// _dummy		(-999),
      _realdummy	(-999.)
{}

// Destructor
inline
HijCrdn::~HijCrdn()
{}

inline float&
HijCrdn::yp	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijcrdn->yp[j-1][i-1];
}

inline float&
HijCrdn::yt	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijcrdn->yt[j-1][i-1];
}

#endif

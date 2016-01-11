
#ifndef HiStrng_h
#define HiStrng_h

extern "C" { void* histrng_address_(void); }
/**
@class HiStrng
@brief       Class definition for HiStrng, which is used
      to modify the Hijing HISTRNG common.
*/
class HiStrng {
public:
    HiStrng();
    ~HiStrng();
    
    int&    	nfp	(int i, int j);
    float&  	pp	(int i, int j);
    int&    	nft	(int i, int j);
    float&  	pt	(int i, int j);
    
    void	init	(void);

    // return common array lengths
    inline int	lenI() const {return _lenI;}
    inline int	lenJ() const {return _lenJ;}

private: 

    // Lengths of array in HiMain2 common
    static const int _lenI	= 300;
    static const int _lenJ	= 15;

    struct HISTRNG;
    friend struct HISTRNG;

    struct HISTRNG
    {
	int	nfp	[_lenJ][_lenI];
	float  	pp	[_lenJ][_lenI];
	int	nft	[_lenJ][_lenI];
	float  	pt	[_lenJ][_lenI];
    };

    int  _dummy;
    float  _realdummy;

    static HISTRNG* _histrng;
};

// set pointer to zero at start
HiStrng::HISTRNG* HiStrng::_histrng =0;

inline void
HiStrng::init(void)
{ if (!_histrng) _histrng = static_cast<HISTRNG*>(histrng_address_()); }

// Constructor
inline
HiStrng::HiStrng()
    : _dummy		(-999),
      _realdummy	(-999.)
{}

// Destructor
inline
HiStrng::~HiStrng()
{}

inline int&
HiStrng::nfp	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _dummy;

    return _histrng->nfp[j-1][i-1];
}

inline float&
HiStrng::pp	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _histrng->pp[j-1][i-1];
}

inline int&
HiStrng::nft	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _dummy;

    return _histrng->nft[j-1][i-1];
}

inline float&
HiStrng::pt	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _histrng->pt[j-1][i-1];
}

#endif

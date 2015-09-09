
#ifndef HiMain2_h
#define HiMain2_h

extern "C" { void* himain2_address_(void); }
/**
@class HiMain2
@brief       Class definition for HiMain2, which is used
      to modify the Hijing HIMAIN2 common.
*/
class HiMain2 {
public:
    HiMain2();
    ~HiMain2();
    
    int&	katt	(int i, int j);
    float&	patt	(int i, int j);
    float&	vatt	(int i, int j);

    void	init	(void);

    // return common array lengths
    inline int	leniKatt() const {return _leniKatt;}
    inline int	lenjKatt() const {return _lenjKatt;}
    inline int	leniPatt() const {return _leniPatt;}
    inline int	lenjPatt() const {return _lenjPatt;}
    inline int	leniVatt() const {return _leniVatt;}
    inline int	lenjVatt() const {return _lenjVatt;}

private: 

    // Lengths of array in HiMain2 common
    static const int _leniKatt	= 130000;
    static const int _lenjKatt	= 4;
    static const int _leniPatt	= 130000;
    static const int _lenjPatt	= 4;
    static const int _leniVatt	= 130000;
    static const int _lenjVatt	= 4;

    struct HIMAIN2;
    friend struct HIMAIN2;

    struct HIMAIN2
    {
	int  	katt[_lenjKatt][_leniKatt];
	float  	patt[_lenjPatt][_leniPatt];
	float  	vatt[_lenjVatt][_leniVatt];
    };

    int  _dummy;
    float  _realdummy;
    
    static HIMAIN2* _himain2;
};

// set pointer to zero at start
HiMain2::HIMAIN2* HiMain2::_himain2 =0;

inline void
HiMain2::init(void)
{ if (!_himain2) _himain2 = static_cast<HIMAIN2*>(himain2_address_()); }

// Constructor
inline
HiMain2::HiMain2()
    : _dummy		(-999),
      _realdummy	(-999.)
{}

// Destructor
inline
HiMain2::~HiMain2()
{}

inline int&
HiMain2::katt	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > leniKatt() ||
	j < 1 || j > lenjKatt() ) return _dummy;

    return _himain2->katt[j-1][i-1];
}

inline float&
HiMain2::patt	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > leniPatt() ||
	j < 1 || j > lenjPatt() ) return _realdummy;

    return _himain2->patt[j-1][i-1];
}

inline float&
HiMain2::vatt	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > leniVatt() ||
	j < 1 || j > lenjVatt() ) return _realdummy;

    return _himain2->vatt[j-1][i-1];
}

#endif

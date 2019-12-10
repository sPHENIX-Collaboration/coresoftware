
#ifndef HiMain1_h
#define HiMain1_h

extern "C" { void* himain1_address_(void); }
/**
@class HiMain1

@brief       Class definition for HiMain1, which is used
      to modify the Hijing HIMAIN1 common.
*/

class HiMain1 {
public:
    HiMain1();
    ~HiMain1();
    
    int&    	natt	(void);
    float&  	eatt	(void);
    int&    	jatt	(void);
    int&    	nt     	(void);
    int&    	np     	(void);
    int&    	n0     	(void);
    int&    	n01    	(void);
    int&    	n10    	(void);
    int&    	n11   	(void);

    //+++BAC 
    // 
    //    Added error status variable to HIMAIN1
    //
    //---BAC
    int&        ierrstat(void);

    void	init	(void);

private: 

    struct HIMAIN1;
    friend struct HIMAIN1;

    struct HIMAIN1
    {
      int    	natt;
      float  	eatt;
      int    	jatt;
      int    	nt;
      int    	np;
      int    	n0;
      int    	n01;
      int    	n10;
      int    	n11;

      //+++BAC 
      // 
      //    Added error status variable to HIMAIN1
      //
      //---BAC

      int       ierrstat;
    };

    static HIMAIN1* _himain1;
};

// set pointer to zero at start
HiMain1::HIMAIN1* HiMain1::_himain1 =0;

inline void
HiMain1::init(void)
{ if (!_himain1) _himain1 = static_cast<HIMAIN1*>(himain1_address_()); }

// Constructor
inline
HiMain1::HiMain1()
{}

// Destructor
inline
HiMain1::~HiMain1()
{}

inline int&
HiMain1::natt	(void)
{
    init();
    return _himain1->natt;
}

inline float&
HiMain1::eatt	(void)
{
    init();
    return _himain1->eatt;
}

inline int&
HiMain1::jatt	(void)
{
    init();
    return _himain1->jatt;
}

inline int&
HiMain1::nt	(void)
{
    init();
    return _himain1->nt;
}

inline int&
HiMain1::np	(void)
{
    init();
    return _himain1->np;
}

inline int&
HiMain1::n0	(void)
{
    init();
    return _himain1->n0;
}

inline int&
HiMain1::n01	(void)
{
    init();
    return _himain1->n01;
}

inline int&
HiMain1::n10	(void)
{
    init();
    return _himain1->n10;
}

inline int&
HiMain1::n11	(void)
{
    init();
    return _himain1->n11;
}

inline int&
HiMain1::ierrstat (void)
{
    init();
    return _himain1->ierrstat;
}

#endif

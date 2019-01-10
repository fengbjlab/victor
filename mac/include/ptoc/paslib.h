#ifndef __PASLIB_H__
#define __PASLIB_H__

#define maxint ((int)(1 << (sizeof(int)*8-1)) - 1)

#ifdef __cplusplus
extern "C" { 
#endif

//double sin(double);
//double cos(double);
//double tan(double);
//double atan(double);

//double exp(double);
//double log(double);

//double sqrt(double);

#ifdef __cplusplus
} 
#endif

/*
 * Macros for some Pascal standard functions
 */


#define ptoc_trunc(x)  ((int)(x))

#define pred(type,x) ((type)((x) - 1))

#define succ(type,x) ((type)((x) + 1))

#define pack(a_l,a_h,a,i,z) memcpy(z, &(a)[(i)-(a_l)], sizeof(z))

#define unpack(z,a_l,a_h,a,i) memcpy(&(a)[(i)-(a_l)], (z), sizeof(z))

#define bitsize(x) (sizeof(x)*8)

#define odd(x) ((x) & 1)

#define chr(n) ((char)(n))

#define ord(c) ((int)(unsigned char)(c))

#ifdef __cplusplus

#ifndef ftoc_abs
#if 0 // abs for ints is defined in stdlib
inline signed char   ftoc_abs(signed char x) { return x < 0 ? -x : x; }
inline short         ftoc_abs(short x)       { return x < 0 ? -x : x; }
inline int           ftoc_abs(int x)         { return x < 0 ? -x : x; }
inline long          ftoc_abs(long x)        { return x < 0 ? -x : x; }
#endif
inline float         ftoc_abs(float x)       { return  x < 0 ? -x : x; }
inline double        ftoc_abs(double x)      { return  x < 0 ? -x : x; }
#endif

inline int           sqr(signed char x)    { return x*x; } 
inline unsigned      sqr(unsigned char x)  { return x*x; } 
inline int           sqr(short x)          { return x*x; } 
inline unsigned      sqr(unsigned short x) { return x*x; } 
inline int           sqr(int x)            { return x*x; } 
inline unsigned      sqr(unsigned x)       { return x*x; } 
inline long          sqr(long x)           { return x*x; } 
inline unsigned long sqr(unsigned long x)  { return x*x; } 
inline float         sqr(float x)          { return x*x; } 
inline double        sqr(double x)         { return x*x; } 

extern "C" void timestamp(int* day,  int* month, int* year, 
			  int* hour, int* min,   int* sec);

extern "C" real get_realtime();

extern "C" int ptoc_round(real);

inline void timestamp(int& day,  int& month, int& year, 
		      int& hour, int& min,   int& sec)
{
    timestamp(&day, &month, &year, &hour, &min, &sec);
}

#else

#ifndef ftoc_abs
#define ftoc_abs(x) ((x) < 0 ? -(x) : (x))
#endif
#define sqr(x) ((x)*(x))

void timestamp(int* day,  int* month, int* year, 
	       int* hour, int* min,   int* sec);

real get_realtime();

int ptoc_round(real);

#endif

#ifdef TURBO_PASCAL

#ifdef USE_CONIO
#include <conio.h>
#define readkey()    _getche()
#define keypressed() _kbhit()
#endif

#define lo(x) ((x) & 0xFF)
#define hi(x) (((x) >> 8) & 0xFF)
#define ptr(seg, ofs) (void*)(((seg) << 16) | (ofs))

char* get_date();
char* get_time();

#define EXTERN extern 

const double pi = 3.14159265359;

#include <ctype.h>

#define ParamCount       paramcount
#define ExitProc         exitproc
#define RandSeed         randseed
#define DirectVideo      directvideo
#define ParamStr(i)      paramstr(i)
#define GetMem(p, size)  getmem(p, size)
#define FreeMem(p, size) freemem(p, size)

#define upcase(c)        ((char)toupper(c))
#define lowcase(c)       ((char)tolower(c))

extern "C" int paramcount;
extern "C" char const* const* param_array;
extern "C" boolean directvideo;
extern "C" pointer exitproc;
extern "C" unsigned randseed;

inline string paramstr(int index) { 
    assert(index <= paramcount);
    return string(param_array[index]);
}

inline void halt(int code = EXIT_FAILURE) { exit(code); }

#define getmem(ptr,size)  ptr = malloc(size)
#define freemem(ptr,size) free(ptr)

extern "C" void     Randomize();
extern "C" unsigned randint(unsigned range);
extern "C" double   randreal();

inline unsigned Random(unsigned range) { return randint(range); }
inline double   Random() { return randreal(); }

#endif

#endif



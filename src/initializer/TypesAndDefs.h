/* Definitions of real and integer variables: 
       if the code should be in double precision, here's 
       the only place where change should be made. 
   Also, some F90 intrinsics are defined here.
 
                      Zarija Lukic, February 2009
                           zarija@lanl.gov
*/

#ifndef TypesDefs_Header_Included
#define TypesDefs_Header_Included

#ifdef DOUBLE_REAL
typedef double real;
#else
typedef float real;
#endif
#ifdef LONG_INTEGER
typedef long integer;
#else
typedef int integer;
#endif

typedef struct{
	double re;
	double im;
} my_fftw_complex;

inline int MOD(int x, int y) { return (x - y*(integer)(x/y));}

#endif

/*
 * Copyright (c) 2003, 2007 Matteo Frigo
 * Copyright (c) 2003, 2007 Massachusetts Institute of Technology
 * Copyright IBM Corp. 2008, 2010
 *
 * The following statement of license applies *only* to header file fftw3.h,
 * fftw2.h, and rfftw.h and *not* to the other program files or any other 
 * files distributed with FFTW or derived therefrom:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef FFTW_H
#define FFTW_H

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#define FFTW_EXTERN extern


/* double precision */


/* If <complex.h> is included, use the C99 complex type.  Otherwise
   define a type bit-compatible with C99 complex */
#if !defined(FFTW_NO_Complex) && defined(_Complex_I) && defined(complex) && defined(I)
  typedef double _Complex fftw_complex;
#else
  typedef double fftw_complex[2];
#endif

typedef struct 
{
  int n;            /* dimension size */
  int is;           /* input stride */
  int os;           /* output stride */
} fftw_iodim;

typedef struct fftw_esv_plan_s *fftw_plan;

FFTW_EXTERN fftw_plan fftw_plan_dft_1d(int n, fftw_complex *in, fftw_complex *out, 
			   int sign,unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_2d(int nx, int ny, 
			   fftw_complex *in, fftw_complex *out, 
			   int sign,unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_3d(int nx, int ny, int nz,
			   fftw_complex *in, fftw_complex *out, 
			   int sign,unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft(int rank, const int *n,
			fftw_complex *in, fftw_complex *out, 
			int sign,unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
			     fftw_complex *in, const int *inembed,
			     int istride, int idist,
			     fftw_complex *out, const int *onembed,
			     int ostride, int odist,  
			     int sign,unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_guru_dft(int rank, const fftw_iodim *dims,
			     int howmany_rank, const fftw_iodim *howmany_dims,
			     fftw_complex *in, fftw_complex *out,
			     int sign,unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_r2c_1d(int n, double *in, 
			       fftw_complex *out, unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_r2c_2d(int nx, int ny, 
			       double *in, fftw_complex *out, 
			       unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_r2c_3d(int nx, int ny, int nz, 
			       double *in, fftw_complex *out, 
			       unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_r2c(int rank, const int *n,
			    double *in, fftw_complex *out, 
			    unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_many_dft_r2c(int rank, const int *n, int howmany,
				 double *in, const int *inembed,
				 int istride, int idist,
				 fftw_complex *out, const int *onembed,
				 int ostride, int odist,  
				 unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_guru_dft_r2c(int rank, const fftw_iodim *dims,
				 int howmany_rank, 
				 const fftw_iodim *howmany_dims,
				 double *in, fftw_complex *out,
				 unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_c2r_1d(int n, fftw_complex *in, 
			       double *out, unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_c2r_2d(int nx, int ny, 
			       fftw_complex *in, double *out, 
			       unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_c2r_3d(int nx, int ny, int nz,
			       fftw_complex *in, double *out, 
			       unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_dft_c2r(int rank, const int *n,
			    fftw_complex *in, double *out, 
			    unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_many_dft_c2r(int rank, const int *n, int howmany,
				 fftw_complex *in, const int *inembed,
				 int istride, int idist,
				 double *out, const int *onembed,
				 int ostride, int odist,  
				 unsigned flags);

FFTW_EXTERN fftw_plan fftw_plan_guru_dft_c2r(int rank, const fftw_iodim *dims,
				 int howmany_rank, 
				 const fftw_iodim *howmany_dims,
				 fftw_complex *in, double *out,
				 unsigned flags);

FFTW_EXTERN void fftw_execute(const fftw_plan p);

FFTW_EXTERN void fftw_execute_dft(const fftw_plan p,
		      fftw_complex *in, fftw_complex *out);

FFTW_EXTERN void fftw_execute_dft_r2c(const fftw_plan p,
			  double *in, fftw_complex *out);

FFTW_EXTERN void fftw_execute_dft_c2r(const fftw_plan p,
			  fftw_complex *in, double *out);

FFTW_EXTERN void fftw_destroy_plan(fftw_plan p);

FFTW_EXTERN void fftw_cleanup(void);

FFTW_EXTERN void *fftw_malloc(size_t size);

FFTW_EXTERN void fftw_free(void *p);


/* single precision */


/* If <complex.h> is included, use the C99 complex type.  Otherwise
   define a type bit-compatible with C99 complex */
#if !defined(FFTW_NO_Complex) && defined(_Complex_I) && defined(complex) && defined(I)
  typedef float _Complex fftwf_complex;
#else
  typedef float fftwf_complex[2];
#endif

typedef struct 
{
  int n;            /* dimension size */
  int is;           /* input stride */
  int os;           /* output stride */
} fftwf_iodim;

typedef struct fftwf_esv_plan_s *fftwf_plan;

FFTW_EXTERN fftwf_plan fftwf_plan_dft_1d(int n, fftwf_complex *in, fftwf_complex *out, 
			   int sign,unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_2d(int nx, int ny, 
			   fftwf_complex *in, fftwf_complex *out, 
			   int sign,unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_3d(int nx, int ny, int nz,
			   fftwf_complex *in, fftwf_complex *out, 
			   int sign,unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft(int rank, const int *n,
			fftwf_complex *in, fftwf_complex *out, 
			int sign,unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_many_dft(int rank, const int *n, int howmany,
			     fftwf_complex *in, const int *inembed,
			     int istride, int idist,
			     fftwf_complex *out, const int *onembed,
			     int ostride, int odist,  
			     int sign,unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_guru_dft(int rank, const fftwf_iodim *dims,
			     int howmany_rank, const fftwf_iodim *howmany_dims,
			     fftwf_complex *in, fftwf_complex *out,
			     int sign,unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_r2c_1d(int n, float *in, 
			       fftwf_complex *out, unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_r2c_2d(int nx, int ny, 
			       float *in, fftwf_complex *out, 
			       unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_r2c_3d(int nx, int ny, int nz, 
			       float *in, fftwf_complex *out, 
			       unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_r2c(int rank, const int *n,
			    float *in, fftwf_complex *out, 
			    unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_many_dft_r2c(int rank, const int *n, int howmany,
				 float *in, const int *inembed,
				 int istride, int idist,
				 fftwf_complex *out, const int *onembed,
				 int ostride, int odist,  
				 unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_guru_dft_r2c(int rank, const fftwf_iodim *dims,
				 int howmany_rank, 
				 const fftwf_iodim *howmany_dims,
				 float *in, fftwf_complex *out,
				 unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_c2r_1d(int n, fftwf_complex *in, 
			       float *out, unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_c2r_2d(int nx, int ny, 
			       fftwf_complex *in, float *out, 
			       unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_c2r_3d(int nx, int ny, int nz,
			       fftwf_complex *in, float *out, 
			       unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_dft_c2r(int rank, const int *n,
			    fftwf_complex *in, float *out, 
			    unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_many_dft_c2r(int rank, const int *n, int howmany,
				 fftwf_complex *in, const int *inembed,
				 int istride, int idist,
				 float *out, const int *onembed,
				 int ostride, int odist,  
				 unsigned flags);

FFTW_EXTERN fftwf_plan fftwf_plan_guru_dft_c2r(int rank, const fftwf_iodim *dims,
				 int howmany_rank, 
				 const fftwf_iodim *howmany_dims,
				 fftwf_complex *in, float *out,
				 unsigned flags);

FFTW_EXTERN void fftwf_execute(const fftwf_plan p);

FFTW_EXTERN void fftwf_execute_dft(const fftwf_plan p,
		      fftwf_complex *in, fftwf_complex *out);

FFTW_EXTERN void fftwf_execute_dft_r2c(const fftwf_plan p,
			  float *in, fftwf_complex *out);

FFTW_EXTERN void fftwf_execute_dft_c2r(const fftwf_plan p,
			  fftwf_complex *in, float *out);

FFTW_EXTERN void fftwf_destroy_plan(fftwf_plan p);

FFTW_EXTERN void fftwf_cleanup(void);

FFTW_EXTERN void *fftwf_malloc(size_t size);

FFTW_EXTERN void fftwf_free(void *p);



#define FFTW_FORWARD (-1)
#define FFTW_BACKWARD (+1)

#define FFTW_NO_TIMELIMIT (-1.0)

/* documented flags */
#define FFTW_MEASURE (0U)
#define FFTW_DESTROY_INPUT (1U << 0)
#define FFTW_UNALIGNED (1U << 1)
#define FFTW_CONSERVE_MEMORY (1U << 2)
#define FFTW_EXHAUSTIVE (1U << 3) /* NO_EXHAUSTIVE is default */
#define FFTW_PRESERVE_INPUT (1U << 4) /* cancels FFTW_DESTROY_INPUT */
#define FFTW_PATIENT (1U << 5) /* IMPATIENT is default */
#define FFTW_ESTIMATE (1U << 6)

/* undocumented beyond-guru flags */
#define FFTW_ESTIMATE_PATIENT (1U << 7)
#define FFTW_BELIEVE_PCOST (1U << 8)
#define FFTW_NO_DFT_R2HC (1U << 9)
#define FFTW_NO_NONTHREADED (1U << 10)
#define FFTW_NO_BUFFERING (1U << 11)
#define FFTW_NO_INDIRECT_OP (1U << 12)
#define FFTW_ALLOW_LARGE_GENERIC (1U << 13) /* NO_LARGE_GENERIC is default */
#define FFTW_NO_RANK_SPLITS (1U << 14)
#define FFTW_NO_VRANK_SPLITS (1U << 15)
#define FFTW_NO_VRECURSE (1U << 16)
#define FFTW_NO_SIMD (1U << 17)
#define FFTW_NO_SLOW (1U << 18)
#define FFTW_NO_FIXED_RADIX_LARGE_N (1U << 19)
#define FFTW_ALLOW_PRUNING (1U << 20)

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTW_H */

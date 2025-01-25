
#ifndef CCOMPLEX
#define CCOMPLEX
#include "SEbasic.h"

typedef struct _complexStruct {
	float r;
	float i;
} se_complex;

/* se_complex a plus se_complex b and return a se_complex structure c=a+b*/
se_complex se_cadd (se_complex a, se_complex b);

/* se_complex a sub se_complex b and return a se_complex structure c=a-b*/
se_complex se_csub (se_complex a, se_complex b);

/* se_complex a multiply se_complex b and return a se_complex structure c=a*b*/
se_complex se_cmul (se_complex a, se_complex b);

/* se_complex a divide se_complex b and return a se_complex structure c=a/b*/
se_complex se_cdiv (se_complex a, se_complex b);


/* *************************************************
cmplx	make a se_complex number from two real numbers
conjg	se_complex conjugate of a se_complex number 
cneg	negate a se_complex number
cinv	invert a se_complex number
csqrt	se_complex square root of a se_complex number
cexp	se_complex exponential of a se_complex number
crmul	multiply a se_complex number by a real number 
rcabs	real magnitude of a se_complex number
************************************************* */

se_complex se_cmplx (float re, float im);

se_complex se_conjg (se_complex z);

se_complex se_cneg (se_complex z);

se_complex se_cinv (se_complex z);

se_complex se_csqrt (se_complex z);

se_complex se_cexp (se_complex z);

se_complex se_crmul (se_complex a, float x);

float se_rcabs (se_complex z);

se_complex *alloc1secomplex (size_t n1);

void free1secomplex (se_complex *p);

se_complex **alloc2secomplex (size_t n1, size_t n2);

void free2secomplex (se_complex **p);

se_complex ***alloc3secomplex (size_t n1, size_t n2, size_t n3);

void free3secomplex (se_complex ***p);

se_complex ****alloc4secomplex (size_t n1, size_t n2, size_t n3, size_t n4);

void free4secomplex (se_complex ****p);

se_complex *****alloc5secomplex (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);

void free5secomplex (se_complex *****p);

#endif





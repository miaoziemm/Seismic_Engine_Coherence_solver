/*****************************************************
Structure:
typedef struct _complexStruct {  se_complex number
	float r,i;
} se_complex;

******************************************************
Function Prototypes:
se_complex cadd (se_complex a, se_complex b);
se_complex csub (se_complex a, se_complex b);
se_complex cmul (se_complex a, se_complex b);
se_complex cdiv (se_complex a, se_complex b);
se_complex cmplx (float re, float im);
se_complex conjg (se_complex z);
se_complex cneg (se_complex z);
se_complex cinv (se_complex z);
se_complex csqrt (se_complex z);
se_complex cexp (se_complex z);
se_complex crmul (se_complex a, float x);
float rcabs (se_complex z);

******************************************************/

/* ***************************************************
Author :     LEON LI and ZHANGJT
Data   :     2005.3.24
		    China University of Petroleum(East China)

All the functions can be found in su software by
cwp, so it can be used freely if you want to use 
these function to program free softare or commercial
softare 
*******************************************************/
#include "../include/SEcomplex.h"
se_complex se_cadd(se_complex a, se_complex b)
{
	se_complex c;
	c.r = a.r+b.r;
	c.i = a.i+b.i;
	return c;
}

se_complex se_csub(se_complex a, se_complex b)
{
	se_complex c;
	c.r = a.r-b.r;
	c.i = a.i-b.i;
	return c;
}

se_complex se_cmul(se_complex a, se_complex b)
{
	se_complex c;
	c.r = a.r*b.r-a.i*b.i;
	c.i = a.i*b.r+a.r*b.i;
	return c;
}

se_complex se_cdiv(se_complex a, se_complex b)
{
	se_complex c;
	float r,den;
	if (fabs(b.r)>=fabs(b.i)) {
		r = b.i/b.r;
		den = b.r+r*b.i;
		c.r = (a.r+r*a.i)/den;
		c.i = (a.i-r*a.r)/den;
	} else {
		r = b.r/b.i;
		den = b.i+r*b.r;
		c.r = (a.r*r+a.i)/den;
		c.i = (a.i*r-a.r)/den;
	}
	return c;
}

se_complex se_cmplx(float re, float im)
{
	se_complex c;
	c.r = re;
	c.i = im;
	return c;
}

se_complex se_conjg(se_complex z)
{
	se_complex c;
	c.r = z.r;
	c.i = -z.i;
	return c;
}

se_complex se_cneg(se_complex z)
{
	se_complex c;
	c.r = -z.r;
	c.i = -z.i;
	return c;
}

se_complex se_cinv(se_complex z)
{
	se_complex c;
	float s;
	s = 1.0/(z.r*z.r+z.i*z.i);
	c.r = z.r*s;
	c.i = -z.i*s;
	return c;
}

se_complex se_csqrt(se_complex z)
{
	se_complex c;
	float x,y,w,r;
	if (z.r==0.0 && z.i==0.0) {
		c.r = c.i = 0.0;
		return c;
	} else {
		x = fabs(z.r);
		y = fabs(z.i);
		if (x>=y) {
			r = y/x;
			w = sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r = x/y;
			w = sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r>=0.0) {
			c.r = w;
			c.i = z.i/(2.0*w);
		} else {
			c.i = (z.i>=0.0) ? w : -w;
			c.r = z.i/(2.0*c.i);
		}
		return c;
	}
}

se_complex se_cexp(se_complex z)
{
	float a;
	se_complex c;
	a = exp(z.r);
	c.r = a*cos(z.i);
	c.i = a*sin(z.i);
	return c;
}

se_complex se_crmul(se_complex a, float x)
{
	se_complex c;
	c.r = x*a.r;
	c.i = x*a.i;
	return c;
}

float se_rcabs(se_complex z)
{
	float x,y,ans,temp;
	x = fabs(z.r);
	y = fabs(z.i);
	if (x==0.0)
		ans = y;
	else if (y==0.0)
		ans = x;
	else if (x>y) {
		temp = y/x;
		ans = x*sqrt(1.0+temp*temp);
	} else {
		temp =x/y;
		ans = y*sqrt(1.0+temp*temp);
	}
	return ans;
}

se_complex *alloc1secomplex (size_t n1){
    return (se_complex*)alloc1(n1,sizeof(se_complex));
}

void free1secomplex (se_complex *p){
    free1(p);
}

se_complex **alloc2secomplex (size_t n1, size_t n2){
    return (se_complex**)alloc2(n1,n2,sizeof(se_complex));
}

void free2secomplex (se_complex **p){
    free2((void**)p);
}

se_complex ***alloc3secomplex (size_t n1, size_t n2, size_t n3){
    return (se_complex***)alloc3(n1,n2,n3,sizeof(se_complex));
}

void free3secomplex (se_complex ***p){
    free3((void***)p);
}

se_complex ****alloc4secomplex (size_t n1, size_t n2, size_t n3, size_t n4){
    return (se_complex****)alloc4(n1,n2,n3,n4,sizeof(se_complex));
}

void free4secomplex (se_complex ****p){
    free4((void****)p);
}

se_complex *****alloc5secomplex (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5){
    return (se_complex*****)alloc5(n1,n2,n3,n4,n5,sizeof(se_complex));
}

void free5secomplex (se_complex *****p){
    free5((void*****)p);
}

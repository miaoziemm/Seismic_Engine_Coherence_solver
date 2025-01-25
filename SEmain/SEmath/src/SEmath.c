#include "../include/SEmath.h"
#include "../include/SEprint.h"

#define CS 362436.0/16777216.0
#define CD 7654321.0/16777216.0
#define CM 16777213.0/16777216.0
#define NBITS 24

/* internal state variables */
static int i=16,j=4;
static float c=CS;
static float u[]={
	0.8668672834288,  0.3697986366357,  0.8008968294805,
	0.4173889774680,  0.8254561579836,  0.9640965269077,
	0.4508667414265,  0.6451309529668,  0.1645456024730,
	0.2787901807898,  0.06761531340295, 0.9663226330820,
	0.01963343943798, 0.02947398211399, 0.1636231515294,
	0.3976343250467,  0.2631008574685
};

void hilbert_by_conv(const float* d, const int n, float*  h) {

    int i;
    for (i=0; i<n; i++) {
        h[i]=0;
        int j;
        for (j=0; j<n; j++) {
            if ( (i-j) & (0x1) ) {
                h[i] += d[j]/(i-j) ;
            }
        }
        h[i]*=(float)(2/PI);
    }
}


// double sinc(double x)
// {
//   static float eps=1e-15;
//   float pix;

//   if(fabs(x)>eps) {
//     pix=PI*x;
//     return sin(pix)/pix;
//   } return 1;

// }

//=============================================================
double bessi0(double x)
/*< Evaluate modified Bessel function In(x) and n=0.  >*/
{
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}

double kaiser_windowed_sinc(double x, double dx, int r)
{
  double xr, w;  
  static float b1[10]={1.24,2.94,4.53,6.31,7.91,9.42,10.95,12.53,14.09,14.18};//monopole
  /* static float b2[10]={0.00,3.33,4.96,6.42,7.77,9.52,11.11,12.52,14.25,16.09};//dipole */
  //float b =b2[r-1];
  //w =(cos(PI*x)-sinc(x))*bessi0(b2_*sx)/(x*i0b2);
  float b = b1[r-1];//r=4
  
  w = 0.;
  xr = fabs(x/(dx*r));
  if(xr<1){
    if(xr>1e-7) w = sinc(x)*bessi0(b*sqrt(1.-xr*xr) )/bessi0(b);
    else w = 1.;
  }
  
  return w;
}

void laplace_filter(grid_cfg *grid, float ***in, float ***out)
{
  int i1, i2 ,i3;
  float diff1, diff2, diff3;
  int i1min = 1;
  int i1max = grid->n1-2;
  int i2min = 1;
  int i2max = grid->n2-2;
  int i3min = (grid->n3>1)?1:0;
  int i3max = (grid->n3>1)?grid->n3-2:0;
  int n123 = grid->n1*grid->n2*grid->n3;

  memset(&out[0][0][0], 0, n123*sizeof(float));
  for(i3=i3min; i3<=i3max; i3++){
    for(i2=i2min; i2<=i2max; i2++){
      for(i1=i1min; i1<=i1max; i1++){
	diff1 = (in[i3][i2][i1+1] - 2.*in[i3][i2][i1] + in[i3][i2][i1-1])/(grid->d1*grid->d1);
	diff2 = (in[i3][i2+1][i1] - 2.*in[i3][i2][i1] + in[i3][i2-1][i1])/(grid->d2*grid->d2);
	diff3 = (grid->n3>1)?(in[i3+1][i2][i1] - 2.*in[i3][i2][i1] + in[i3-1][i2][i1])/(grid->d3*grid->d3):0.;
	out[i3][i2][i1] = diff1 + diff2 + diff3;
      }
    }
  }


}



float franuni (void)
/*****************************************************************************
return a pseudo-random float between 0.0 (inclusive) and 1.0 (exclusive)
******************************************************************************
Returned:	pseudo-random float
*****************************************************************************/
{
	float uni;

	/* basic generator is Fibonacci */
	uni = u[i]-u[j];
	if (uni<0.0) uni += 1.0;
	u[i] = uni;
	i--;
	if (i<0) i = 16;
	j--;
	if (j<0) j = 16;
	
	/* second generator is congruential */
	c -= CD;
	if (c<0.0) c += CM;
	
	/* combination generator */
	uni -= c;
	if (uni<0.0) uni += 1.0;
	return uni;
}

float franuni2(float x1, float x2)
/*****************************************************************************
return a pseudo-random float between x1 (inclusive) and x2 (exclusive)
******************************************************************************
Returned:	pseudo-random float
*****************************************************************************/
{
  return x1 + (x2-x1)*franuni();
}

void sranuni (int seed)
/*****************************************************************************
seed random number generator
******************************************************************************
Input:
seed		different seeds yield different sequences of random numbers.
*****************************************************************************/
{
	int ii,jj,i1,j01,k1,l1,m1;
	float s,t;
	
	/* convert seed to four smallish positive integers */
	i1 = (ABS(seed)%177)+1;
	j01 = (ABS(seed)%167)+1;
	k1 = (ABS(seed)%157)+1;
	l1 = (ABS(seed)%147)+1;
	
	/* generate random bit pattern in array based on given seed */
	for (ii=0; ii<17; ii++) {
		s = 0.0;
		t = 0.5;
		
		/* loop over bits in the float mantissa */
		for (jj=0; jj<NBITS; jj++) {
			m1 = (((i1*j01)%179)*k1)%179;
			i1 = j01;
			j01 = k1;
			k1 = m1;
			l1 = (53*l1+1)%169;
			if (((l1*m1)%64)>=32) s += t;
			t *= 0.5;
		}
		u[ii] = s;
	}
	
	/* initialize generators */
	i = 16;
	j = 4;
	c = CS;
}

float distance_2_dim(int sx, int sy, int gx, int gy, int isx, int isy, int igx, int igy)
{
  return sqrt((float)((sx-isx)*(sx-isx) + (sy-isy)*(sy-isy) + (gx-igx)*(gx-igx) + (gy-igy)*(gy-igy)));
}



/*质因数分解*/
/* The trial divisor increment wheel.  Use it to skip over divisors that
   are composites of 2, 3, 5, 7, or 11.  The part from WHEEL_START up to
   WHEEL_END is reused periodically, while the "lead in" is used to test
   for those primes and to jump onto the wheel.
   The first 4 elements correspond to the incremental offsets of the first 5
   primes (2 3 5 7 11).  The 5(th) element is the difference between that
   last prime and the next largest integer that is not a multiple of
   those primes.  The remaining numbers define the wheel.
   For more information, see

拭除增量轮盘.用它来跳过2 3 5 7 11的合数的除数.当"lead in"用来测试这些质数并在wheel上
跳跃，从WHEEL_START道WHEEL_END的部分被周期性重复使用.
头4个元素对应头五个质数（2 3 5 7 11）的增量。第5个元素是最后一个质数和不是这些质数的
倍数挨着的最大整数的差.剩下的数字定义了这个wheel.
   http://www.utm.edu/research/primes/glossary/WheelFactorization.html */
#define WHEEL_SIZE 5 
static const unsigned int wheel_tab[] = {
    1, 2, 2, 4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4,
    6, 8, 4, 2, 4, 2, 4, 14, 4, 6, 2, 10, 2, 6, 6, 4, 2, 4, 6, 2, 10, 2,
    4, 2, 12, 10, 2, 4, 2, 4, 6, 2, 6, 4, 6, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8,
    4, 2, 4, 6, 8, 6, 10, 2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 6, 10,
    2, 10, 2, 4, 2, 4, 6, 8, 4, 2, 4, 12, 2, 6, 4, 2, 6, 4, 6, 12, 2, 4,
    2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 10, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
    10, 2, 4, 6, 6, 2, 6, 6, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 6, 4,
    8, 6, 4, 6, 2, 4, 6, 8, 6, 4, 2, 10, 2, 6, 4, 2, 4, 2, 10, 2, 10, 2,
    4, 2, 4, 8, 6, 4, 2, 4, 6, 6, 2, 6, 4, 8, 4, 6, 8, 4, 2, 4, 2, 4, 8,
    6, 4, 6, 6, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2, 10, 2, 6,
    4, 6, 2, 6, 4, 2, 4, 6, 6, 8, 4, 2, 6, 10, 8, 4, 2, 4, 2, 4, 8, 10, 6,
    2, 4, 8, 6, 6, 4, 2, 4, 6, 2, 6, 4, 6, 2, 10, 2, 10, 2, 4, 2, 4, 6, 2,
    6, 4, 2, 4, 6, 6, 2, 6, 6, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 8, 4,
    6, 2, 6, 6, 4, 2, 4, 6, 8, 4, 2, 4, 2, 10, 2, 10, 2, 4, 2, 4, 6, 2,
    10, 2, 4, 6, 8, 6, 4, 2, 6, 4, 6, 8, 4, 6, 2, 4, 8, 6, 4, 6, 2, 4, 6,
    2, 6, 6, 4, 6, 6, 2, 6, 6, 4, 2, 10, 2, 10, 2, 4, 2, 4, 6, 2, 6, 4, 2,
    10, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 12, 6, 4, 6, 2, 4, 6, 2,
    12, 4, 2, 4, 8, 6, 4, 2, 4, 2, 10, 2, 10, 6, 2, 4, 6, 2, 6, 4, 2, 4,
    6, 6, 2, 6, 4, 2, 10, 6, 8, 6, 4, 2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6,
    6, 4, 6, 2, 6, 4, 2, 4, 2, 10, 12, 2, 4, 2, 10, 2, 6, 4, 2, 4, 6, 6,
    2, 10, 2, 6, 4, 14, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6, 4, 2,
    4, 6, 2, 6, 4, 2, 4, 12, 2, 12 
};

#define WHEEL_START (wheel_tab + WHEEL_SIZE)
#define WHEEL_END (wheel_tab + (sizeof wheel_tab / sizeof wheel_tab[0]))


int factor(uint64_t n0, uint64_t *factors)
{
    uint64_t n = n0, d, q;
    int n_factors = 0;
    unsigned int const *w = wheel_tab;

    if (n < 1)
        return n_factors;

    /* The exit condition in the following loop is correct because
       any time it is tested one of these 3 conditions holds:
       (1) d divides n
       (2) n is prime
       (3) n is composite but has no factors less than d.
       If (1) or (2) obviously the right thing happens.
       If (3), then since n is composite it is >= d^2. */
    d = 2;
    do {
        q = n / d;
        while (n == q * d) {
            factors[n_factors++] = d;    //后加
            n = q;
            q = n / d;
        }
        d += *(w++);
        if (w == WHEEL_END)       //循环
            w = WHEEL_START;
    } while (d <= q);  //<sqrt(n)

    if (n != 1 || n0 == 1) {
        factors[n_factors++] = n;
    }

    return n_factors;
}

//将number分成nfct个近似相等的因子
//数组 fct  个数 nfct
void split_factors(uint64_t number, uint64_t* fct, size_t nfct)/** Splits number in nfct aproximatively equal factors */
{
    uint64_t factors[MAX_N_FACTORS];
    int i, ifct, n = factor(number, factors);
//初始化 nfct 已经固定大小
    for(ifct = 0; ifct < nfct; ++ifct) {
        fct[ifct] = 1;
    }
//n 是 factors的个数
    for(i = n-1; i >= 0; --i) {
        int ifctmin = 0;      //ifctmin局部变量
        uint64_t min = fct[0]*factors[i];    //局部变量 最小为1*factors[i]
        for(ifct = 1; ifct < nfct; ++ifct) {
            uint64_t k = fct[ifct]*factors[i];
            if(k < min) {
                min = k;
                ifctmin = ifct;  //ifctmin=ifct
            }
        }
        fct[ifctmin] *= factors[i];
    }
}
/** Returns the smallest integer bigger or equals to number that has a
    max prime factor smaller or equal than max_factor. If growth_limit
    is > 0, and no integer smaller than number+growth_limit is found,
    then the one with the smallest bigger prime factor between number
    and number+growth_limit is returned. */
uint64_t next_int_with_max_factor(uint64_t number, uint64_t max_factor, int64_t growth_limit)
{
    uint64_t result, candidate;
    uint64_t min = max_factor;  //max_factor=31

    for(result = candidate = number; 
        growth_limit <= 0 || candidate < number+growth_limit; ++candidate) {

        uint64_t factors[MAX_N_FACTORS];
        int n_factors = factor(candidate, factors);

        if( factors[n_factors-1] <= max_factor ) {  //最大因子小于max_factor 返回最小的大于或者等于number的数，
//且这个数有一个小于或者等于max_factor的质因子 如果growth_limit>0 ，不会有比比number+growth_limit小的数
//被选取。然后，如果number太小 ok 如果 质因子>max_factor 


//返回挨着的 最小的整数 大于或者等于number 且这个整数有一个最大质因数 小于或者等于 max_factor
            result = candidate;
            break;
        }

        if( min > factors[n_factors-1] ) {
            min = factors[n_factors-1];
            result = candidate;
        }
    }

    return result;
}





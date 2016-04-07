/* -------------------------------------------------------------- */
/* C code for Matsumoto & Nishimura - Mersene Twister             */
/*                                                                */
/* This version is compact with a limited API and 3 new functions */
/* 2 functions to save and restore statuses                       */
/* 1 function mtRand optimized for speed (MT is already one if    */
/*   the fastest high quality generator, if not the fastest with  */
/*   MLFG from Michael Mascagni et al                             */
/*                                                                */
/* Updates by D. Hill - March 2016                                */
/* -------------------------------------------------------------- */

extern "C" {

#include "MT.h"

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* ----------------------------------------------------------- */
/* saveStatus         Saves the MT status in a text File       */
/*                                                             */
/* Input: inFileName  The file name                    D. Hill */
/* ----------------------------------------------------------- */

void saveStatus(char * inFileName)
{
    FILE * fp  = fopen(inFileName, "w");
    int    i   = 0;

    fprintf(fp, "%d\n", mti);

    for(i = 0 ; i < N ; i++)
    {
       fprintf(fp, "%ld ", mt[i]);
    }

    fclose(fp);
}

/* ----------------------------------------------------------- */
/* restoreStatus     Restores the MT status from a text File   */
/*                                                             */
/* Input: inFileName  The file name                    D. Hill */
/* ----------------------------------------------------------- */

void restoreStatus(char * inFileName)
{
    FILE * fp  = fopen(inFileName, "r");
    int    i   = 0;

    fscanf(fp, "%d", &mti);

    for(i = 0 ; i < N ; i ++)
    {
       fscanf(fp, "%ld", &mt[i]);
    }

    fclose(fp);
}

/* -------------------------------------------------------------- */
/* mtRand     Fast generation of a double precision pseudo random */
/*            number with the Mersenne Twister algorithm          */
/*            Uniform on [0,1)-real-interval                      */
/*                                                                */
/* This version is embedding the code of genrand_int32 for speed  */
/* and it also avoids the local allocation of y - interesting     */
/* for large simulations when the generator is called very often. */
/* D. Hill - March 2016                                           */
/* -------------------------------------------------------------- */

double mtRand(void)
{
    static unsigned long y;
    /* y is set static for speed - avoids reallocation at each call */
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk = 0; kk < N-M; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < N-1; kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y *(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* --------------------------------------------------- */
/* init_genrand                                        */
/*                                                     */
/* Input: initializes mt[N] with a simple int seed - s */
/*                                                     */
/* This can be confusing since the MT status is huge   */
/* --------------------------------------------------- */

void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* --------------------------------------------------- */
/* init_by_array     initializes MT by an array        */
/*                                                     */
/* Input:  init_key is the array for initializing keys */
/*         key_length is its length                    */
/*         slight changes for C++, 2004/2/26           */
/* --------------------------------------------------- */

void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* -------------------------------------------------------------- */
/* genrand_int32        generates an integer random number        */
/*                                                                */
/* Output: a random number on the [0,0xffffffff]-interval         */
/* -------------------------------------------------------------- */

unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* -------------------------------------------------------------- */
/* genrand_real1   generates an random number in double precision */
/*                                                                */
/* Output: a random number on the [0,1]-real-interval             */
/* -------------------------------------------------------------- */

double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* -------------------------------------------------------------- */
/* genrand_real2   generates an random number in double precision */
/*                                                                */
/* Output: a random number on the [0,1)-real-interval             */
/* -------------------------------------------------------------- */

double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* -------------------------------------------------------------- */

//int main(void)
//{
//    int i;
//    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
//    init_by_array(init, length);


//   saveStatus("status0");
//    printf("\n10^6 outputs of mtRand()\n");

//    for (i = 0; i < 1000000; i++) {
//      mtRand();
//    }
//    saveStatus("status1");

//    printf("\n10 outputs of mtRand()\n");

//    for (i = 0; i < 10; i++) {
//      printf("%10.8f ", mtRand());
//      if (i%5 == 4) printf("\n");
//    }

//    restoreStatus("status0");
//    printf("\n10 outputs of mtRand()\n");

//    for (i = 0; i < 10; i++) {
//      printf("%10.8f ", mtRand());
//      if (i%5 == 4) printf("\n");
//    }

//    restoreStatus("status1");
//    printf("\n10 outputs of mtRand()\n");

//    for (i = 0; i < 10; i++) {
//      printf("%10.8f ", mtRand());
//      if (i%5 == 4) printf("\n");
//    }

// * D. Hill: Remove the comments to checking for numerical
// * reproducibility and portability.
// * The following code produces the output which
// * has to be compared with mt19937ar.out
// * given by Makoto Matsumoto

//    printf("1000 outputs of genrand_int32()\n");
//    for (i=0; i<1000; i++) {
//      printf("%10lu ", genrand_int32());
//      if (i%5==4) printf("\n");
//    }

//    printf("\n1000 outputs of genrand_real2()\n");
//    for (i=0; i<1000; i++) {
//      printf("%10.8f ", genrand_real2());
//      if (i%5==4) printf("\n");
//    }
//*/

//    for(i = 0; i < 1000000000; i++) 
//    {
//       mtRand();
//    }

//    return 0;
//}

/* -------------------------------------------------------------- */

}

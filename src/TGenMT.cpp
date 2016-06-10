

//#define N 624
//#define M 397
//#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
//#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
//#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

#include "TGenMT.hpp"


TGenMT::TGenMT() {
	mti=N+1;
	init_genrand(5489UL);
}

TGenMT::TGenMT(ulong s) {
	mti=N+1;
	init_genrand(s);
}

TGenMT::TGenMT(ulong init_key[], int key_length) {
	mti=N+1;
	init_by_array(init_key, key_length);
}

TGenMT::TGenMT(std::string const filename) {
	restoreStatus( filename.c_str() );
}

/* ----------------------------------------------------------- */
/* saveStatus         Saves the MT status in a text File       */
/*                                                             */
/* Input: inFileName  The file name                    D. Hill */
/* ----------------------------------------------------------- */
void TGenMT::saveStatus(char * inFileName) {
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
void TGenMT::restoreStatus(const char * inFileName) {
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
double TGenMT::mtRand() {
    static ulong y;
    /* y is set static for speed - avoids reallocation at each call */
    static ulong mag01[2]={0x0UL, MATRIX_A};
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
void TGenMT::init_genrand(ulong s) {
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
void TGenMT::init_by_array(ulong init_key[], int key_length) {
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
ulong TGenMT::genrand_int32() {
    ulong y;
    static ulong mag01[2]={0x0UL, MATRIX_A};
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
double TGenMT::genrand_real1() {
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* -------------------------------------------------------------- */
/* genrand_real2   generates an random number in double precision */
/*                                                                */
/* Output: a random number on the [0,1)-real-interval             */
/* -------------------------------------------------------------- */
double TGenMT::genrand_real2(void) {
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}
